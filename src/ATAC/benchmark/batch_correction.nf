nextflow.enable.dsl=2

params.outdir = "results"

include { dim_reduct_snapatac as snapatac; } from '../software/snapatac.nf'
include { dim_reduct_cosine as snapatac2 } from '../software/snapatac2.nf'
include { dim_reduct_signac_1 as signac } from '../software/signac.nf'
include { dim_reduct_archr_2 as archr } from '../software/archr.nf'
include { dim_reduct_peakvi as peakvi } from '../software/peakvi.nf'
include { dim_reduct_pycistopic as cisTopic } from '../software/pycistopic.nf'
include { dim_reduct_episcanpy as epiScanpy } from '../software/episcanpy.nf'
include { dim_reduct_scbasset as scBasset } from '../software/scbasset.nf'

include { eval_embedding; output_metrics; plot_metrics
        } from '../../common/metrics.nf' params(resultDir: "${params.outdir}/ATAC/batch_correction")
include { json; genBenchId } from '../../common/utils.gvy'
include { download_genome } from '../../common/download.nf'

process mnc_correct {
    container 'kaizhang/snapatac2:2.3.1'
    tag "${json(metadata).data_name}"
    cpus 4
    errorStrategy 'ignore'

    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple stdout, path("corrected_embed.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import json

    metadata = json.loads('$metadata')
    adata = snap.read("data.h5ad", backed=None)
    batch_key = metadata['batch_key']
    adata.obsm['X_embed'] = np.loadtxt("reduced_dim.tsv", delimiter="\t")

    corrected = snap.pp.mnc_correct(adata, use_rep='X_embed', batch=batch_key, inplace=False)
    np.savetxt("corrected_embed.tsv", corrected, delimiter="\t")

    metadata['method'] += '/mnc'
    print(json.dumps(metadata), end='')
    """
}

process harmony_correct {
    container 'kaizhang/snapatac2:2.3.1'
    tag "${json(metadata).data_name}"
    cpus 4
    errorStrategy 'ignore'

    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple stdout, path("corrected_embed.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import json

    metadata = json.loads('$metadata')
    adata = snap.read("data.h5ad", backed=None)
    batch_key = metadata['batch_key']
    adata.obsm['X_embed'] = np.loadtxt("reduced_dim.tsv", delimiter="\t")

    corrected = snap.pp.harmony(adata, use_rep='X_embed', batch=batch_key, inplace=False)
    np.savetxt("corrected_embed.tsv", corrected, delimiter="\t")

    metadata['method'] += '/harmony'
    print(json.dumps(metadata), end='')
    """
}

workflow bench {
    take: datasets
    main:
        genomes = Channel.of("hg38") | download_genome

        embedding = snapatac2(datasets).concat(
            snapatac(datasets),
            cisTopic(datasets),
            archr(datasets),
            signac(datasets),
            epiScanpy(datasets),
        )
            | map({ [json(it[0]).data_name, it] })
            | combine(
                datasets | map({ [json(it[0]).data_name, it[1]] }),
                by: 0
            )
            | map ({ [it[1][0], it[1][1], it[2]] })

        corrected_embedding = peakvi(datasets).concat(
            scBasset(
                datasets
                    | map ({ [json(it[0]).genome, it[0], it[1]] })
                    | combine(genomes, by: 0)
                    | map({ [it[1], it[2], it[3]] })
            ),
 
            mnc_correct(embedding),
            harmony_correct(embedding),
        )

        metrics = corrected_embedding
            | map({ [json(it[0]).data_name, it] })
            | combine(
                datasets | map({ [json(it[0]).data_name, it[1]] }),
                by: 0
            )
            | map ({ [genBenchId(it[1][0]), it[1][1], it[2]] })
            | eval_embedding

        metrics | map { "'" + it + "'" } | toSortedList | output_metrics | plot_metrics

}