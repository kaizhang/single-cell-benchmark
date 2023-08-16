nextflow.enable.dsl=2

params.outdir = "results"

include { dim_reduct_snapatac as snapatac; } from '../software/snapatac.nf'
include { dim_reduct_jaccard as snapatac2_jaccard;
          dim_reduct_cosine as snapatac2_cosine;
          dim_reduct_svd as snapatac2_svd;
        } from '../software/snapatac2.nf'
include { dim_reduct_signac_1 as signac_1;
          dim_reduct_signac_2 as signac_2;
          dim_reduct_signac_3 as signac_3;
          dim_reduct_signac_4 as signac_4;
        } from '../software/signac.nf'
include { dim_reduct_archr_1 as archr_1;
          dim_reduct_archr_2 as archr_2;
          dim_reduct_archr_3 as archr_3;
        } from '../software/archr.nf'
include { dim_reduct_peakvi as peakvi } from '../software/peakvi.nf'
include { dim_reduct_scale as scale } from '../software/scale.nf'
include { dim_reduct_pycistopic as cisTopic } from '../software/pycistopic.nf'
include { dim_reduct_scbasset as scBasset } from '../software/scbasset.nf'
include { dim_reduct_episcanpy as epiScanpy } from '../software/episcanpy.nf'

include { bench_embedding
        } from '../../common/benchmark.nf' params(resultDir: "${params.outdir}/ATAC/dim_reduct")
include { json; genBenchId } from '../../common/utils.gvy'

process download_genome {
    container 'kaizhang/scatac-bench:0.2.0'
    input:
      val(genome)

    output:
      tuple val(genome), path("*.decomp")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import os
    os.environ["SNAP_DATA_DIR"] = "./"

    if "${genome}" == "mm10":
        snap.genome.mm10.fetch_fasta()
    elif "${genome}" == "hg38":
        snap.genome.hg38.fetch_fasta()
    elif "${genome}" == "hg19":
        snap.genome.hg19.fetch_fasta()
    else:
        raise ValueError("Unknown genome: ${genome}")
    """
}

workflow bench {
    take: datasets
    main:
        genomes = Channel.of("mm10", "hg38", "hg19") | download_genome

        embedding = snapatac2_cosine(datasets) | concat(
            snapatac2_jaccard(datasets),
            snapatac2_svd(datasets),

            snapatac(datasets),

            epiScanpy(datasets),

            signac_1(datasets),
            signac_2(datasets),
            signac_3(datasets),
            signac_4(datasets),

            archr_1(datasets),
            archr_2(datasets),
            archr_3(datasets),

            cisTopic(datasets),

            scale(datasets),

            peakvi(datasets),

            scBasset(
                datasets
                    | map ({ [json(it[0]).genome, it[0], it[1]] })
                    | combine(genomes, by: 0)
                    | map({ [it[1], it[2], it[3]] })
            )
        )

        embedding
            | map({ [json(it[0]).data_name, it] })
            | combine(
                datasets | map({ [json(it[0]).data_name, it[1]] }),
                by: 0
            )
            | map ({ [genBenchId(it[1][0]), it[1][1], it[2]] })
            | bench_embedding
}