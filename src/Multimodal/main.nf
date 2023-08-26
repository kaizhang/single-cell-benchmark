nextflow.enable.dsl=2

params.outdir = "results"

include { highly_variable_genes; scale_features; subset_atac_data } from './preprocess'
include { download_dataset } from '../common/download'

include { dim_reduct as snapatac2 } from './software/snapatac2.nf'
include { dim_reduct as multivi } from './software/multivi.nf'
include { dim_reduct as mofa2 } from './software/mofapy2.nf'
include { dim_reduct as cobolt } from './software/cobolt.nf'
include { dim_reduct as mira } from './software/mira.nf'

include { eval_embedding; output_metrics; plot_metrics; umap; plot_umap
        } from '../common/metrics.nf' params(resultDir: "${params.outdir}/Multiome")
include { json; genBenchId } from '../common/utils.gvy'

workflow bench_multiome {
    take: metadata

    main:
        data = metadata
            | filter {it.assay == "Multiome:RNA" || it.assay == "Multiome:ATAC"}
            | download_dataset
            | branch {
                rna: json(it[0]).assay == "Multiome:RNA" 
                atac: json(it[0]).assay == "Multiome:ATAC" 
            }

        data_rna_hvg = highly_variable_genes(data.rna, 3000)

        multiome_data = data.atac
            | map { [json(it[0]).data_name, it] }
            | combine(
                data_rna_hvg | map({ [json(it[0]).data_name, it] }),
                by: 0
            )
            | map ({ [it[1][0], [it[1][1], it[2][1]]] })

        multiome_data_small = subset_atac_data(data.atac, 200000)
            | map { [json(it[0]).data_name, it] }
            | combine(
                data_rna_hvg | map({ [json(it[0]).data_name, it] }),
                by: 0
            )
            | map ({ [it[1][0], [it[1][1], it[2][1]]] })

        //data_scaled = scale_features(data_hvg)

        embedding = snapatac2(multiome_data).concat(
            mofa2(multiome_data_small),
            cobolt(multiome_data),
            mira(multiome_data),
        )
            | map({ [json(it[0]).data_name, it] })
            | combine(
                data.rna | map({ [json(it[0]).data_name, it[1]] }),
                by: 0
            )
            | map ({ [genBenchId(it[1][0]), it[1][1], it[2]] })

        embedding
            | eval_embedding
            | map { "'" + it + "'" } | toSortedList
            | output_metrics
            | plot_metrics

        embedding
            | umap
            | map { [json(it[0]).data_name, json(it[0]).batch_key, it[1]] }
            | groupTuple(by: [0,1], sort: true)
            | plot_umap
}