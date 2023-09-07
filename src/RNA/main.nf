nextflow.enable.dsl=2

params.outdir = "results"

include { highly_variable_genes; scale_features } from './preprocess'
include { download_dataset } from '../common/download'
include { dim_reduct as scanpy } from './software/scanpy.nf'
include { dim_reduct as snapatac2 } from './software/snapatac2.nf'
include { dim_reduct as scvi } from './software/scvi.nf'

include { eval_embedding; output_metrics; plot_metrics; umap; plot_umap
        } from '../common/metrics.nf' params(resultDir: "${params.outdir}/RNA")

include { json; genBenchId } from '../common/utils.gvy'

workflow bench_rna {
    take: metadata

    main:
        data_full = metadata | filter {it.assay == "RNA" } | download_dataset

        data_hvg = highly_variable_genes(data_full, params.num_hvg)
        data_scaled = scale_features(data_hvg)

        embedding = snapatac2(data_hvg) | concat(
            scvi(data_hvg),
            scanpy(data_hvg | concat(data_scaled)),
        )
            | map({ [json(it[0]).data_name, it] })
            | combine(
                data_full | map({ [json(it[0]).data_name, it[1]] }),
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