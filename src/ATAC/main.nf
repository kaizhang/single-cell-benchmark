nextflow.enable.dsl=2

params.outdir = 'results'

include { download_dataset } from '../common/download'
include { bench as dim_reduct } from './benchmark/dim_reduct.nf'
include { bench as subsample } from './benchmark/subsample.nf'
include { bench as clustering } from './benchmark/clustering.nf'
include { bench as leiden } from './benchmark/leiden_resolution.nf'
include { bench as batch_correction } from './benchmark/batch_correction.nf'
include { bench as scalability } from './benchmark/scalability/main.nf'

include { json } from '../common/utils.gvy'

workflow bench_atac {
    take: metadata

    main:
        data = metadata
            | filter {it.assay == "ATAC" }
            | download_dataset
            | branch {
                without_batch: !json(it[0]).containsKey('batch_key')
                with_batch: json(it[0]).containsKey('batch_key')
            }

        if (params.components == null || params.components.contains('dim_reduct')) {
            dim_reduct(data.without_batch)
        }

        if (params.components == null || params.components.contains('batch_correction')) {
            batch_correction(data.with_batch)
        }

        if (params.components == null || params.components.contains('clustering')) {
            clustering(data.without_batch)
        }

        if (params.components == null || params.components.contains('scalability')) {
            scalability()
        }

        //bench_subsample(datasets | map { [it[0], it[1]] })
        //bench_leiden(data_full)
}

workflow bench_atac_simulated {
    datasets = download_simulated_dataset()
    bench_dim_reduct(datasets)
}