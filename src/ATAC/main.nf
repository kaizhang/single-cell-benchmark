nextflow.enable.dsl=2

params.outdir = 'results'

include { download_dataset } from '../common/download'
include { bench as dim_reduct } from './benchmark/dim_reduct.nf'
include { bench as subsample } from './benchmark/subsample.nf'
include { bench as clustering } from './benchmark/clustering.nf'
include { bench as leiden } from './benchmark/leiden_resolution.nf'
include { bench as batch_correction } from './benchmark/batch_correction.nf'

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

        dim_reduct(data.without_batch)
        clustering(data.without_batch)
        //bench_subsample(datasets | map { [it[0], it[1]] })
        //bench_leiden(data_full)
        batch_correction(data.with_batch)
}

workflow bench_atac_simulated {
    datasets = download_simulated_dataset()
    bench_dim_reduct(datasets)
}