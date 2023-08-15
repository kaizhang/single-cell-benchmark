nextflow.enable.dsl=2

params.outdir = 'results'

include { download_dataset } from '../common/download'
include { bench_dim_reduct } from './benchmark/dim_reduct.nf'
include { bench_subsample } from './benchmark/subsample.nf'
include { bench_clustering } from './benchmark/clustering.nf'
include { bench_leiden } from './benchmark/leiden.nf'

workflow bench_atac {
    take: metadata

    main:
        data_full = metadata | filter {it.type == "ATAC" } | download_dataset

        bench_dim_reduct(data_full)
        //bench_subsample(datasets | map { [it[0], it[1]] })
        //bench_clustering(datasets)
        //bench_leiden(datasets)
}

workflow bench_atac_simulated {
    datasets = download_simulated_dataset()
    bench_dim_reduct(datasets)
}