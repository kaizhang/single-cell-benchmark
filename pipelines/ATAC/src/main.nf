nextflow.enable.dsl=2

include { download_dataset } from './dataset.nf'
include { bench_dim_reduct } from './benchmark/dim_reduct.nf'
include { bench_subsample } from './benchmark/subsample.nf'

workflow {
    datasets = download_dataset()
    bench_dim_reduct(datasets)
    bench_subsample(datasets | map { [it[0], it[1]] })
}