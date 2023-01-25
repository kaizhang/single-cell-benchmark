nextflow.enable.dsl=2

include { import_dataset } from './dataset.nf'
include { bench_dim_reduct } from './benchmark/dim_reduct.nf'

workflow {
    datasets = Channel.fromList([
        tuple("Real", "dataset/Real"),
        tuple("Synthetic", "dataset/Synthetic"),
    ]) | import_dataset | flatten | collate(3)

    bench_dim_reduct(datasets)
}