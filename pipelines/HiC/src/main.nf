nextflow.enable.dsl=2

include { download_dataset } from './dataset'
include { dim_reduct as snapatac2 } from './software/snapatac2.nf'
include { dim_reduct_higashi as higashi } from './software/higashi.nf'

include { benchmark } from '../../common/benchmark.nf'

workflow {
    data = download_dataset()
    reports = snapatac2(data) | concat(higashi(data))
        | combine(data | map {[it[0], it[1]]}, by: 0)
        | benchmark
}
