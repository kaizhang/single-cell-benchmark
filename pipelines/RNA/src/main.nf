nextflow.enable.dsl=2

include { download_dataset } from './dataset'
include { dim_reduct as scanpy;
          dim_reduct_scaled as scanpy_scaled;
        } from './software/scanpy.nf'
include { dim_reduct as snapatac2
          dim_reduct_svd as snapatac2_svd
        } from './software/snapatac2.nf'

include { benchmark } from '../../common/benchmark.nf'

workflow {
    data = download_dataset()
    reports = snapatac2(data) | concat(
        snapatac2_svd(data),
        scanpy(data),
        scanpy_scaled(data),
    )
        | combine(data, by: 0)
        | benchmark
}
