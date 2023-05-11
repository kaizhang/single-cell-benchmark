nextflow.enable.dsl=2

params.outdir = "results"

include { download_dataset } from './dataset'
include { dim_reduct as scanpy;
          dim_reduct_scaled as scanpy_scaled;
        } from './software/scanpy.nf'
include { dim_reduct as snapatac2
          dim_reduct_svd as snapatac2_svd
        } from './software/snapatac2.nf'

include { bench_embedding
        } from '../common/benchmark.nf' params(resultDir: "${params.outdir}/RNA")

workflow bench_RNA {
    data = download_dataset()
    reports = snapatac2(data) | concat(
        snapatac2_svd(data),
        scanpy(data),
        scanpy_scaled(data),
    )
        | combine(data, by: 0)
        | bench_embedding
}
