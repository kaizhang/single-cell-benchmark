nextflow.enable.dsl=2

include { preproc_snapatac2;
          preproc_archr;
          preproc_snapatac_1;

          dim_reduct_snapatac2;
          dim_reduct_snapatac;
          dim_reduct_archr;
          dim_reduct_pycistopic;
          dim_reduct_peakvi;
          dim_reduct_signac;
          dim_reduct_scale;

          clust_snapatac;
          clust_snapatac2;
          clust_archr;
        } from './software'
        
include { download_dataset;
          merge;
          subsample;
          prep_h5ad;
        } from './dataset'

workflow {
    data = download_dataset() | merge

    bench_data = subsample(
        data,
        Channel.fromList([5000, 10000, 20000, 40000, 80000, 160000]),
    )
        | concat(data)

    //preproc_snapatac2(bench_data)
    //preproc_archr(bench_data)

    h5ads = bench_data | prep_h5ad
    dim_reduct_snapatac2(h5ads)
    dim_reduct_archr(h5ads)
    dim_reduct_signac(h5ads)
    dim_reduct_snapatac(h5ads | filter { it[0] <= 20000 })
    dim_reduct_pycistopic(h5ads | filter { it[0] <= 20000 })
    dim_reduct_peakvi(h5ads | filter {it[0] <= 20000 })
    dim_reduct_scale(h5ads | filter {it[0] <= 20000 })
}
