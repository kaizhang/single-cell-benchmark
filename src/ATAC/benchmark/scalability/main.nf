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
          dim_reduct_episcanpy;
          dim_reduct_scbasset;

          clust_snapatac;
          clust_snapatac2;
          clust_archr;
        } from './software'
        
include { download_dataset;
          merge;
          subsample;
          prep_h5ad;
        } from './dataset'
include { download_genome } from '../../../common/download.nf'

process subset_atac_data {
    container 'kaizhang/snapatac2:2.3.1'
    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), path("out.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    adata = snap.read("data.h5ad", backed=None)
    snap.pp.select_features(adata, n_features=500000, filter_lower_quantile=0, filter_upper_quantile=0)
    adata[:, adata.var['selected']].write("out.h5ad", compression="gzip")
    """
}

workflow bench {
    genome = Channel.of("hg38") | download_genome | map { it[1] }
    data = download_dataset() | merge

    bench_data = subsample(
        data,
        Channel.fromList([5000, 10000, 20000, 30000, 40000, 80000, 160000, 200000]),
    )


    preproc_snapatac2(bench_data | filter { it[0] != 30000 && it[0] < 160000 })
    preproc_archr(bench_data | filter { it[0] != 30000 && it[0] < 160000})

    h5ads = bench_data | prep_h5ad
    h5ads_subset = h5ads | subset_atac_data
    dim_reduct_snapatac2(h5ads)
    dim_reduct_archr(h5ads)
    dim_reduct_signac(h5ads)
    dim_reduct_episcanpy(h5ads)
    dim_reduct_snapatac(h5ads | filter { it[0] <= 40000 })
    dim_reduct_pycistopic(h5ads | filter { it[0] <= 80000 })

    dim_reduct_peakvi(h5ads_subset)
    dim_reduct_scale(h5ads_subset)
    dim_reduct_scbasset(h5ads_subset | combine(genome))
}
