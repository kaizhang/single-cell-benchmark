nextflow.enable.dsl=2

include { preproc_snapatac_2;
          preproc_archr;
          preproc_snapatac_1;

          dim_reduct_snapatac_2;
          dim_reduct_snapatac_1;
          dim_reduct_archr;

          clust_snapatac_1;
          clust_snapatac_2;
          clust_archr;
        } from './software'
        
include { import_dataset; } from './dataset'

workflow {
    datasets = import_dataset(Channel.of(
        "dataset/fixed_number"
    )).flatten()

    //clust_snapatac_1(dim_reduct_snapatac_1(preproc_snapatac_1(datasets)))
    clust_snapatac_2(dim_reduct_snapatac_2(preproc_snapatac_2(datasets)))
    clust_archr(dim_reduct_archr(preproc_archr(datasets)))
}