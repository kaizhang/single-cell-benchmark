nextflow.enable.dsl=2

include { preproc_snapatac2;
          preproc_archr;
          preproc_snapatac_1;

          dim_reduct_snapatac2;
          dim_reduct_snapatac;
          dim_reduct_archr;

          clust_snapatac;
          clust_snapatac2;
          clust_archr;
        } from './software'
        
include { import_dataset; } from './dataset'

workflow {
    datasets = import_dataset(Channel.of(
        "dataset/fixed_number"
    )) | flatten() | collate(4)

    dim_reduct_snapatac2(preproc_snapatac2(datasets))

    //clust_snapatac(dim_reduct_snapatac(preproc_snapatac_1(datasets)))
    //clust_snapatac2(dim_reduct_snapatac2(preproc_snapatac2(datasets)))
    //clust_archr(dim_reduct_archr(preproc_archr(datasets)))
}
