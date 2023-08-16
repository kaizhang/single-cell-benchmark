nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct {
    container 'kaizhang/snapatac2:2.3.0'
    tag "${json(metadata).data_name}"
    cpus 4

    when: is_included("pca", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'PCA')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import scanpy as sc
    import numpy as np
    adata = sc.read("data.h5ad")
    sc.tl.pca(adata, n_comps=30)
    np.savetxt("reduced_dim.tsv", adata.obsm['X_pca'], delimiter="\t")
    """
}