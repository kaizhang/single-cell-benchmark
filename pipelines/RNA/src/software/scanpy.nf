nextflow.enable.dsl=2

process dim_reduct {
    //container 'kaizhang/snapatac2:1.99.99.7'
    cpus 4
    errorStrategy 'ignore'
    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val("PCA (scanpy)"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import scanpy as sc
    import numpy as np

    adata = sc.read("data.h5ad")
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000) 
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)
    np.savetxt("reduced_dim.tsv", adata.obsm['X_pca'], delimiter="\t")
    """
}