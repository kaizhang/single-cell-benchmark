nextflow.enable.dsl=2

process dim_reduct_episcanpy {
    container 'kaizhang/episcanpy:0.4.0'
    cpus 4
    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val("epiScanpy (PCA)"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python
    import anndata as ad
    import episcanpy as epi
    import numpy as np
    data = ad.read("data.h5ad")
    data.X.data = data.X.data.astype(np.float64)
    epi.pp.normalize_per_cell(data)
    epi.pp.log1p(data)
    epi.pp.pca(data, n_comps=30)
    np.savetxt("reduced_dim.tsv", data.obsm['X_pca'], delimiter="\t")
    """
}