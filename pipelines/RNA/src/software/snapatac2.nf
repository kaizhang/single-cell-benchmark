nextflow.enable.dsl=2

process dim_reduct {
    container 'kaizhang/snapatac2:2.3.0'
    cpus 4
    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val("SnapATAC2"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import scanpy as sc
    import snapatac2 as snap
    import numpy as np

    adata = sc.read("data.h5ad")
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=5000)
    adata = adata[:, adata.var.highly_variable]
    snap.tl.spectral(adata, features=None)
    np.savetxt("reduced_dim.tsv", adata.obsm['X_spectral'], delimiter="\t")
    """
}

process dim_reduct_svd {
    container 'kaizhang/snapatac2:2.3.0'
    cpus 4
    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val("SnapATAC2 (svd)"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import scanpy as sc
    import snapatac2 as snap
    import scipy as sp
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=5000)
    adata = adata[:, adata.var.highly_variable]
    feature_weights = snap.tl._embedding.idf(adata)
    X = adata.X @ sp.sparse.diags(feature_weights)
    s = 1 / np.sqrt(np.ravel(sp.sparse.csr_matrix.power(X, 2).sum(axis = 1)))
    X = sp.sparse.diags(s) @ X
    D = np.ravel(X @ X.sum(axis = 0).T)
    X = sp.sparse.diags(1 / np.sqrt(D)) @ X
    evecs,evals,_ = sp.sparse.linalg.svds(X, k=30)
    ix = evals.argsort()[::-1]
    evals = evals[ix]
    evecs = evecs[:, ix]
    idx = [i for i in range(evals.shape[0]) if evals[i] > 0]
    evals = evals[idx]
    evecs = evecs[:, idx] * evals
    np.savetxt("reduced_dim.tsv", evecs, delimiter="\t")
    """
}