nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct {
    container 'kaizhang/snapatac2:2.3.1'
    tag "${json(metadata).data_name}"
    cpus 4

    when: is_included("snapatac2", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'SnapATAC2')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np

    adata = snap.read("data.h5ad", backed=None)
    snap.tl.spectral(adata, features=None)
    np.savetxt("reduced_dim.tsv", adata.obsm['X_spectral'], delimiter="\t")
    """
}

process dim_reduct_svd {
    container 'kaizhang/snapatac2:2.3.1'
    tag "${json(metadata).data_name}"
    cpus 4

    when: is_included("snapatac2 (svd)", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'SnapATAC2 (svd)')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import scipy as sp
    import numpy as np

    adata = snap.read("data.h5ad", backed=None)
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