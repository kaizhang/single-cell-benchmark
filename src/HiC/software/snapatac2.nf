nextflow.enable.dsl=2

process dim_reduct {
    container 'kaizhang/snapatac2:2.3.0'
    tag "$name"
    cpus 4
    input:
      tuple val(name), path("data.h5ad"), val(_), val(_), val(_), val(_)
    output:
      tuple val(name), val("SnapATAC2"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np

    adata = snap.read("data.h5ad", backed=None)
    snap.pp.select_features(adata, n_features=500000)
    snap.tl.spectral(adata, distance_metric="cosine")
    if 'batch' in adata.obs:
        snap.pp.mnc_correct(adata, batch='batch')
        adata.obsm['X_spectral'] = adata.obsm['X_spectral_mnn']
    np.savetxt("reduced_dim.tsv", adata.obsm['X_spectral'], delimiter="\t")
    """
}
