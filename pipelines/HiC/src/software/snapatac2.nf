nextflow.enable.dsl=2

process dim_reduct {
    //container 'kaizhang/snapatac2:1.99.99.7'
    cpus 4
    errorStrategy 'ignore'
    input:
      tuple val(name), path("data.h5ad"), val(_), val(_), val(_)
    output:
      tuple val(name), val("SnapATAC2"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)
    _, evecs = snap.tl.spectral(adata, features=None, inplace=False, distance_metric="cosine")
    np.savetxt("reduced_dim.tsv", evecs, delimiter="\t")
    """
}