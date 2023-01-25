nextflow.enable.dsl=2

process dim_reduct_peakvi {
    container 'kaizhang/scvi-tools:0.19.0'
    errorStrategy 'ignore'
    containerOptions '--nv'

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val('PeakVI'), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import numpy as np
    import scvi
    import math
    data = ad.read("data.h5ad")
    scvi.model.PEAKVI.setup_anndata(data)
    pvi = scvi.model.PEAKVI(data, n_latent=30)
    pvi.train()
    latent = pvi.get_latent_representation()
    np.savetxt("reduced_dim.tsv", latent, delimiter="\t")
    """
}
