nextflow.enable.dsl=2

process dim_reduct_scale {
    container 'kaizhang/scale:1.1.2'
    errorStrategy 'ignore'
    containerOptions '--nv'
    cpus 16

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val('SCALE'), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import numpy as np
    import subprocess
    data = ad.read("data.h5ad")
    k = np.unique(data.obs["cell_annotation"]).size
    subprocess.run([
      "SCALE.py", "-d", "data.h5ad", "-k", str(k), "--min_peaks", "0",
      "--min_cells", "0", "-l", "30",
    ], check = True)
    data = ad.read("output/adata.h5ad")
    np.savetxt("reduced_dim.tsv", data.obsm['latent'], delimiter="\t")
    """
}

