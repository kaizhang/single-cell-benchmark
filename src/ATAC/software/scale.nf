nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct_scale {
    container 'kaizhang/scale:1.1.2'
    tag "${json(metadata).data_name}"
    errorStrategy 'ignore'
    containerOptions '--nv'
    cpus 8

    when: is_included("scale", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'SCALE')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import numpy as np
    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory(dir='./') as temp_dir:
        data = ad.read("data.h5ad")
        k = np.unique(data.obs["cell_annotation"]).size
        input_file = f"{temp_dir}/data.h5ad"
        data.X = data.X.astype(np.float64)
        data.write(input_file)
        subprocess.run([
          "SCALE.py", "-d", input_file, "-k", str(k), "--min_peaks", "0",
          "--min_cells", "0", "-l", "30",
        ], check = True)
        data = ad.read("output/adata.h5ad")
        np.savetxt("reduced_dim.tsv", data.obsm['latent'], delimiter="\t")
    """
}

