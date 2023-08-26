nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct_scale {
    container 'kaizhang/scale:1.1.2'
    tag "${json(metadata).data_name}"
    errorStrategy 'ignore'
    containerOptions '--nv'
    label "gpu"

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

    def get_free_gpus(threshold_vram_usage=200):
        import subprocess
        import os

        # Get the list of GPUs via nvidia-smi
        smi_query_result = subprocess.check_output(
            "nvidia-smi -q -d Memory | grep -A4 GPU", shell=True
        )
        # Extract the usage information
        gpu_info = smi_query_result.decode("utf-8").split("\\n")
        gpu_info = list(filter(lambda info: "Used" in info, gpu_info))
        gpu_info = [
            int(x.split(":")[1].replace("MiB", "").strip()) for x in gpu_info
        ]  # Remove garbage
        # Keep gpus under threshold only
        free_gpus = [
            str(i) for i, mem in enumerate(gpu_info) if mem < threshold_vram_usage
        ]
        return free_gpus

    with tempfile.TemporaryDirectory(dir='./') as temp_dir:
        data = ad.read("data.h5ad")
        k = np.unique(data.obs["cell_annotation"]).size
        input_file = f"{temp_dir}/data.h5ad"
        data.X = data.X.astype(np.float64)
        data.write(input_file)

        gpus = get_free_gpus()
        if len(gpus) > 0:
            subprocess.run([
              "SCALE.py", "-d", input_file, "-k", str(k), "--min_peaks", "0",
              "--min_cells", "0", "-l", "30", "--gpu", gpus[0],
            ], check = True)
        else:
            raise Exception("No free GPU available")
        data = ad.read("output/adata.h5ad")
        np.savetxt("reduced_dim.tsv", data.obsm['latent'], delimiter="\t")
    """
}

