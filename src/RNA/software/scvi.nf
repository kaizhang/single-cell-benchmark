nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct {
    container 'kaizhang/scvi-tools:1.0.3'
    tag "${json(metadata).data_name}"
    errorStrategy 'ignore'
    cpus 4
    containerOptions '--nv'
    label "gpu"

    when: is_included("scvi", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'scVI')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import scvi
    import numpy as np
    scvi.settings.seed = 0
    adata = ad.read("data.h5ad")
    scvi.model.SCVI.setup_anndata(adata, layer="counts")
    model = scvi.model.SCVI(adata, n_latent=30)
    model.train()
    latent = model.get_latent_representation()

    np.savetxt("reduced_dim.tsv", latent, delimiter="\t")
    """
}