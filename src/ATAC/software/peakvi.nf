nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct_peakvi {
    container 'kaizhang/scvi-tools:0.19.0'
    tag "${json(metadata).data_name}"
    errorStrategy 'ignore'
    containerOptions '--nv'
    label "gpu"

    when: is_included("peakvi", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'PeakVI')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import numpy as np
    import scvi
    import math
    import json

    data = ad.read("data.h5ad")
    metadata = json.loads('${metadata}')

    if 'batch_key' in metadata:
        batch_key = metadata['batch_key']
        scvi.model.PEAKVI.setup_anndata(data, batch_key=batch_key)
    else:
        scvi.model.PEAKVI.setup_anndata(data)

    pvi = scvi.model.PEAKVI(data, n_latent=30)
    pvi.train()
    latent = pvi.get_latent_representation()
    np.savetxt("reduced_dim.tsv", latent, delimiter="\t")
    """
}
