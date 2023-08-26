nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct {
    container 'kaizhang/cobolt:1.0.1'
    tag "${json(metadata).data_name}"
    errorStrategy 'ignore'
    containerOptions '--nv'
    label "gpu"

    when: is_included("cobolt", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path(h5ads, stageAs: "?.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'Cobolt')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    from cobolt.utils import SingleData, MultiomicDataset
    from cobolt.model import Cobolt
    import anndata as ad
    import numpy as np
    import pandas as pd

    atac_data = ad.read('1.h5ad')
    rna_data = ad.read('2.h5ad')
    cobolt_rna = SingleData(
        "GeneExpr", "Multiome", rna_data.var_names,
        rna_data.layers['counts'].astype(np.float64), rna_data.obs_names
    )
    cobolt_atac = SingleData(
        "ChromAccess", "Multiome", atac_data.var_names,
        atac_data.X.astype(np.float64), atac_data.obs_names
    )

    multi_dt = MultiomicDataset.from_singledata(cobolt_rna, cobolt_atac)

    model = Cobolt(dataset=multi_dt, lr=0.002, n_latent=30)
    model.train()
    model.calc_all_latent()
    latent, barcodes = model.get_all_latent()
    barcodes = pd.DataFrame(index=[x.split("Multiome~")[1] for x in barcodes])
    latent = ad.AnnData(X=latent, obs=barcodes)[rna_data.obs_names].copy().X

    np.savetxt("reduced_dim.tsv", latent, delimiter="\t")
    """
}
