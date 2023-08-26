nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct {
    container 'kaizhang/scvi-tools:1.0.3'
    tag "${json(metadata).data_name}"
    errorStrategy 'ignore'
    containerOptions '--nv'
    label "gpu"

    when: is_included("multivi", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path(h5ads, stageAs: "?.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'MultiVI')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import numpy as np
    import scvi
    import torch

    scvi.settings.seed = 0
    torch.set_float32_matmul_precision("high")

    adata_atac = ad.read("1.h5ad")
    adata_atac.var['modality'] = "Peaks"
    adata_rna = ad.read("2.h5ad")
    adata_rna.var['modality'] = "Gene Expression"
    adata_rna.X = adata_rna.layers["counts"]
    adata_paired = ad.concat([adata_rna, adata_atac], axis=1, join='inner', merge='same')

    adata_mvi = scvi.data.organize_multiome_anndatas(adata_paired)
    adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
    scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="modality")
    model = scvi.model.MULTIVI(
        adata_mvi,
        n_genes=(adata_mvi.var["modality"] == "Gene Expression").sum(),
        n_regions=(adata_mvi.var["modality"] == "Peaks").sum(),
        n_latent=30,
    )
    model.train()
    latent = model.get_latent_representation()
    np.savetxt("reduced_dim.tsv", latent, delimiter="\t")
    """
}
