nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../common/utils.gvy'

process highly_variable_genes {
    container 'kaizhang/snapatac2:2.3.0'
    tag "${json(metadata).data_name}"
    input:
      tuple val(metadata), path("data.h5ad")
      val(n_top_genes)
    output:
      tuple val("${add_meta(metadata, 'hvg', n_top_genes)}"), path("out.h5ad")

    """
    #!/usr/bin/env python3
    import scanpy as sc
    adata = sc.read("data.h5ad")
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=${n_top_genes})
    adata = adata[:, adata.var.highly_variable]
    adata.write("out.h5ad", compression="gzip")
    """
}

process scale_features {
    container 'kaizhang/snapatac2:2.3.0'
    tag "${json(metadata).data_name}"
    input:
      tuple val(metadata), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'scaled', true)}"), path("out.h5ad")

    """
    #!/usr/bin/env python3
    import scanpy as sc
    adata = sc.read("data.h5ad")
    sc.pp.scale(adata, max_value=10)
    adata.write("out.h5ad", compression="gzip")
    """
}

process subset_atac_data {
    container 'kaizhang/snapatac2:2.3.1'
    tag "${json(metadata).data_name}"
    input:
      tuple val(metadata), path("data.h5ad")
      val(n_features)
    output:
      tuple val("${add_meta(metadata, 'atac_feats', n_features)}"), path("out.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    adata = snap.read("data.h5ad", backed=None)
    snap.pp.select_features(adata, n_features=${n_features}, filter_lower_quantile=0, filter_upper_quantile=0)
    adata[:, adata.var['selected']].write("out.h5ad", compression="gzip")
    """
}