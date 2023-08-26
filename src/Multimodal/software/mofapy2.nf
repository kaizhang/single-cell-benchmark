nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct {
    container 'kaizhang/mofapy2:0.7.0'
    tag "${json(metadata).data_name}"
    errorStrategy 'ignore'
    cpus 4

    when: is_included("mofa+", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path(h5ads, stageAs: "?.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'MOFA+')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import numpy as np
    import scanpy as sc
    import muon as mu
    from muon import atac as ac

    mdata = mu.MuData({'rna': sc.read('2.h5ad'), 'atac': sc.read('1.h5ad')})
    atac = mdata.mod['atac']
    ac.pp.tfidf(atac, scale_factor=1e4)
    mu.tl.mofa(mdata, n_factors=30, seed=2023)
    embedding = mdata.obsm['X_mofa']
    np.savetxt("reduced_dim.tsv", embedding, delimiter="\t")
    """
}