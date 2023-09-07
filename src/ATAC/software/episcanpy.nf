nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct_episcanpy {
    container 'kaizhang/episcanpy:0.4.0'
    tag "${json(metadata).data_name}"
    cpus 4
    errorStrategy 'ignore'

    when: is_included("pca", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'epiScanpy (PCA)')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python
    import anndata as ad
    import episcanpy as epi
    import scanpy as sc
    import numpy as np
    def find_elbow(x):
        n = len(x)
        marks = []
        saturation = 0.01
        accum_gap = 0
        for i in range(1, n):
            gap = x[i-1] - x[i]
            accum_gap += gap
            if gap > saturation * accum_gap:
                marks.append(i)
        return min(n - 1, max(marks) + 1)
    data = ad.read("data.h5ad")
    data.X.data = data.X.data.astype(np.float64)
    epi.pp.normalize_per_cell(data)
    epi.pp.log1p(data)
    epi.pp.pca(data, n_comps=30)
    sc.pl.pca_variance_ratio(data, log=True, save=".png")
    n_pc = find_elbow(data.uns['pca']['variance'])
    print(f"Elbow point: {n_pc}")
    embedding = data.obsm['X_pca'][:, :n_pc]
    np.savetxt("reduced_dim.tsv", embedding, delimiter="\t")
    """
}