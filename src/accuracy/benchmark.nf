nextflow.enable.dsl=2

process benchmark_dim_reduct {
    //container 'kaizhang/snapatac2:1.99.99.7'
    input:
      tuple val(method), val(data), path(reduced_dim)
    output:
      tuple val(method), val(data), path("benchmark.txt"), path("umap.txt.gz")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.metrics import adjusted_rand_score
    import numpy as np
    adata = snap.read("${data.anndata}", mode = "r")
    try:
        low_dim = np.loadtxt("${reduced_dim}")
        score = -1
        best_dim = -1
        for n_dim in [5, 15, 30]:
            knn = snap.pp.knn(low_dim, use_dims = n_dim, inplace = False)
            for i in np.arange(0.1, 1.8, 0.1):
                clusters = snap.tl.leiden(knn, resolution = i, inplace = False)
                sc = adjusted_rand_score(clusters, adata.obs["cell_annotation"])
                if sc > score:
                    score = sc
                    best_dim = n_dim
        with open("benchmark.txt", "w") as f: print(score, file=f)

        np.savetxt(
            "umap.txt.gz",
            snap.tl.umap(low_dim, use_dims = best_dim, inplace = False)
        )
        
    except ValueError:
        with open("benchmark.txt", "w") as f: print("nan", file=f)
        with open("umap.txt.gz", "w") as f: print("nan", file=f)
    """
}

process benchmark_end_to_end {
    //container 'kaizhang/snapatac2:1.99.99.7'
    publishDir 'result'
    input:
      tuple val(method), val(data), val(dim_reduct), val(clusters)

    output:
      tuple val(method), val(data), path("${data.name}-${method}.pdf")

    """
    #!/usr/bin/env python3
    from sklearn.metrics import adjusted_rand_score
    import numpy as np
    import pandas as pd
    import snapatac2 as snap
    import matplotlib.pyplot as plt
    import seaborn as sns
    import colorcet
    meta = pd.read_csv("${data.metadata}", index_col = 0, sep = '\t')
    data = zip(
        ${dim_reduct.collect { '"' + it + '"' }},
        ${clusters.collect { '"' + it + '"' }},
    )

    reduced_dim = None
    cell_cluster = None
    cell_labels = None
    max_score = -1
    for dim_file, cluster in data:
        result = pd.read_csv(cluster, header=None, index_col = 0, sep = '\t', na_filter = False)
        expected = meta.loc[result.index]["subclass"]

        score = -1
        for (_, clusters) in result.iteritems():
            sc = adjusted_rand_score(clusters, expected)
            if sc > score: score = sc

        if score > max_score:
            max_score = score
            reduced_dims = pd.read_csv(
                dim_file, header=None, index_col = 0, sep = '\t', na_filter = False
            )
            idx = [result.index.get_loc(i) for i in reduced_dims.index]
            cell_cluster = clusters[idx]
            cell_labels = expected[idx]

    umap = snap.tl.umap(reduced_dims.to_numpy(), inplace=False)
    n = umap.shape[0]
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharey=True, figsize=(15,6))
    df = pd.DataFrame({
        "UMAP1": umap[:, 0],
        "UMAP2": umap[:, 1],
        "cluster": cell_cluster,
        "cell type": cell_labels,
    })
    size = 0.8 if n > 20000 else 2
    sns.scatterplot(
        x="UMAP1", y="UMAP2", data=df,
        hue = 'cluster', edgecolor="none", s=size, alpha=0.5,
        palette = colorcet.glasbey[:np.unique(cell_cluster).size],
        ax = ax1,
    )
    sns.scatterplot(
        x="UMAP1", y="UMAP2", data=df,
        hue = 'cell type', edgecolor="none", s=size, alpha=0.5,
        palette = colorcet.glasbey[:np.unique(cell_labels).size],
        ax = ax2,
    )
    ax1.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=2)
    ax2.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.subplots_adjust(left = 0.05, right = 0.88, wspace = 0.4)
    plt.title("ARI: " + str(max_score))
    plt.savefig("${data.name}-${method}.pdf")
    """
}


process benchmark_clustering {
    //container 'kaizhang/snapatac2:1.99.99.7'
    input:
      tuple val(method), val(data), path(reduced_dim)
    output:
      tuple val(method), val(data), path("benchmark.txt")

    """
    #!/usr/bin/env python3
    from tempfile import TemporaryDirectory
    import shutil
    import snapatac2 as snap
    from sklearn.metrics import adjusted_rand_score
    import numpy as np
    with TemporaryDirectory() as tmp_dir:
        shutil.copyfile("${data.anndata}", tmp_dir + "/tmp.h5ad")
        adata = snap.read(tmp_dir + "/tmp.h5ad")
        low_dim = np.loadtxt("${reduced_dim}")
        result = {}
        score = -1
        best_dim = -1
        n_clusters = len(np.unique(adata.obs["cell_annotation"]))

        # Standard KNN + leiden clustering
        for n_dim in [5, 15, 30]:
            adata.obsm["X_spectral"] = low_dim[:, 0:n_dim]
            snap.pp.knn(adata, use_approximate_search=False)
            for i in np.arange(0.1, 1.8, 0.1):
                snap.tl.leiden(adata, resolution = i)
                sc = adjusted_rand_score(adata.obs["leiden"], adata.obs["cell_annotation"])
                if sc > score:
                    score = sc
                    best_dim = n_dim
        result["knn-leiden"] = score

        # approximate KNN
        score = -1
        adata.obsm["X_spectral"] = low_dim[:, 0:best_dim]
        snap.pp.knn(adata, use_approximate_search=True)
        for i in np.arange(0.1, 1.8, 0.1):
            snap.tl.leiden(adata, resolution = i)
            sc = adjusted_rand_score(adata.obs["leiden"], adata.obs["cell_annotation"])
            if sc > score: score = sc
        result["approximate-knn-leiden"] = score

        # K-means
        adata.obsm["X_spectral"] = low_dim[:, 0:best_dim]
        snap.tl.kmeans(adata, n_clusters)
        result["kmeans"] = adjusted_rand_score(adata.obs["kmeans"], adata.obs["cell_annotation"])

        # HDBSCAN
        adata.obsm["X_spectral"] = low_dim[:, 0:best_dim]
        snap.tl.hdbscan(adata)
        result["hdbscan"] = adjusted_rand_score(adata.obs["hdbscan"], adata.obs["cell_annotation"])

        with open("benchmark.txt", "w") as f:
            result = "\\n".join([x + '\t' + str(y) for x, y in result.items()])
            print(result, file=f)
        """
}
