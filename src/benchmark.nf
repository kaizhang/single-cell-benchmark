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
                clusters = snap.tl.leiden(adata, adjacency = knn, resolution = i, inplace = False)
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
