nextflow.enable.dsl=2

include { add_meta; json } from '../../common/utils.gvy'

process knn_leiden_exact {
    container 'kaizhang/snapatac2:2.3.1'
    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'bench_id', 'knn+leiden (exact)')}"), path("clusters.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)
    embedding = np.genfromtxt("reduced_dim.tsv")
    n_cluster = np.unique(adata.obs["cell_annotation"]).size

    knn = snap.pp.knn(embedding, method="exact", inplace=False)
    prev_n = -100
    prev_clusters = None
    for i in np.arange(0.1, 3, 0.1):
        clusters = snap.tl.leiden(knn, resolution = i, inplace = False)
        n = np.unique(clusters).size
        if n == n_cluster:
            break
        elif n > n_cluster:
            if n - n_cluster < n_cluster - prev_n:
                break
            else:
                clusters = prev_clusters
                break
        else:
            prev_clusters = clusters
            prev_n = n
    with open("clusters.txt", "w") as f:
        print('\\n'.join(clusters), file=f)
    """
}

process knn_leiden_exact_weighted {
    container 'kaizhang/snapatac2:2.3.1'
    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'bench_id', 'weighted_knn+leiden (exact)')}"), path("clusters.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)
    embedding = np.genfromtxt("reduced_dim.tsv")
    n_cluster = np.unique(adata.obs["cell_annotation"]).size

    knn = snap.pp.knn(embedding, method="exact", inplace=False)
    prev_n = -100
    prev_clusters = None
    for i in np.arange(0.1, 3, 0.1):
        clusters = snap.tl.leiden(knn, resolution=i, weighted=True, inplace=False)
        n = np.unique(clusters).size
        if n == n_cluster:
            break
        elif n > n_cluster:
            if n - n_cluster < n_cluster - prev_n:
                break
            else:
                clusters = prev_clusters
                break
        else:
            prev_clusters = clusters
            prev_n = n
    with open("clusters.txt", "w") as f:
        print('\\n'.join(clusters), file=f)
    """
}

process knn_leiden_hora {
    container 'kaizhang/snapatac2:2.3.1'
    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'bench_id', 'knn+leiden (hora)')}"), path("clusters.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)
    embedding = np.genfromtxt("reduced_dim.tsv")
    n_cluster = np.unique(adata.obs["cell_annotation"]).size

    knn = snap.pp.knn(embedding, method="hora", inplace=False)
    prev_n = -100
    prev_clusters = None
    for i in np.arange(0.1, 3, 0.1):
        clusters = snap.tl.leiden(knn, resolution = i, inplace = False)
        n = np.unique(clusters).size
        if n == n_cluster:
            break
        elif n > n_cluster:
            if n - n_cluster < n_cluster - prev_n:
                break
            else:
                clusters = prev_clusters
                break
        else:
            prev_clusters = clusters
            prev_n = n
    with open("clusters.txt", "w") as f:
        print('\\n'.join(clusters), file=f)
    """
}

process kmeans {
    container 'kaizhang/snapatac2:2.3.1'
    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'bench_id', 'Kmeans')}"), path("clusters.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)
    embedding = np.genfromtxt("reduced_dim.tsv")
    n_cluster = np.unique(adata.obs["cell_annotation"]).size
    clusters = snap.tl.kmeans(embedding, n_clusters=n_cluster, inplace=False)
    with open("clusters.txt", "w") as f:
        print('\\n'.join(clusters), file=f)
    """
}