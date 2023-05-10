nextflow.enable.dsl=2

include { dim_reduct_cosine as snapatac2 } from '../software/snapatac2.nf'

include { compute_accuracy_metrics } from '../../../common/benchmark.nf'

process knn_exact {
    container 'kaizhang/snapatac2:2.3.0'
    input:
      tuple val(name), val(method), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple val(name), val("exact"), path("clusters.txt")

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

process knn_hora {
    container 'kaizhang/snapatac2:2.3.0'
    input:
      tuple val(name), val(method), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple val(name), val(method), path("clusters.txt")

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

process report {
    publishDir 'result'

    input:
        val result
    output:
        path "subsample_benchmark.tsv"

    exec:
    outputFile = task.workDir.resolve("subsample_benchmark.tsv")
    outputFile.text = "dataset\tfraction\talgorithm\tARI\n"
    for (x in result) {
        outputFile.append([
            x[0],
            x[1],
            x[2],
            x[3].text
        ].join('\t'))
    }
}

workflow bench_knn {
    take: datasets
    main:
        data = datasets | map { [it[0], it[1]] }
        embeddings = snapatac2(data) | combine(data, by: 0)

        //knn_exact(embeddings) | combine(data, by: 0) | compute_accuracy_metrics
}
