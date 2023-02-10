nextflow.enable.dsl=2

include { dim_reduct_nystrom as snapatac2 } from '../software/snapatac2.nf'

include { dim_reduct_archr_subsample as archr } from '../software/archr.nf'

process count_cells {
    //container 'kaizhang/snapatac2:1.99.99.7'
    input:
      tuple val(type), val(name), path(h5ad)
    output:
      tuple val(name), path(h5ad), path("n_cells.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    adata = snap.read("${h5ad}", backed='r')
    with open("n_cells.txt", "w") as f:
        f.write(str(adata.n_obs))
    """
}

process accuracy {
    //container 'kaizhang/snapatac2:1.99.99.7'
    input:
      tuple val(name), val(fraction), val(method), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple val(name), val(fraction), val(method), path("benchmark.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.metrics import adjusted_rand_score, silhouette_score
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)
    embedding = np.genfromtxt(
        "reduced_dim.tsv",
        missing_values="NA",
        filling_values = np.nan,
    )
    n_cluster = np.unique(adata.obs["cell_annotation"]).size

    if np.isnan(embedding).any():
        with open("benchmark.txt", "w") as f: print("nan", file=f)
    else:
        score = -1
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
        sc = adjusted_rand_score(clusters, adata.obs["cell_annotation"])
        if sc > score:
            score = sc
        with open("benchmark.txt", "w") as f:
            print(str(score), file=f)
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

workflow bench_subsample {
    take: datasets
    main:
        data = datasets
            | count_cells
            | filter { it[2].text.toInteger() > 20000 }
            | map { [it[0], it[1]] }
        bench_data = data | combine([0.1, 0.2, 0.5, 1.0])

        snapatac2(bench_data)
            | concat(archr(bench_data))
            | combine(data, by: 0)
            | accuracy
            | toSortedList
            | report
}
