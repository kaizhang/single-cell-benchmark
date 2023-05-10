nextflow.enable.dsl=2

include { dim_reduct_nystrom as snapatac2_subsample;
          dim_reduct_nystrom2 as snapatac2_subsample_degree;
          dim_reduct_cosine as snapatac2;
        } from '../software/snapatac2.nf'

include { dim_reduct_archr_subsample as archr_subsample;
          dim_reduct_archr_2 as archr;
        } from '../software/archr.nf'

process count_cells {
    container 'kaizhang/snapatac2:2.3.0'
    input:
      tuple val(name), path(h5ad)
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

process clustering {
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



process accuracy {
    container 'kaizhang/snapatac2:2.3.0'
    input:
      tuple val(name), val(_), val(fraction), val(method), path("reduced_dim.tsv"), path("clusters.txt")
    output:
      tuple val(name), val(fraction), val(method), path("benchmark.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.metrics import adjusted_rand_score, silhouette_score
    import numpy as np
    ground_truth = np.genfromtxt("clusters.txt", dtype=str)
    embedding = np.genfromtxt(
        "reduced_dim.tsv",
        missing_values="NA",
        filling_values = np.nan,
    )
    n_cluster = np.unique(ground_truth).size

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
        sc = adjusted_rand_score(clusters, ground_truth)
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

process plot {
    publishDir 'result'
    container 'kaizhang/scatac-bench:0.1.0'
    input:
        path "benchmark.tsv"
    output:
        path("*.pdf")

    """
    #!/usr/bin/env python3
    import pandas as pd
    from plotnine import *
    import glasbey
    import numpy as np

    data = pd.read_csv('benchmark.tsv', sep="\t", index_col=0).sort_values(['dataset', 'algorithm', 'fraction'])
    data['%total cells'] = data['fraction'] * 100
    df = pd.DataFrame([
        pd.Series(data.groupby(['algorithm', '%total cells']).sem()['ARI'], name="SE"),
        pd.Series(data.groupby(['algorithm', '%total cells']).mean()['ARI'], name="Mean"),
    ]).T.reset_index()

    ( ggplot(df, aes(x="%total cells", y="Mean", color="factor(algorithm)"))
    + geom_errorbar(aes(x="%total cells", ymin="Mean-SE", ymax="Mean+SE"), colour="black", width=2)
    + geom_point(size=0.6)
    + geom_line()
    + theme_light(base_size=7)
    + scale_color_manual(glasbey.create_palette(palette_size=2))
    + theme(
        figure_size=(3, 2),
        legend_key=element_rect(color = "white"),
    )
    ).save("subsample.pdf")
    """
}

workflow bench_subsample {
    take: datasets
    main:
        data = datasets
            | count_cells
            | filter { it[2].text.toInteger() > 5000 }
            | map { [it[0], it[1]] }
        subsampled_data = data | combine([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

        clusters = snapatac2(data) | map { [it[0], "SnapATAC2", it[2]] } | combine(data, by: 0)
            | concat(archr(data) | map { [it[0], "ArchR", it[2]] } | combine(data, by: 0))
            | clustering

        snapatac2_subsample(subsampled_data) | map { [it[0], "SnapATAC2", it[1], it[2], it[3]] }
            | concat(archr_subsample(subsampled_data) | map { [it[0], "ArchR", it[1], it[2], it[3]] })
            | combine(clusters, by: [0, 1])
            | accuracy
            | toSortedList
            | report
            | plot
}
