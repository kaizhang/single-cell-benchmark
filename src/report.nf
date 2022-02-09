nextflow.enable.dsl=2

process benchmark_dim_reduct {
    input:
      tuple val(method), val(data), path(reduced_dim)
    output:
      tuple val(method), val(data), path("benchmark.txt"), path("umap.txt.gz")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import scanpy
    from sklearn.metrics import adjusted_rand_score
    import numpy as np
    adata = snap.read("${data.anndata}")
    try:
        low_dim = np.loadtxt("${reduced_dim}")
        score = -1
        best_dim = -1
        for n_dim in [5, 15, 30]:
            adata.obsm["X_spectral"] = low_dim[:, 0:n_dim]
            snap.pp.knn(adata)
            for i in np.arange(0.1, 1.8, 0.1):
                snap.tl.leiden(adata, resolution = i, use_weights=False)
                sc = adjusted_rand_score(adata.obs["leiden"], adata.obs["cell_annotation"])
                if sc > score:
                    score = sc
                    best_dim = n_dim
        with open("benchmark.txt", "w") as f: print(score, file=f)

        adata.obsm["X_spectral"] = low_dim[:, 0:best_dim]
        snap.tl.umap(adata)
        np.savetxt("umap.txt.gz", adata.obsm["X_umap"])
        
    except ValueError:
        with open("benchmark.txt", "w") as f: print("nan", file=f)
        with open("umap.txt.gz", "w") as f: print("nan", file=f)
    """
}


process generate_report {
    publishDir 'result'

    input:
        val result
    output:
        file "benchmark.tsv"

    exec:
    outputFile = task.workDir.resolve("benchmark.tsv")
    outputFile.text = "dataset\ttype\talgorithm\tsample_fraction\trandom_seed\tARI\n"
    for (x in result.sort { it.datasetName + it.methodName }) {
        outputFile.append([
            x.datasetName,
            x.datasetType,
            x.methodName,
            x.samplingFraction,
            x.randomSeed,
            x.resultOutput.text,
        ].join('\t'))
    }
}

process plot_report {
    publishDir 'result'

    input:
        val report
    output:
        file "*.pdf"

    """
    #!/usr/bin/env python3
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    data = pd.read_csv("${report}", sep="\t")

    plot = sns.FacetGrid(
        data[data["sample_fraction"] == 1],
        col="dataset",
        gridspec_kws={"bottom": 0.4}
    )
    plot.map(sns.barplot, "algorithm", "ARI")

    plot.set_xticklabels(rotation=45, horizontalalignment='right')

    plt.savefig('report.pdf')

    methods = set(data[data["sample_fraction"] < 1]["algorithm"])
    plot = sns.FacetGrid(
        data.loc[data["algorithm"].isin(methods)],
        col="dataset",
        hue="algorithm",
    )
    plot.map(sns.lineplot, "sample_fraction", "ARI",
        err_style="bars",
        ci=95,
        markers=True,
        marker="o",
        dashes=False,
    )
    plot.add_legend()
    plt.savefig('subsample.pdf')
    """
}

process plot_umap {
    publishDir 'result'

    input:
        val result
    output:
        file "*.pdf"

    """
    #!/usr/bin/env python3
    """
}

