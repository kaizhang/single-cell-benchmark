nextflow.enable.dsl=2

process gen_dim_reduct_report {
    publishDir 'result'

    input:
        val result
    output:
        file "benchmark.tsv"

    exec:
    outputFile = task.workDir.resolve("benchmark.tsv")
    outputFile.text = "dataset\ttype\talgorithm\tsample_fraction\trandom_seed\tanndata\tUMAP\tARI\n"
    for (x in result.sort { it.datasetName + it.methodName }) {
        outputFile.append([
            x.datasetName,
            x.datasetType,
            x.methodName,
            x.samplingFraction,
            x.randomSeed,
            x.datasetMatrix,
            x.umap,
            x.resultOutput.text,
        ].join('\t'))
    }
}

process plot_umap {
    publishDir 'result/umap'
    input:
        val report
    output:
        file "*.pdf"

    """
    #!/usr/bin/env python3
    import seaborn as sns
    import pandas as pd
    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    data = pd.read_csv("${report}", sep="\t")
    data = data[data["sample_fraction"] == 1]
    dfs = []
    for index, row in data.iterrows():
        annot = ad.read(row['anndata']).obs[['cell_annotation']]
        umap = np.loadtxt(row['UMAP'])
        annot["UMAP1"] = umap[:, 0]
        annot["UMAP2"] = umap[:, 1]
        annot["dataset"] = row["dataset"]
        annot["method"] = row["algorithm"]
        dfs.append(annot)
    df = pd.concat(dfs, ignore_index=False)
    for name, group in df.groupby("dataset"):
        plot = sns.FacetGrid(
            group,
            col="method",
            col_wrap=4,
            hue="cell_annotation",
            sharex=False,
            sharey=False,
        )
        plot.map(sns.scatterplot, "UMAP1", "UMAP2", s=2)
        plot.add_legend(markerscale=3)
        plt.savefig(name + '.pdf')
    """
}

process plot_dim_reduct_report {
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

process gen_clust_report {
    publishDir 'result'

    input:
        val result
    output:
        file "benchmark_clustering.tsv"

    exec:
    outputFile = task.workDir.resolve("benchmark_clustering.tsv")
    outputFile.text = "dataset\ttype\talgorithm\tARI\n"
    for (item in result) {
        for (line in item.resultOutput.text.split('\n')) {
            outputFile.append([
                item.datasetName,
                item.datasetType,
                line + "\n",
            ].join('\t'))
        }
    }
}

process plot_clust_report {
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
        data,
        col="dataset",
        gridspec_kws={"bottom": 0.4}
    )
    plot.map(sns.barplot, "algorithm", "ARI")

    plot.set_xticklabels(rotation=45, horizontalalignment='right')

    plt.savefig('report_clust.pdf')
    """
}

