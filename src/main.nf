nextflow.enable.dsl=2

include { dim_reduct_snapatac_2;
          dim_reduct_snapatac_2_svd;
          dim_reduct_snapatac_2_cosine;
          dim_reduct_snapatac_2_nystrom;
          dim_reduct_snapatac_2_v2_nystrom_full;
          dim_reduct_snapatac_1;
          dim_reduct_snapatac_1_nystrom;
          dim_reduct_archr_1_log_tf_idf;
          dim_reduct_archr_1_logtf_logidf;
          dim_reduct_archr_1_tf_logidf;
          dim_reduct_archr_1_subsample;
        } from './software'
include { import_dataset } from './dataset'

process benchmark_dim_reduct {
    input:
      tuple val(method), val(data), path(reduced_dim)
    output:
      tuple val(method), val(data), path("benchmark.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.metrics import adjusted_rand_score
    import numpy as np
    adata = snap.read("${data.anndata}")
    try:
        low_dim = np.loadtxt("${reduced_dim}")
        score = -1
        for n_dim in [5, 15, 30]:
            adata.obsm["X_spectral"] = low_dim[:, 0:n_dim]
            snap.pp.knn(adata)
            for i in np.arange(0.1, 1.8, 0.1):
                snap.tl.leiden(adata, resolution = i, use_weights=False)
                sc = adjusted_rand_score(adata.obs["leiden"], adata.obs["cell_annotation"])
                if sc > score: score = sc
        with open("benchmark.txt", "w") as f: print(score, file=f)
    except ValueError:
        with open("benchmark.txt", "w") as f: print("nan", file=f)
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
        col="dataset"
    )
    plot.map(sns.barplot, "algorithm", "ARI")

    plot.set_xticklabels(rotation=45, horizontalalignment='right')
    #for axes in plot.axes.flat:
    #    _ = axes.set_xticklabels(axes.get_xticklabels(), rotation=45, horizontalalignment='right')

    plt.savefig('report.pdf')

    methods = set(data[data["sample_fraction"] < 1]["algorithm"])
    plot = sns.FacetGrid(
        data.loc[data["algorithm"].isin(methods)],
        col="dataset",
        hue="algorithm",
    )
    plot.map(sns.lineplot, "sample_fraction", "ARI",
        err_style="bars",
        ci=68,
        markers=True,
        dashes=False,
    )
    plot.add_legend()
    plt.savefig('subsample.pdf')
    """
}

workflow {
    datasets = import_dataset(Channel.fromList([
        tuple("Real", "dataset/Real"),
        tuple("Synthetic", "dataset/Synthetic"),
    ])).flatten()
    
    benchData = datasets.map { tuple(it, 50) }

    nystromBenchData = datasets.flatMap { data -> [
            [0.1, 0.2, 0.5],
            [1, 2, 3, 4, 5],
        ].combinations().collect { x -> 
            d = data.clone()
            d["samplingFraction"] = x[0]
            d["randomSeed"] = x[1]
            tuple(d, 50)
        }
    }

    runResult = dim_reduct_snapatac_2(benchData).concat(
        dim_reduct_archr_1_log_tf_idf(benchData),
        dim_reduct_archr_1_logtf_logidf(benchData),
        dim_reduct_archr_1_tf_logidf(benchData),

        dim_reduct_snapatac_1(benchData),
        dim_reduct_snapatac_2_cosine(benchData),
        dim_reduct_snapatac_2_svd(benchData),

        dim_reduct_snapatac_1_nystrom(nystromBenchData),
        dim_reduct_snapatac_2_nystrom(nystromBenchData),
        dim_reduct_archr_1_subsample(nystromBenchData),
    )

    benchResult = benchmark_dim_reduct(runResult)
        .map { method, data, result -> [
            "datasetType": data.dataType,
            "datasetName": data.name,
            "methodName": method,
            "samplingFraction": data.get("samplingFraction", 1),
            "randomSeed": data.randomSeed,
            "resultOutput": result
        ]}.collect()

    plot_report(generate_report(benchResult))
}
