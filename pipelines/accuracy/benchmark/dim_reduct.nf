nextflow.enable.dsl=2

include { download_hg38; download_hg19 } from '../dataset.nf'

include { dim_reduct_snapatac as snapatac; } from '../software/snapatac.nf'

include { dim_reduct_jaccard as snapatac2_jaccard;
          dim_reduct_cosine as snapatac2_cosine;
          dim_reduct_svd as snapatac2_svd;
        } from '../software/snapatac2.nf'

include { dim_reduct_signac_1 as signac_1;
          dim_reduct_signac_2 as signac_2;
          dim_reduct_signac_3 as signac_3;
          dim_reduct_signac_4 as signac_4;
        } from '../software/signac.nf'

include { dim_reduct_archr_1 as archr_1;
          dim_reduct_archr_2 as archr_2;
          dim_reduct_archr_3 as archr_3;
        } from '../software/archr.nf'

include { dim_reduct_peakvi as peakvi } from '../software/peakvi.nf'
include { dim_reduct_scale as scale } from '../software/scale.nf'
include { dim_reduct_pycistopic as cisTopic } from '../software/pycistopic.nf'
include { dim_reduct_scbasset as scBasset } from '../software/scbasset.nf'
include { dim_reduct_episcanpy as epiScanpy } from '../software/episcanpy.nf'

process accuracy {
    //container 'kaizhang/snapatac2:1.99.99.7'
    errorStrategy 'ignore'
    input:
      tuple val(name), val(method), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple val(name), val(method), path("benchmark.txt"), path("umap.txt.gz"), path("data.h5ad")

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
        with open("umap.txt.gz", "w") as f: print("nan", file=f)
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
        umap = snap.tl.umap(embedding, inplace = False)
        np.savetxt("umap.txt.gz", umap)
        with open("benchmark.txt", "w") as f:
            ss = silhouette_score(umap, adata.obs["cell_annotation"], sample_size = 10000)
            print("\t".join([str(score), str(ss)]), file=f)
    """
}

process report {
    publishDir 'result'

    input:
        val result
    output:
        path "benchmark.tsv"

    exec:
    outputFile = task.workDir.resolve("benchmark.tsv")
    outputFile.text = "dataset\talgorithm\tUMAP\tanndata\tARI\tsilhouette_score\n"
    for (x in result) {
        outputFile.append([
            x[0],
            x[1],
            x[3],
            x[4],
            x[2].text
        ].join('\t'))
    }
}

process plot {
    publishDir 'result'

    input:
        path "benchmark.tsv"
    output:
        file "*.pdf"

    """
    #!/usr/bin/env python3
    import pandas as pd
    from plotnine import *
    import numpy as np
    import colorcet

    data = pd.read_csv("benchmark.tsv", sep="\t")
    df = data[["algorithm", "dataset", "ARI", "silhouette_score"]]
    my_order = df.groupby(by=["algorithm"])["ARI"].mean().sort_values().index[::-1]
    df['algorithm'] = pd.Categorical(df['algorithm'], ordered=True, categories=my_order)


    ( ggplot(df, aes(x="algorithm", y="ARI"))
        + geom_violin(df)
        + theme(axis_text_x=element_text(angle=90, hjust=1))
        + geom_jitter()
    ).save(filename='benchmark1.pdf')

    n_label = np.unique(df['algorithm'].to_numpy()).size
    n = -(np.unique(df['dataset'].to_numpy()).size // -3)
    ( ggplot(df, aes(x='silhouette_score', y='ARI', color='factor(algorithm)'))
        + geom_point(aes(shape='factor(algorithm)'))
        + facet_wrap('dataset', scales="free", ncol=3)
        + scale_fill_manual(colorcet.glasbey[:n_label])
        + scale_color_manual(colorcet.glasbey[:n_label])
        + theme_light(base_size=7, base_family="Arial")
        + theme(
            figure_size=(2 * 3, 2 * n),
            subplots_adjust={'hspace':0.2, 'wspace':0.2},
            legend_key=element_rect(color = "white"),
        )
    ).save(filename='benchmark2.pdf')
    """
}

process umap {
    publishDir 'result/umap'
    input:
        path report
    output:
        file "*.pdf"

    """
    #!/usr/bin/env python3
    import seaborn as sns
    import pandas as pd
    import anndata as ad
    import matplotlib.pyplot as plt
    import colorcet
    import numpy as np
    data = pd.read_csv("${report}", sep="\t")
    dfs = []
    for index, row in data.iterrows():
        try:
            umap = np.loadtxt(row['UMAP'])
            annot = ad.read(row['anndata']).obs[['cell_annotation']]
            n = umap.shape[0]
            annot["UMAP1"] = umap[:, 0]
            annot["UMAP2"] = umap[:, 1]
            annot["dataset"] = row["dataset"]
            annot["method"] = row["algorithm"]
            dfs.append(annot)
        except:
            pass
    df = pd.concat(dfs, ignore_index=False)
    size = 0.8 if n > 20000 else 2
    n_label = np.unique(np.array(list(map(str, df['cell_annotation'])))).size
    for name, group in df.groupby("dataset"):
        plot = sns.FacetGrid(
            group,
            col="method",
            col_wrap=4,
            hue="cell_annotation",
            palette = colorcet.glasbey[:n_label],
            sharex=False,
            sharey=False,
        )
        plot.map(
            sns.scatterplot, "UMAP1", "UMAP2",
            edgecolor="none", s=size, alpha=0.5,
        )
        plot.add_legend()
        for lh in plot._legend.legendHandles: 
            lh.set_alpha(1)
            lh._sizes = [30]
        plt.savefig(name + '.pdf')
    """
}

workflow bench_dim_reduct {
    take: datasets
    main:
        dr_bench_data = datasets | map { [it[1], it[2]] }

        hg19_dataset = ["BoneMarrow_Chen_2019", "Buenrostro_2018"]
        dr_bench_data_hg19 = dr_bench_data
            | filter { hg19_dataset.contains(it[0]) }
            | combine(download_hg19())

        reports = snapatac2_jaccard(dr_bench_data) | concat(
            snapatac2_cosine(dr_bench_data),
            snapatac2_svd(dr_bench_data),

            snapatac(dr_bench_data),

            epiScanpy(dr_bench_data),

            signac_1(dr_bench_data),
            signac_2(dr_bench_data),
            signac_3(dr_bench_data),
            signac_4(dr_bench_data),

            archr_1(dr_bench_data),
            archr_2(dr_bench_data),
            archr_3(dr_bench_data),

            peakvi(dr_bench_data),
            scale(dr_bench_data),
            scBasset(dr_bench_data_hg19),

            cisTopic(dr_bench_data),
        )
            | combine(dr_bench_data, by: 0)
            | accuracy
            | toSortedList
            | report

        reports | plot
        reports | umap
}
