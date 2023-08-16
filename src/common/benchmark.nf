nextflow.enable.dsl=2

include { json } from './utils.gvy'

params.resultDir = 'results'

process eval_embedding {
    container 'kaizhang/scatac-bench:0.2.0'
    tag "${json(metadata).data_name}"
    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple stdout, path("umap.tsv.gz")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.metrics import adjusted_rand_score, silhouette_score, adjusted_mutual_info_score
    import pandas as pd
    import numpy as np
    import json

    metadata = json.loads('${metadata}')

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
        ari = adjusted_rand_score(clusters, adata.obs["cell_annotation"])
        ami = adjusted_mutual_info_score(clusters, adata.obs["cell_annotation"])
        umap = snap.tl.umap(embedding, inplace = False)
        pd.DataFrame({
            "UMAP-1": umap[:, 0],
            "UMAP-2": umap[:, 1],
            "cell_annotation": adata.obs["cell_annotation"],
            "method": metadata['bench_id'],
        }).to_csv("umap.tsv.gz", sep="\t", index=False, header=True)

        metadata.update({
            "ARI": float(ari),
            "AMI": float(ami),
            "silhouette_score": float(silhouette_score(umap, adata.obs["cell_annotation"], sample_size = 10000)),
        })
        print(json.dumps(metadata), end='')
    """
}

process eval_cluster {
    container 'kaizhang/scatac-bench:0.2.0'
    input:
      tuple val(metadata), path("clusters.tsv"), path("data.h5ad")
    output:
      stdout

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.metrics import adjusted_rand_score, silhouette_score, adjusted_mutual_info_score
    import pandas as pd
    import numpy as np
    import json
    metadata = json.loads('${metadata}')
    adata = snap.read("data.h5ad", backed=None)
    clusters = np.genfromtxt("clusters.tsv")
    metadata.update({
        "ARI": float(adjusted_rand_score(clusters, adata.obs["cell_annotation"])),
        "AMI": float(adjusted_mutual_info_score(clusters, adata.obs["cell_annotation"])),
        "silhouette_score": float(silhouette_score(adata.X, adata.obs["cell_annotation"], sample_size = 10000)),
    })
    print(json.dumps(metadata), end='')
    """
}

process output_metrics {
    container 'kaizhang/scatac-bench:0.2.0'
    publishDir "${params.resultDir}/metrics", mode: 'copy'

    input:
        val result
    output:
        path "benchmark.tsv"

    """
    #!/usr/bin/env python3
    import json
    import pandas as pd

    dicts = []
    for item in ${result}:
        metadata = json.loads(item)
        dicts.append({
            "dataset": metadata['data_name'],
            "algorithm": metadata['bench_id'],
            "ARI": metadata['ARI'],
            "AMI": metadata['AMI'],
            "silhouette_score": metadata['silhouette_score'],
        })
    pd.DataFrame(dicts).to_csv("benchmark.tsv", sep="\t", index=False, header=True)
    """
}

process plot_metrics {
    publishDir "${params.resultDir}/plots", mode: 'copy'
    container 'kaizhang/scatac-bench:0.2.0'

    input:
        path "benchmark.tsv"
    output:
        file "*.pdf"

    """
    #!/usr/bin/env python3
    import pandas as pd
    from plotnine import *
    import numpy as np
    import itertools
    import glasbey

    data = pd.read_csv('benchmark.tsv', sep="\t")
    data = data[['algorithm', 'dataset', 'ARI', 'AMI', 'silhouette_score']]
    my_order = data.groupby(by=["algorithm"])["ARI"].mean().sort_values().index[::-1]
    data['algorithm'] = pd.Categorical(data['algorithm'], ordered=True, categories=my_order)

    n_label = np.unique(data['algorithm'].to_numpy()).size
    n = -(np.unique(data['dataset'].to_numpy()).size // -3)
    shapes = ['o', 's', '^', 'v', 'D', 'X', '*', 'p', 'h', 'P', 'd', '<', '>']
    shapes = list(itertools.islice(itertools.cycle(shapes), n_label))
    colors = glasbey.create_palette(palette_size=n_label)

    ( ggplot(data, aes(x='algorithm', y='ARI', fill='factor(algorithm)'))
        + geom_col()
        + facet_wrap('dataset', scales="free_y", ncol=3)
        + scale_fill_manual(colors)
        + scale_color_manual(colors)
        + theme_light(base_size=7, base_family="Arial")
        + theme(
            #axis_text_x=element_text(rotation=90, hjust=1),
            axis_text_x=element_text(rotation=90),
            figure_size=(2 * 3, 2 * n),
            subplots_adjust={'hspace':0.2, 'wspace':0.2},
            legend_key=element_rect(color = "white"),
        )
    ).save(filename='ARI.pdf')

    ( ggplot(data, aes(x='algorithm', y='AMI', fill='factor(algorithm)'))
        + geom_col()
        + facet_wrap('dataset', scales="free_y", ncol=3)
        + scale_fill_manual(colors)
        + scale_color_manual(colors)
        + theme_light(base_size=7, base_family="Arial")
        + theme(
            #axis_text_x=element_text(rotation=90, hjust=1),
            axis_text_x=element_text(rotation=90),
            figure_size=(2 * 3, 2 * n),
            subplots_adjust={'hspace':0.2, 'wspace':0.2},
            legend_key=element_rect(color = "white"),
        )
    ).save(filename='AMI.pdf')

    ( ggplot(data, aes(x='AMI', y='ARI', color='factor(algorithm)'))
        + geom_point(aes(shape='factor(algorithm)'))
        + facet_wrap('dataset', scales="free", ncol=3)
        + scale_fill_manual(colors)
        + scale_color_manual(colors)
        + scale_shape_manual(values=shapes)
        + theme_light(base_size=7, base_family="Arial")
        + theme(
            figure_size=(2 * 3, 2 * n),
            subplots_adjust={'hspace':0.2, 'wspace':0.2},
            legend_key=element_rect(color = "white"),
        )
    ).save(filename='ARI_vs_AMI.pdf')
    """
}

process plot_umap {
    publishDir "${params.resultDir}/plots/UMAP/", mode: 'copy'
    container 'kaizhang/scatac-bench:0.2.0'
    input:
        tuple val(name), path(umap, stageAs: "?.tsv.gz")
    output:
        file "*.pdf"

    """
    #!/usr/bin/env python3
    from plotnine import *
    import pandas as pd
    import anndata as ad
    import glasbey
    import numpy as np

    dfs = []
    for file in "${umap}".split():
        df = pd.read_csv(file, sep="\t")
        n_points = df.shape[0]
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=False)

    n_label = np.unique(np.array(list(map(str, df['cell_annotation'])))).size
    n = -(np.unique(df['method'].to_numpy()).size // -4)
    size = 0.1 if n_points > 20000 else 0.25

    ( ggplot(df, aes(x='UMAP-1', y='UMAP-2', fill='factor(cell_annotation)'))
        + geom_point(size=size, raster=True, stroke=0)
        + facet_wrap('method', scales="free", ncol=4)
        + scale_fill_manual(glasbey.create_palette(palette_size=n_label))
        + theme_light(base_size=7)
        + theme(
            axis_ticks=element_blank(),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            figure_size=(1.7 * 4, 1.5 * n),
            subplots_adjust={'hspace':0.4, 'wspace':0.2},
            legend_key=element_rect(color = "white"),
        )
        + guides(fill = guide_legend(override_aes = {'size': 3}))
    ).save("${name}.pdf", dpi=300)
    """
}

workflow bench_embedding {
    take: datasets
    main:
        metrics = datasets | eval_embedding
        metrics | map { "'" + it[0] + "'" } | toSortedList | output_metrics | plot_metrics
        metrics
            | map { [json(it[0]).data_name, it[1]] }
            | groupTuple(sort: true)
            | plot_umap
}

workflow bench_clustering {
    take: datasets
    main:
        metrics = datasets | eval_cluster
        metrics | map { "'" + it + "'" } | toSortedList | output_metrics | plot_metrics
}