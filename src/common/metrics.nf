nextflow.enable.dsl=2

include { json } from './utils.gvy'

params.resultDir = 'results'

process eval_embedding {
    container 'kaizhang/scatac-bench:0.2.1'
    tag "${json(metadata).data_name}"
    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      stdout

    """
    #!/usr/bin/env python
    import snapatac2 as snap
    import numpy as np
    import json
    import scib_metrics
    from sklearn.metrics import adjusted_rand_score, silhouette_score, adjusted_mutual_info_score
    import sys

    # Redirect stdout to stderr
    original_stdout = sys.stdout
    sys.stdout = sys.stderr

    metadata = json.loads('${metadata}')
    metrics = {}
    adata = snap.read("data.h5ad", backed=None)

    embedding = np.genfromtxt(
        "reduced_dim.tsv",
        missing_values="NA",
        filling_values = np.nan,
    )
    umap = snap.tl.umap(embedding, n_comps=3, inplace=False)

    knn_50 = snap.pp.knn(embedding, n_neighbors=50, method="exact", inplace=False)
    knn_90 = snap.pp.knn(embedding, n_neighbors=90, method="exact", inplace=False)

    n_cluster = np.unique(adata.obs["cell_annotation"]).size
    prev_n = -100
    prev_clusters = None
    for i in np.arange(0.1, 3, 0.1):
        clusters = snap.tl.leiden(knn_50, resolution=i, inplace=False)
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

    cell_anno = adata.obs["cell_annotation"]
    if knn_50.data.min() == 0:
        knn_50.data = knn_50.data + 1e-6
    if knn_90.data.min() == 0:
        knn_90.data = knn_90.data + 1e-6

    metrics["ARI"] = float(adjusted_rand_score(clusters, cell_anno))
    metrics["AMI"] = float(adjusted_mutual_info_score(clusters, cell_anno))
    metrics["Cell_type_ASW"] = scib_metrics.silhouette_label(umap, cell_anno)
    metrics['cLISI'] = float(np.median(scib_metrics.clisi_knn(knn_90, cell_anno)))
    if 'batch_key' in metadata:
        batch_key = metadata['batch_key']
        batch = adata.obs[batch_key]
        metrics['isolated_label_silhouette'] = float(scib_metrics.isolated_labels(
            umap, cell_anno, batch))
        metrics["kBET"] = float(scib_metrics.kbet(knn_50, batch)[0])
        metrics['iLISI'] = float(np.median(scib_metrics.ilisi_knn(knn_90, batch)))
        metrics['Graph_conn'] = float(scib_metrics.graph_connectivity(knn_50, cell_anno))
        metrics['Batch_ASW'] = float(scib_metrics.silhouette_batch(
            umap, cell_anno, batch
        ))

    metadata["metrics"] = metrics
    sys.stdout = original_stdout
    print(json.dumps(metadata), end='')
    """
}

process rare_label_asw {
    container 'kaizhang/scatac-bench:0.2.1'
    tag "${json(metadata).data_name}"
    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      stdout

    """
    #!/usr/bin/env python
    import snapatac2 as snap
    import numpy as np
    import json
    import scib_metrics
    from sklearn.metrics import silhouette_samples
    import sys

    # Redirect stdout to stderr
    original_stdout = sys.stdout
    sys.stdout = sys.stderr

    metadata = json.loads('${metadata}')
    adata = snap.read("data.h5ad", backed=None)
    rare_label = metadata['rare_label']
    embedding = np.genfromtxt(
        "reduced_dim.tsv",
        missing_values="NA",
        filling_values = np.nan,
    )
    umap = snap.tl.umap(embedding, n_comps=3, inplace=False)

    cell_anno = adata.obs["cell_annotation"]
    sw = silhouette_samples(umap, cell_anno)
    rare_label_ASW = np.mean(sw[cell_anno == rare_label])

    metadata["metrics"] = {"Rare_label_ASW": float(rare_label_ASW)}
    sys.stdout = original_stdout
    print(json.dumps(metadata), end='')
    """
}



process umap {
    container 'kaizhang/scatac-bench:0.2.1'
    tag "${json(metadata).data_name}"
    input:
      tuple val(metadata), path("reduced_dim.tsv"), path("data.h5ad")
    output:
      tuple val(metadata), path("umap.tsv.gz")

    """
    #!/usr/bin/env python
    import snapatac2 as snap
    import numpy as np
    import pandas as pd
    import json

    metadata = json.loads('${metadata}')
    metrics = {}
    adata = snap.read("data.h5ad", backed=None)

    embedding = np.genfromtxt(
        "reduced_dim.tsv",
        missing_values="NA",
        filling_values = np.nan,
    )
    umap = snap.tl.umap(embedding, inplace=False)

    if 'batch_key' in metadata:
        batch_key = metadata['batch_key']
        umap = pd.DataFrame({
            "UMAP-1": umap[:, 0],
            "UMAP-2": umap[:, 1],
            "cell_type": adata.obs["cell_annotation"],
            "batch": adata.obs[batch_key],
            "method": metadata['bench_id'],
        })
    else:
        umap = pd.DataFrame({
            "UMAP-1": umap[:, 0],
            "UMAP-2": umap[:, 1],
            "cell_type": adata.obs["cell_annotation"],
            "method": metadata['bench_id'],
        })
    umap.to_csv("umap.tsv.gz", sep="\t", index=False, header=True)
    """
}

process output_metrics {
    container 'kaizhang/scatac-bench:0.2.1'
    publishDir "${params.resultDir}", mode: 'copy'

    input:
        val result
    output:
        path "metrics.tsv"

    """
    #!/usr/bin/env python3
    import json
    import pandas as pd

    dicts = []
    for item in ${result}:
        metadata = json.loads(item)
        metric = metadata['metrics']
        metric['dataset'] = metadata['data_name']
        metric['algorithm'] = metadata['bench_id']
        metric['group'] = metadata['group'] if 'group' in metadata else 'main'
        dicts.append(metric)
    pd.DataFrame(dicts).to_csv("metrics.tsv", sep="\t", index=False, header=True)
    """
}

process output_rare_label_asw {
    container 'kaizhang/scatac-bench:0.2.1'
    publishDir "${params.resultDir}", mode: 'copy'

    input:
        val result
    output:
        path "rare_label_asw.tsv"

    """
    #!/usr/bin/env python3
    import json
    import pandas as pd

    dicts = []
    for item in ${result}:
        metadata = json.loads(item)
        metric = metadata['metrics']
        metric['dataset'] = metadata['data_name']
        metric['algorithm'] = metadata['bench_id']
        metric['group'] = metadata['group'] if 'group' in metadata else 'main'
        dicts.append(metric)
    pd.DataFrame(dicts).to_csv("rare_label_asw.tsv", sep="\t", index=False, header=True)
    """
}

process plot_metrics {
    publishDir "${params.resultDir}/figures", mode: 'copy'
    container 'kaizhang/scatac-bench:0.2.1'

    input:
        path "metrics.tsv"
    output:
        file "*.pdf"

    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from plottable import ColumnDefinition, ColDef, Table
    from plottable.plots import bar
    from plottable.cmap import normed_cmap
    from matplotlib.colors import LinearSegmentedColormap
    from sklearn.preprocessing import MinMaxScaler, StandardScaler

    df = pd.read_csv("metrics.tsv", sep="\t")
    datasets = np.unique(df['dataset'])

    cmap_fn = lambda col_data: normed_cmap(col_data, cmap=matplotlib.cm.PRGn, num_stds=2.5)
    summary = []
    for name in datasets:
        plot_df = df[df['dataset'] == name]
        # remove columns that contain only NA values
        plot_df = plot_df.loc[:, plot_df.notna().any(axis=0)]

        # Perform scaling column-wise
        numeric_cols = plot_df.select_dtypes(include=['float64', 'int64']).columns
        min_max_scaler = MinMaxScaler()
        plot_df[numeric_cols] = min_max_scaler.fit_transform(plot_df[numeric_cols])

        aggregation_cols = []
        bio_conservation_cols = {
            "ARI": "ARI",
            "AMI": "AMI",
            "Cell_type_ASW": "Cell type ASW",
            "cLISI": "Graph cLISI",
            "isolated_label_silhouette": "Isolated label silhouette",
        }
        bio_conservation_cols = {k: v for k, v in bio_conservation_cols.items() if k in plot_df.columns}
        if len(bio_conservation_cols) > 0:
            plot_df["Bio conservation"] = plot_df[list(bio_conservation_cols.keys())].mean(axis=1)
            aggregation_cols.append("Bio conservation")

        batch_correction_cols = {
            "Batch_ASW": "Batch ASW",
            "Graph_conn": "Graph connectivity",
            "kBET": "kBET",
            "iLISI": "Graph iLISI",
        }
        batch_correction_cols = {k: v for k, v in batch_correction_cols.items() if k in plot_df.columns}
        if len(batch_correction_cols) > 0:
            plot_df["Batch correction"] = plot_df[list(batch_correction_cols.keys())].mean(axis=1)
            aggregation_cols.append("Batch correction")

        if len(aggregation_cols) >= 2:
            plot_df["Overall"] = 0.6 * plot_df["Bio conservation"] + 0.4 * plot_df["Batch correction"]
            aggregation_cols.append("Overall")
        else:
            plot_df = plot_df.rename(columns={aggregation_cols[0]: "Overall"})
            aggregation_cols = ["Overall"]

        column_definitions = [
            ColumnDefinition(
                "algorithm",
                title="Method",
                width=plot_df["algorithm"].str.len().max() * 0.1,
                textprops={"ha": "left", "weight": "bold"}),
        ]

        column_definitions += [
            ColumnDefinition(
                col,
                width=1,
                title=title.replace(" ", "\\n", 1),
                textprops={
                    "ha": "center",
                    "bbox": {"boxstyle": "circle", "pad": 0.25},
                },
                cmap=cmap_fn(plot_df[col]),
                group="Bio conservation",
                formatter="{:.2f}",
            )
            for col, title in bio_conservation_cols.items()
        ]
        column_definitions += [
            ColumnDefinition(
                col,
                width=1,
                title=title.replace(" ", "\\n", 1),
                textprops={
                    "ha": "center",
                    "bbox": {"boxstyle": "circle", "pad": 0.25},
                },
                cmap=cmap_fn(plot_df[col]),
                group="Batch correction",
                formatter="{:.2f}",
            )
            for col, title in batch_correction_cols.items()
        ]
        column_definitions += [
            ColumnDefinition(
                col,
                width=1,
                title=col.replace(" ", "\\n", 1),
                plot_fn=bar,
                plot_kw={
                    "cmap": LinearSegmentedColormap.from_list(
                        name="bugw", colors=["#ffffff", "#f2fbd2", "#c9ecb4", "#93d3ab", "#35b0ab"], N=256
                    ),
                    "plot_bg_bar": False,
                    "annotate": True,
                    "height": 0.9,
                    "formatter": "{:.2f}",
                },
                group="Aggregate score",
                border="left" if i == 0 else None,
            )
            for i, col in enumerate(aggregation_cols)
        ]

        summary.append(plot_df[["group", "algorithm", "dataset", "Overall"]].rename(columns={"Overall": "score"}))
        plot_df = ( plot_df[["algorithm"] +
                    list(bio_conservation_cols.keys()) +
                    list(batch_correction_cols.keys()) +
                    aggregation_cols])
        plot_df = plot_df.sort_values(by="Overall", ascending=False)

        with matplotlib.rc_context({"svg.fonttype": "none"}):
            fig, ax = plt.subplots(figsize=(len(plot_df.columns) * 1.25, 3 + 0.5 * plot_df.shape[0]))
            tab = Table(
                    plot_df,
                    cell_kw={
                        "linewidth": 0,
                        "edgecolor": "k",
                    },
                    column_definitions=column_definitions,
                    ax=ax,
                    row_dividers=True,
                    footer_divider=True,
                    textprops={"fontsize": 10, "ha": "center"},
                    row_divider_kw={"linewidth": 1, "linestyle": (0, (1, 5))},
                    col_label_divider_kw={"linewidth": 1, "linestyle": "-"},
                    column_border_kw={"linewidth": 1, "linestyle": "-"},
                    index_col="algorithm",
                ).autoset_fontcolors(colnames=plot_df.columns)
            fig.savefig(name + ".pdf", facecolor=ax.get_facecolor(), bbox_inches='tight', pad_inches=0, dpi=300)
    summary = pd.concat(summary)
    for group in np.unique(summary['group']):
        plot_df = summary[summary['group'] == group]
        plot_df = plot_df.drop(columns=['group'])
        datasets = np.unique(plot_df['dataset'])
        plot_df = plot_df.pivot(index='algorithm', columns='dataset', values='score')
        plot_df['algorithm'] = plot_df.index
        plot_df['Overall'] = plot_df[datasets].mean(axis=1)
        plot_df = plot_df.sort_values(by="Overall", ascending=False)

        column_definitions = [
            ColumnDefinition(
                "algorithm",
                title="Method",
                width=plot_df["algorithm"].str.len().max() * 0.1,
                textprops={"ha": "left", "weight": "bold"}),
        ]

        column_definitions += [
            ColumnDefinition(
                col,
                width=1,
                title=col.replace(" ", "\\n", 1),
                textprops={
                    "rotation": 45,
                    "ha":"left", "va":"bottom",
                },
                plot_fn=bar,
                plot_kw={
                    "cmap": LinearSegmentedColormap.from_list(
                        name="bugw", colors=["#ffffff", "#f2fbd2", "#c9ecb4", "#93d3ab", "#35b0ab"], N=256
                    ),
                    "plot_bg_bar": False,
                    "annotate": True,
                    "height": 0.9,
                    "formatter": "{:.2f}",
                },
                group="Dataset",
                border="left" if i == 0 else None,
            )
            for i, col in enumerate(datasets)
        ]

        column_definitions += [
            ColumnDefinition(
                "Overall",
                width=1,
                plot_fn=bar,
                textprops={"weight": "bold"},
                plot_kw={
                    "cmap": LinearSegmentedColormap.from_list(
                        name="bugw", colors=["#FFE5B4", "#FFF3B0", "#FFD6C1", "#FFC4A1", "#FFD1A1"], N=256
                    ),
                    "plot_bg_bar": False,
                    "annotate": True,
                    "height": 0.9,
                    "formatter": "{:.2f}",
                },
                border="left",
            )
        ]

        with matplotlib.rc_context({"svg.fonttype": "none"}):
            fig, ax = plt.subplots(figsize=(len(plot_df.columns) * 1.25, 3 + 0.5 * plot_df.shape[0]))
            tab = Table(
                    plot_df,
                    cell_kw={
                        "linewidth": 0,
                        "edgecolor": "k",
                    },
                    column_definitions=column_definitions,
                    ax=ax,
                    row_dividers=True,
                    footer_divider=True,
                    textprops={"fontsize": 10, "ha": "center"},
                    row_divider_kw={"linewidth": 1, "linestyle": (0, (1, 5))},
                    col_label_divider_kw={"linewidth": 1, "linestyle": "-"},
                    column_border_kw={"linewidth": 1, "linestyle": "-"},
                    index_col="algorithm",
                ).autoset_fontcolors(colnames=plot_df.columns)
            fig.savefig(group + "_summary.pdf", facecolor=ax.get_facecolor(), bbox_inches='tight', pad_inches=0, dpi=300)
    """
}

process plot_umap {
    publishDir "${params.resultDir}/figures/UMAP/", mode: 'copy'
    container 'kaizhang/scatac-bench:0.2.1'
    input:
        tuple val(name), val(batch), path(umap, stageAs: "?.tsv.gz")
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

    def output_umap(df, label, file):
        n_label = np.unique(np.array(list(map(str, df[label])))).size
        n = -(np.unique(df['method'].to_numpy()).size // -4)
        size = 0.1 if n_points > 20000 else 0.25
        ( ggplot(df, aes(x='UMAP-1', y='UMAP-2', fill=label))
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
        ).save(file, dpi=300)

    output_umap(df, "cell_type", "${name}.pdf")

    if "${batch == null}" == 'false':
        output_umap(df, "batch", "${name}_batch.pdf")
    """
}

