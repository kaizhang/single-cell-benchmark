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
    metrics["ARI"] = float(adjusted_rand_score(clusters, cell_anno))
    metrics["AMI"] = float(adjusted_mutual_info_score(clusters, cell_anno))
    metrics["Cell_type_ASW"] = scib_metrics.silhouette_label(embedding, cell_anno)
    if knn_90.data.min() == 0:
        knn_90.data = knn_90.data + 1e-6
    metrics['cLISI'] = float(np.median(scib_metrics.clisi_knn(knn_90, cell_anno)))

    if 'batch_key' in metadata:
        batch_key = metadata['batch_key']
        batch = adata.obs[batch_key]
        metrics['isolated_label_silhouette'] = float(scib_metrics.isolated_labels(
            embedding, cell_anno, batch))
        metrics["kBET"] = float(scib_metrics.kbet(knn_50, batch)[0])
        metrics['iLISI'] = float(np.median(scib_metrics.ilisi_knn(knn_90, batch)))
        metrics['Graph_conn'] = float(scib_metrics.graph_connectivity(knn_50, cell_anno))
        metrics['Batch_ASW'] = float(scib_metrics.silhouette_batch(
            embedding, cell_anno, batch
        ))

    metadata["metrics"] = metrics
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

process plot_metrics {
    publishDir "${params.resultDir}/figures", mode: 'copy'
    container 'kaizhang/scatac-bench:0.2.1'

    input:
        path "metrics.tsv"
    output:
        file "*.pdf"

    """
    #!/usr/bin/env Rscript
    library(funkyheatmap)
    library(dplyr, warn.conflicts = FALSE)
    library(tibble, warn.conflicts = FALSE)
    library(ggplot2)

    min_max_scaling <- function(x) {
        min_x <- min(x, na.rm = TRUE)
        max_x <- max(x, na.rm = TRUE)
        (x - min_x) / (max_x - min_x)
    }

    bio_conservation <- function(df, column_info) {
        labels <- (column_info %>% filter(group == "Bio conservation"))['id']\$id
        df <- df[, labels]
        apply(df, 1, mean, na.rm=T)
    }

    batch_correction <- function(df, column_info) {
        labels <- (column_info %>% filter(group == "Batch correction"))['id']\$id
        df <- df[, labels]
        apply(df, 1, mean, na.rm=T)
    }

    data <- read.table('metrics.tsv', sep='\t',header=T) %>% rename(id = algorithm)

    for (name in unlist(unique(data['dataset']))) {
        df <- data[data['dataset'] == name,]

        # remove columns that contain only NA values
        df <- df %>% select(where(~ !all(is.na(.))))

        # Perform min-max scaling column-wise
        df <- df %>% mutate(across(where(is.numeric), min_max_scaling))

        contain_batch <- "kBET" %in% colnames(df)

        column_info <- tribble(
            ~id,                         ~group,             ~name,                       ~geom,       ~palette,    ~options,
            "id",                        NA,                 "",                          "text",      NA,          list(hjust = 0, width = 6),
        )

        if (contain_batch) {
            column_info <- bind_rows(column_info, tribble(
                ~id,                         ~group,             ~name,                       ~geom,       ~palette,    ~options,
                "overall",                   "",                 "Overall score",             "bar",       "palette3",  list(width = 3, legend = FALSE),
            ))
        }

        column_info <- bind_rows(column_info, tribble(
            ~id,                         ~group,             ~name,                       ~geom,       ~palette,    ~options,
            "bio_conservation",          "",                 "Bio conservation",          "bar",       "palette2",  list(width = 3, legend = FALSE),
            "ARI",         "Bio conservation", "ARI cluster label",         "funkyrect", "palette2",  lst(),
            "AMI",         "Bio conservation", "NMI cluster label",         "funkyrect", "palette2",  lst(),
            "Cell_type_ASW",                 "Bio conservation", "Cell type ASW",             "funkyrect", "palette2",  lst(),
            "cLISI",                     "Bio conservation", "Graph cLISI",               "funkyrect", "palette2",  lst(),
        ))

        if (contain_batch) {
            column_info <- bind_rows(column_info, tribble(
                ~id,                         ~group,             ~name,                       ~geom,       ~palette,    ~options,
                "isolated_label_silhouette", "Bio conservation", "Isolated label silhouette", "funkyrect", "palette2",  lst(),
                "batch_correction",          "",                 "Batch correction",          "bar",       "palette1",  list(width = 3, legend = FALSE),
                "Batch_ASW",                 "Batch correction", "Batch ASW",                 "funkyrect", "palette1",  lst(),
                "Graph_conn",                "Batch correction", "Graph connectivity",        "funkyrect", "palette1",  lst(),
                "kBET",                      "Batch correction", "kBET",                      "funkyrect", "palette1",  lst(),
                "iLISI",                     "Batch correction", "Graph iLISI",               "funkyrect", "palette1",  lst(),
            ))
        }


        df['bio_conservation'] <- bio_conservation(df, column_info)
        df <- df %>% arrange(desc(bio_conservation))

        if (contain_batch) {
            df['batch_correction'] <- batch_correction(df, column_info)
            df['overall'] <- 0.6* df['bio_conservation'] + 0.4 * df['batch_correction']
            df <- df %>% arrange(desc(overall))
        }

        g <- funky_heatmap(df, column_info = as.data.frame(column_info), expand = list(xmax = 4))
        ggsave(paste0(name, ".pdf"), g, device = cairo_pdf, width = g\$width, height = g\$height)
    }
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

