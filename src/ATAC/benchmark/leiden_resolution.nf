nextflow.enable.dsl=2

include { dim_reduct_cosine as snapatac2 } from '../software/snapatac2.nf'
include { json } from '../../common/utils.gvy'

process leiden {
    container 'kaizhang/snapatac2:2.3.1'
    input:
      tuple val(name), path("reduced_dim.tsv"), path('data.h5ad')
    output:
      path("result.csv.gz")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import snapatac2_contrib as snap_contrib
    import numpy as np
    import pandas as pd
    from sklearn.metrics import (adjusted_rand_score, silhouette_score,
        calinski_harabasz_score, davies_bouldin_score, adjusted_mutual_info_score)
    adata = snap.read("data.h5ad", backed=None)
    embedding = np.genfromtxt("reduced_dim.tsv")
    knn = snap.pp.knn(embedding, method="exact", inplace=False)

    resolutions = list(np.arange(0.1, 3, 0.1))
    aris = []
    amis = []
    silhouettes = []
    calinskis = []
    davies_bouldins = []
    disp = []
    pac = []
    for i in resolutions:
        clusters = snap.tl.leiden(knn, resolution=i, inplace=False)
        aris.append(adjusted_rand_score(clusters, adata.obs["cell_annotation"]))
        amis.append(adjusted_mutual_info_score(clusters, adata.obs["cell_annotation"]))
        silhouettes.append(silhouette_score(embedding, clusters, sample_size = 10000))
        calinskis.append(calinski_harabasz_score(embedding, clusters))
        davies_bouldins.append(davies_bouldin_score(embedding, clusters))
        d, p = snap_contrib.metrics.cemba(knn, resolution=i)
        disp.append(d)
        pac.append(p)
    df = pd.DataFrame({
        "resolution": resolutions,
        "ARI": aris,
        "AMI": amis,
        "silhouette": silhouettes,
        "calinski": calinskis,
        "davies_bouldin": davies_bouldins,
        "dispersion": disp,
        "purity adjusted correlation": pac,
    })
    df.insert(0, 'dataset', "${name}")
    df.to_csv("result.csv.gz", index=False)
    """
}

process report {
    publishDir "${params.outdir}/ATAC/leiden", mode: 'copy'
    container 'kaizhang/scatac-bench:0.2.0'
    input:
      path(files, stageAs: '?.tsv.gz') 
    output:
      path "benchmark.tsv"

    """
    #!/usr/bin/env python3
    import pandas as pd

    files = ${files.collect { "\"${it}\"" }}

    df = pd.concat([pd.read_csv(f) for f in files])
    df.to_csv("benchmark.tsv", index=False, sep="\t")
    """
}

process plot {
    publishDir "${params.outdir}/ATAC/leiden", mode: 'copy'
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

    def scale(lst):
        min_value = min(lst)
        max_value = max(lst)
        # Handle the case where all values are the same
        if max_value == min_value:
            return [0.5 for _ in lst]
        return [(i - min_value) / (max_value - min_value) for i in lst]

    data = pd.read_csv('benchmark.tsv', sep="\t")
    for c in data.columns[2:]:
        data[c] = scale(data[c])
    data = data.melt(
        id_vars=['dataset', 'resolution'],
        value_vars=data.columns[2:],
        var_name='Metric',
        value_name='Scaled value',
    )

    n_label = np.unique(data['Metric'].to_numpy()).size
    n = -(np.unique(data['dataset'].to_numpy()).size // -3)
    shapes = ['o', 's', '^', 'v', 'D', 'X', '*', 'p', 'h', 'P', 'd', '<', '>']
    shapes = list(itertools.islice(itertools.cycle(shapes), n_label))
    colors = glasbey.create_palette(palette_size=n_label)

    ( ggplot(data)
        + aes('resolution', 'Scaled value', color='Metric', group='Metric')
        + geom_line()
        + facet_wrap('dataset', scales="free_y", ncol=3)
        + scale_fill_manual(colors)
        + scale_color_manual(colors)
        + theme_light(base_size=7, base_family="Arial")
        + theme(
            #axis_text_x=element_text(rotation=90, hjust=1),
            axis_text_x=element_text(rotation=90),
            figure_size=(2.5 * 3, 2.0 * n),
            subplots_adjust={'hspace':0.2, 'wspace':0.2},
            legend_key=element_rect(color = "white"),
        )
    ).save(filename='benchmark.pdf')
    """
}

workflow bench {
    take: datasets
    main:
        embedding = snapatac2(datasets)

        embedding
            | map({ [json(it[0]).data_name, it[1]] })
            | combine(
                datasets | map({ [json(it[0]).data_name, it[1]] }),
                by: 0
            )
            | leiden
            | toSortedList
            | report
            | plot
}
