nextflow.enable.dsl=2

include { dim_reduct_cosine as snapatac2 } from '../software/snapatac2.nf'

include { bench_clustering as bench
        } from '../../common/benchmark.nf' params(resultDir: "${params.outdir}/ATAC/dim_reduct")
include { knn_leiden_exact;
          knn_leiden_hora;
          kmeans;
        } from '../../common/algorithms/clustering.nf'

process leiden {
    container 'kaizhang/snapatac2:2.3.0'
    input:
      tuple val(name), val(method), path("reduced_dim.tsv"), path('data.h5ad')
    output:
      tuple val(name), val("leiden"), path("result.csv.gz")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    import pands as pd
    from sklearn.metrics import adjusted_rand_score, silhouette_score, adjusted_mutual_info_score
    embedding = np.genfromtxt("reduced_dim.tsv")
    knn = snap.pp.knn(embedding, method="exact", inplace=False)

    resolutions = list(np.arange(0.1, 3, 0.1))
    aris = []
    silhouettes = []
    for i in resolutions:
        clusters = snap.tl.leiden(knn, resolution=i, inplace=False)
        aris.append(adjusted_rand_score(clusters, adata.obs["cell_annotation"]))
        silhouettes.append(silhouette_score(embedding, adata.obs["cell_annotation"], sample_size = 10000))
    df = pd.DataFrame({"resolution": resolutions, "ARI": aris, "silhouette": silhouettes})
    df['dataset'] = "${name}"
    df.to_csv("result.csv.gz", index=False)
    """
}

workflow bench_leiden {
    take: datasets
    main:
        data = datasets | map { [it[0], it[1]] }
        snapatac2(data) | combine(data, by: 0) | leiden

}
