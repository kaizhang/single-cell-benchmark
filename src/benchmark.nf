nextflow.enable.dsl=2

include { dim_reduct_snapatac_2;
          dim_reduct_snapatac_2_nystrom;
          dim_reduct_snapatac_2_kmeans_nystrom;
          dim_reduct_snapatac_1;
          dim_reduct_archr_1;
        } from './software'
include { import_dataset } from './dataset'

process benchmark_dim_reduct {
    input:
      tuple val(name), val(data), path(reduced_dim)
    output:
      tuple val(name), val(data), path("benchmark.txt")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from sklearn.metrics import adjusted_rand_score
    import numpy as np

    adata = snap.read("${data.anndata}")
    adata.obsm["X_spectral"] = np.loadtxt("${reduced_dim}")
    snap.tl.umap(adata)
    snap.pp.knn(adata)
    score = -1
    for i in np.arange(0, 2, 0.1):
        snap.tl.leiden(adata, resolution = i) 
        sc = adjusted_rand_score(adata.obs["leiden"], adata.obs["cell_annotation"])
        if sc > score: score = sc
    with open("benchmark.txt", "w") as f: print(score, file=f)
    """
}

/*
process generate_report {
    publishDir 'result'

    input:
      tuple val(name), val(data), path(reduced_dim)
    output:
      tuple val(name), val(data), path("benchmark.txt")

}
*/

workflow {
    datasets = import_dataset(Channel.fromList([
        "dataset/Real",
        "dataset/Synthetic",
    ])).flatten()
    
    benchData = datasets.map { data -> tuple(data, 50) }

    nystromBenchData = datasets.map { data -> tuple(data, 50, 0.01) }.concat(
        datasets.map { data -> tuple(data, 50, 0.1) },
        datasets.map { data -> tuple(data, 50, 0.2) },
        datasets.map { data -> tuple(data, 50, 0.5) },
        datasets.map { data -> tuple(data, 50, 1.0) }
    )

    runResult = dim_reduct_snapatac_2(benchData).concat(
        dim_reduct_archr_1(benchData),
        dim_reduct_snapatac_1(benchData),
        dim_reduct_snapatac_2_nystrom(nystromBenchData),
    )

    benchmark_dim_reduct(runResult)
        .map { tuple(name, data, result) -> tuple(data.name, name, result) }
        .collect()
}