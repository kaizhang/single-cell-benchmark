nextflow.enable.dsl=2

process dim_reduct_jaccard {
    container 'kaizhang/snapatac2:2.3.0'
    cpus 4
    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val("SnapATAC2 (jaccard)"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)
    _, evecs = snap.tl.spectral(adata, features=None, inplace=False, distance_metric="jaccard")
    np.savetxt("reduced_dim.tsv", evecs, delimiter="\t")
    """
}

process dim_reduct_cosine {
    container 'kaizhang/snapatac2:2.3.0'
    cpus 4
    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val("SnapATAC2 (cosine)"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    seed = 1
    adata = snap.read("data.h5ad", backed=None)
    _, evecs = snap.tl.spectral(adata, features=None, inplace=False, distance_metric="cosine")
    np.savetxt("reduced_dim.tsv", evecs, delimiter="\t")
    """
}

process dim_reduct_svd {
    container 'kaizhang/snapatac2:2.3.0'
    cpus 4
    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val("SnapATAC2 (svd)"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import scipy as sp
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)

    feature_weights = snap.tl._embedding.idf(adata)
    X = adata.X @ sp.sparse.diags(feature_weights)
    s = 1 / np.sqrt(np.ravel(sp.sparse.csr_matrix.power(X, 2).sum(axis = 1)))
    X = sp.sparse.diags(s) @ X
    D = np.ravel(X @ X.sum(axis = 0).T)
    X = sp.sparse.diags(1 / np.sqrt(D)) @ X
    evecs,evals,_ = sp.sparse.linalg.svds(X, k=30)
    ix = evals.argsort()[::-1]
    evals = evals[ix]
    evecs = evecs[:, ix]
    idx = [i for i in range(evals.shape[0]) if evals[i] > 0]
    evals = evals[idx]
    evecs = evecs[:, idx] * evals
    np.savetxt("reduced_dim.tsv", evecs, delimiter="\t")
    """
}


process dim_reduct_nystrom {
    container 'kaizhang/snapatac2:2.3.0'
    input:
      tuple val(name), path("data.h5ad"), val(fraction)
    output:
      tuple val(name), val(fraction), val("SnapATAC2"), path("reduced_dim.tsv")
    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    seed = 1
    adata = snap.read("data.h5ad", backed='r')
    result = snap.tl.spectral(adata, features=None, chunk_size=500,
        distance_metric="cosine",
        sample_size=${fraction}, random_state=0,
        inplace=False
    )
    np.savetxt("reduced_dim.tsv", result[1], delimiter="\t")
    """
}

process dim_reduct_nystrom2 {
    container 'kaizhang/snapatac2:2.3.0'
    input:
      tuple val(name), path("data.h5ad"), val(fraction)
    output:
      tuple val(name), val(fraction), val("SnapATAC2 (degree)"), path("reduced_dim.tsv")
    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    seed = 1
    adata = snap.read("data.h5ad", backed='r')
    result = snap.tl.spectral(adata, features=None, chunk_size=500,
        distance_metric="cosine",
        sample_size=${fraction}, random_state=0,
        sample_method="degree",
        inplace=False
    )
    np.savetxt("reduced_dim.tsv", result[1], delimiter="\t")
    """
}

process end_to_end_snapatac_2 {
    container 'kaizhang/snapatac2:2.3.0'
    input:
      val(data)
    output:
      tuple val("SnapATAC2"), val(data), path("*.dim"), path("*.cluster")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from tempfile import TemporaryDirectory
    import numpy as np
    import pandas as pd
    if "${data.name}".startswith("human"):
        genome = snap.genome.hg38
    elif "${data.name}".startswith("mouse"):
        genome = snap.genome.mm10
    else:
        genome = None
    data = snap.pp.import_data(
        fragment_file = "${data.fragments_name_sorted}",
        gff_file = "${data.gene_annotations}",
        chrom_size = genome,
        min_tsse = 0,
        min_num_fragments = 0,
    )
    snap.pp.make_tile_matrix(data)
    snap.pp.select_features(data)
    snap.tl.spectral(data)

    for d in [5, 15, 30]:
        pd.DataFrame(
            data.obsm["X_spectral"][:, 0:d], index = data.obs_names
        ).to_csv(str(d) + ".dim", header=False, sep="\t")
        snap.pp.knn(data, use_dims = d)
        clusters = []
        for i in np.arange(0.1, 1.501, 0.1):
            snap.tl.leiden(data, resolution=i)
            clusters.append(data.obs['leiden'])
        df = pd.DataFrame(np.array(clusters).T, index = data.obs_names)
        df.to_csv(str(d) + ".cluster", header=False, sep="\t")
    """
}

process end_to_end_snapatac_2_cosine {
    container 'kaizhang/snapatac2:2.3.0'
    input:
      val(data)
    output:
      tuple val("SnapATAC2 (cosine)"), val(data), path("*.dim"), path("*.cluster")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    from tempfile import TemporaryDirectory
    import numpy as np
    import pandas as pd
    if "${data.name}".startswith("human"):
        genome = snap.genome.hg38
    elif "${data.name}".startswith("mouse"):
        genome = snap.genome.mm10
    else:
        genome = None
    data = snap.pp.import_data(
        fragment_file = "${data.fragments_name_sorted}",
        gff_file = "${data.gene_annotations}",
        chrom_size = genome,
        min_tsse = 0,
        min_num_fragments = 0,
    )
    snap.pp.make_tile_matrix(data)
    snap.pp.select_features(data)
    snap.tl.spectral(data, distance_metric="cosine")

    for d in [5, 15, 30]:
        pd.DataFrame(
            data.obsm["X_spectral"][:, 0:d], index = data.obs_names
        ).to_csv(str(d) + ".dim", header=False, sep="\t")
        snap.pp.knn(data, use_dims = d)
        clusters = []
        for i in np.arange(0.1, 1.501, 0.1):
            snap.tl.leiden(data, resolution=i)
            clusters.append(data.obs['leiden'])
        df = pd.DataFrame(np.array(clusters).T, index = data.obs_names)
        df.to_csv(str(d) + ".cluster", header=False, sep="\t")
    """
}
