nextflow.enable.dsl=2

process dim_reduct_jaccard {
    //container 'kaizhang/snapatac2:1.99.99.7'
    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val("SnapATAC2 (jaccard)"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("data.h5ad", backed=None)
    _, evecs = snap.tl.spectral(adata, features=None, inplace=False)
    np.savetxt("reduced_dim.tsv", evecs, delimiter="\t")
    """
}

process dim_reduct_cosine {
    //container 'kaizhang/snapatac2:1.99.99.7'
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

process dim_reduct_nystrom {
    //container 'kaizhang/snapatac2:1.99.99.7'
    input:
      val(data)
    output:
      tuple val("SnapATAC2"), val(data), path("${data.name}_snapatac2_nystrom_${data.samplingFraction}_${data.randomSeed}_reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("${data.anndata}", backed=None)
    result = snap.tl.spectral(adata, features=None, chunk_size=500,
        sample_size=${data.samplingFraction}, random_state=${data.randomSeed},
        inplace=False
    )
    output = "${data.name}_snapatac2_nystrom_${data.samplingFraction}_${data.randomSeed}_reduced_dim.tsv"
    np.savetxt(output, result[1], delimiter="\t")
    """
}

process dim_reduct_cosine_nystrom {
    //container 'kaizhang/snapatac2:1.99.99.7'
    input:
      val(data)
    output:
      tuple val("SnapATAC2_cosine"), val(data), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("${data.anndata}", backed=None)
    result = snap.tl.spectral(
        adata,
        distance_metric="cosine",
        features=None,
        chunk_size=500,
        sample_size=${data.samplingFraction},
        random_state=${data.randomSeed},
        inplace=False,
    )
    np.savetxt("reduced_dim.tsv", result[1], delimiter="\t")
    """
}

process end_to_end_snapatac_2 {
    //container 'kaizhang/snapatac2:1.99.99.7'
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
    //container 'kaizhang/snapatac2:1.99.99.7'
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
