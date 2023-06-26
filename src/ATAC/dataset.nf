nextflow.enable.dsl=2

process download_mm10 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      path("*.decomp")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import os
    os.environ["SNAP_DATA_DIR"] = "./"
    snap.genome.mm10.fetch_fasta()
    """
}

process download_hg38 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      path("*.decomp")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import os
    os.environ["SNAP_DATA_DIR"] = "./"
    snap.genome.hg38.fetch_fasta()
    """
}

process download_hg19 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      path("*.decomp")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import os
    os.environ["SNAP_DATA_DIR"] = "./"
    snap.genome.hg19.fetch_fasta()
    """
}

process data_BoneMarrow_Chen_2019 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("BoneMarrow_Chen_2019"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/r3xuy/",
        "sha256:b6a559e7be058bdf174f00f1a6f10749ebf025f0995bbcca2288b9ed0fd51b0e",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_Buenrostro_2018 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("Buenrostro_2018"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/z52xh/",
        "sha256:03561b4de5cb9e84e188e544ae7c4836e0e6b5e880c22e9e8670e1f38a27a6cc",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_10x_Brain5k {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("10x_Brain5k"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/r2gj3/",
        "sha256:2495bf7c168348daf69b62b38042305712e047a626e4733478bad6b3ad96aecf",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_10x_PBMC10k {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("10x_PBMC10k"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/bvuk7/",
        "sha256:540e8e3bf98b5e0ac547165a699fc8a3c121458dcef2b4eec48d6a1d17e4f79c",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_Chen_NBT_2019 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("Chen_NBT_2019"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/j2vg4/",
        "sha256:458bd9d1d2648d82c1935b4578ddb5f3b16856a6c2365b7f94aa1152d5fbe745",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_GSE194122 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("GSE194122"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/dk5zf/",
        "sha256:ee235d584390e65a5f4484020e544d5bae0d1634257ec9f49db8a2cb025c9e4d",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_Ma_Cell_2020 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("Ma_Cell_2020"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/xv3m9/",
        "sha256:5bb29c3bd20ffc738c0e041d6dac6f29617afd3d0a146600959524ca3a2325ec",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_Trevino_Cell_2021 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("Trevino_Cell_2021"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/dfsp7/",
        "sha256:927fff9d06a7dcf790f9c52a528db9d042bcd19d0112a9a40b55b5c922717aeb",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_Yao_Nature_2021 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("Yao_Nature_2021"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/86cx9/",
        "sha256:36db6624adf9ae3ac49c41b4ab9d5abd7c3d1ca127abae2ae36a2813452ef5e3",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_Zemke_2023_human {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("Zemke_2023_human"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/qpx2f/",
        "sha256:f6f3914b22c5f76785c4cac6e491f1a3dda208814769092c33be8d99c98e3454",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_Zemke_2023_mouse {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("Zemke_2023_mouse"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/3pexa/",
        "sha256:38e2b2daee63dfe318ff49a63480e78cf5e37a54ee62b6ef47be42a92221ef99",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_Zhang_Cell_2021_GI {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("Zhang_Cell_2021"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/k6tys/",
        "sha256:046a33f87f573f2730df107ae36bcb7ef78bcd7295e5a4bcaa50b14df73e3bf7",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

workflow download_dataset {
    main:
        hg19 = data_BoneMarrow_Chen_2019() | concat(
            data_Buenrostro_2018(),
        ) | combine(download_hg19())

        hg38 = data_Trevino_Cell_2021() | concat(
            //data_Zhang_Cell_2021_GI(),
            data_10x_PBMC10k(),
            data_GSE194122(),
            data_Zemke_2023_human(),
        ) | combine(download_hg38())

        mm10 = data_10x_Brain5k() | concat(
            data_Chen_NBT_2019(),
            data_Ma_Cell_2020(),
            data_Yao_Nature_2021(),
            data_Zemke_2023_mouse(),
        ) | combine(download_mm10())

        data = hg19.concat(hg38, mm10)
    emit:
        data
}