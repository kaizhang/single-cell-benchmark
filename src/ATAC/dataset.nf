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

process data_BoneMarrow_Chen_2019_clean {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("BoneMarrow_Chen_2019_clean"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/a5nh7/",
        "sha256:aefb5b88798757057dd55ff89f2c580cb67b16ef4b12c091aea6373da717d6d6",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_BoneMarrow_Chen_2019_cov250 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("BoneMarrow_Chen_2019_cov250"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/6yms7/",
        "sha256:8125ecf60a6c2e3f69f8c0d3434aacd8f8cd4e5c33d9e49a3156895c06c5712f",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_BoneMarrow_Chen_2019_cov500 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("BoneMarrow_Chen_2019_cov500"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/n7djf/",
        "sha256:1ba87eca49b87c032aae85c0fa7e6800d64268694abea075773ebd5623f9ca35",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_BoneMarrow_Chen_2019_cov1000 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("BoneMarrow_Chen_2019_cov1000"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/367pn/",
        "sha256:871d837c0405e46d31e6821e5a55c8dd508de246d5f5d81ef6a652159c1300d4",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_BoneMarrow_Chen_2019_cov2500 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("BoneMarrow_Chen_2019_cov2500"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/28wby/",
        "sha256:d872f0a161fd838e9d1bb4c2830e8741e229da0a6b39820b53b73c254ab17d1d",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_BoneMarrow_Chen_2019_cov5000 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("BoneMarrow_Chen_2019_cov5000"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/dkzp9/",
        "sha256:3a631d0bea8d683c0ecbeede3028a6628c112d9904a3f42fb27daf7dca317dae",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_BoneMarrow_Chen_2019_noise_p2 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("BoneMarrow_Chen_2019_noise_p2"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/xtwh9/",
        "sha256:99769f835a4eb1daa1140b246079b5b53e7b507dc7f8f9796b2b820e6b21b852",
        fname="matrix.h5ad",
        path="./",
    )
    """
}

process data_BoneMarrow_Chen_2019_noise_p4 {
    container 'kaizhang/scatac-bench:0.2.0'
    output:
      tuple val("BoneMarrow_Chen_2019_noise_p4"), path("matrix.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/5g9jk/",
        "sha256:5aac8d5b6c730d18ced677bc8773e10c9cbb36b5eab888bb93b980d122d7ae73",
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

workflow download_simulated_dataset {
    main:
        data = data_BoneMarrow_Chen_2019_clean() | concat(
            data_BoneMarrow_Chen_2019_cov250(),
            data_BoneMarrow_Chen_2019_cov500(),
            data_BoneMarrow_Chen_2019_cov1000(),
            data_BoneMarrow_Chen_2019_cov2500(),
            data_BoneMarrow_Chen_2019_cov5000(),
            data_BoneMarrow_Chen_2019_noise_p2(),
            data_BoneMarrow_Chen_2019_noise_p4(),
        ) | combine(download_hg19())
    emit:
        data
}