nextflow.enable.dsl=2

process download_dataset {
    container 'kaizhang/scatac-bench:0.2.1'
    output:
      path("*.tsv.gz")

    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/y3r5c/",
        "sha256:7faeb0c1c11c2fd594789e07177a8ce350a473bd6d421b2cd82c14b09bf13aa0",
        fname="0.tsv.gz",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/63dafaa16946a002797a5805/",
        "sha256:49e596398f77b2131096b241238f8121fc6ba132875b96cb69e31c011a9505cb",
        fname="1.tsv.gz",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/yv2zw/",
        "sha256:f285adc654a7a5b44f60c8235753f4ad470b172e42a7bf15aefa054e41c14e32",
        fname="2.tsv.gz",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/63daff0a1e96860298b25cae/",
        "sha256:2b3e6ec5f0154452f9c9509e8b8e722fbe3e6f80f0b51849b243e94bff3d09b9",
        fname="3.tsv.gz",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/63db02291e96860297b25ce8/",
        "sha256:83483a0c8961626529c8f86b3732546c31b92738e624ef86337d0801980b25bf",
        fname="4.tsv.gz",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/63db04ed78a623028999e07d/",
        "sha256:655170234e2ef3c2ade29bc0a9d83d10acf73762a8f3b233bc67d1786c6ebab7",
        fname="5.tsv.gz",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/63db06786946a0027d7a5503/",
        "sha256:585d7695abb63deb236692c93b2408b90718367656125716b58bc38808690b15",
        fname="6.tsv.gz",
        path="./",
    )
    """
}

process merge {
    container 'kaizhang/scatac-bench:0.2.1'
    input:
      path(files)
    output:
      tuple val(320000), path("nsrt.tsv.gz"), path("srt.tsv.gz")

    """
    zcat ${files.join(' ')} | gzip > nsrt.tsv.gz
    zcat nsrt.tsv.gz | sort -T ./ -S 4G -k1,1 -k2,2n -k3,3n | bgzip > srt.tsv.gz
    """
}


process subsample {
    container 'kaizhang/scatac-bench:0.2.1'
    input:
      tuple val(n), path("fragment.tsv.gz"), path("_")
      each m
    output:
      tuple val(m), path("nsrt.tsv.gz"), path("srt.tsv.gz")

    """
    subsample fragment.tsv.gz nsrt.tsv.gz $m
    zcat nsrt.tsv.gz | sort -T ./ -S 4G -k1,1 -k2,2n -k3,3n | bgzip > srt.tsv.gz
    """
}

process prep_h5ad {
    container 'kaizhang/scatac-bench:0.2.1'
    tag "$name"
    input:
      tuple val(name), path("nsrt.tsv.gz"), path("_")
    output:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env python

    import snapatac2 as snap
    data = snap.pp.import_data(
        fragment_file="nsrt.tsv.gz",
        genome=snap.genome.hg38,
        min_tsse=0,
        min_num_fragments=0,
    )
    data = snap.pp.add_tile_matrix(data, inplace=False)
    snap.pp.select_features(
        data,
        n_features=500000000,
        filter_lower_quantile=0,
        filter_upper_quantile=0,
    )
    data[:, data.var["selected"]].write("data.h5ad", compression="gzip")
    """
}