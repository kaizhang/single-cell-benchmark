nextflow.enable.dsl=2

process data_4dn_kim {
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("4DN_Kim"), path("4DN_data.h5ad"), path("4DN_data.txt"), path("4DN_config.JSON"), path("4DN_hg19.chrom.sizes")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/x7fcz/",
        "sha256:d4e099b17cea4172bcfe980fb8b9fb6d2aa7b61a45566ed94ba23f40c2dd0640",
        fname="4DN_data.h5ad",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/k583f/",
        "sha256:9f851dd8ce9f8a1b23d1c91b76d8a7a1698d821f3e5504b8e6ddec356f34e487",
        fname="_4DN_data.txt.gz",
        processor=pooch.Decompress(method='gzip', name="4DN_data.txt"),
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/apm6w/",
        "sha256:f9d8a9be76aabbb17fddb4047773291fabdf4a0a96d6ccea5c6233e36dc0873d",
        fname="4DN_config.JSON",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/j7uke/",
        "sha256:7b2a4e4dcd483ddd3043bfdc685309c497bc4454877f89a04f6bbd41f1b42765",
        fname="4DN_hg19.chrom.sizes",
        path="./",
    )
    """
}

workflow download_dataset {
    main:
        data_4dn_kim()
    emit:
        data_4dn_kim.out
}
