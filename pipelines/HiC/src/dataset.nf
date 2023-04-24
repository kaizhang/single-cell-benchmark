nextflow.enable.dsl=2

process data_Lee {
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("Lee"), path("Lee.h5ad"), path("Lee.txt"), path("Lee.JSON"), path("Lee.chrom.sizes"), path("label_info.pickle")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/q7xyk/",
        "sha256:86b0d7934386443ce5ec363ba0058678b24aedc86af884261f619b985bcf514a",
        fname="Lee.h5ad",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/k3m8e/",
        "sha256:f1a4703c696a8f143b99dd9491602a0ffeebddfb698925b6ad4888edb87e9020",
        fname="_data.txt.xz",
        processor=pooch.Decompress(method='xz', name="Lee.txt"),
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/eaj4g/",
        "sha256:cbc342c70a85d735c6831950182c358f91c60cf98dc4745ca314fe062173ed04",
        fname="Lee.JSON",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/j7uke/",
        "sha256:7b2a4e4dcd483ddd3043bfdc685309c497bc4454877f89a04f6bbd41f1b42765",
        fname="Lee.chrom.sizes",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/h8wjn/",
        "sha256:baed8588fec74c0fc41112e3dead54d137dc35264eb4465aced61f39c07e82e3",
        fname="label_info.pickle",
        path="./",
    )
    """
}

process data_4dn_kim {
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("4DN_Kim"), path("4DN_data.h5ad"), path("4DN_data.txt"), path("4DN_config.JSON"), path("4DN_hg19.chrom.sizes"), path("label_info.pickle")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/s73be/",
        "sha256:a7c0d3dde2790ea1132ef292f4f0455514c9dc9c87645863cd0fab77ac758b19",
        fname="4DN_data.h5ad",
        path="./",
    )
    pooch.retrieve(
        "https://osf.io/download/eg3bk/",
        "sha256:d3e9c78e7992602d36593e1430de0bc9b0ac6b4ba3d168005d6177ec249bf91d",
        fname="_4DN_data.txt.xz",
        processor=pooch.Decompress(method='xz', name="4DN_data.txt"),
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
    pooch.retrieve(
        "https://osf.io/download/eaw8r/",
        "sha256:a9f11135925bd5871604a2865cc6d49afe415158bfaf16f4d254fa9c170a6ad3",
        fname="label_info.pickle",
        path="./",
    )
    """
}

workflow download_dataset {
    main:
        data = data_4dn_kim() | concat(
            data_Lee(),
        )
    emit:
        data
}
