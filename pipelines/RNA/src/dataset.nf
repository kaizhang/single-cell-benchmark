nextflow.enable.dsl=2

process data_Koh { 
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("Koh"), path("Koh.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/97pfy/",
        "sha256:184b2f045fd4067b9101084abc5488da20fee77848b696a6bf129c5ec904cbed",
        fname="Koh.h5ad",
        path="./",
    )
    """
}

process data_Kumar { 
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("Kumar"), path("Kumar.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/gu6qv/",
        "sha256:67ade87f04955bf7f78e54786ec3ffb25035c4d1b832e90de6f19556ffc4d103",
        fname="Kumar.h5ad",
        path="./",
    )
    """
}

process data_Zhengmix4eq {
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("Zhengmix4eq"), path("Zhengmix4eq.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/k8zsu/",
        "sha256:160e41a2a360b2a98c373085f1719c535fa594b9f521ad64e2de1084c7ba4ead",
        fname="Zhengmix4eq.h5ad",
        path="./",
    )
    """
}

process data_Zhengmix4uneq {
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("Zhengmix4uneq"), path("Zhengmix4uneq.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/xm7h6/",
        "sha256:28747c0dc5b01deeb8d83441705eb5110913adcb0e4a29541a535d9ce911c8bf",
        fname="Zhengmix4uneq.h5ad",
        path="./",
    )
    """
}

process data_Zhengmix8eq {
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("Zhengmix8eq"), path("Zhengmix8eq.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/mbs6t/",
        "sha256:bf27df9727bca093984c40a3b9730b6c803ca954f5fa9eef6a40511050cd8536",
        fname="Zhengmix8eq.h5ad",
        path="./",
    )
    """
}

workflow download_dataset {
    main:
        data = data_Zhengmix8eq() | concat(
            data_Zhengmix4eq(),
            data_Zhengmix4uneq(),
            data_Koh(),
            data_Kumar(),
        )

    emit:
        data
}
