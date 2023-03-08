nextflow.enable.dsl=2

process data_Trapnell { 
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("Trapnell"), path("Trapnell.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/j5w29/",
        "sha256:25b6e3614e5dc2409b4351ea7e88db24ed9d02e902730f8ffcbb1fe1a85885c1",
        fname="Trapnell.h5ad",
        path="./",
    )
    """
}

process data_TrapnellTCC { 
    container 'kaizhang/scatac-bench:0.1.0'
    output:
      tuple val("TrapnellTCC"), path("TrapnellTCC.h5ad")
    """
    #!/usr/bin/env python
    import pooch
    pooch.retrieve(
        "https://osf.io/download/y54zb/",
        "sha256:94e090f5092fadb9f7e6c9bac25b554e4ce002ea90f0691502a864796cb65163",
        fname="TrapnellTCC.h5ad",
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
            data_Trapnell(),
            data_TrapnellTCC(),
        )

    emit:
        data
}
