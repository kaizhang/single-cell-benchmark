nextflow.enable.dsl=2

process download_hg38 {
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

process import_dataset {
    cache false
    input:
      tuple val(type), val(dir)
    output:
      val list

    exec:
    list = []
    file("../../" + dir).eachDir { item -> 
        list << tuple(type, item.getBaseName(), item + "/matrix.h5ad")
    }
}
