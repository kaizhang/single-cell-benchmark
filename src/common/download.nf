nextflow.enable.dsl=2

include { json_string } from './utils.gvy'

process download_dataset {
    container 'kaizhang/scatac-bench:0.2.0'
    input:
      val(metadata)
    output:
      tuple val("${json_string(metadata)}"), path("data.h5ad")

    """
    #!/usr/bin/env python
    import pooch
    import json
    metadata = json.loads('${json_string(metadata)}')
    hash = metadata['data_hash'] if 'data_hash' in metadata else None
    pooch.retrieve(metadata['data_url'], hash, fname="data.h5ad", path="./")
    """
}

process download_genome {
    container 'kaizhang/scatac-bench:0.2.0'
    input:
      val(genome)

    output:
      tuple val(genome), path("*.decomp")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import os
    os.environ["SNAP_DATA_DIR"] = "./"

    if "${genome}" == "mm10":
        snap.genome.mm10.fetch_fasta()
    elif "${genome}" == "hg38":
        snap.genome.hg38.fetch_fasta()
    elif "${genome}" == "hg19":
        snap.genome.hg19.fetch_fasta()
    else:
        raise ValueError("Unknown genome: ${genome}")
    """
}

