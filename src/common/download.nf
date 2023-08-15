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