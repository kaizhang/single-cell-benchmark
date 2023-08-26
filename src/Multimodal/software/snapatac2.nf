nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct {
    container 'kaizhang/snapatac2:2.3.1'
    tag "${json(metadata).data_name}"
    errorStrategy 'ignore'
    cpus 4

    when: is_included("snapatac2", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path(h5ads, stageAs: "?.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'SnapATAC2')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np

    files = "${h5ads}".split()
    adatas = [snap.read(f, backed=None) for f in files]
    embedding = snap.tl.multi_spectral(adatas, features=None)[1]
    np.savetxt("reduced_dim.tsv", embedding, delimiter="\t")
    """
}