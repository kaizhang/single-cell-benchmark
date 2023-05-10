nextflow.enable.dsl=2

process dim_reduct_higashi {
    container 'kaizhang/higashi:latest'
    tag "$name"
    errorStrategy 'ignore'
    containerOptions '--nv'
    label "gpu"

    input:
      tuple val(name), val(_), path("data.txt"), path("config.JSON"), path("chrom.sizes"), path("label_info.pickle")
    output:
      tuple val(name), val('Higashi'), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    from higashi.Higashi_wrapper import *
    import pickle
    import numpy as np

    higashi_model = Higashi("config.JSON")
    higashi_model.config['genome_reference_path'] = "chrom.sizes"
    higashi_model.config['dimensions'] = 30
    higashi_model.config['input_format'] = 'higashi_v1'
    higashi_model.process_data()
    higashi_model.prep_model()
    higashi_model.train_for_embeddings()
    latent = higashi_model.fetch_cell_embeddings()
    np.savetxt("reduced_dim.tsv", latent, delimiter="\t")
    """
}
