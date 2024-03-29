nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct_pycistopic {
    container 'kaizhang/pycistopic:latest'
    tag "${json(metadata).data_name}"
    cpus 8
    errorStrategy 'ignore'

    when: is_included("cistopic", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("data.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'cisTopic')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import numpy as np
    from pycisTopic.cistopic_class import *

    data = ad.read("data.h5ad")
    obj = create_cistopic_object(
        data.X.T.tocsr(),
        cell_names=list(data.obs_names),
        region_names=list(data.var_names),
    )
    models = run_cgs_models(
        obj,
        n_topics=[5,10,15,20,25,30],
        n_iter=300,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        n_cpu=6,
        save_path=None,
        _temp_dir="/tmp",
    )
    model = evaluate_models(
        models,
        return_model=True,
        plot=False,
        plot_metrics=False,
    )
    np.savetxt("reduced_dim.tsv", model.cell_topic.T.to_numpy(), delimiter="\t")
    """
}