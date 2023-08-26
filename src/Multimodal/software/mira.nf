nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct {
    container 'kaizhang/mira:2.1.0'
    tag "${json(metadata).data_name}"
    cpus 4
    errorStrategy 'ignore'
    containerOptions '--nv'
    label "gpu"

    when: is_included("mira", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path(h5ads, stageAs: "?.h5ad")
    output:
      tuple val("${add_meta(metadata, 'method', 'MIRA')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import mira
    import anndata as ad
    import scanpy as sc
    import numpy as np

    def find_elbow(x):
        n = len(x)
        marks = []
        saturation = 0.01
        accum_gap = 0
        for i in range(1, n):
            gap = x[i-1] - x[i]
            accum_gap += gap
            if gap > saturation * accum_gap:
                marks.append(i)
        return min(n - 1, max(marks) + 1)

    rna_data = sc.read('2.h5ad')

    rna_model = mira.topics.make_model(
        rna_data.n_obs, rna_data.n_vars,
        feature_type = 'expression',
        counts_layer='counts',
    )
    rna_model.get_learning_rate_bounds(rna_data)
    rna_model.set_learning_rates(1e-3, 0.1)

    topic_contributions = mira.topics.gradient_tune(rna_model, rna_data)
    rna_model.set_params(num_topics=find_elbow(topic_contributions)).fit(rna_data)

    atac_data = sc.read('1.h5ad')

    atac_model = mira.topics.make_model(
        atac_data.n_obs, atac_data.n_vars,
        feature_type = 'accessibility',
    )
    atac_model.get_learning_rate_bounds(atac_data)
    atac_model.set_learning_rates(1e-3, 0.1)

    topic_contributions = mira.topics.gradient_tune(atac_model, atac_data)
    atac_model.set_params(num_topics=find_elbow(topic_contributions)).fit(atac_data)

    atac_model.predict(atac_data)
    rna_model.predict(rna_data)
    rna_model.get_umap_features(rna_data, box_cox=0.25)
    atac_model.get_umap_features(atac_data, box_cox=0.25)
    rna_data_new, _ = mira.utils.make_joint_representation(rna_data, atac_data)
    rna_data_new = rna_data_new[rna_data.obs_names].copy()

    embedding = rna_data_new.obsm['X_joint_umap_features']
    np.savetxt("reduced_dim.tsv", embedding, delimiter="\t")
    """
}