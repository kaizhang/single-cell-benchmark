nextflow.enable.dsl=2

include { is_included; add_meta; json; } from '../../common/utils.gvy'

process dim_reduct_scbasset {
    container 'kaizhang/scbasset:latest'
    tag "${json(metadata).data_name}"
    errorStrategy 'ignore'
    containerOptions '--nv'
    label "gpu"

    when: is_included("scbasset", params.method_include, params.method_exclude)

    input:
      tuple val(metadata), path("data.h5ad"), path("fasta.fa")
    output:
      tuple val("${add_meta(metadata, 'method', 'scBasset')}"), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import pandas as pd
    import numpy as np
    import h5py
    import gc
    import json
    from scbasset.utils import *

    def assign_free_gpus(threshold_vram_usage=1500):
        import subprocess
        import os
        def _check():
            # Get the list of GPUs via nvidia-smi
            smi_query_result = subprocess.check_output(
                "nvidia-smi -q -d Memory | grep -A4 GPU", shell=True
            )
            # Extract the usage information
            gpu_info = smi_query_result.decode("utf-8").split("\\n")
            gpu_info = list(filter(lambda info: "Used" in info, gpu_info))
            gpu_info = [
                int(x.split(":")[1].replace("MiB", "").strip()) for x in gpu_info
            ]  # Remove garbage
            # Keep gpus under threshold only
            free_gpus = [
                str(i) for i, mem in enumerate(gpu_info) if mem < threshold_vram_usage
            ]
            gpus_to_use = ",".join(free_gpus)
            return gpus_to_use

        gpus_to_use = _check()
        if not gpus_to_use:
            raise RuntimeError("No free GPUs found")
        os.environ["CUDA_VISIBLE_DEVICES"] = gpus_to_use
    

    metadata = json.loads('${metadata}')
    adata = ad.read_h5ad("data.h5ad")
    adata.var['chr'] = [x.split(':')[0] for x in adata.var_names]
    adata.var['start'] = [int(x.split(':')[1].split('-')[0]) for x in adata.var_names]
    adata.var['end'] = [int(x.split(':')[1].split('-')[1]) for x in adata.var_names]

    # Preprocess data
    train_ids, test_ids, val_ids = split_train_test_val(np.arange(adata.shape[1]))
    ad_train = adata[:,train_ids]
    ad_test = adata[:,test_ids]
    ad_val = adata[:,val_ids]
    make_h5_sparse(ad_train, 'train_seqs.h5', "fasta.fa")
    make_h5_sparse(ad_test, 'test_seqs.h5', "fasta.fa")
    make_h5_sparse(ad_val, 'val_seqs.h5', "fasta.fa")

    # Train model
    bottleneck_size = 30
    batch_size = 128
    lr = 0.01
    epochs = 500
    
    train_data = 'train_seqs.h5'
    val_data = 'val_seqs.h5'
    n_cells = adata.shape[0]
    
    # convert to csr matrix
    m = adata.X.tocoo().transpose().tocsr()
    print_memory()
    
    m_train = m[train_ids,:]
    m_val = m[val_ids,:]
    del m
    gc.collect()

    assign_free_gpus()
    
    # generate tf.datasets
    train_ds = tf.data.Dataset.from_generator(
        generator(train_data, m_train), 
        output_signature=(
             tf.TensorSpec(shape=(1344,4), dtype=tf.int8),
             tf.TensorSpec(shape=(n_cells), dtype=tf.int8),
        )
    ).shuffle(2000, reshuffle_each_iteration=True).batch(128).prefetch(tf.data.AUTOTUNE)
    
    val_ds = tf.data.Dataset.from_generator(
        generator(val_data, m_val), 
        output_signature=(
             tf.TensorSpec(shape=(1344,4), dtype=tf.int8),
             tf.TensorSpec(shape=(n_cells), dtype=tf.int8),
        )
    ).batch(128).prefetch(tf.data.AUTOTUNE)
    
    # build model
    if 'batch_key' in metadata:
        batch_key = metadata['batch_key']
        batch_m = pd.get_dummies(adata.obs[batch_key])
        model = make_model_bc(bottleneck_size, n_cells, batch_m, l2_1=0, l2_2=0)
    else:
        model = make_model(bottleneck_size, n_cells)

    # compile model
    loss_fn = tf.keras.losses.BinaryCrossentropy()
    optimizer = tf.keras.optimizers.Adam(learning_rate=lr,beta_1=0.95,beta_2=0.9995)
    model.compile(loss=loss_fn, optimizer=optimizer,
                  metrics=[tf.keras.metrics.AUC(curve='ROC', multi_label=True),
                           tf.keras.metrics.AUC(curve='PR', multi_label=True)])

    callbacks = [
        tf.keras.callbacks.TensorBoard('./'),
        tf.keras.callbacks.ModelCheckpoint('best_model.h5', save_best_only=True, 
                                           save_weights_only=True, monitor='auc', mode='max'),
        tf.keras.callbacks.EarlyStopping(monitor='auc', min_delta=1e-6, 
                                         mode='max', patience=50, verbose=1),
    ]
    
    # train the model
    model.fit(
        train_ds,
        epochs=epochs,
        callbacks=callbacks,
        validation_data=val_ds)

    proj = get_cell_embedding(model) # get_cell_embedding function
    np.savetxt("reduced_dim.tsv", proj, delimiter="\t")
    """
}
