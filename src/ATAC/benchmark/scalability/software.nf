nextflow.enable.dsl=2

include { is_included } from '../../../common/utils.gvy'

////////////////////////////////////////////////////////////////////////////////
// SnapATAC2                                                                  //
////////////////////////////////////////////////////////////////////////////////

process preproc_snapatac2 {
    container 'kaizhang/snapatac2:2.3.1'
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'

    when: is_included("snapatac2", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("nsrt.tsv.gz"), path("srt.tsv.gz")
    output:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    data = snap.pp.import_data(
        fragment_file="nsrt.tsv.gz",
        genome=snap.genome.hg38,
        file="data.h5ad",
        min_tsse=0,
        min_num_fragments=0,
    )
    snap.pp.add_tile_matrix(data)
    """
}

process dim_reduct_snapatac2 {
    container 'kaizhang/snapatac2:2.3.1'
    stageInMode "copy"
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'

    when: is_included("snapatac2", params.method_include, params.method_exclude)

    input:
      tuple val(name), path(data)

    """
    #!/usr/bin/env python3
    import snapatac2 as snap

    data = snap.read("$data")
    snap.tl.spectral(data, features=None, distance_metric="cosine")
    """
}

process clust_snapatac2 {
    container 'kaizhang/snapatac2:2.3.1'
    stageInMode "copy"
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'

    when: is_included("snapatac2", params.method_include, params.method_exclude)

    input:
      tuple val(name), path(data)
    output:
      tuple val(name), path(data)

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    data = snap.read("$data")
    snap.pp.knn(data)
    snap.tl.leiden(data)
    """
}

////////////////////////////////////////////////////////////////////////////////
// SnapATAC R version                                                         //
////////////////////////////////////////////////////////////////////////////////

process preproc_snapatac_1 {
    container 'kaizhang/snapatac:1.0'
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    errorStrategy 'ignore'

    when: is_included("snapatac", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("name_srt.bed.gz"), path("srt.bed.gz"), path("anno.gff3.gz")
    output:
      tuple val(name), path("data.rds")

    """
    #!/usr/bin/env Rscript
    library(SnapATAC);
    utils::download.file(
        "https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/mm10.chrom.sizes",
        "chrom_size.txt"
    )
    cmd <- paste(
        "snaptools",
        "snap-pre",
        "--input-file=name_srt.bed.gz",
        "--output-snap=data.snap",
        "--genome-name=mm10",
        "--genome-size=chrom_size.txt",
        "--min-mapq=30",
        "--min-flen=0",
        "--max-flen=2000",
        "--keep-chrm=FALSE",
        "--keep-single=TRUE",
        "--keep-secondary=True",
        "--overwrite=True",
        "--max-num=200000000",
        "--min-cov=0"
    )
    system(cmd)

    cmd <- paste(
        "snaptools",
        "snap-add-bmat",
        "--snap-file=data.snap",
        "--bin-size-list 500"
    )
    system(cmd)

    x.sp <- createSnap(
        file="data.snap",
        sample="sample",
        num.cores=4,
    )
    x.sp <- addBmatToSnap(x.sp, bin.size=500, num.cores=4)
    x.sp <- makeBinary(x.sp, mat="bmat")
    saveRDS(x.sp, file="data.rds")
    """
}

process dim_reduct_snapatac {
    container 'kaizhang/snapatac:1.0'
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    errorStrategy 'ignore'

    when: is_included("snapatac", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env Rscript
    library("Matrix")
    library("rhdf5")
    set.seed(2022)

    file = H5Fopen("data.h5ad", flags = "H5F_ACC_RDONLY")
    x <- H5Gopen(file, "X")
    data <- sparseMatrix(j = x\$indices, p = x\$indptr,
        x = c(x\$data), index1=F, repr = "C", dims=H5Aread(H5Aopen(x, "shape"))
    )

    x.sp <- SnapATAC::newSnap()
    x.sp@bmat <- data
    x.sp <- SnapATAC::makeBinary(x.sp, mat="bmat")
    SnapATAC::runDiffusionMaps(
            obj=x.sp,
            input.mat="bmat", 
            num.eigs=30
    )
    """
}

process clust_snapatac {
    container 'kaizhang/snapatac:1.0'
    stageInMode "copy"
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    errorStrategy 'ignore'

    when: is_included("snapatac", params.method_include, params.method_exclude)

    input:
      tuple val(name), path(data)
    output:
      tuple val(name), path("data.rds")

    """
    #!/usr/bin/env Rscript
    library(SnapATAC)
    x.sp <- readRDS("${data}")
    x.sp <- runKNN(
      obj=x.sp,
      eigs.dims=1:50,
      k=15,
    )
    x.sp <- runCluster(
      obj=x.sp,
      tmp.folder=tempdir(),
      louvain.lib="R-igraph",
      seed.use=10
    )
    x.sp@metaData\$cluster <- x.sp@cluster
    save(x.sp, file="data.rds")
    """
}


////////////////////////////////////////////////////////////////////////////////
// ArchR                                                                      //
////////////////////////////////////////////////////////////////////////////////
process  preproc_archr {
    container 'kaizhang/archr:1.0.1'
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    errorStrategy 'ignore'

    when: is_included("archr", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("nsrt.tsv.gz"), path("srt.tsv.gz")
    output:
      tuple val(name), path("archr")

    """
    #!/usr/bin/env Rscript
    library("ArchR")
    library("parallel")
    set.seed(1)
    addArchRGenome("hg38")

    ArrowFiles <- createArrowFiles(
        inputFiles = c("srt.tsv.gz"),
        sampleNames = c("sample"),
        minTSS = 0,
        minFrags = 0, 
        maxFrags = 1e+10,
        excludeChr = c("chrM"),
        addTileMat = T,
        addGeneScoreMat = F,
    )
    proj <- ArchRProject(
        ArrowFiles = ArrowFiles, 
        outputDirectory = "archr",
        copyArrows = T,
    )
    saveArchRProject(ArchRProj = proj, outputDirectory = "archr", load = F)
    """
}

process dim_reduct_archr {
    container 'kaizhang/archr:1.0.1'
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    errorStrategy 'ignore'

    when: is_included("archr", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env Rscript
    library("ArchR")
    library("Matrix")
    library("rhdf5")
    set.seed(2022)
    file = H5Fopen("data.h5ad", flags = "H5F_ACC_RDONLY")
    x <- H5Gopen(file, "X")
    data <- sparseMatrix(i = x\$indices, p = x\$indptr,
        x = c(x\$data), index1=F, repr = "C", dims=rev(H5Aread(H5Aopen(x, "shape")))
    )
    ArchR:::.computeLSI(mat = data,
        LSIMethod = 2,
        scaleTo = 10^4,
        nDimensions = 30,
        binarize = TRUE, 
        outlierQuantiles = NULL,
        seed = 1
    )
    """
}

process clust_archr {
    container 'kaizhang/archr:1.0.1'
    stageInMode "copy"
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    errorStrategy 'ignore'

    when: is_included("archr", params.method_include, params.method_exclude)

    input:
      tuple val(name), path(archr)
    output:
      tuple val(name), path("archr")

    """
    #!/usr/bin/env Rscript
    library("ArchR")
    library("parallel")
    set.seed(1)
    proj <- loadArchRProject(path = "$archr", force = FALSE, showLogo = F)
    proj <- addClusters(input = proj, reducedDims = "IterativeLSI", force=T)
    saveArchRProject(ArchRProj = proj, outputDirectory = "archr", load = F)
    """
}


////////////////////////////////////////////////////////////////////////////////
// CisTopic                                                                   //
////////////////////////////////////////////////////////////////////////////////
process dim_reduct_pycistopic {
    container 'kaizhang/pycistopic:latest'
    tag "$name"
    cpus 4
    memory '120 GB'
    errorStrategy 'ignore'

    when: is_included("cistopic", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("data.h5ad")

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
    run_cgs_models(
        obj,
        n_topics=[30],
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
    """
}


////////////////////////////////////////////////////////////////////////////////
// PeakVI                                                                     //
////////////////////////////////////////////////////////////////////////////////
process dim_reduct_peakvi {
    container 'kaizhang/scvi-tools:0.19.0'
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    containerOptions '--nv'
    label "gpu"
    errorStrategy 'ignore'

    when: is_included("peakvi", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import scvi
    scvi.settings.seed = 0
    data = ad.read("data.h5ad")
    scvi.model.PEAKVI.setup_anndata(data)
    pvi = scvi.model.PEAKVI(data, n_latent=30)
    pvi.train(max_epochs=10, early_stopping=False)
    """
}

////////////////////////////////////////////////////////////////////////////////
// Signac                                                                     //
////////////////////////////////////////////////////////////////////////////////
process dim_reduct_signac {
    container 'kaizhang/signac:1.6'
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    errorStrategy 'ignore'

    when: is_included("signac", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env Rscript
    library("rhdf5")
    set.seed(2022)

    file = H5Fopen("data.h5ad", flags = "H5F_ACC_RDONLY")
    x <- H5Gopen(file, "X")
    data <- Matrix::sparseMatrix(i = x\$indices, p = x\$indptr,
        x = c(x\$data), index1=F, repr = "C", dims=rev(H5Aread(H5Aopen(x, "shape")))
    )
    result <- Signac:::RunTFIDF.default(data,method = 1)
    Signac:::RunSVD.default(result, n = 30)
    """
}

////////////////////////////////////////////////////////////////////////////////
// Scale                                                                      //
////////////////////////////////////////////////////////////////////////////////
process dim_reduct_scale {
    container 'kaizhang/scale:1.1.2'
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    containerOptions '--nv'
    label "gpu"
    errorStrategy 'ignore'

    when: is_included("scale", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import numpy as np
    import subprocess
    import tempfile

    def get_free_gpus(threshold_vram_usage=200):
        import subprocess
        import os

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
        return free_gpus

    with tempfile.TemporaryDirectory(dir='./') as temp_dir:
        data = ad.read("data.h5ad")
        k = 15
        input_file = f"{temp_dir}/data.h5ad"
        data.X = data.X.astype(np.float64)
        data.write(input_file)

        gpus = get_free_gpus()
        if len(gpus) > 0:
            subprocess.run([
              "SCALE.py", "-d", input_file, "-k", str(k), "--n_feature", "-1",
              "-i", "10000",
              "--min_cells", "0", "-l", "30", "--gpu", gpus[0],
            ], check = True)
        else:
            raise Exception("No free GPU available")
    """
}

////////////////////////////////////////////////////////////////////////////////
// epiScanpy
////////////////////////////////////////////////////////////////////////////////
process dim_reduct_episcanpy {
    container 'kaizhang/episcanpy:0.4.0'
    tag "$name"
    cpus 4
    time '3h'
    memory '120 GB'
    errorStrategy 'ignore'

    when: is_included("pca", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env python
    import anndata as ad
    import episcanpy as epi
    import numpy as np
    data = ad.read("data.h5ad")
    data.X.data = data.X.data.astype(np.float64)
    epi.pp.normalize_per_cell(data)
    epi.pp.log1p(data)
    epi.pp.pca(data, n_comps=30)
    """
}

process dim_reduct_scbasset {
    container 'kaizhang/scbasset:latest'
    tag "$name"
    errorStrategy 'ignore'
    containerOptions '--nv'
    label "gpu"
    cpus 4
    time '3h'
    memory '120 GB'

    when: is_included("scbasset", params.method_include, params.method_exclude)

    input:
      tuple val(name), path("data.h5ad"), path("fasta.fa")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import pandas as pd
    import numpy as np
    import h5py
    import gc
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
    epochs = 10
    
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
    """
}
