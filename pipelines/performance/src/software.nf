nextflow.enable.dsl=2

////////////////////////////////////////////////////////////////////////////////
// SnapATAC2                                                                  //
////////////////////////////////////////////////////////////////////////////////

process preproc_snapatac2 {
    container 'kaizhang/snapatac2:2.3.0'
    tag "$name"
    cpus 4
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
    container 'kaizhang/snapatac2:2.3.0'
    stageInMode "copy"
    tag "$name"
    cpus 4
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
    container 'kaizhang/snapatac2:2.3.0'
    stageInMode "copy"
    tag "$name"
    cpus 4
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

    input:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import numpy as np
    import scvi
    import math
    data = ad.read("data.h5ad")
    scvi.model.PEAKVI.setup_anndata(data)
    pvi = scvi.model.PEAKVI(data, n_latent=30)
    pvi.train()
    """
}

////////////////////////////////////////////////////////////////////////////////
// Signac                                                                     //
////////////////////////////////////////////////////////////////////////////////
process dim_reduct_signac {
    container 'kaizhang/signac:1.6'
    tag "$name"
    cpus 4

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

    input:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env python3
    import anndata as ad
    import numpy as np
    import subprocess
    data = ad.read("data.h5ad")
    subprocess.run([
      "SCALE.py", "-d", "data.h5ad", "-k", "10", "--min_peaks", "0",
      "--min_cells", "0", "-l", "30",
    ], check = True)
    """
}

////////////////////////////////////////////////////////////////////////////////
// epiScanpy
////////////////////////////////////////////////////////////////////////////////
process dim_reduct_episcanpy {
    container 'kaizhang/episcanpy:0.4.0'
    tag "$name"
    cpus 4

    input:
      tuple val(name), path("data.h5ad")

    """
    #!/usr/bin/env python
    import anndata as ad
    import episcanpy.api as epi
    import scanpy as sc
    import numpy as np

    data = ad.read("data.h5ad")
    epi.pp.normalize_per_cell(data)
    epi.pp.log1p(data)
    epi.pp.pca(data, n_comps=30)
    """
}