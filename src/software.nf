nextflow.enable.dsl=2

/*******************************************************************************
// SnapATAC2
*******************************************************************************/

process dim_reduct_snapatac_2 {
    publishDir 'result/reduced_dim/'
    input:
      tuple val(data), val(nDims)
    output:
      tuple val("SnapATAC2"), val(data), path("${data.name}_snapatac2_reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("${data.anndata}")
    snap.tl.spectral(adata)
    np.savetxt("${data.name}_snapatac2_reduced_dim.tsv", adata.obsm["X_spectral"], delimiter="\t")
    """
}

process dim_reduct_snapatac_2_cosine {
    input:
      tuple val(data), val(nDims)
    output:
      tuple val("SnapATAC2_cosine"), val(data), path("reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("${data.anndata}")
    snap.tl.spectral(adata, distance_metric="cosine")
    np.savetxt("reduced_dim.tsv", adata.obsm["X_spectral"], delimiter="\t")
    """
}


process dim_reduct_snapatac_2_svd {
    publishDir 'result/reduced_dim/'
    input:
      tuple val(data), val(nDims)
    output:
      tuple val("SnapATAC2_SVD"), val(data), path("${data.name}_snapatac2_svd_reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np

    adata = snap.read("${data.anndata}")
    snap.tl.laplacian(adata)
    output = "${data.name}_snapatac2_svd_reduced_dim.tsv"
    np.savetxt(output, adata.obsm["X_spectral"], delimiter="\t")
    """
}

process dim_reduct_snapatac_2_nystrom {
    publishDir 'result/reduced_dim/'
    input:
      tuple val(data), val(nDims)
    output:
      tuple val("SnapATAC2"), val(data), path("${data.name}_snapatac2_nystrom_${data.samplingFraction}_${data.randomSeed}_reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np

    adata = snap.read("${data.anndata}")

    snap.tl.spectral(adata, chunk_size=500, sample_size=${data.samplingFraction}, random_state=${data.randomSeed})
    output = "${data.name}_snapatac2_nystrom_${data.samplingFraction}_${data.randomSeed}_reduced_dim.tsv"
    np.savetxt(output, adata.obsm["X_spectral"], delimiter="\t")
    """
}

process dim_reduct_snapatac_2_v2_nystrom_full {
    publishDir 'result/reduced_dim/'
    input:
      tuple val(data), val(nDims)
    output:
      tuple val("SnapATAC2_v2_full"), val(data), path("${data.name}_snapatac2_v2_nystrom_full_${data.samplingFraction}_${data.randomSeed}_reduced_dim.tsv")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("${data.anndata}")
    snap.tl.spectral2(adata, sample_size=${data.samplingFraction}, random_state=${data.randomSeed}, full_nystrom=True)
    output = "${data.name}_snapatac2_v2_nystrom_full_${data.samplingFraction}_${data.randomSeed}_reduced_dim.tsv"
    np.savetxt(output, adata.obsm["X_spectral"], delimiter="\t")
    """
}


/*******************************************************************************
// SnapATAC R version
*******************************************************************************/

process dim_reduct_snapatac_1 {
    container 'kaizhang/snapatac:1.0'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val("SnapATAC-v1.0"), val(data), path("reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("Matrix")
    set.seed(2022)
    data <- as(readMM("${data.matrixMarket}"), "dgCMatrix")
    x.sp <- SnapATAC::newSnap()
    x.sp@bmat <- data
    x.sp <- SnapATAC::makeBinary(x.sp, mat="bmat")
    x.sp <- SnapATAC::runDiffusionMaps(
            obj=x.sp,
            input.mat="bmat", 
            num.eigs=$nDims
        )
    result <- SnapATAC:::weightDimReduct(x.sp@smat, 1:$nDims, weight.by.sd=F)
    output <- "reduced_dim.tsv"
    write.table(result, file=output, row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_snapatac_1_nystrom {
    container 'kaizhang/snapatac:1.0'
    publishDir 'result/reduced_dim/'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val("SnapATAC-v1.0"), val(data), path("${data.name}_snapatac_nystrom_${data.samplingFraction}_${data.randomSeed}_reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("Matrix")
    set.seed($data.randomSeed)
    data <- as(readMM("${data.matrixMarket}"), "dgCMatrix")
    sample_size <- trunc(nrow(data) * $data.samplingFraction)
    n_dim <- min(sample_size - 2, $nDims)
    reference <- data[sample(nrow(data),size=sample_size, replace=F),]

    x.ref <- SnapATAC::newSnap()
    x.ref@bmat <- reference
    x.ref <- SnapATAC::makeBinary(x.ref, mat="bmat")
    x.ref <- SnapATAC::runDiffusionMaps(
        obj=x.ref,
        input.mat="bmat", 
        num.eigs=n_dim
    )

    x.sp <- SnapATAC::newSnap()
    x.sp@bmat <- data
    x.sp <- SnapATAC::makeBinary(x.sp, mat="bmat")
    x.sp <- SnapATAC::runDiffusionMapsExtension(
	      obj1 = x.ref,
        obj2 = x.sp,
        input.mat="bmat"
    )

    result <- SnapATAC:::weightDimReduct(x.sp@smat, 1:n_dim, weight.by.sd=F)
    output <- "${data.name}_snapatac_nystrom_${data.samplingFraction}_${data.randomSeed}_reduced_dim.tsv"
    write.table(result, file=output, row.names=F, col.names=F, sep="\t")
    """
}


/*******************************************************************************
// ArchR
*******************************************************************************/

process dim_reduct_archr_1_tf_logidf {
    container 'kaizhang/archr:1.0.1'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val('ArchR (tf-logidf)'), val(data), path("reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("ArchR")
    set.seed(2022)
    data <- as(t(readMM("${data.matrixMarket}")), "dgCMatrix")
    result <- ArchR:::.computeLSI(mat = data,
        LSIMethod = 1,
        scaleTo = 10^4,
        nDimensions = $nDims,
        binarize = TRUE, 
        outlierQuantiles = NULL,
        seed = 1
    )
    output <- "reduced_dim.tsv"
    write.table(result\$matSVD, file=output, row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_archr_1_log_tf_idf {
    container 'kaizhang/archr:1.0.1'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val('ArchR (log(tf-idf))'), val(data), path("reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("ArchR")
    set.seed(2022)
    data <- as(t(readMM("${data.matrixMarket}")), "dgCMatrix")
    result <- ArchR:::.computeLSI(mat = data,
        LSIMethod = 2,
        scaleTo = 10^4,
        nDimensions = $nDims,
        binarize = TRUE, 
        outlierQuantiles = NULL,
        seed = 1
    )
    output <- "reduced_dim.tsv"
    write.table(result\$matSVD, file=output, row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_archr_1_logtf_logidf {
    container 'kaizhang/archr:1.0.1'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val('ArchR (logtf-logidf)'), val(data), path("reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("ArchR")
    set.seed(2022)
    data <- as(t(readMM("${data.matrixMarket}")), "dgCMatrix")
    result <- ArchR:::.computeLSI(mat = data,
        LSIMethod = 3,
        scaleTo = 10^4,
        nDimensions = $nDims,
        binarize = TRUE, 
        outlierQuantiles = NULL,
        seed = 1
    )
    output <- "reduced_dim.tsv"
    write.table(result\$matSVD, file=output, row.names=F, col.names=F, sep="\t")
    """
}


process dim_reduct_archr_1_subsample {
    container 'kaizhang/archr:1.0.1'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val('ArchR-v1.0.1'), val(data), path("reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("ArchR")
    set.seed($data.randomSeed)
    data <- as(t(readMM("${data.matrixMarket}")), "dgCMatrix")
    sample_size <- trunc(ncol(data) * $data.samplingFraction)
    n_dim <- min(sample_size - 2, $nDims)
    reference_data <- data[, sample(ncol(data),size=sample_size, replace=F)]

    ref <- ArchR:::.computeLSI(mat = reference_data,
        LSIMethod = 2,
        scaleTo = 10^4,
        nDimensions = n_dim,
        binarize = TRUE, 
        outlierQuantiles = NULL,
        seed = 1
    )
    result <- ArchR:::.projectLSI(
        mat = data,
        LSI = ref,
        returnModel = FALSE, 
        verbose = FALSE, 
        tstart = NULL,
        logFile = NULL
    )
    output <- "reduced_dim.tsv"
    write.table(result, file=output, row.names=F, col.names=F, sep="\t")
    """
}


/*******************************************************************************
// cisTopic
*******************************************************************************/

process dim_reduct_cistopic_3 {
    container 'kaizhang/cistopic:3.0'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val('cisTopic-v3.0'), val(data), path('reduced_dim.tsv')

    """
    #!/usr/bin/env Rscript
    library("cisTopic")
    library("anndata")
    library("stringr")
    set.seed(0)
    mk_chr_name <- function(x) {
        r <- str_split(x, "_|:|-")[[1]]
        return(paste(r[1], paste(r[2], r[3], sep="-"), sep=":"))
    }
    adata <- read_h5ad("${data.anndata}")\$T
    mat <- adata\$X
    rownames(mat) <- sapply(adata\$obs_names, mk_chr_name)
    colnames(mat) <- adata\$var_names
    cisTopicObject <- createcisTopicObject(
        mat,
        project.name='cisTopic',
        min.cells = 0,
        min.regions = 0,
        is.acc = 0,
        keepCountsMatrix=F,
    )

    cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=50,
        seed=2022, nCores=10, addModels=FALSE
    )

      cisTopicObject <- selectModel(cisTopicObject)
      fm_cisTopic <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
      return(fm_cisTopic)
}



    data <- createcisTopicObject(
        data,
        project.name = "cisTopicProject",
        min.cells = 0,
        min.regions = 0,
        is.acc = 0,
        keepCountsMatrix=F,
    )

    write.table(result\$matSVD, file="reduced_dim.tsv", row.names=F, col.names=F, sep="\t")
    """
}