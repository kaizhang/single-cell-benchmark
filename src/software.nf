nextflow.enable.dsl=2

process dim_reduct_snapatac_2 {
    input:
      tuple val(data), val(nDims)
    output:
      tuple val("SnapATAC2"), val(data), path('reduced_dim.tsv')

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("${data.anndata}")
    snap.tl.spectral(adata)
    np.savetxt("reduced_dim.tsv", adata.obsm["X_spectral"], delimiter="\t")
    """
}

process dim_reduct_snapatac_2_nystrom {
    input:
      tuple val(data), val(nDims), val(samplingRate)
    output:
      tuple val("SnapATAC2(nystrom-${samplingRate})"), val(data), path('reduced_dim.tsv')

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("${data.anndata}")
    snap.tl.spectral(adata, sample_size=$samplingRate)
    np.savetxt("reduced_dim.tsv", adata.obsm["X_spectral"], delimiter="\t")
    """
}

process dim_reduct_snapatac_2_kmeans_nystrom {
    input:
      tuple val(data), val(nDims), val(samplingRate)
    output:
      tuple val("SnapATAC2(kmeans-nystrom-${samplingRate})"), val(data), path('reduced_dim.tsv')

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    import numpy as np
    adata = snap.read("${data.anndata}")
    snap.tl.spectral(adata, sample_size=$samplingRate, sampling_method="kmeans")
    np.savetxt("reduced_dim.tsv", adata.obsm["X_spectral"], delimiter="\t")
    """
}

process dim_reduct_snapatac_1 {
    container 'kaizhang/snapatac:1.0'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val("SnapATAC-v1.0"), val(data), path('reduced_dim.tsv')

    """
    #!/usr/bin/env Rscript
    library("Matrix")
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
    write.table(result, file="reduced_dim.tsv", row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_archr_1 {
    container 'kaizhang/archr:1.0.1'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val('ArchR-v1.0.1'), val(data), path('reduced_dim.tsv')

    """
    #!/usr/bin/env Rscript
    library("ArchR")
    data <- as(t(readMM("${data.matrixMarket}")), "dgCMatrix")
    result <- ArchR:::.computeLSI(mat = data,
        LSIMethod = 1,
        scaleTo = 10^4,
        nDimensions = $nDims,
        binarize = TRUE, 
        outlierQuantiles = NULL,
        seed = 1
    )
    write.table(result\$matSVD, file="reduced_dim.tsv", row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_cistopic_3 {
    container 'kaizhang/cistopic:3.0'

    input:
      tuple val(data), val(nDims)
    output:
      tuple val('cisTopic-v3.0'), val(data), path('reduced_dim.tsv')

    """
    #!/usr/bin/env Rscript
    library("cisTopic")
    data <- as(t(readMM("${data.matrixMarket}")), "dgCMatrix")

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