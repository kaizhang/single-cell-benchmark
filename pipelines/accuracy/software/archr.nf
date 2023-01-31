nextflow.enable.dsl=2

process dim_reduct_archr_1 {
    container 'kaizhang/archr:1.0.1'
    cpus 4

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val('ArchR (TF-logIDF)'), path("reduced_dim.tsv")

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
    result <- ArchR:::.computeLSI(mat = data,
        LSIMethod = 1,
        scaleTo = 10^4,
        nDimensions = 30,
        binarize = TRUE, 
        outlierQuantiles = NULL,
        seed = 1
    )
    output <- "reduced_dim.tsv"
    write.table(result\$matSVD, file=output, row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_archr_2 {
    container 'kaizhang/archr:1.0.1'
    cpus 4

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val('ArchR (log(TF-IDF))'), path("reduced_dim.tsv")

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
    result <- ArchR:::.computeLSI(mat = data,
        LSIMethod = 2,
        scaleTo = 10^4,
        nDimensions = 30,
        binarize = TRUE, 
        outlierQuantiles = NULL,
        seed = 1
    )
    output <- "reduced_dim.tsv"
    write.table(result\$matSVD, file=output, row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_archr_3 {
    container 'kaizhang/archr:1.0.1'
    cpus 4

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val('ArchR (logTF-logIDF)'), path("reduced_dim.tsv")

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
    result <- ArchR:::.computeLSI(mat = data,
        LSIMethod = 3,
        scaleTo = 10^4,
        nDimensions = 30,
        binarize = TRUE, 
        outlierQuantiles = NULL,
        seed = 1
    )
    output <- "reduced_dim.tsv"
    write.table(result\$matSVD, file=output, row.names=F, col.names=F, sep="\t")
    """
}