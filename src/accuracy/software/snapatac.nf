nextflow.enable.dsl=2

process dim_reduct_snapatac {
    container 'kaizhang/snapatac:1.0'

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val("SnapATAC"), path("reduced_dim.tsv")

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

    if (nrow(data) <= 20000) {
        x.sp <- SnapATAC::newSnap()
        x.sp@bmat <- data
        x.sp <- SnapATAC::makeBinary(x.sp, mat="bmat")
        x.sp <- SnapATAC::runDiffusionMaps(
                obj=x.sp,
                input.mat="bmat", 
                num.eigs=30
            )
        result <- SnapATAC:::weightDimReduct(x.sp@smat, 1:30, weight.by.sd=T)
    } else {
        sample_size <- 20000
        n_dim <- min(sample_size - 2, 30)
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

        result <- SnapATAC:::weightDimReduct(x.sp@smat, 1:30, weight.by.sd=T)
    }

    write.table(result, file="reduced_dim.tsv", row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_snapatac_nystrom {
    container 'kaizhang/snapatac:1.0'

    input:
      val(data)
    output:
      tuple val("SnapATAC"), val(data), path("${data.name}_snapatac_nystrom_${data.samplingFraction}_${data.randomSeed}_reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("Matrix")
    library("rhdf5")
    set.seed($data.randomSeed)

    file = H5Fopen("${data.anndata}", flags = "H5F_ACC_RDONLY")
    x <- H5Gopen(file, "X")
    data <- sparseMatrix(j = x\$indices, p = x\$indptr,
        x = c(x\$data), index1=F, repr = "C", dims=H5Aread(H5Aopen(x, "shape"))
    )
    sample_size <- trunc(nrow(data) * $data.samplingFraction)
    n_dim <- min(sample_size - 2, 50)
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