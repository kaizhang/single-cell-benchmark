nextflow.enable.dsl=2

process dim_reduct_signac_1 {
    container 'kaizhang/signac:1.6'
    tag "$name"
    cpus 4
    errorStrategy 'ignore'

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val('Signac (log(TF-IDF))'), path("reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("rhdf5")
    set.seed(2022)

    find_elbow <- function(x) {
        saturation <- 0.01
        accum_gap <- 0
        for (i in 2:length(x)) {
            gap <- x[i-1] - x[i]
            accum_gap <- accum_gap + gap
            if (gap < saturation * accum_gap) {
                return(i)
            }
        }
        return(i)
    }

    file = H5Fopen("data.h5ad", flags = "H5F_ACC_RDONLY")
    x <- H5Gopen(file, "X")
    data <- Matrix::sparseMatrix(i = x\$indices, p = x\$indptr,
        x = c(x\$data), index1=F, repr = "C", dims=rev(H5Aread(H5Aopen(x, "shape")))
    )
    result <- Signac:::RunTFIDF.default(data,method = 1)
    result <- Signac:::RunSVD.default(result, n = 30)

    output <- "reduced_dim.tsv"
    i <- find_elbow(result@stdev)
    write.table(result@cell.embeddings[, 1:i], file=output, row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_signac_2 {
    container 'kaizhang/signac:1.6'
    tag "$name"
    cpus 4
    errorStrategy 'ignore'

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val('Signac (TF-logIDF)'), path("reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("rhdf5")
    set.seed(2022)

    find_elbow <- function(x) {
        saturation <- 0.01
        accum_gap <- 0
        for (i in 2:length(x)) {
            gap <- x[i-1] - x[i]
            accum_gap <- accum_gap + gap
            if (gap < saturation * accum_gap) {
                return(i)
            }
        }
        return(i)
    }

    file = H5Fopen("data.h5ad", flags = "H5F_ACC_RDONLY")
    x <- H5Gopen(file, "X")
    data <- Matrix::sparseMatrix(i = x\$indices, p = x\$indptr,
        x = c(x\$data), index1=F, repr = "C", dims=rev(H5Aread(H5Aopen(x, "shape")))
    )
    result <- Signac:::RunTFIDF.default(data,method = 2)
    result <- Signac:::RunSVD.default(result, n = 30)

    output <- "reduced_dim.tsv"
    i <- find_elbow(result@stdev)
    write.table(result@cell.embeddings[, 1:i], file=output, row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_signac_3 {
    container 'kaizhang/signac:1.6'
    tag "$name"
    cpus 4
    errorStrategy 'ignore'

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val('Signac (logTF-logIDF)'), path("reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("rhdf5")
    set.seed(2022)

    find_elbow <- function(x) {
        saturation <- 0.01
        accum_gap <- 0
        for (i in 2:length(x)) {
            gap <- x[i-1] - x[i]
            accum_gap <- accum_gap + gap
            if (gap < saturation * accum_gap) {
                return(i)
            }
        }
        return(i)
    }

    file = H5Fopen("data.h5ad", flags = "H5F_ACC_RDONLY")
    x <- H5Gopen(file, "X")
    data <- Matrix::sparseMatrix(i = x\$indices, p = x\$indptr,
        x = c(x\$data), index1=F, repr = "C", dims=rev(H5Aread(H5Aopen(x, "shape")))
    )
    result <- Signac:::RunTFIDF.default(data,method = 3)
    result <- Signac:::RunSVD.default(result, n = 30)

    output <- "reduced_dim.tsv"
    i <- find_elbow(result@stdev)
    write.table(result@cell.embeddings[, 1:i], file=output, row.names=F, col.names=F, sep="\t")
    """
}

process dim_reduct_signac_4 {
    container 'kaizhang/signac:1.6'
    tag "$name"
    cpus 4
    errorStrategy 'ignore'

    input:
      tuple val(name), path("data.h5ad")
    output:
      tuple val(name), val('Signac (IDF)'), path("reduced_dim.tsv")

    """
    #!/usr/bin/env Rscript
    library("rhdf5")
    set.seed(2022)

    find_elbow <- function(x) {
        saturation <- 0.01
        accum_gap <- 0
        for (i in 2:length(x)) {
            gap <- x[i-1] - x[i]
            accum_gap <- accum_gap + gap
            if (gap < saturation * accum_gap) {
                return(i)
            }
        }
        return(i)
    }

    file = H5Fopen("data.h5ad", flags = "H5F_ACC_RDONLY")
    x <- H5Gopen(file, "X")
    data <- Matrix::sparseMatrix(i = x\$indices, p = x\$indptr,
        x = c(x\$data), index1=F, repr = "C", dims=rev(H5Aread(H5Aopen(x, "shape")))
    )
    result <- Signac:::RunTFIDF.default(data,method = 4)
    result <- Signac:::RunSVD.default(result, n = 30)
    output <- "reduced_dim.tsv"
    i <- find_elbow(result@stdev)
    write.table(result@cell.embeddings[, 1:i], file=output, row.names=F, col.names=F, sep="\t")
    """
}

process end_to_end_signac {
    container 'kaizhang/signac:1.6'
    errorStrategy 'ignore'

    input:
      val(data)
    output:
      tuple val('Signac'), val(data), path("clusters.tsv")

    """
    #!/usr/bin/env Rscript
    library(Signac)
    library(Seurat)
    library(GenomeInfoDb)
    library(EnsDb.Hsapiens.v75)
    library(ggplot2)
    library(patchwork)
    set.seed(1234)

    frags <- CreateFragmentObject(path = "${data.fragments_name_sorted}")

    GenomeBinMatrix(
        fragments = fragments,
        genome = genome,
        binsize = 1000
    )

    clusters <- c()
    for (i in seq(from = 0.1, to = 1.5, by = 0.1)) {
        proj <- addClusters(
            input = proj,
            reducedDims = "IterativeLSI",
            maxClusters = 100,
            filterBias = F,
            force=T,
            resolution = i,
        )
        clusters <- cbind(clusters, proj\$Clusters)
    }
    rownames(clusters) <- sapply(proj\$cellNames, function(x) strsplit(x, '#')[[1]][2])
    write.table(
        clusters, file = "clusters.tsv", quote = F, sep = "\t", row.names = T,
        col.names = F, 
    )
    """
}