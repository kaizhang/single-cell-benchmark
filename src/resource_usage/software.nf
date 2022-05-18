nextflow.enable.dsl=2

/*******************************************************************************
// SnapATAC2
*******************************************************************************/

process preproc_snapatac_2 {
    //container 'kaizhang/snapatac2:1.99.99.7'
    tag "$data.name"
    input:
      val(data)
    output:
      tuple val("$data.name"), path("data.h5ad")

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    data = snap.pp.import_data(
        fragment_file = "${data.fragments_name_sorted}",
        gff_file = "${data.gene_annotations}",
        chrom_size = snap.genome.mm10,
        file = "data.h5ad",
        min_tsse = 0,
        min_num_fragments = 0,
        chunk_size = 100,
    )
    snap.pp.make_tile_matrix(data, chunk_size = 100)
    """
}

process dim_reduct_snapatac_2 {
    cache false
    //container 'kaizhang/snapatac2:1.99.99.7'
    tag "$name"
    input:
      tuple val(name), path(data)
    output:
      tuple val(name), path(data)

    """
    #!/usr/bin/env python3
    import snapatac2 as snap
    data = snap.read("$data")
    snap.pp.select_features(data)
    snap.tl.spectral(data)
    """
}

process clust_snapatac_2 {
    //container 'kaizhang/snapatac2:1.99.99.7'
    tag "$name"
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

/*******************************************************************************
// SnapATAC R version
*******************************************************************************/

process preproc_snapatac_1 {
    container 'kaizhang/snapatac:1.0'
    tag "$data.name"

    input:
      val(data)
    output:
      tuple val("$data.name"), path("data.rds")

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
        "--input-file=${data.fragments_name_sorted}",
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

process dim_reduct_snapatac_1 {
    container 'kaizhang/snapatac:1.0'
    tag "$name"
    input:
      tuple val(name), path(data)
    output:
      tuple val(name), path("data.rds")

    """
    #!/usr/bin/env Rscript
    library(SnapATAC)
    x.sp <- readRDS("${data}")
    bin.cov <- log10(Matrix::colSums(x.sp@bmat)+1)
    bin.cutoff <- quantile(bin.cov[bin.cov > 0], 0.95)
    idy <- which(bin.cov <= bin.cutoff & bin.cov > 0)
    x.sp <- x.sp[, idy, mat="bmat"]

    x.sp = runDiffusionMaps(
      obj=x.sp,
      input.mat="bmat", 
      num.eigs=50
    )
    saveRDS(x.sp, file="data.rds")
    """
}

process clust_snapatac_1 {
    container 'kaizhang/snapatac:1.0'
    tag "$name"
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

/*******************************************************************************
// ArchR
*******************************************************************************/

process  preproc_archr {
    container 'kaizhang/archr:1.0.1'
    tag "$data.name"

    input:
      val(data)
    output:
      tuple val("$data.name"), path("archr")

    """
    #!/usr/bin/env Rscript
    library("ArchR")
    library("parallel")
    set.seed(1)
    addArchRGenome("mm10")

    ArrowFiles <- createArrowFiles(
        inputFiles = c("${data.fragments_sorted}"),
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

process  dim_reduct_archr {
    container 'kaizhang/archr:1.0.1'
    tag "$name"

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
    proj <- addIterativeLSI(
        ArchRProj = proj,
        useMatrix = "TileMatrix",
        name = "IterativeLSI",
        force=T,
    )
    saveArchRProject(ArchRProj = proj, outputDirectory = "archr", load = F)
    """
}

process clust_archr {
    container 'kaizhang/archr:1.0.1'
    tag "$name"

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