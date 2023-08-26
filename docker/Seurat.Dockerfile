FROM bioconductor/bioconductor_docker
RUN Rscript -e "install.packages( \
    c('rhdf5', 'remotes', 'httr', 'png', 'leiden', 'GenomeInfoDb', 'GenomicRanges', 'IRanges', \
    'Rsamtools', 'S4Vectors', 'BiocGenerics', 'Signac', 'Seurat'), \
    repos = BiocManager::repositories())"
RUN Rscript -e "remotes::install_github('mojaveazure/seurat-disk', repos = BiocManager::repositories())"