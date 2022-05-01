FROM bioconductor/bioconductor_docker
RUN Rscript -e "install.packages( \
    c('rhdf5', 'devtools', 'httr', 'png', 'leiden', 'GenomeInfoDb', 'GenomicRanges', 'IRanges', \
    'Rsamtools', 'S4Vectors', 'BiocGenerics'), \
    repos = BiocManager::repositories())"
RUN Rscript -e "install.packages('Signac', repos = BiocManager::repositories())"