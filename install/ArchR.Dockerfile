FROM bioconductor/bioconductor_docker
RUN Rscript -e "install.packages( \
    c('rhdf5', 'devtools', 'httr', 'png', 'leiden', 'GenomeInfoDb', 'GenomicRanges', 'IRanges', \
    'Rsamtools', 'S4Vectors', 'BiocGenerics'), \
    repos = BiocManager::repositories())"
RUN Rscript -e "devtools::install_github('GreenleafLab/ArchR', ref='v1.0.1', repos = BiocManager::repositories())"