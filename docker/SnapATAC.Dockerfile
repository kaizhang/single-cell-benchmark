FROM bioconductor/bioconductor_docker
RUN Rscript -e "install.packages( \
    c('rhdf5', 'devtools', 'httr', 'png', 'leiden', 'GenomeInfoDb', 'GenomicRanges', 'IRanges', \
    'Rsamtools', 'S4Vectors', 'BiocGenerics'), \
    repos = BiocManager::repositories())"
RUN Rscript -e "install.packages(c('doSNOW', 'plot3D', 'Rtsne'), repos = BiocManager::repositories())"
RUN Rscript -e "devtools::install_github('r3fang/SnapATAC', repos = BiocManager::repositories())"
RUN pip install snaptools