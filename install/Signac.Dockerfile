FROM rstudio/r-base:4.1.2-opensuse153
RUN zypper refresh && zypper -n update
RUN zypper -n install libpng-devel openssl-devel libxml2-devel lapack-devel arpack-ng-devel python3-pip
RUN pip3 install anndata
RUN Rscript -e "install.packages('BiocManager', repos = 'https://cran.us.r-project.org')"
RUN Rscript -e "install.packages( \
    c('anndata', 'devtools', 'httr', 'png', 'leiden', 'GenomeInfoDb', 'GenomicRanges', 'IRanges', \
    'Rsamtools', 'S4Vectors', 'BiocGenerics'), \
    repos = BiocManager::repositories())"
RUN Rscript -e "install.packages('Signac', repos = BiocManager::repositories())"