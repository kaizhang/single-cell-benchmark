FROM rstudio/r-base:4.1.2-opensuse153
RUN zypper refresh && zypper -n update
RUN zypper -n install libpng-devel openssl-devel libxml2-devel lapack-devel arpack-ng-devel python3-pip
RUN pip3 install anndata
RUN Rscript -e "install.packages('BiocManager', repos = 'https://cran.us.r-project.org')"
RUN Rscript -e "install.packages( \
    c('anndata', 'devtools', 'httr', 'png', 'leiden', 'GenomeInfoDb', 'GenomicRanges', 'IRanges', \
    'Rsamtools', 'S4Vectors', 'BiocGenerics'), \
    repos = BiocManager::repositories())"

RUN zypper -n install --oldpackage gdal-devel proj-devel proj sqlite3-devel libsqlite3-0-3.28.0 geos-devel
RUN Rscript -e "install.packages(c('doSNOW', 'plot3D', 'Rtsne'), repos = BiocManager::repositories())"
RUN Rscript -e "devtools::install_github('r3fang/SnapATAC', repos = BiocManager::repositories())"