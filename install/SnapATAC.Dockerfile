FROM ubuntu:20.04

RUN apt-get update

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    dirmngr gnupg apt-transport-https ca-certificates \
    software-properties-common build-essential gfortran \
    libcurl4-gnutls-dev libxml2-dev libssl-dev libgsl-dev libgdal-dev r-base
RUN Rscript -e "install.packages(c('doSNOW', 'devtools', 'plot3D', 'BiocManager', 'irlba', 'Rtsne', 'edgeR', 'igraph'))"
RUN Rscript -e "BiocManager::install('edgeR')"
RUN Rscript -e "devtools::install_github('r3fang/SnapATAC')"

ENTRYPOINT ["Rscript", "--vanilla"]