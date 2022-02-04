FROM ubuntu:20.04

RUN apt-get update

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    dirmngr gnupg apt-transport-https ca-certificates \
    software-properties-common build-essential \
    libcurl4-gnutls-dev libxml2-dev libssl-dev libgsl-dev
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

RUN apt-get install -y r-base
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "devtools::install_github('GreenleafLab/ArchR', ref='v1.0.1', repos = BiocManager::repositories())"

ENTRYPOINT ["Rscript", "--vanilla"]