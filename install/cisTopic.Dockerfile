FROM rstudio/r-base:3.5.3-bionic
RUN apt-get update -y
RUN apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev

RUN Rscript -e "install.packages('devtools', repos = 'https://cran.us.r-project.org')"
RUN Rscript -e "devtools::install_github('aertslab/AUCell')"
RUN Rscript -e "devtools::install_github('aertslab/RcisTarget')"

RUN Rscript -e "install.packages(c('RcppParallel', 'mlapi', 'sparsepp', 'foreach'), repos = 'https://cran.us.r-project.org')"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/irlba/irlba_2.3.3.tar.gz', repo=NULL, type='source')"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/text2vec/text2vec_0.5.1.tar.gz', repo=NULL, type='source')"
RUN Rscript -e "devtools::install_github('aertslab/cisTopic')"

RUN apt-get install -y libpng-dev
RUN Rscript -e "install.packages('anndata', repos = 'https://cran.us.r-project.org')"
RUN Rscript -e "reticulate::install_miniconda()"
RUN Rscript -e "anndata::install_anndata(method='conda')"