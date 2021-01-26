FROM r-base:4.0.3

RUN apt-get update && apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev pkg-config
RUN R -e "install.packages('data.table', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('NOISeq')"
RUN R -e "BiocManager::install('GSVA')"
RUN R -e "install.packages('remotes', repos = 'http://cran.us.r-project.org')"
RUN R -e "remotes::install_github('icbi-lab/immunedeconv')"
RUN R -e "install.packages('survival', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('glmnet', repos = 'http://cran.us.r-project.org')"

COPY szabo_inflammation_signature.R /szabo_inflammation_signature.R

ENTRYPOINT ["Rscript", "/szabo_inflammation_signature.R","/data/GRCh37ERCC_refseq105_genes_count.csv","/output/predictions.csv"]
