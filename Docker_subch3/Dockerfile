FROM r-base:4.0.3

RUN apt-get update && apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev pkg-config
RUN R -e "install.packages('data.table', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('NOISeq')"
RUN R -e "BiocManager::install('GSEABase')"
RUN R -e "BiocManager::install('GSVA')"
RUN R -e "install.packages('remotes', repos = 'http://cran.us.r-project.org')"
RUN R -e "remotes::install_github('icbi-lab/immunedeconv')"
RUN R -e "install.packages('survival', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('glmnet', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('XML', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('lpSolve', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('RSQLite', repos = 'http://cran.us.r-project.org')"


COPY FICAN-OSCAR-subchallenge3_5.R /FICAN-OSCAR-subchallenge3_5.R
COPY c6.all.v7.2.symbols.gmt /c6.all.v7.2.symbols.gmt
COPY c7.all.v7.2.symbols.gmt /c7.all.v7.2.symbols.gmt
COPY h.all.v7.2.symbols.gmt /h.all.v7.2.symbols.gmt
COPY c2.all.v7.2.symbols.gmt /c2.all.v7.2.symbols.gmt
COPY selfmade.gmt /selfmade.gmt


ENTRYPOINT ["Rscript", "/FICAN-OSCAR-subchallenge3_5.R","/data/GRCh37ERCC_refseq105_genes_tpm.csv","/data/clinical_data.csv","/output/predictions.csv"]
