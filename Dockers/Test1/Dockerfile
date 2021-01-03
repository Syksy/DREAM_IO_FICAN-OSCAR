FROM r-base:4.0.3

RUN R -e "install.packages('hamlet', repos = 'http://cran.us.r-project.org')"

COPY test.R /test.R

CMD ["Rscript", "/test.R", "/data/GRCh37ERCC_refseq105_genes_count.csv", "/output/predictions.csv"]
