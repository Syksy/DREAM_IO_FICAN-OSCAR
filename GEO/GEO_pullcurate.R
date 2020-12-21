###
#
# Pull & curate suitable datasets from GEO
#
###

# A lot of the may be utilizable best from: https://github.com/JasonACarter/IMPRES_Correspondence

# Studies suggested by the DREAM organizers available in GEO:

# GSE115821 - Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115821
# Auslander et al.
# Cancer(s): Melanoma/Neuroblastoma (?)
# Treatment(s): N=31 treated with anti-PD-1, N=10 treated with anti-CTLA-4
# N(s): 37, for some patients multiple samples at different time points or different biopsies at same time point
# Response types: Binary (?)
# Year: 2018
# Notes: 
# - Identified 15 pairwise transcriptomics relations between immune checkpoint genes; AUC=0.83; validated on 6 published existing cohorts at the time (2018) (n total 297)
# - Criticized at https://www.nature.com/articles/s41591-019-0671-4

# GSE121810 - Neoadjuvant anti-PD-1 immunotherapy promotes a survival benefit with intratumoral and systemic immune responses in recurrent glioblastoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121810
# Cloughesy et al.
# Cancer(s): Glioblastoma (recurrent, surgically resectable)
# Treatment(s): Neoadjuvant Pembrolizumab vs. neoadjuvant post-surgical PD-1 blockade alone
# N(s): 29 (?)
# Response types: Survival
# Year: 2019
# Notes: 
# - Recurrent disease
# - "This dataset contains the transcriptomes of recurrent glioblastoma with either neoadjuvant (1 dose) or adjuvant pembrolizumab treatment"
# - "This trial was registered with ClinicalTrials.gov under the identifier NCT02852655 (https://clinicaltrials.gov/ct2/show/NCT02852655)."

# GSE78220 - mRNA expressions in pre-treatment melanomas undergoing anti-PD-1 checkpoint inhibition therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220
# Hugo et al.
# Cancer(s): Melanoma
# Treatment(s): All with anti-PD-1, some responded others didn't
# N(s): 28 (?)
# Response types: Response vs. no response (binary)
# Year: 2016
# Notes:
# - "Mutations in cell adhesion genes and the DNA repair gene BRCA2 were enriched in responding tumors, and a high mutational load associated with improved survival."
# - "Innately resistant tumors displayed frequent transcriptomic up-expression of genes that enriched for mesenchymal transition, cell adhesion, ECM organization, wound-healing and angiogenesis."
# - "he transcriptomes of innate resistance also enriched for signatures indicating up-regulation of these processes. Notably, MAPK-targeted therapy (MAPKi) induced similar signatures in melanoma, suggesting that a form of MAPKi resistance mediates cross-resistance to anti-PD-1 therapy. Co-enrichment of IPRIM (Innate anti-PD-1 Resistance Induced by MAPKi) signatures defined a transcriptomic subset across advanced cancers, suggesting that attenuating processes underlying these signatures may augment anti-PD1 responses."

# GSE52562 - Gene expression profiling of tumor biopsies before and after pidilizumab therapy in patients with relapsed follicular lymphoma grade 1 or grade 2.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52562
# Westin et al.
# Cancer(s): Follicular lymphoma (grade 1 or 2)
# Treatment(s): Pidilizumab
# N(s): N=8 pairs for pre vs. post-treatment, N=10 additional pre-treatment samples
# Response types: Progression free survival
# Year: 2014
# Notes: 
# - Comparison before and after pidilizumab

# GSE67501 - Expression data in human renal cell carcinoma samples from patients who did or did not respond to anti-PD-1 (Nivolumab) immunotherapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67501
# Ascierto et al. 2016
# Cancer(s): Renal cell carcinoma
# Treatment(s): Anti-PD-1 (Nivolumab)
# N(s): N=4 responders, N=7 non-responders
# Response types: Binary (Response vs. no response)
# Year: 2016
# Notes:
# - Focus on positive tumor expression of PD-L1
# - Microarray expression data

# GSE79691 - Transcriptional mechanisms of resistance to anti-PD-1 therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79691
# Ascierto et al. 2017
# Cancer(s): Melanoma, skin metastases
# Treatment(s): Nivolumab
# N(s): response (N=6), no response (N=4)
# Response types: Binary
# Year: 2017
# Notes: 
# - All samples from a single patient
# - Patient included in Ascierto et al. 2016?

## Other identified relevant studies

# GSE91061 - Molecular portraits of tumor mutational and micro-environmental sculpting by immune checkpoint blockade therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
# Riaz et al. 
# Cancer(s): Advanced melanoma
# Treatment(s): Ipilimumab, Nivolumab
# N(s): N=58 on-treatment and N=51 pre-treatment, 64 patients
# Response types:
# Year:
# Notes
# - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5685550/
# - Ipilimumab-naive
# - Abstract: "The mechanisms by which immune checkpoint blockade modulates tumor evolution during therapy are unclear. We assessed genomic changes in tumors from 68 patients with advanced melanoma, who progressed on ipilimumab or were ipilimumab-naive, before and after nivolumab initiation (CA209-038 study). Tumors were analyzed by whole-exome, transcriptome, and/or T cell receptor (TCR) sequencing. In responding patients, mutation and neoantigen load were reduced from baseline, and analysis of intratumoral heterogeneity during therapy demonstrated differential clonal evolution within tumors and putative selection against neoantigenic mutations on-therapy. Transcriptome analyses before and during nivolumab therapy revealed increases in distinct immune cell subsets, activation of specific transcriptional networks, and upregulation of immune checkpoint genes that were more pronounced in patients with response. Temporal changes in intratumoral TCR repertoire revealed expansion of T cell clones in the setting of neoantigen loss. Comprehensive genomic profiling data in this study provide insight into nivolumab's mechanism of action."

# GSE93157 - Programmed death 1 receptor blockade and immune-related gene expression profiling in non-small cell lung carcinoma, head and neck squamous cell carcinoma and melanoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93157
# Prat et al.
# Cancer(s): Non-small cell lung carcinoma, head and neck squamous cell carcinoma, melanoma
# Treatment(s): Anti-PD1 (pembrolizumab or nivolumab)
# N(s): 
# Response types:
# Year:
# Notes
# - Has non-small cell lung carcinoma
# - https://cancerres.aacrjournals.org/content/77/13/3540.long
# - Focus heavily on immune cell decomposition

## Other

# Template
# Name
# URL
# Cite
# Cancer(s):
# Treatment(s):
# N(s):
# Response types:
# Year:
# Notes
# - 1
# - 2
# - 3

###
#
# Download and process
#
###

setwd("D:\\Gits\\DREAM_2020_IO\\download")
# GEOquery
library(GEOquery)
# Read excel
library(readxl)
# RNA-seq normalization
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")

# Generate a human genome mapping helper data frame
library("biomaRt")
# Generate gene names
#> grep("entrez", listAttributes(biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))[,1], value=TRUE)
#[1] "entrezgene_trans_name"  "entrezgene_description" "entrezgene_accession"   "entrezgene_id"
# Fetch gene names for various annotations
genes <- biomaRt::getBM(
	attributes = 
		c(
			# ENSEMBL
			'ensembl_gene_id', 'ensembl_transcript_id',
			# entrez
			'entrezgene_id',
			# Hugo
			'hgnc_symbol',
			# RefSeq
			'refseq_mrna',
			# Chromosomal information
			'chromosome_name','start_position','end_position',
			# Description
			'description'
		),
	mart = biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
)
# Sort using chromosomes and then bp locations
genes <- genes[order(genes$chromosome_name, genes$start_position, genes$end_position),]
# Omit row names (wrong order indices)
rownames(genes) <- NULL



# Follow guidelines e.g. in http://genomicsclass.github.io/book/pages/GEOquery.html
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Regarding gene length normalization in DESeq2: https://www.biostars.org/p/140090/

# Extract GPLs https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/

# GSE115821 - Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115821
gse_auslander <- GEOquery::getGEO("GSE115821", GSEMatrix = TRUE)
sup_auslander <- GEOquery::getGEOSuppFiles("GSE115821")
GEOquery::gunzip(rownames(sup_auslander)[1])
sup_auslander <- read.csv("GSE115821/GSE115821_MGH_counts.csv")

# GSE121810 - Neoadjuvant anti-PD-1 immunotherapy promotes a survival benefit with intratumoral and systemic immune responses in recurrent glioblastoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121810
gse_cloughesy <- GEOquery::getGEO("GSE121810", GSEMatrix = TRUE)
sup_cloughesy <- GEOquery::getGEOSuppFiles("GSE121810")
sup_cloughesy <- readxl::read_excel(rownames(sup_cloughesy)[1])

# GSE78220 - mRNA expressions in pre-treatment melanomas undergoing anti-PD-1 checkpoint inhibition therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220
gse_hugo <-  GEOquery::getGEO("GSE78220", GSEMatrix = TRUE)
sup_hugo <- GEOquery::getGEOSuppFiles("GSE78220")
sup_hugo <- readxl::read_excel(rownames(sup_hugo)[1])
# sup_hugo -> tibble table, cast to matrix preferably

# GSE52562 - Gene expression profiling of tumor biopsies before and after pidilizumab therapy in patients with relapsed follicular lymphoma grade 1 or grade 2.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52562
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558
gse_westin <-  GEOquery::getGEO("GSE52562", GSEMatrix = TRUE, getGPL = TRUE)
# Mapping/collapsing:
# gse_westin[[1]]@featureData@data[1:2,]
# column "ID" <-> "ILMN_Gene"

# GSE79691 - Transcriptional mechanisms of resistance to anti-PD-1 therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67501
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL18281
gse_ascierto <- GEOquery::getGEO("GSE67501", GSEMatrix = TRUE, getGPL = TRUE)
# Mapping/collapsing:
# gse_ascierto[[1]]@featureData@data[1:2,]
# column "ID" <-> "ILMN_Gene"

# GSE79691 - Transcriptional mechanisms of resistance to anti-PD-1 therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79691
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14951
gse_ascierto2 <- GEOquery::getGEO("GSE79691", GSEMatrix = TRUE, getGPL = TRUE)
# Mapping/collapsing:
# gse_ascierto2[[1]]@featureData@data[1:2,]
# column "ID" <-> "ILMN_Gene"

# GSE91061 - Molecular portraits of tumor mutational and micro-environmental sculpting by immune checkpoint blockade therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
# Note: 'Gene IDs in the processed data files are NCBI Entrez Gene IDs.'
gse_riaz <- GEOquery::getGEO("GSE91061", GSEMatrix = TRUE, getGPL = TRUE)
sup_riaz <- GEOquery::getGEOSuppFiles("GSE91061")
lapply(rownames(sup_riaz), FUN=GEOquery::gunzip)
sup_riaz1 <- read.csv(gsub(".gz", "", rownames(sup_riaz)[1]))
sup_riaz2 <- read.csv(gsub(".gz", "", rownames(sup_riaz)[2]))
sup_riaz3 <- read.csv(gsub(".gz", "", rownames(sup_riaz)[3]))
sup_riaz4 <- read.csv(gsub(".gz", "", rownames(sup_riaz)[4]))
# sup_riaz2 are raw counts, using those with DESeq2
gpl_riaz <- getGEO('GPL9052', destdir=".")

# GSE93157 - Programmed death 1 receptor blockade and immune-related gene expression profiling in non-small cell lung carcinoma, head and neck squamous cell carcinoma and melanoma
gse_prat <- GEOquery::getGEO("GSE93157", GSEMatrix = TRUE)

# Examine main characteristics of the reported phenodata
head(pData(gse_auslander[[1]]))[1:2,]
head(pData(gse_cloughesy[[1]]))[1:2,]
head(pData(gse_hugo[[1]]))[1:2,]
head(pData(gse_westin[[1]]))[1:2,]
head(pData(gse_ascierto[[1]]))[1:2,]
head(pData(gse_ascierto2[[1]]))[1:2,]
head(pData(gse_riaz[[1]]))[1:2,]
head(pData(gse_prat[[1]]))[1:2,]
# Response abbreviations:
#
# From e.g. Hugo et al.
# PD: Progressive disease
# CR: Complete response
# PR: Partial response


# Create gene expression matrices

##
## Auslander et al.
##
exprs(gse_auslander[[1]])
# Botched
# 0 features, 23 + 14 samples? Raw counts in supplementary files
colnames(sup_auslander1 <- gsub("X|.bam", "", colnames(sup_auslander))
# First batch
# Substrata
str_auslander1 <- pData(gse_auslander[[1]])$`treatment state:ch1`
# DESeq2
gex_auslander1 <- DESeq2::DESeqDataSetFromMatrix(countData=sup_auslander[,gsub(".bam", "", gsub("-", ".", pData(gse_auslander[[1]])$title))], 
	colData = data.frame("str_auslander" = pData(gse_auslander[[1]])$`treatment state:ch1`), 
	design = ~ str_auslander)
gex_auslander1 <- estimateSizeFactors(gex_auslander1)
gex_auslander1 <- counts(gex_auslander1, normalized=TRUE)
# Substrata
str_auslander2 <- pData(gse_auslander[[2]])$`treatment state:ch1`
# DESeq2
gex_auslander2 <- DESeq2::DESeqDataSetFromMatrix(countData=sup_auslander[,gsub(".bam", "", gsub("-", ".", pData(gse_auslander[[2]])$title))], 
	colData = data.frame("str_auslander" = pData(gse_auslander[[2]])$`treatment state:ch1`), 
	design = ~ str_auslander)
gex_auslander2 <- estimateSizeFactors(gex_auslander2)
gex_auslander2 <- counts(gex_auslander2, normalized=TRUE)
# Combine the two batches
gex_auslander <- cbind(gex_auslander1, gex_auslander2)
# Collapse over unique gene IDs, use mean
gex_auslander <- do.call("rbind", by(gex_auslander, INDICES=sup_auslander[,1], FUN=function(z) { apply(z, MARGIN=2, FUN=mean) }))
# Spurious gene names such as '1-Dec', '1-Mar', '1-Sep', '10-Mar', '10-Sep', ... ??
# Gene names shouldn't be repeated
# Include only ones included in genes from biomaRt
gex_auslander <- gex_auslander[which(rownames(gex_auslander) %in% genes$hgnc_symbol),]
#> dim(gex_auslander)
#[1] 20245    37
# Much more reasonable row count - previously over 70k?!

##
## Westin et al.
##
#str(gse_westin[[1]]@featureData)
#head(gse_westin[[1]]@featureData@data)


##
## Ascierto et al.
##

##
## Ascierto et al. 2
##


##
## Riaz et al.
##
#rownames(sup_riaz2) <- genes[match(sup_riaz2[,1], genes[,"entrezgene_id"]),"hgnc_symbol"]
#> table(table(genes[match(sup_riaz2[,1], genes[,"entrezgene_id"]),"hgnc_symbol"]))
#
#    1     2     3   217 
#20584    25     1     1
#
# 20584 had unique hugo symbol mapping, filtering the rest ("" occurred 217 times)
sup_riaz2[,1] <- genes[match(sup_riaz2[,1], genes[,"entrezgene_id"]),"hgnc_symbol"]
sup_riaz2 <- sup_riaz2[-which(sup_riaz2[,1] %in% names(table(sup_riaz2[,1])>1)[which(table(sup_riaz2[,1])>1)]),]
sup_riaz2 <- sup_riaz2[-which(is.na(sup_riaz2[,1])),]
rownames(sup_riaz2) <- sup_riaz2[,1]
sup_riaz2 <- sup_riaz2[,-1]
# Substrata
str_riaz <- pData(gse_riaz[[1]])$`characteristics_ch1`
# DESeq2
gex_riaz <- DESeq2::DESeqDataSetFromMatrix(countData=sup_riaz2,
	colData = data.frame("str_riaz" = str_riaz),
	design = ~ str_riaz)
gex_riaz <- estimateSizeFactors(gex_riaz)
gex_riaz <- counts(gex_riaz, normalized=TRUE)

##
## Prat et al.
##
gex_prat <- exprs(gse_prat[[1]])
#> dim(gex_prat)
#[1] 765  65
gex_prat <- as.matrix(gex_prat)
# Removing NA-rows with no values observed
gex_prat <- gex_prat[-which(apply(gex_prat, MARGIN=1, FUN=function(z) { all(is.na(z)) })),]
#> dim(gex_prat)
#[1] 725  65
#> sum(is.na(gex_prat))
#[1] 0
