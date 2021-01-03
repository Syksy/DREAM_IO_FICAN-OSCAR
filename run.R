###
#
# DREAM 2020; Anti-PD1 Response Prediction DREAM Challenge
#
###

# First setwd to correct project root!
# setwd("...")
setwd("D:\\Gits\\DREAM_2020_IO\\")

# Just in case ...
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
update.packages(ask=FALSE, checkBuilt=TRUE)

# Generate biomaRt gene mapping data.frame for convenience
#BiocManager::install("biomaRt")
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

### Synthetic dataset from DREAM IO 2020 (no signal, must right composition for data)

# Read in synthetic data and examine fields
cli_synthetic <- read.csv(".\\Synthetic\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
gex_synthetic_ensembl75_genes_count <- read.csv(".\\Synthetic\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_ensembl75_genes_count.csv", row.names=1)
gex_synthetic_ensembl75_genes_tpm <- read.csv(".\\Synthetic\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_ensembl75_genes_tpm.csv", row.names=1)
gex_synthetic_ensembl75_isoforms_count <- read.csv(".\\Synthetic\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_ensembl75_isoforms_count.csv", row.names=1)
gex_synthetic_ensembl75_isoforms_tpm <- read.csv(".\\Synthetic\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_ensembl75_isoforms_tpm.csv", row.names=1)
gex_synthetic_refseq105_genes_count <- read.csv(".\\Synthetic\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_count.csv", row.names=1)
gex_synthetic_refseq105_genes_tpm <- read.csv(".\\Synthetic\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
gex_synthetic_refseq105_isoforms_count <- read.csv(".\\Synthetic\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_isoforms_count.csv", row.names=1)
gex_synthetic_refseq105_isoforms_tpm <- read.csv(".\\Synthetic\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_isoforms_tpm.csv", row.names=1)

# Sanity check of dimensions
#> dim(gex_synthetic_ensembl75_genes_count)
#[1] 57997    56
#> dim(gex_synthetic_refseq105_genes_count)
#[1] 29182    56
#> dim(gex_synthetic_ensembl75_genes_tpm)
#[1] 57997    56
#> dim(gex_synthetic_refseq105_genes_tpm)
#[1] 29182    56
#> dim(gex_synthetic_ensembl75_isoforms_count)
#[1] 196593     56
#> dim(gex_synthetic_refseq105_isoforms_count)
#[1] 74375    56

## Non-comforming dimensions, except when normalized count -> tpm




### Corresponding datasets from TCGA

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("RTCGA")
library("RTCGA")

RTCGA::checkTCGA(what="Dates", c("LUAD", "LUSC"))
RTCGA::checkTCGA(what="DataSets", "LUAD")
RTCGA::checkTCGA(what="DataSets", "LUSC")

# Latest "2016-01-28"

# Clinical annotations
RTCGA::downloadTCGA(c("LUAD", "LUSC"), "Merge_Clinical.Level_1", "TCGA", date="2016-01-28", untarFile=TRUE, removeTar=TRUE)

#> list.files("TCGA")
#[1] "gdac.broadinstitute.org_LUAD.Merge_Clinical.Level_1.2016012800.0.0"
#[2] "gdac.broadinstitute.org_LUSC.Merge_Clinical.Level_1.2016012800.0.0"

RTCGA::downloadTCGA(c("LUAD", "LUSC"), "mRNAseq_Preprocess.Level_3", "TCGA", date="2016-01-28", untarFile=TRUE, removeTar=TRUE)

## use cBioPortal and its package cgdsr instead

# install.packages("cgdsr")
library("cgdsr")

mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
#getCancerStudies(mycgds)

#> getCancerStudies(mycgds)[grep("luad", getCancerStudies(mycgds)$cancer_study_id),-3]
#                 cancer_study_id                                         name
#136                   luad_broad       Lung Adenocarcinoma (Broad, Cell 2012)
#137              luad_mskcc_2015    Lung Adenocarcinoma (MSKCC, Science 2015)
#138             luad_oncosg_2020 Lung Adenocarcinoma (OncoSG, Nat Genet 2020)
#139                    luad_tcga  Lung Adenocarcinoma (TCGA, Firehose Legacy)
#140                luad_tcga_pub      Lung Adenocarcinoma (TCGA, Nature 2014)
#141 luad_tcga_pan_can_atlas_2018  Lung Adenocarcinoma (TCGA, PanCancer Atlas)
#142                     luad_tsp       Lung Adenocarcinoma (TSP, Nature 2008)
#> getCancerStudies(mycgds)[grep("lusc", getCancerStudies(mycgds)$cancer_study_id),-3]
#                 cancer_study_id                                                 name
#143                    lusc_tcga Lung Squamous Cell Carcinoma (TCGA, Firehose Legacy)
#144                lusc_tcga_pub     Lung Squamous Cell Carcinoma (TCGA, Nature 2012)
#145 lusc_tcga_pan_can_atlas_2018 Lung Squamous Cell Carcinoma (TCGA, PanCancer Atlas)

#> getCaseLists(mycgds, "luad_tcga")[,-5]
#                  case_list_id                              case_list_name                                        case_list_description cancer_study_id
#1                luad_tcga_all                                 All samples                                    All samples (586 samples)            2715
#2      luad_tcga_3way_complete                            Complete samples Samples with mutation, CNA and expression data (230 samples)            2715
#3                luad_tcga_cna                       Samples with CNA data                          Samples with CNA data (516 samples)            2715
#4    luad_tcga_methylation_all               Samples with methylation data                  Samples with methylation data (580 samples)            2715
#5   luad_tcga_methylation_hm27        Samples with methylation data (HM27)           Samples with methylation data (HM27) (150 samples)            2715
#6  luad_tcga_methylation_hm450       Samples with methylation data (HM450)          Samples with methylation data (HM450) (492 samples)            2715
#7               luad_tcga_mrna Samples with mRNA data (Agilent microarray)               Samples with mRNA expression data (32 samples)            2715
#8    luad_tcga_rna_seq_v2_mrna         Samples with mRNA data (RNA Seq V2)              Samples with mRNA expression data (517 samples)            2715
#9             luad_tcga_cnaseq          Samples with mutation and CNA data             Samples with mutation and CNA data (230 samples)            2715
#10         luad_tcga_sequenced                  Samples with mutation data                     Samples with mutation data (230 samples)            2715
#11              luad_tcga_rppa            Samples with protein data (RPPA)                    Samples protein data (RPPA) (365 samples)            2715

#cli_luad <- cgdsr::getClinicalData(mycgds, caseList=cgdsr::getCaseLists(mycgds, "luad_tcga")[1,"case_ids"])
cli_luad <- cgdsr::getClinicalData(mycgds, caseList="luad_tcga_all")
#cli_lusc <- cgdsr::getClinicalData(mycgds, caseList=cgdsr::getCaseLists(mycgds, "lusc_tcga")[1,"case_ids"])
cli_lusc <- cgdsr::getClinicalData(mycgds, caseList="lusc_tcga_all")

## Response binarization:
#> table(cli_luad$TREATMENT_OUTCOME_FIRST_COURSE)
#
#                            Complete Remission/Response  Partial Remission/Response         Progressive Disease              Stable Disease 
#                        435                         130                           2                           7                          12 
## Samples should be preoperative ("" = ???):
#> table(cli_luad$PERFORMANCE_STATUS_TIMING)
#
#                                      Other Post-Adjuvant Therapy  Pre-Adjuvant Therapy          Preoperative 
#                  355                    23                    11                    20                   177
## Discretized notation of tobacco smoking
#> table(cli_luad$TOBACCO_SMOKING_HISTORY_INDICATOR)
#
#  1   2   3   4   5 
# 76 122 137 171   4
## level 1 -> never smoked

#> cgdsr::getGeneticProfiles(mycgds, "luad_tcga")[,-c(3,4,6)]
#                                    genetic_profile_id                                                   genetic_profile_name genetic_alteration_type
#1                                       luad_tcga_rppa                                              Protein expression (RPPA)           PROTEIN_LEVEL
#2                               luad_tcga_rppa_Zscores                                     Protein expression z-scores (RPPA)           PROTEIN_LEVEL
#3                                     luad_tcga_gistic                           Putative copy-number alterations from GISTIC  COPY_NUMBER_ALTERATION
#4                                       luad_tcga_mrna                                           mRNA expression (microarray)         MRNA_EXPRESSION
#5                        luad_tcga_mrna_median_Zscores      mRNA expression z-scores relative to diploid samples (microarray)         MRNA_EXPRESSION
#6                            luad_tcga_rna_seq_v2_mrna                                      mRNA expression (RNA Seq V2 RSEM)         MRNA_EXPRESSION
#7             luad_tcga_rna_seq_v2_mrna_median_Zscores mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM)         MRNA_EXPRESSION
#8                                 luad_tcga_linear_CNA                              Capped relative linear copy-number values  COPY_NUMBER_ALTERATION
#9                           luad_tcga_methylation_hm27                                                     Methylation (HM27)             METHYLATION
#10                         luad_tcga_methylation_hm450                                                    Methylation (HM450)             METHYLATION
#11                                 luad_tcga_mutations                                                              Mutations       MUTATION_EXTENDED
#12 luad_tcga_rna_seq_v2_mrna_median_all_sample_Zscores mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM)         MRNA_EXPRESSION
#13            luad_tcga_mrna_median_all_sample_Zscores      mRNA expression z-scores relative to all samples (log microarray)         MRNA_EXPRESSION

#> cgdsr::getGeneticProfiles(mycgds, "lusc_tcga")[,-c(3,4,6)]
#                                    genetic_profile_id                                                        genetic_profile_name genetic_alteration_type
#1                                       lusc_tcga_rppa                                                   Protein expression (RPPA)           PROTEIN_LEVEL
#2                               lusc_tcga_rppa_Zscores                                          Protein expression z-scores (RPPA)           PROTEIN_LEVEL
#3                                     lusc_tcga_gistic                                Putative copy-number alterations from GISTIC  COPY_NUMBER_ALTERATION
#4                                  lusc_tcga_mrna_U133                                      mRNA expression (U133 microarray only)         MRNA_EXPRESSION
#5                          lusc_tcga_mrna_U133_Zscores mRNA expression z-scores relative to diploid samples (U133 microarray only)         MRNA_EXPRESSION
#6                                       lusc_tcga_mrna                                                mRNA expression (microarray)         MRNA_EXPRESSION
#7                        lusc_tcga_mrna_median_Zscores           mRNA expression z-scores relative to diploid samples (microarray)         MRNA_EXPRESSION
#8                            lusc_tcga_rna_seq_v2_mrna                                           mRNA expression (RNA Seq V2 RSEM)         MRNA_EXPRESSION
#9             lusc_tcga_rna_seq_v2_mrna_median_Zscores      mRNA expression z-scores relative to diploid samples (RNA Seq V2 RSEM)         MRNA_EXPRESSION
#10                                lusc_tcga_linear_CNA                                   Capped relative linear copy-number values  COPY_NUMBER_ALTERATION
#11                          lusc_tcga_methylation_hm27                                                          Methylation (HM27)             METHYLATION
#12                         lusc_tcga_methylation_hm450                                                         Methylation (HM450)             METHYLATION
#13                                 lusc_tcga_mutations                                                                   Mutations       MUTATION_EXTENDED
#14 lusc_tcga_rna_seq_v2_mrna_median_all_sample_Zscores      mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM)         MRNA_EXPRESSION
#15            lusc_tcga_mrna_median_all_sample_Zscores           mRNA expression z-scores relative to all samples (log microarray)         MRNA_EXPRESSION
#16              lusc_tcga_mrna_U133_all_sample_Zscores     mRNA expression z-scores relative to all samples (U133 microarray only)         MRNA_EXPRESSION


## MUTATION_COUNT --> Number of genes with relevant mutations (proxy for TMB?)

generate_cbioportal <- function(
  genes = sort(unique(curatedPCaData_genes$hgnc_symbol)), 
  geneticProfiles = c("luad_tcga_rna_seq_v2_mrna_median_all_sample_Zscores", # Lung adenocarcinoma gex 
                      "lusc_tcga_rna_seq_v2_mrna_median_all_sample_Zscores", # Lung squamous gex
                      ), # for cgdsr calls, platform and dataset specific string
  caseList = c("luad_tcga_all", # Lung adenocarcinoma
               "lusc_tcga_all"  # Lung squamous
               ), # for cgdsr calls, platform and dataset specific string
  delay = 0.05, 
  splitsize = 100, 
  verb = TRUE
){
  # If given genes is a list (with slots for various annotation types), try to extract hugo gene symbols
  if(class(genes)=="list"){
    genes <- genes$hgnc_symbol
  }
  # Establisigh connection to cBioPortal
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
  # Split gene name vector into suitable lengths
  genesplit <- rep(1:ceiling(length(genes)/splitsize), 
                   each = splitsize)[1:length(genes)]
  splitgenes <- split(genes, f = genesplit)
  # Bind the API calls as per columns
  gex <- as.matrix(do.call("cbind", lapply(1:length(splitgenes), FUN = function(z){
      if(verb == TRUE) cat(paste(z, "of", length(splitgenes), "\n"))
      # Sleep if necessary to avoid API call overflow
      Sys.sleep(delay)
      # Fetch each split gene name list from the URL-based API, essentially a wrapper for cgdsr's own function
      cgdsr::getProfileData(mycgds, genes = splitgenes[[z]], 
                            geneticProfiles = geneticProfiles, caseList = caseList)
    })))  
  
  gex <- t(gex)
  gex
}  

## -> 
# luad_tcga_rna_seq_v2_mrna_median_all_sample_Zscores
# lusc_tcga_rna_seq_v2_mrna_median_all_sample_Zscores

empty.omit <- function(x) { if("" %in% x){ x[-which(x=="")] } else { x }}

## Download TCGAs for lung for LUAD and LUSC

gex_luad <- generate_cbioportal(
	genes = empty.omit(sort(unique(genes$hgnc_symbol))),
	geneticProfiles = "luad_tcga_rna_seq_v2_mrna_median_all_sample_Zscores",
	caseList = "luad_tcga_all"
)
gex_lusc <- generate_cbioportal(
	genes = empty.omit(sort(unique(genes$hgnc_symbol))),
	geneticProfiles = "lusc_tcga_rna_seq_v2_mrna_median_all_sample_Zscores",
	caseList = "lusc_tcga_all"
)

# Clean up clinical information from cli_luad and cli_lusc
cli_luad <- data.frame(
	PFS = cli_luad$OS_MONTHS,
	PFS.Event = cli_luad$,
	OS = ,
	OS.Event = 
	Response =
	Age = ,
	Tobacco,
	T = ,
	N = ,
	M = ,
	RACE = cli_luad$RACE
	SEX = cli_luad$SEX
)
	





### Datasets downloaded from TIDE

# Lauss et al., Nat Commun 2017
# PFS, OS and 0/1 Response available
gex_lauss <- read.table(".\\TIDE\\Lauss2017_ACT_Melanoma_RNASeq\\ICB.Lauss2017_ACT_Melanoma.self_subtract", header=TRUE)
cli_lauss <- read.table(".\\TIDE\\Lauss2017_ACT_Melanoma_RNASeq\\ICB.Lauss2017_ACT_Melanoma.clinical", header=TRUE)
# Map to ensembl/hugo/etc
gex_lauss[,1] <- genes[match(gex_lauss[,1], genes$entrezgene_id),"hgnc_symbol"]
gex_lauss <- do.call("rbind", by(gex_lauss, INDICES=gex_lauss[,1], FUN=function(z){ apply(z[,-1], MARGIN=2, FUN=mean) }))
# Take intersection of patients present in both
pat_lauss <- intersect(rownames(cli_lauss), colnames(gex_lauss))
gex_lauss <- gex_lauss[,pat_lauss]
cli_lauss <- cli_lauss[pat_lauss,]

# Kim et al., Nat Medicine 2018
# (Pembrolizumab)
# Response 0/1 available

gex_kim <- read.table(".\\TIDE\\Kim2018_PD1_Gastric_RNASeq\\ICB.Kim2018_Pembrolizumab_Gastric.self_subtract", header=TRUE)
cli_kim <- read.table(".\\TIDE\\Kim2018_PD1_Gastric_RNASeq\\ICB.Kim2018_Pembrolizumab_Gastric.clinical", header=TRUE, nrows=57)
# Map to ensembl/hugo/etc
gex_kim[,1] <- genes[match(gex_kim[,1], genes$entrezgene_id),"hgnc_symbol"]
gex_kim <- do.call("rbind", by(gex_kim, INDICES=gex_kim[,1], FUN=function(z){ apply(z[,-1], MARGIN=2, FUN=mean) }))
# cli_kim missing some responses in the end
# '-' symbols transform to '.' in column names
cli_kim$patient <- gsub('-', '.', cli_kim$patient)
rownames(cli_kim) <- cli_kim[,1]
cli_kim <- cli_kim[,-1,drop=FALSE]
# Take intersection of patients present in both
pat_kim <- intersect(rownames(cli_kim), colnames(gex_kim))
gex_kim <- gex_kim[,pat_kim]
cli_kim <- cli_kim[pat_kim,,drop=FALSE]


# Chen et al., Cancer Discov 2016
# Response 0/1 available
gex_chen <- read.table(".\\TIDE\\Chen2016_PD1_Melanoma_Nanostring_Ipi.Prog\\ICB.Chen2016_PD1_Melanoma_Ipi.Prog.self_subtract", header=TRUE)
cli_chen <- read.table(".\\TIDE\\Chen2016_PD1_Melanoma_Nanostring_Ipi.Prog\\ICB.Chen2016_PD1_Melanoma_Ipi.Prog.clinical", header=TRUE)
# P6 duplicate omitted
gex_chen <- gex_chen[,-6]
cli_chen <- cli_chen[-5,]
# Map to ensembl/hugo/etc
gex_chen[,1] <- genes[match(gex_chen[,1], genes$entrezgene_id),"hgnc_symbol"]
gex_chen <- do.call("rbind", by(gex_chen, INDICES=gex_chen[,1], FUN=function(z){ apply(z[,-1], MARGIN=2, FUN=mean) }))
rownames(cli_chen) <- cli_chen[,1]
cli_chen <- cli_chen[,-1,drop=FALSE]
# Take intersection of patients present in both
pat_chen <- intersect(rownames(cli_chen), colnames(gex_chen))
gex_chen <- gex_chen[,pat_chen]
cli_chen <- cli_chen[pat_chen,,drop=FALSE]













###
#
# Derive simple sample GSVA results to refine the genes to pathway-level information
#
###

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
library(GSVA)

# Read chosen set of GMTs; namely: hallmarks, oncogenic, and immunologic pathways
# Downloaded from e.g. https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C7
gmt_h <- GSEABase::getGmt(".\\MSigDB\\h.all.v7.2.symbols.gmt")
gmt_c6 <- GSEABase::getGmt(".\\MSigDB\\c6.all.v7.2.symbols.gmt")
gmt_c7 <- GSEABase::getGmt(".\\MSigDB\\c7.all.v7.2.symbols.gmt")

gmt_lauss <- rbind(
	GSVA::gsva(gex_lauss, gmt_h),	# Hallmarks
	GSVA::gsva(gex_lauss, gmt_c6),	# Oncogenic
	GSVA::gsva(gex_lauss, gmt_c7)	# Immunology
)
gmt_kim <- rbind(
	GSVA::gsva(gex_kim, gmt_h),	# Hallmarks
	GSVA::gsva(gex_kim, gmt_c6),	# Oncogenic
	GSVA::gsva(gex_kim, gmt_c7)	# Immunology
)
gmt_chen <- rbind(
	GSVA::gsva(gex_chen, gmt_h),	# Hallmarks
	GSVA::gsva(gex_chen, gmt_c6),	# Oncogenic
	GSVA::gsva(gex_chen, gmt_c7)	# Immunology
)



###
#
# Immune cell deconvolution, using immunedeconv-package
#
###

#install.packages("remotes")
#remotes::install_github("icbi-lab/immunedeconv")
library(immunedeconv)

## Lauss
idc_lauss <- rbind(
	immunedeconv::deconvolute(gex_lauss, method="xcell"),
	immunedeconv::deconvolute(gex_lauss, method="mcp_counter")
)
tmp <- idc_lauss[,1]
idc_lauss <- as.matrix(idc_lauss[,-1])
rownames(idc_lauss) <- c(tmp)[[1]]
## Kim
idc_kim <- rbind(
	immunedeconv::deconvolute(gex_kim, method="xcell"),
	immunedeconv::deconvolute(gex_kim, method="mcp_counter")
)
tmp <- idc_kim[,1]
idc_kim <- as.matrix(idc_kim[,-1])
rownames(idc_kim) <- c(tmp)[[1]]
## Chen
# - unable to run, key genes missing -


rm(tmp)

###
#
# DeMixT
#
###

# https://github.com/wwylab/DeMixT
# devtools::install_github("wwylab/DeMixT")



###
#
# ESTIMATE
#
###

# https://www.nature.com/articles/ncomms3612
# https://bioinformatics.mdanderson.org/public-software/estimate/
# http://r-forge.r-project.org/R/?group_id=2237
# install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)
library(estimate)

# Nevermind, requires .gct input

# Preliminary testing whether any of the gmt / gene sets are found by lasso
if(FALSE){
	library(glmnet)
	library(survival)
	# LASSO testing

	### Lauss

	set.seed(1)
	# PFS
	lasso_sub1_lauss <- glmnet(x=t(gmt_lauss), y=survival::Surv(time=cli_lauss[,"PFS"], event=cli_lauss[,"PFS.Event"]), family="cox")
	lassocv_sub1_lauss <- cv.glmnet(x=t(gmt_lauss), y=survival::Surv(time=cli_lauss[,"PFS"], event=cli_lauss[,"PFS.Event"]), family="cox", nfolds=5)
	plot(lassocv_sub1_lauss)
	# No non-zeros...
	# OS
	lasso_sub2_lauss <- glmnet(x=t(gmt_lauss), y=survival::Surv(time=cli_lauss[,"OS"], event=cli_lauss[,"OS.Event"]), family="cox")
	lassocv_sub2_lauss <- cv.glmnet(x=t(gmt_lauss), y=survival::Surv(time=cli_lauss[,"OS"], event=cli_lauss[,"OS.Event"]), family="cox", nfolds=5)
	plot(lassocv_sub2_lauss)
	# No non-zeros...
	# Response
	lasso_sub3_lauss <- glmnet(x=t(gmt_lauss), y=cli_lauss[,"Response"], family="binomial")
	lassocv_sub3_lauss <- cv.glmnet(x=t(gmt_lauss), y=cli_lauss[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_lauss)

	### Kim

	set.seed(1)
	# Response
	lasso_sub3_kim <- glmnet(x=t(gmt_kim), y=cli_kim[,"Response"], family="binomial")
	lassocv_sub3_kim <- cv.glmnet(x=t(gmt_kim), y=cli_kim[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_kim)
	#> rownames(gmt_kim)[predict(lasso_sub3_kim, s=lassocv_sub3_kim$lambda.min, type="nonzero")[,1]]
	# [1] "IL15_UP.V1_DN"                                                            
	# [2] "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN"                                      
	# [3] "GSE17721_4H_VS_24H_POLYIC_BMDC_DN"                                        
	# [4] "GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_UP"                                    
	# [5] "GSE29615_CTRL_VS_DAY7_LAIV_FLU_VACCINE_PBMC_DN"                           
	# [6] "GSE36476_YOUNG_VS_OLD_DONOR_MEMORY_CD4_TCELL_40H_TSST_ACT_UP"             
	# [7] "GSE8515_IL1_VS_IL6_4H_STIM_MAC_DN"                                        
	# [8] "GSE22601_DOUBLE_NEGATIVE_VS_IMMATURE_CD4_SP_THYMOCYTE_UP"                 
	# [9] "GSE2585_CD80_HIGH_VS_LOW_MTEC_UP"                                         
	#[10] "GSE3920_IFNA_VS_IFNG_TREATED_FIBROBLAST_DN"                               
	#[11] "GSE6875_WT_VS_FOXP3_KO_TREG_DN"                                           
	#[12] "GSE7831_CPG_VS_INFLUENZA_STIM_PDC_4H_UP"                                  
	#[13] "GSE12484_HEALTHY_VS_PERIDONTITIS_NEUTROPHILS_DN"                          
	#[14] "GSE7459_UNTREATED_VS_IL6_TREATED_ACT_CD4_TCELL_DN"                        
	#[15] "GSE19374_UNINF_VS_LISTERIA_INFECTED_MACROPHAGE_UP"                        
	#[16] "GSE22432_CONVENTIONAL_CDC_VS_PLASMACYTOID_PDC_DN"                         
	#[17] "GSE34156_NOD2_LIGAND_VS_NOD2_AND_TLR1_TLR2_LIGAND_24H_TREATED_MONOCYTE_UP"
	#[18] "GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_DN"                     
	#[19] "GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN"                      
	#[20] "GSE42724_MEMORY_VS_B1_BCELL_UP"                                           
	#[21] "GSE46606_IRF4_KO_VS_WT_CD40L_IL2_IL5_1DAY_STIMULATED_BCELL_UP"
	#> rownames(gmt_kim)[predict(lasso_sub3_kim, s=lassocv_sub3_kim$lambda.1se, type="nonzero")[,1]]
	#[1] "IL21_UP.V1_DN"                                                "GSE10325_MYELOID_VS_LUPUS_MYELOID_UP"                        
	#[3] "GSE36476_YOUNG_VS_OLD_DONOR_MEMORY_CD4_TCELL_40H_TSST_ACT_UP" "GSE2585_CD80_HIGH_VS_LOW_MTEC_UP"                            
	#[5] "GSE41176_UNSTIM_VS_ANTI_IGM_STIM_TAK1_KO_BCELL_3H_DN"

	### Chen

	set.seed(1)
	# Response
	lasso_sub3_chen <- glmnet(x=t(gmt_chen), y=cli_chen[,"Response"], family="binomial")
	## Of note:
	#Warning message:
	#In lognet(x, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
	#  one multinomial or binomial class has fewer than 8  observations; dangerous ground
	lassocv_sub3_chen <- cv.glmnet(x=t(gmt_chen), y=cli_chen[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_chen)
	#rownames(gmt_chen)[predict(lasso_sub3_chen, s=lassocv_sub3_chen$lambda.min, type="nonzero")[,1]]
	#> rownames(gmt_chen)[predict(lasso_sub3_chen, s=lassocv_sub3_chen$lambda.min, type="nonzero")[,1]]
	#[1] "GSE22601_IMMATURE_CD4_SINGLE_POSITIVE_VS_DOUBLE_POSITIVE_THYMOCYTE_UP"
	#[2] "GSE37301_HEMATOPOIETIC_STEM_CELL_VS_CD4_TCELL_DN"


}






# Preliminary testing of whether glmnet grabs anything
if(FALSE){

	### Lauss

	set.seed(1)
	# PFS
	lasso_sub1_lauss <- glmnet(x=t(idc_lauss), y=survival::Surv(time=cli_lauss[,"PFS"], event=cli_lauss[,"PFS.Event"]), family="cox")
	lassocv_sub1_lauss <- cv.glmnet(x=t(idc_lauss), y=survival::Surv(time=cli_lauss[,"PFS"], event=cli_lauss[,"PFS.Event"]), family="cox", nfolds=5)
	plot(lassocv_sub1_lauss)
	#> rownames(idc_lauss)[predict(lasso_sub1_lauss, s=lassocv_sub1_lauss$lambda.min, type="nonzero")[,1]]
	# [1] "Myeloid dendritic cell activated" "T cell CD4+ (non-regulatory)"     "T cell CD4+ central memory"      
	# [4] "T cell CD4+ effector memory"      "T cell CD8+ effector memory"      "Common lymphoid progenitor"      
	# [7] "Common myeloid progenitor"        "Cancer associated fibroblast"     "Hematopoietic stem cell"         
	#[10] "Macrophage M2"                    "Mast cell"                        "Neutrophil"                      
	#[13] "NK cell"                          "T cell NK"                        "B cell plasma"                   
	#[16] "T cell CD4+ Th1"                  "T cell CD4+ Th2"                  "T cell regulatory (Tregs)"       
	#[19] "stroma score"                     "cytotoxicity score"               "NK cell"                         
	#[22] "Myeloid dendritic cell"           "Neutrophil"
	# --> All cell types except couple...
	# OS
	lasso_sub2_lauss <- glmnet(x=t(idc_lauss), y=survival::Surv(time=cli_lauss[,"OS"], event=cli_lauss[,"OS.Event"]), family="cox")
	lassocv_sub2_lauss <- cv.glmnet(x=t(idc_lauss), y=survival::Surv(time=cli_lauss[,"OS"], event=cli_lauss[,"OS.Event"]), family="cox", nfolds=5)
	plot(lassocv_sub2_lauss)
	#> rownames(idc_lauss)[predict(lasso_sub2_lauss, s=lassocv_sub2_lauss$lambda.min, type="nonzero")[,1]]
	# [1] "Myeloid dendritic cell activated" "T cell CD4+ (non-regulatory)"     "T cell CD4+ effector memory"     
	# [4] "T cell CD8+ naive"                "T cell CD8+ effector memory"      "Class-switched memory B cell"    
	# [7] "Common lymphoid progenitor"       "Cancer associated fibroblast"     "Hematopoietic stem cell"         
	#[10] "Monocyte"                         "Neutrophil"                       "NK cell"                         
	#[13] "T cell NK"                        "T cell gamma delta"               "T cell CD4+ Th1"                 
	#[16] "T cell regulatory (Tregs)"        "stroma score"                     "T cell"                          
	#[19] "NK cell"                          "Myeloid dendritic cell"           "Neutrophil"                      
	#[22] "Endothelial cell"                 "Cancer associated fibroblast"    
	#> rownames(idc_lauss)[predict(lasso_sub2_lauss, s=lassocv_sub2_lauss$lambda.1se, type="nonzero")[,1]]
	# [1] "Myeloid dendritic cell activated" "T cell CD4+ effector memory"      "T cell CD8+ naive"               
	# [4] "T cell CD8+ effector memory"      "Common lymphoid progenitor"       "Granulocyte-monocyte progenitor" 
	# [7] "Hematopoietic stem cell"          "Macrophage M2"                    "T cell CD4+ Th1"                 
	#[10] "T cell regulatory (Tregs)"        "cytotoxicity score"               "NK cell"
	# --> Highly informative for both PFS/OS
	# Response
	lasso_sub3_lauss <- glmnet(x=t(idc_lauss), y=cli_lauss[,"Response"], family="binomial")
	lassocv_sub3_lauss <- cv.glmnet(x=t(idc_lauss), y=cli_lauss[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_lauss)
	# No non-zeros...

	### Kim

	set.seed(1)
	# Response
	lasso_sub3_kim <- glmnet(x=t(idc_kim), y=cli_kim[,"Response"], family="binomial")
	lassocv_sub3_kim <- cv.glmnet(x=t(idc_kim), y=cli_kim[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_kim)
	#> rownames(idc_kim)[predict(lasso_sub3_kim, s=lassocv_sub3_kim$lambda.min, type="nonzero")[,1]]
	#[1] "T cell CD4+ central memory"   "Macrophage"                   "T cell CD4+ Th1"              "Endothelial cell"            
	#[5] "Cancer associated fibroblast"
	#> rownames(idc_kim)[predict(lasso_sub3_kim, s=lassocv_sub3_kim$lambda.1se, type="nonzero")[,1]]
	#[1] "Macrophage"                   "T cell CD4+ Th1"              "Endothelial cell"             "Cancer associated fibroblast"

}

