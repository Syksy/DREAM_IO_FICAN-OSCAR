###
#
# DREAM 2020; Anti-PD1 Response Prediction DREAM Challenge
#
###

# First setwd to correct project root!
# setwd("...")
setwd("C:\\Users\\Daniel\\DREAM_2020_IO\\")

# Just in case ...
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

library(survival)

# Datasets downloaded from TIDE

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




###
#
# Immune cell deconvolution
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

