## FICAN-OSCAR ##
## DREAM challenge:  subchallenge 1 ##

library(data.table)
library(EPIC)
library(parallel)
library(BiocGenerics)
library(rlang)
library(immunedeconv)
library(survival)
library(stats4)
library(S4Vectors)
library(IRanges)
library(Biobase)
library(XML)
library(AnnotationDbi)
library(annotate)
library(graph)
library(GSEABase)
library(GSVA)
library(dplyr)


###
###############
### data in ###
###############
# Getting command line args: 
# args[1] input rna-seq gene tpm data
# args[2] clinical data
# args[3] output file: two columns (patient id, prediction) .csv file
args        <- commandArgs(trailingOnly=TRUE)
gex_input      <- fread(args[1],data.table = F)  ## Gene tpm
dat_input      <- fread(args[2],data.table = F)  ## clinical data

rownames(gex_input) <- gex_input[,1]
gex_input  <- gex_input[,-1]
gex_input <- gex_input[order(rownames(gex_input)),]
gex_input <- as.matrix(gex_input)
rownames(dat_input) <- dat_input[,1]
dat_input <- dat_input[,-1]
print("Done reading in tpm and clinical data")

# DREAM IO: !PD response

# Transformation for gene expression data
logz <- function(x) { 
	tmp <- scale(log(x+1)) 
	if(any(!is.finite(tmp))){
		tmp[!is.finite(tmp)] <- 0
	}
	tmp
}

# Aggregation
aggregateX <- function(
	gex,
	dat,
	normalize = TRUE
){
	library(GSVA)
	library(immunedeconv)
	
	# Whether whole gene matrix should be normalized
	if(normalize){
		nam <- colnames(gex)
		gex <- apply(gex, MARGIN=1, FUN=logz)
		rownames(gex) <- nam
		gex <- t(gex) # Transpose needed for GSVA etc
	}
	
	# 242 mutations threshold for 'high' mutatation; impute NAs as if being not highly mutated
	X <- data.frame(isTMBhigh = as.integer(dat[,"TMB"] > 242))
	rownames(X) <- rownames(dat)
	# Impute zero indicators if there are NA values
	if(any(is.na(X[,"isTMBhigh"]))) X[is.na(X[,"isTMBhigh"]),"isTMBhigh"] <- 0
	
	X <- cbind(X, isMale = as.integer(dat[,"SEX"] == "M"))

	X <- cbind(X, isSquamous = as.integer(dat[,"CRFHIST"] == "SQUAMOUS"))

	X <- cbind(X, isSquamous_above5PDL1 = ifelse(dat[,"PDL1"]>=5, X[,"isSquamous"], 0))

	X <- cbind(X, isEversmoker = as.integer(dat[,"TOBACUSE"] %in% c("FORMER", "CURRENT")))

	X <- cbind(X, isECOG0 = as.integer(dat[,"ECOGPS"] == 0))
	
	# GSVA 
	gmt_custom <- GSEABase::getGmt("selfmade.gmt")
	res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE, mx.diff=TRUE)) # Custom GMTs, mx.diff = TRUE as the genes in the FICAN-OSCAR panel are concordantly in one direction
	X <- cbind(X, CUSTOM_FOPANEL = res_gsva[,"CUSTOM_FOPANEL"])
	
	# Cell expression as reported by xCell
	#tmp <- immunedeconv::deconvolute(gex, method="xcell")
	#tmp <- as.matrix(tmp)
	#rownames(tmp) <- paste("xce_", gsub(" ", "_", tmp[,1]), sep="")
	#tmp <- tmp[,-1]
	#class(tmp) <- "numeric"
	#X <- cbind(X, xce_Endothelial_cell = t(tmp)[,"xce_Endothelial_cell"])	
	
	# Return X as a matrix
	as.matrix(X)
}
# Xb
predictX <- function(
	X,	# curateX-generated X data matrix
	b	# Vector of beta coefficient values; these will be name-matched in combining the correct fields from data matrix X and produce final prediction value
){
	# Make sure X is handled as a numeric matrix
	X <- as.matrix(X)
	class(X) <- "numeric"
	# Subset to make the dimension fit for beta
	if(any(!names(b) %in% colnames(X))){
		stop(paste("Warning! Variables not found in X but present in beta (will be omitted):", b[!names(b) %in% colnames(X)]), collapse=", ")
		b <- b[which(names(b) %in% colnames(X))]
	}	
	Xhat <- X[,names(b)]
	# Patient IDs ought to be present as rownames in X
	IDs <- rownames(Xhat)
	# Return a data frame with patientID and prediction (matrix product) columns
	ret <- data.frame(
		patientID = IDs,
		prediction = Xhat %*% b
	)
	# Flip sign
	ret[,"prediction"] <- -ret[,"prediction"]
	rownames(ret) <- NULL
	ret
}

# Beta vector
## Subchallenge 3 submission 1
#b_sub = c(	
#	"CUSTOM_MCP_ENDOTHELIAL.CELLS" = -1.0566, # Presence of endothelial cells according to GSVA reduces response likelihood
#	"CUSTOM_IFNG3" = 1.3909 # Presence of Interferon-gamma signalling increases response likelihood
#)
## -> DSS -0.0436, Nivo 0.5149, Chemo 0.4052

## Subchallenge 3 submission 2
#b_sub = c(	
#	"CD274" = 0.03635064, # mRNA expression of PD-L1; increases likelihood of response; aggregate estimate from all RESP data
#	"log10TMB" = 1, # log10 of tumor mutational burden; increases likelihood of response; literature based estimate
#	"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = -0.5 # EMT, conservative estimate derived from Hugo et al.
#)
## -> DSS -0.0404, Nivo 0.4319, Chemo 0.6026

# To test: CD8A and isMale came out from Prat et al, others borderline (CD274 p=0.11, isEversmoker p=0.2); TMB known from literature to predict and does not correlate with PD-L1 (or CD274)
# Coefs from Prat et al. in PD prediction (based on available limited sample size and variable data):
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.23187    0.21616   1.073   0.2908    
#CD274         0.12319    0.07513   1.640   0.1100    
#CD8A          0.15857    0.07597   2.087   0.0442 *  
#isEversmoker -0.35441    0.27158  -1.305   0.2004    
#isMale        0.85753    0.18241   4.701 3.95e-05 ***
#
## TCGA control didn't identify coefs as significant


## Subchallenge 3 submission 3
#b_sub = c(
#	"CD274" = 0.12319, # Estimated from Prat et al.
#	"CD8A" = 0.15857, # Estimated from Prat et al.
#	"isPDL1high" = 1, # a priori
#	"isTMBhigh" = 1, # a priori
#	"isPDL1high_isTMBhigh" = 2, # a priori
#	"isMale" = 0.25 # Being male had a significant effect in predicting !PD in Prat et al., coef estimate a priori
#)
## -> DSS -0.068, Nivo 0.4564, Chemo 0.3961

## Subchallenge 3 submission 4
#b_sub = c(
#	"CUSTOM_FOPANEL" = 1, # The GSVA is bimodal around -0.5, and 0.5 with individual variation; giving equal weight to the gene set is done if both have the same coefficient
#	"isTMBhigh" = 1 # a priori equal importance, albeit binary indicator
#)
## -> DSS -0.0029 [-0.003, 0.0035], Nivo 0.4434, Chemo 0.4813

## Subchallenge 3 submission 5 - using correlation between OS and ORR as basis for reusing OS model (Shukuya et al. 2016)
# also Gyawali et al. which included CheckMates: "The treatment effect sizes between PFS and OS for RCTs of non–small cell lung cancer 
# trials were similar (rHR, 1.14; 95% CI, 0.99-1.31), whereas for the other cancer types, there was a greater effect on OS than on PFS 
# (rHR, 1.23; 95% CI, 1.05-1.44)" 
b_sub = c(
	"CUSTOM_FOPANEL" = -log(0.5), # Custom trained panel estimated using GSVA; estimate from aggregated estimates from training data (e.g. Prat, Gide)
	"isTMBhigh" = -log(0.7), # In chemo low TMB advantage but not high or medium, in nivo high TMB advantage; pick high separately, estimate from Cristescu et al. omitting GEP
	"isMale" = -log(0.9), # Brahmer HR = 0.57 (SQ),  Borghaei HR = 0.73 (non-SQ); conservative estimate; significant effect in Prat et al.
	"isSquamous" = -log(0.82), # Carbone et al. for OS; HR = 0.82 if PD-L1 any
	"isSquamous_above5PDL1" = -log(0.95), # Carbone et al. for OS; HR = 0.77 if PD-L1 >=5% alternatively; adding a small additive effect for PDL1 >= 5% IHC
	"isEversmoker" = -log(0.8), # A manually curated average from Brahmer et al. and Borghaei et al., seemed a consistent OS univariate effect between SQ vs. non-SQ
	"isECOG0" = -log(0.9) # A manually curated average from Brahmer et al. and Borghaei et al., seemed a consistent OS univariate effect between SQ vs. non-SQ
)
## -> DSS x

## REPLACE
#setwd("..")
#gex_input <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
#gex_input <- gex_input[order(rownames(gex_input)),]
#gex_input <- as.matrix(gex_input)
#dat_input <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
#setwd("Docker_subch3")
## REPLACE ENDS

X_sub <- aggregateX(
	gex = gex_input,	# Gene expression matrix of actual data
	dat = dat_input	# Clinical variable matrix of actual data
)

# Xb
ret <- predictX(
	X = X_sub, # Data matrix
	b = b_sub # Beta coefficients
)
# Sanity checking
#head(cbind(ret, dat_input[,c("TMB")], X_sub[,c("isTMBhigh", "isMale", "isSquamous", "isSquamous_above5PDL1", "CUSTOM_FOPANEL")]))

# Write output
write.csv(ret, file = args[3], quote=F, row.names=FALSE)
print("Done writing out result")

