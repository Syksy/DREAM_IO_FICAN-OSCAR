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

# Submission 3, PD or !PD

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
	dat
){
	library(GSVA)
	library(immunedeconv)
	
	# 242 mutations threshold for 'high' mutatation; impute NAs as if being not highly mutated
	X <- data.frame(isTMBhigh = dat[,"TMB"]>=242)
	rownames(X) <- rownames(dat)
	# Impute zero indicators if there are NA values
	if(any(is.na(X[,"isTMBhigh"]))) X[is.na(X[,"isTMBhigh"]),"isTMBhigh"] <- 0
	
	X <- cbind(X, isPDL1high = as.integer(dat[,"PDL1"]>=50)
	if(any(is.na(X[,"isPDL1high"]))) X[is.na(X[,"isPDL1high"]),"isPDL1high"] <- 0
	
	X <- cbind(X, isPDL1high_isTMBhigh = X[,"isTMBhigh"]*X[,"isPDL1high"])
	# No NAs since imputed in previous steps
	
	X <- cbind(X, isMale = as.integer(dat[,"SEX"] == "M"))
	
	# CD274 expression level modelled as a surrogate for PD-L1 IHC
	# Normalized expressions between various platforms and their respective distributional characteristics
	#X <- cbind(X, CD274 = logz(gex["CD274",]))
	
	# GSVA 
	#gmt_custom <- GSEABase::getGmt("selfmade.gmt")
	#gmt_h <- GSEABase::getGmt("h.all.v7.2.symbols.gmt")
	#res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE)) # Custom GMTs
	#X <- cbind(X, res_gsva)
	#res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_h, verbose=FALSE)) # Hallmark GMTs
	#X <- cbind(X, res_gsva)
	#
	# Epithelial cell expression as reported by xCell
	#tmp <- immunedeconv::deconvolute(gex, method="xcell")
	#tmp <- as.matrix(tmp)
	#rownames(tmp) <- paste("xce_", gsub(" ", "_", tmp[,1]), sep="")
	#tmp <- tmp[,-1]
	#class(tmp) <- "numeric"
	#X <- cbind(X, t(tmp))	
	
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
## Subchallenge 3 submission 2
#b_sub = c(	
#	"CD274" = 0.03635064, # mRNA expression of PD-L1; increases likelihood of response; aggregate estimate from all RESP data
#	"log10TMB" = 1, # log10 of tumor mutational burden; increases likelihood of response; literature based estimate
#	"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = -0.5 # EMT, conservative estimate derived from Hugo et al.
#)
## Subchallenge 3 submission 3
b_sub = c(	
	"isPDL1high" = 1, # 
	"isTMBhigh" = 1, # 
	"isPDL1high_isTMBhigh" = 1, #
	"isMale" = 0.5 # Being male had a significant effect in predicting PD in Prat et al.
)
# To test: CD8A and isMale came out from Prat et al, others borderline (CD274 p=0.11, isEversmoker p=0.2); TMB known from literature to predict and does not correlate with PD-L1 (or CD274)
# Coefs from Prat et al. in PD prediction (based on available limited sample size and variable data):
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.23187    0.21616   1.073   0.2908    
#CD274         0.12319    0.07513   1.640   0.1100    
#CD8A          0.15857    0.07597   2.087   0.0442 *  
#isEversmoker -0.35441    0.27158  -1.305   0.2004    
#isMale        0.85753    0.18241   4.701 3.95e-05 ***


## REPLACE
#setwd("..")
#gex_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
#gex_synthetic <- gex_synthetic[order(rownames(gex_synthetic)),]
#gex_synthetic <- as.matrix(gex_synthetic)
#dat_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
#setwd("Docker_subch3")
## REPLACE ENDS

X_sub <- aggregateX(
	gex = gex_input,	# Gene expression matrix of actual data
	dat = dat_input	# Clinical variable matrix of actual data
	#gex = gex_synthetic,
	#dat = dat_synthetic
)

# Xb
ret <- predictX(
	X = X_sub, # Data matrix
	b = b_sub # Beta coefficients
)
# Write output
write.csv(ret, file = args[3], quote=F, row.names=FALSE)
print("Done writing out result")

