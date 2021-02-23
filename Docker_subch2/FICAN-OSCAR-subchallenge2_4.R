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

# DREAM IO: OS

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
	# Flip the sign
	ret[,"prediction"] <- -ret[,"prediction"]
	rownames(ret) <- NULL
	ret
}

# Beta vector
## Subchallenge 2 submission 1
#b_sub = c(	
#	"TMB.PDL1" = log(0.8), # Educated estimate of 0.8 HR transformed into linear coefficient, as interaction is exponentially increasing
#	"TMB.GEP" = log(0.7) # HR estimate 0.7 based on Cristescu et al. in Science, of TMB & GEP sign interaction
#)
## -> DSS -0.0143, Nivo 0.5359, Chemo 0.5636

## Subchallenge 2 submission 2
#b_sub = c(
#	"CUSTOM_FOPANEL" = -0.5, # Custom trained panel estimated using GSVA; estimate from aggregated estimates from training data (e.g. Prat, Gide)
#	"isTMBhigh" = log(0.7), # In chemo low TMB advantage but not high or medium, in nivo high TMB advantage; pick high separately, estimate from Cristescu et al. omitting GEP
#	"isMale" = log(0.9), # Brahmer HR = 0.57 (SQ),  Borghaei HR = 0.73 (non-SQ); conservative estimate
#	"isSquamous" = log(0.82) # Carbone et al. supplementary effect for OS; 0.72 if PD-L1 >5% alternatively
#)
## -> DSS 0.0295 [0.0237, 0.0386], Nivo 0.5883, Chemo 0.5025

## Subchallenge 2 submission 3
#b_sub = c(
#	"CUSTOM_FOPANEL" = -0.5, # Custom trained panel estimated using GSVA; estimate from aggregated estimates from training data (e.g. Prat, Gide)
#	"isTMBhigh" = log(0.7), # In chemo low TMB advantage but not high or medium, in nivo high TMB advantage; pick high separately, estimate from Cristescu et al. omitting GEP
#	"isMale" = log(0.9), # Brahmer HR = 0.57 (SQ),  Borghaei HR = 0.73 (non-SQ); conservative estimate; significant effect in Prat et al.
#	"isSquamous" = log(0.82), # Carbone et al. for OS; HR = 0.82 if PD-L1 any
#	"isSquamous_above5PDL1" = log(0.95), # Carbone et al. for OS; HR = 0.77 if PD-L1 >=5% alternatively; adding a small additive effect for PDL1 >= 5% IHC
#	"isEversmoker" = log(0.8), # A manually curated average from Brahmer et al. and Borghaei et al., seemed a consistent OS univariate effect between SQ vs. non-SQ
#	"isECOG0" = log(0.7) # A manually curated average from Brahmer et al. and Borghaei et al., seemed a consistent OS univariate effect between SQ vs. non-SQ
#)
## -> DSS 0.0313 [0.0142, 0.0416], Nivo 0.5941, Chemo 0.5341

## Subchallenge 2 submission 4
b_sub = c(
	"CUSTOM_FOPANEL" = log(0.5), # Custom trained panel estimated using GSVA; estimate from aggregated estimates from training data (e.g. Prat, Gide)
	"isTMBhigh" = log(0.7), # In chemo low TMB advantage but not high or medium, in nivo high TMB advantage; pick high separately, estimate from Cristescu et al. omitting GEP
	"isMale" = log(0.9), # Brahmer HR = 0.57 (SQ),  Borghaei HR = 0.73 (non-SQ); conservative estimate; significant effect in Prat et al.
	"isSquamous" = log(0.82), # Carbone et al. for OS; HR = 0.82 if PD-L1 any
	"isSquamous_above5PDL1" = log(0.95), # Carbone et al. for OS; HR = 0.77 if PD-L1 >=5% alternatively; adding a small additive effect for PDL1 >= 5% IHC
	"isEversmoker" = log(0.8), # A manually curated average from Brahmer et al. and Borghaei et al., seemed a consistent OS univariate effect between SQ vs. non-SQ
	"isECOG0" = log(0.9) # A manually curated average from Brahmer et al. and Borghaei et al., seemed a consistent OS univariate effect between SQ vs. non-SQ
)
## -> DSS 0.0447 [0.0266, 0.064], Nivo 0.6074, Chemo 0.4895


## REPLACE
#setwd("..")
#gex_input <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
#gex_input <- gex_input[order(rownames(gex_input)),]
#gex_input <- as.matrix(gex_input)
#dat_input <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
#setwd("Docker_subch2")
## REPLACE ENDS

# https://onlinelibrary.wiley.com/doi/epdf/10.1002/cam4.2807 -> never smoker is a marker for poor prognosis (or vice versa?)
X_sub <- aggregateX(
	gex = gex_input,	# Gene expression matrix of actual data
	dat = dat_input	# Clinical variable matrix of actual data
)
# Xb
ret <- predictX(
	X = X_sub, # Data matrix
	b = b_sub # Beta coefficients
)
# cbind(ret, TMB = dat_input[,"TMB"], FOPANEL = X_sub[,"CUSTOM_FOPANEL"], SEX = dat_input[,"SEX"], CRFHIST = dat_input[,"CRFHIST"], 
# 	isECOG0 = X_sub[,"isECOG0"], isTMBhigh = X_sub[,"isTMBhigh"], isEversmoker = X_sub[,"isEversmoker"], isSquamous = X_sub[,"isSquamous"], isSquamous_above5PDL1 = X_sub[,"isSquamous_above5PDL1"], PDL1 = dat_input[,"PDL1"], 
# Base panel
# 	CD274 = gex_input["CD274",], PDCD1 = gex_input["PDCD1",], TIGIT = gex_input["TIGIT",], CXCL9 = gex_input["CXCL9",], CXCR6 = gex_input["CXCR6",], 
# Part of wider panel
# 	CD8A = gex_input["CD8A",], CCL5 = gex_input["CCL5",])
# Write output
write.csv(ret, file = args[3], quote=F, row.names=FALSE)
print("Done writing out result")

