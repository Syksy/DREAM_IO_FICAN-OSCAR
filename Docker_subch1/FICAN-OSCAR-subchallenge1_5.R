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

# DREAM IO, PFS

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
	# TMB as a continuous value
	X <- cbind(X, TMB = dat[,"TMB"])
	if(any(is.na(X[,"TMB"]))) X[is.na(X[,"TMB"]),"TMB"] <- median(dat[,"TMB"], na.rm=TRUE)
	
	# Very often used characteristics
	X <- cbind(X, isMale = as.integer(dat[,"SEX"] == "M"))
	X <- cbind(X, isSquamous = as.integer(dat[,"CRFHIST"] == "SQUAMOUS"))
	X <- cbind(X, isSquamous_above5PDL1 = ifelse(dat[,"PDL1"]>=5, X[,"isSquamous"], 0))	
	X <- cbind(X, isEversmoker = as.integer(dat[,"TOBACUSE"] %in% c("FORMER", "CURRENT")))
	X <- cbind(X, isECOG0 = as.integer(dat[,"ECOGPS"] == 0))

	# GSVA, gene panel developed for DREAM IO having the biggest effect on prediction
	gmt_custom <- GSEABase::getGmt("selfmade.gmt")
	res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE, mx.diff=TRUE)) # Custom GMTs, assuming identified gene effects are unidirectional with mx.diff
	X <- cbind(X, CUSTOM_FOPANEL = res_gsva[,"CUSTOM_FOPANEL"])
	
	# xCell
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
## Subchallenge 1 submission 1
#b_sub = c(	
#	"CD274" = -0.5326, # Estimated from aggregated PFS data
#	"isTMBhigh" = log(0.62) # Carbone et al. cutoff and HR into a coef
#)
## -> DSS 0.0353, Nivo 0.4064 (flip 1-0.4064), Chemo 0.5118

## Subchallenge 1 submission 2 
# missing exact coefs, last submitted coefs buried in version control
#b_sub = c(	
#	"CD274" = ,
#	"isTMBhigh" =, # In chemo low TMB advantage but not high or medium, in nivo high TMB advantage; pick high separately
#	"isMale" =, # Brahmer HR = 0.57 (SQ),  Borghaei HR = 0.73 (non-SQ)
#	"isSquamous" = log(0.77) # Carbone et al. supplementary effect for OS
#	"isSquamousEversmoker" = log(0.59), # Brahmer et al. for SQ eversmoker HR
#	"isNonsquamousEversmoker" = log(0.7), # Borhaei et al for non-SQ eversmoker HR
#	"isEversmoker" = log(0.7) # Huant et al. OncoImmunology 2018 pooled meta-analy
#)
## -> DSS 0.0291, Nivo 0.5897, Chemo 0.4638

## Subchallenge 1 submission 3
#b_sub = c(
#	"TMB" = -0.0025, # Estimate a continuous value for the TMB effect, which seems to focus towards hypermutated samples
#	"CUSTOM_FOPANELv2" = -1, # Custom trained panel estimated using GSVA; estimate from aggregated estimates from training data (e.g. Prat, Gide)
#	"isMale" = log(0.9), # Brahmer HR = 0.57 (SQ),  Borghaei HR = 0.73 (non-SQ); lesser effect expected in this CheckMate trial
#	"isSquamous" = log(0.83) # Carbone et al. PFS effect HR
#)
## -> DSS 0.0125, Nivo 0.5591, Chemo 0.5007

## Subchallenge 1 submission 4
#b_sub = c(
#	"isTMBhigh" = log(0.62), # Carbone et al. cutoff and HR into a coef
#	"CUSTOM_FOPANELv2" = log(0.5), # Custom trained panel estimated using GSVA; estimate from aggregated estimates from training data (e.g. Prat, Gide)
#	"isSquamous" = log(0.83) # Carbone et al. PFS effect HR; in Brahmer et al (SQ) PFS don't cross, in Borghaei et al. (non-SQ) the PFS vs Chemo crosses
#)
## -> DSS 0.0144 [0.0027, 0.0214], Nivo 0.56, Chemo 0.5074

## Subchallenge 1 submission 5
b_sub = c(
	"CUSTOM_FOPANEL" = log(0.5), # Custom trained panel estimated using GSVA; estimate from aggregated estimates from training data (e.g. Prat, Gide)
	"isTMBhigh" = log(0.7), # In chemo low TMB advantage but not high or medium, in nivo high TMB advantage; pick high separately, estimate from Cristescu et al. omitting GEP, supported by literature reviews for Nivo-specificity
	"isMale" = log(0.9), # Brahmer HR = 0.57 (SQ),  Borghaei HR = 0.73 (non-SQ); conservative estimate; significant effect in Prat et al., supported by literature reviews for Nivo-specificity
	"isSquamous" = log(0.82), # Derived from Carbone et al.; HR = 0.82 if PD-L1 any, supported by literature reviews for Nivo-specificity
	"isSquamous_above5PDL1" = log(0.95), # Derived from Carbone et al.; HR = 0.77 if PD-L1 >=5% alternatively; adding a small additive effect for PDL1 >= 5% IHC
	"isEversmoker" = log(0.8), # A manually curated average from Brahmer et al. and Borghaei et al., seemed a consistent PFS univariate effect between SQ vs. non-SQ, supported by literature reviews for Nivo-specificity
	"isECOG0" = log(0.9) # A manually curated average from Brahmer et al. and Borghaei et al., seemed a consistent PFS univariate effect between SQ vs. non-SQ, supported by literature reviews for Nivo-specificity
)
## -> DSS 0.0143 [0.0094, 0.0152], Nivo 0.5676, Chemo 0.4849

## REPLACE
#setwd("..")
#gex_input <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
#gex_input <- gex_input[order(rownames(gex_input)),]
#gex_input <- as.matrix(gex_input)
#dat_input <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
#setwd("Docker_subch2")
## REPLACE ENDS

X_sub <- aggregateX(
	gex = gex_input,	# Gene expression matrix of actual data
	dat = dat_input	# Clinical variable matrix of actual data
)
## For sanity checking in synthetic data input
# cbind(ret, TMB = dat_input[,"TMB"], FOPANEL = X_sub[,"CUSTOM_FOPANEL"], SEX = dat_input[,"SEX"], CRFHIST = dat_input[,"CRFHIST"], CD274 = gex_input["CD274",], PDCD1 = gex_input["PDCD1",], TIGIT = gex_input["TIGIT",], CXCL9 = gex_input["CXCL9",], CXCR6 = gex_input["CXCR6",], CD8A = gex_input["CD8A",], CCL5 = gex_input["CCL5",])
# head(cbind(ret, dat_input[,c("TMB")], X_sub[,c("isTMBhigh", "isMale", "isSquamous", "isSquamous_above5PDL1", "CUSTOM_FOPANEL")]))

# Xb
ret <- predictX(
	X = X_sub, # Data matrix
	b = b_sub # Beta coefficients
)
# Write output
write.csv(ret, file = args[3], quote=F, row.names=FALSE)
print("Done writing out result")

