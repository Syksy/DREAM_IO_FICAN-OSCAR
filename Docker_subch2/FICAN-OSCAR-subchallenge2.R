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

# Submission 1, PFS

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
	
	# Hard-coded Carbone et al. thresholds
	# 242 mutations threshold for 'high' mutatation; impute NAs as if being not highly mutated
	X <- data.frame(isTMBhigh = as.integer(dat[,"TMB"] > 242))
	rownames(X) <- rownames(dat)
	# Impute zero indicators if there are NA values
	if(any(is.na(X[,"isTMBhigh"]))) X[is.na(X[,"isTMBhigh"]),"isTMBhigh"] <- 0
	
	# There was no significant association between tumor-mutation burden and 
	# PD-L1 expression level (Pearson’s correlation coefficient = 0.059).	
	# --> Complementary information
	
	# However, in the nivolumab group, patients with both a high tumor-mutation burden and a 
	# PD-L1 expression level of 50% or more had a higher response rate (75%) than those with 
	# only one of these factors (32% among patients with a high tumor-mutation burden only and 
	# 34% among those with a PD-L1 expression level of >50% only) or neither factor (16%).
	
	# 75% -> 3/4
	# 32%, 34 % ~ 1/3
	# 16% ~ 1/7 - 1/8
	# looks like a non-linear trend?
	
	# Thirds for quantiles
	TMBq3s <- quantile(dat[,"TMB"], probs=c(0, 1/3, 2/3, 1), na.rm=TRUE)
	dat[is.na(dat[,"TMB"]),"TMB"] <- median(dat[,"TMB"], na.rm=TRUE)
	PDL1q3s <- quantile(dat[,"PDL1"], probs=c(0, 1/3, 2/3, 1), na.rm=TRUE)
	X <- cbind(X, TMBq3s = findInterval(x=dat[,"TMB"], vec=TMBq3s, rightmost.closed = TRUE))
	X <- cbind(X, PDL1q3s = findInterval(x=dat[,"PDL1"], vec=PDL1q3s, rightmost.closed = TRUE))
	X <- cbind(X, TMB.PDL1 = X[,"TMBq3s"] * X[,"PDL1q3s"])
	# CD274 expression level modelled as a surrogate for PD-L1 IHC
	# Normalized expressions between various platforms and their respective distributional characteristics
	X <- cbind(X, CD274 = gex["CD274",])
	
	# GSVA for CUSTOM_IFNG3
	gmt_custom <- GSEABase::getGmt(".\\selfmade.gmt")
	res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE)) # Custom GMTs
	#X <- cbind(X, IFNG = res_gsva[,grep("CUSTOM_IFNG3", colnames(res_gsva))])
	X <- cbind(X, res_gsva)
	
	# Interaction with T-cell inflammatory signal Gene Expression Profile (GEP) and high tumor mutational burden
	GEPq2s <- median(res_gsva[,"CUSTOM_GEP"])
	X <- cbind(X, GEPq2s = as.integer(res_gsva[,"CUSTOM_GEP"] > GEPq2s))
	# TMB x GEP interaction, upper tertile of TMB vs. upper median of GEP signature
	X <- cbind(X, TMB.GEP = as.integer(X[,"TMBq3s"] == 2) * X[,"GEPq2s"])
	
	# GSVA for Hallmarks
	gmt_hallmarks <- GSEABase::getGmt(".\\h.all.v7.2.symbols.gmt")
	res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_hallmarks, verbose=FALSE)) # Custom GMTs
	#X <- cbind(X, IFNG = res_gsva[,grep("CUSTOM_IFNG3", colnames(res_gsva))])
	X <- cbind(X, res_gsva)
	
	# Epithelial cell expression as reported by xCell
	tmp <- immunedeconv::deconvolute(gex, method="xcell")
	tmp <- as.matrix(tmp)
	rownames(tmp) <- paste("xce_", gsub(" ", "_", tmp[,1]), sep="")
	tmp <- tmp[,-1]
	class(tmp) <- "numeric"
	X <- cbind(X, t(tmp))	
	
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
	# Flip sign of prediction
	ret[,"prediction"] <- -ret[,"prediction"]
	rownames(ret) <- NULL
	ret
}

# Beta vector
b_sub = c(	
	"TMB.PDL1" = log(0.8), # Educated estimate of 0.8 HR transformed into linear coefficient, as interaction is exponentially increasing
	"TMB.GEP" = log(0.7) # HR estimate 0.7 based on Cristescu et al. in Science, of TMB & GEP sign interaction
)

X_sub <- aggregateX(
	gex = gex_input,	# Gene expression matrix of actual data
	dat = dat_input	# Clinical variable matrix of actual data
)

# Xb
ret <- predictX(
	X = X_sub, # Data matrix
	b = b_sub # Beta coefficients
)
# Write output
write.csv(ret, file = args[3], quote=F, row.names=FALSE)
print("Done writing out result")

