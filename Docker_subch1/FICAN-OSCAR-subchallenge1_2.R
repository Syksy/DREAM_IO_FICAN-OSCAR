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

# Submission 2, PFS

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
	# Carbone et al.
	# For efficacy analyses, patients were grouped in thirds according to tumor-mutation burden. 
	# The boundaries for these three
	# groups were a tumor-mutation burden of 0 to less than 100 mutations (low burden), 100 to
	# 242 mutations (medium burden), and 243 or more mutations (high burden)
	# ...
	# Among the patients with a high tumor-mutation burden, the response rate was higher in the
	# nivolumab group than in the chemotherapy group (47% vs. 28%), and progression-free
	# survival was longer (median, 9.7 vs. 5.8 months; hazard ratio for disease progression or
	# death, 0.62; 95% CI, 0.38 to 1.00) (Fig. 2C).
	# ...
	# However, in the nivolumab group, patients with both a high
	# tumor-mutation burden and a PD-L1 expression level of 50% or more had a higher response
	# rate (75%) than those with only one of these factors
	# ...
	# Patients with both a high tumor-mutation
	# burden and a PD-L1 expression level of 50% or more may have a greater likelihood of re
	# sponse to nivolumab than those with only one or neither of these factors. Overall, the current
	# findings are consistent with the hypothesis that immunotherapy may have enhanced activity
	# in patients with a high tumor-mutation burden.14 However, because this was an exploratory
	#
	# ---> Fig 2
	library(GSVA)
	library(immunedeconv)
	
	# 242 mutations threshold for 'high' mutatation; impute NAs as if being not highly mutated
	X <- data.frame(isTMBhigh = as.integer(dat[,"TMB"] > 242))
	rownames(X) <- rownames(dat)
	# Impute zero indicators if there are NA values
	if(any(is.na(X[,"isTMBhigh"]))) X[is.na(X[,"isTMBhigh"]),"isTMBhigh"] <- 0
	
	X <- cbind(X, isMale = as.integer(dat[,"SEX"] == "M"))

	X <- cbind(X, isSquamous = as.integer(dat[,"CRFHIST"] == "SQUAMOUS"))

	X <- cbind(X, isEversmoker = as.integer(dat[,"TOBACCOUSE"] %in% c("FORMER", "CURRENT")))
	
	# CD274 expression level modelled as a surrogate for PD-L1 IHC
	# Normalized expressions between various platforms and their respective distributional characteristics
	X <- cbind(X, CD274 = logz(gex["CD274",]))
	
	# GSVA for CUSTOM_IFNG3
	#gmt_custom <- GSEABase::getGmt("selfmade.gmt")
	#res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE)) # Custom GMTs
	#X <- cbind(X, IFNG = res_gsva[,grep("CUSTOM_IFNG3", colnames(res_gsva))])
	
	# Epithelial cell expression as reported by xCell
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
	rownames(ret) <- NULL
	ret
}

# Beta vector
## Subchallenge 1 submission 1
#b_sub = c(	
#	"CD274" = -0.5326, # Estimated from aggregated PFS data
#	"isTMBhigh" = log(0.62) # Carbone et al. cutoff and HR into a coef
#)
## Subchallenge 1 submission 2
b_sub = c(
	"CD274" = -0.5326, # = same as HR 0.5870766; has the most effect per unit
	"isTMBhigh" = log(0.62), # In chemo low TMB advantage but not high or medium, in nivo high TMB advantage; pick high separately
	"isMale" = log(0.8), # Brahmer HR = 0.57 (SQ),  Borghaei HR = 0.73 (non-SQ)
	"isSquamous" = log(0.77) # Carbone et al. supplementary effect for OS
	"isEversmoker" = log(0.7) # Huant et al. OncoImmunology 2018 pooled meta-analysis
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

