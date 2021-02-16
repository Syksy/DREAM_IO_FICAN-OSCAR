# Submission 1, PFS
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
	
	# CD274 expression level modelled as a surrogate for PD-L1 IHC
	# Normalized expressions between various platforms and their respective distributional characteristics
	X <- cbind(X, CD274 = logz(gex["CD274",]))
	
	# GSVA for CUSTOM_IFNG3
	gmt_custom <- GSEABase::getGmt(".\\selfmade.gmt")
	res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE)) # Custom GMTs
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
	rownames(ret) <- NULL
	ret
}
# Transformation for gene expression data
logz <- function(x) { 
	tmp <- scale(log(x+1)) 
	if(any(!is.finite(tmp))){
		tmp[!is.finite(tmp)] <- 0
	}
	tmp
}
# Beta bector
b_sub = c(	
	"CUSTOM_MCP_ENDOTHELIAL.CELLS" = -1.0566, # Presence of endothelial cells according to GSVA reduces response likelihood
	"CUSTOM_IFNG3" = 1.3909 # Presence of Interferon-gamma signalling increases response likelihood
)
# Modified X data matrix
## REPLACE
setwd("..")
gex_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
gex_synthetic <- gex_synthetic[order(rownames(gex_synthetic)),]
gex_synthetic <- as.matrix(gex_synthetic)
dat_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
setwd("Docker_subch3")
## REPLACE ENDS
# REPLACE THESE WITH THE ACTUAL DATA INPUT
X_sub <- aggregateX(
	gex = gex_synthetic,	# Gene expression matrix of actual data
	dat = dat_synthetic	# Clinical variable matrix of actual data
)
## REPLACE PARAMETERS ABOVE ENDS
# Xb
ret <- predictX(
	X = X_sub, # Data matrix
	b = b_sub # Beta coefficients
)
# Write output
#write.csv(ret, file = "predictions.csv", row.names=FALSE)




## Safety checking


