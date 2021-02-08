##
#
# Source code (functions) for performing predict in the Docker container for the DREAM IO
# A preliminary example for subchallenge 1 (PFS)
#
##

predictX <- function(
	#gex,	# Gene expression matrix
	#dat,	# Data matrix containing clinical variables (in training cohorts also end-points)
	X,	# curateX-generated X data matrix
	b,	# Vector of beta coefficient values; these will be name-matched in combining the correct fields from data matrix X and produce final prediction value
	verb = FALSE	# Should there be verbosity, i.e. debugging information etc
){
	# Double checking input b
	if(verb){
		print("Beta vector:")
		print(b)
		print(class(b))
	}
	# Make sure X is handled as a numeric matrix
	X <- as.matrix(X)
	class(X) <- "numeric"
	# Subset to make the dimension fit for beta
	if(any(!names(b) %in% colnames(X))){
		warning(paste("Warning! Variables not found in X but present in beta (will be omitted):", b[!names(b) %in% colnames(X)]), collapse=", ")
		# NOTE! 
		# INSIDE DOCKER FOR DREAM IO
		# THESE ERROR MESSAGES ARE NOT RETURNED
		# AND MAY RESULT IN A SPURIOUS MODEL PREDICTION!!!
		# By removing this if-based sanity checking model is in these cases invalid, which may be better for sending in predictions
		b <- b[which(names(b) %in% colnames(X))]
	}	
	Xhat <- X[,names(b)]
	if(verb){
		print(head(Xhat))
	}
	# Patient IDs ought to be present as rownames in X
	IDs <- rownames(Xhat)
	if(verb){
		print(head(IDs))
	}
	# Return a data frame with patientID and prediction (matrix product) columns
	ret <- data.frame(
		patientID = IDs,
		prediction = Xhat %*% b
	)
	rownames(ret) <- NULL
	ret
}

# Double-check paths; this assumes they are in subfolders
# NOTE! Replace with actual data here instead of synthetic example
gex_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
gex_synthetic <- gex_synthetic[order(rownames(gex_synthetic)),]
gex_synthetic <- as.matrix(gex_synthetic)
dat_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
# Curate the derived variables and hand-picked gene expressions
X_sub1 <- curateX(
	gex = gex_synthetic,	# Gene expression matrix
	dat = dat_synthetic	# Clinical variable matrix
)

# Named vector for linear model Xb product
# Since we only care of the "ordering of the patients, it doesn't matter if we use the link function after Xb product or not
# NOTE!! 
# Important to double check that the final model coefficient directions are correct and of correct magnitude!
# Especially if hard-coded manually from the final derived model
b_sub1 = c(	
	"CD274" = -0.12345, 
	"HALLMARK_DNA_REPAIR" = 123.45,
	"KRAS.LUNG_UP.V1_UP" = 1.23,
	"xce_T_cell_CD8plus_effector_memory" = -0.987,
	"xce_Cancer_associated_fibroblast" = 0.456
)


sub1 <- predictX(
	X = X_sub1,
	b = b_sub1,
	verb = TRUE
)

#> head(sub1)
#  patientID  prediction
#1      p227   0.7237105
#2      p359  -2.6058215
#3      p533  -8.3748367
#4      p149 -15.4823535
#5      p160 -19.3129264
#6        p2   0.3425207
#> head(X_sub1[,names(b_sub1)])
#     CD274 HALLMARK_DNA_REPAIR KRAS.LUNG_UP.V1_UP xce_T_cell_CD8plus_effector_memory xce_Cancer_associated_fibroblast
#p227 36.45          0.04122774         0.08759159                       3.910469e-19                     5.737034e-02
#p359  4.18         -0.01544484        -0.14889064                       1.705507e-17                     0.000000e+00
#p533 25.34         -0.04123325        -0.12712911                       0.000000e+00                     1.300881e-18
#p149  8.89         -0.11686467         0.02592495                       2.411474e-18                     2.230791e-02
#p160 28.69         -0.12694975        -0.09972296                       1.035231e-17                     5.144833e-02
#p2   15.80          0.02014395        -0.15751219                       0.000000e+00                     0.000000e+00


# Saving the output "predictions.csv" for subchallenge 1 in this case
# NOTE! Different subchallenges need different output here
# NOTE! May require parameter quote=FALSE
write.csv(sub1, file = "predictions.csv", row.names=FALSE)
