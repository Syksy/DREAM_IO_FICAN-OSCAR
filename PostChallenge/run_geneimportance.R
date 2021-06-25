####
#
# Test gene importances by omitting them one-by-one from GSVA gene sets and then re-running predictions for them
# Output is ranking of genes, added to the CSV file provided by organizers
#
###


# Commented out to ease debugging below
if(FALSE){

	setwd("D:\\Gits\\DREAM_2020_IO\\PostChallenge\\")
	load("fits.RData")

	## GSVA
	library(GSVA)	

	# Hallmark gene sets
	gmt_all <- GSEABase::getGmt("FOPANEL.gmt")
	gmt_noCD274 <- GSEABase::getGmt("FOPANEL_noCD274.gmt")
	gmt_noCXCL9 <- GSEABase::getGmt("FOPANEL_noCXCL9.gmt")
	gmt_noCXCR6 <- GSEABase::getGmt("FOPANEL_noCXCR6.gmt")
	gmt_noPDCD1 <- GSEABase::getGmt("FOPANEL_noPDCD1.gmt")
	gmt_noTIGIT <- GSEABase::getGmt("FOPANEL_noTIGIT.gmt")

	genes <- c("CD274", "CXCL9", "CXCR6", "PDCD1", "TIGIT")

	# Genes as part of GSVA pathway inference, leaving out one gene at a time
	gsva_prat <- do.call("cbind", list(
		gsva_prat_all = t(GSVA::gsva(as.matrix(X_prat), gmt_all, verbose=FALSE))[,1],	
		gsva_prat_noCD274 = t(GSVA::gsva(as.matrix(X_prat), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_prat_noCXCL9 = t(GSVA::gsva(as.matrix(X_prat), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_prat_noCXCR6 = t(GSVA::gsva(as.matrix(X_prat), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_prat_noPDCD1 = t(GSVA::gsva(as.matrix(X_prat), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_prat_noTIGIT = t(GSVA::gsva(as.matrix(X_prat), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_gide <- do.call("cbind", list(
		gsva_gide_all = t(GSVA::gsva(as.matrix(X_gide), gmt_all, verbose=FALSE))[,1],	
		gsva_gide_noCD274 = t(GSVA::gsva(as.matrix(X_gide), gmt_noCD274, verbose=FALSE))[,1],
		gsva_gide_noCXCL9 = t(GSVA::gsva(as.matrix(X_gide), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_gide_noCXCR6 = t(GSVA::gsva(as.matrix(X_gide), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_gide_noPDCD1 = t(GSVA::gsva(as.matrix(X_gide), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_gide_noTIGIT = t(GSVA::gsva(as.matrix(X_gide), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_hugo <- do.call("cbind", list(
		gsva_hugo_all = t(GSVA::gsva(as.matrix(X_hugo), gmt_all, verbose=FALSE))[,1],	
		gsva_hugo_noCD274 = t(GSVA::gsva(as.matrix(X_hugo), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_hugo_noCXCL9 = t(GSVA::gsva(as.matrix(X_hugo), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_hugo_noCXCR6 = t(GSVA::gsva(as.matrix(X_hugo), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_hugo_noPDCD1 = t(GSVA::gsva(as.matrix(X_hugo), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_hugo_noTIGIT = t(GSVA::gsva(as.matrix(X_hugo), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_lauss <- do.call("cbind", list(
		gsva_lauss_all = t(GSVA::gsva(as.matrix(X_lauss), gmt_all, verbose=FALSE))[,1],	
		gsva_lauss_noCD274 = t(GSVA::gsva(as.matrix(X_lauss), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_lauss_noCXCL9 = t(GSVA::gsva(as.matrix(X_lauss), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_lauss_noCXCR6 = t(GSVA::gsva(as.matrix(X_lauss), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_lauss_noPDCD1 = t(GSVA::gsva(as.matrix(X_lauss), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_lauss_noTIGIT = t(GSVA::gsva(as.matrix(X_lauss), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_kim <- do.call("cbind", list(
		gsva_kim_all = t(GSVA::gsva(as.matrix(X_kim), gmt_all, verbose=FALSE))[,1],	
		gsva_kim_noCD274 = t(GSVA::gsva(as.matrix(X_kim), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_kim_noCXCL9 = t(GSVA::gsva(as.matrix(X_kim), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_kim_noCXCR6 = t(GSVA::gsva(as.matrix(X_kim), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_kim_noPDCD1 = t(GSVA::gsva(as.matrix(X_kim), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_kim_noTIGIT = t(GSVA::gsva(as.matrix(X_kim), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_chen <- do.call("cbind", list(
		gsva_chen_all = t(GSVA::gsva(as.matrix(X_chen), gmt_all, verbose=FALSE))[,1],	
		gsva_chen_noCD274 = t(GSVA::gsva(as.matrix(X_chen), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_chen_noCXCL9 = t(GSVA::gsva(as.matrix(X_chen), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_chen_noCXCR6 = t(GSVA::gsva(as.matrix(X_chen), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_chen_noPDCD1 = t(GSVA::gsva(as.matrix(X_chen), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_chen_noTIGIT = t(GSVA::gsva(as.matrix(X_chen), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_westin <- do.call("cbind", list(
		gsva_westin_all = t(GSVA::gsva(as.matrix(X_westin), gmt_all, verbose=FALSE))[,1],	
		gsva_westin_noCD274 = t(GSVA::gsva(as.matrix(X_westin), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_westin_noCXCL9 = t(GSVA::gsva(as.matrix(X_westin), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_westin_noCXCR6 = t(GSVA::gsva(as.matrix(X_westin), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_westin_noPDCD1 = t(GSVA::gsva(as.matrix(X_westin), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_westin_noTIGIT = t(GSVA::gsva(as.matrix(X_westin), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_riaz <- do.call("cbind", list(
		gsva_riaz_all = t(GSVA::gsva(as.matrix(X_riaz), gmt_all, verbose=FALSE))[,1],	
		gsva_riaz_noCD274 = t(GSVA::gsva(as.matrix(X_riaz), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_riaz_noCXCL9 = t(GSVA::gsva(as.matrix(X_riaz), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_riaz_noCXCR6 = t(GSVA::gsva(as.matrix(X_riaz), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_riaz_noPDCD1 = t(GSVA::gsva(as.matrix(X_riaz), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_riaz_noTIGIT = t(GSVA::gsva(as.matrix(X_riaz), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_braun_nivo <- do.call("cbind", list(
		gsva_braun_nivo_all = t(GSVA::gsva(as.matrix(X_braun_nivo), gmt_all, verbose=FALSE))[,1],	
		gsva_braun_nivo_noCD274 = t(GSVA::gsva(as.matrix(X_braun_nivo), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_braun_nivo_noCXCL9 = t(GSVA::gsva(as.matrix(X_braun_nivo), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_braun_nivo_noCXCR6 = t(GSVA::gsva(as.matrix(X_braun_nivo), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_braun_nivo_noPDCD1 = t(GSVA::gsva(as.matrix(X_braun_nivo), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_braun_nivo_noTIGIT = t(GSVA::gsva(as.matrix(X_braun_nivo), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_braun_ever <- do.call("cbind", list(
		gsva_braun_ever_all = t(GSVA::gsva(as.matrix(X_braun_ever), gmt_all, verbose=FALSE))[,1],	
		gsva_braun_ever_noCD274 = t(GSVA::gsva(as.matrix(X_braun_ever), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_braun_ever_noCXCL9 = t(GSVA::gsva(as.matrix(X_braun_ever), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_braun_ever_noCXCR6 = t(GSVA::gsva(as.matrix(X_braun_ever), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_braun_ever_noPDCD1 = t(GSVA::gsva(as.matrix(X_braun_ever), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_braun_ever_noTIGIT = t(GSVA::gsva(as.matrix(X_braun_ever), gmt_noTIGIT, verbose=FALSE))[,1]	
	))
	gsva_tcga <- do.call("cbind", list(
		gsva_tcga_all = t(GSVA::gsva(as.matrix(X_tcga), gmt_all, verbose=FALSE))[,1],	
		gsva_tcga_noCD274 = t(GSVA::gsva(as.matrix(X_tcga), gmt_noCD274, verbose=FALSE))[,1],	
		gsva_tcga_noCXCL9 = t(GSVA::gsva(as.matrix(X_tcga), gmt_noCXCL9, verbose=FALSE))[,1],	
		gsva_tcga_noCXCR6 = t(GSVA::gsva(as.matrix(X_tcga), gmt_noCXCR6, verbose=FALSE))[,1],	
		gsva_tcga_noPDCD1 = t(GSVA::gsva(as.matrix(X_tcga), gmt_noPDCD1, verbose=FALSE))[,1],	
		gsva_tcga_noTIGIT = t(GSVA::gsva(as.matrix(X_tcga), gmt_noTIGIT, verbose=FALSE))[,1]	
	))

	# Genes as univariate predictors without GSVA pathway inference
	unigene_prat <- t(X_prat[genes,])
	unigene_gide <- t(X_gide[genes,])
	unigene_hugo <- t(X_hugo[genes,])
	unigene_lauss <- t(X_lauss[genes,])
	unigene_kim <- t(X_kim[genes,])
	#> genes[!genes %in% rownames(X_chen)]
	#[1] "TIGIT"
	## Chen et al. is missing TIGIT
	unigene_chen <- t(X_chen[genes[-which(genes=="TIGIT")],])
	unigene_westin <- t(X_westin[genes,])
	unigene_riaz <- t(X_riaz[genes,])
	unigene_braun_nivo <- t(X_braun_nivo[genes,])
	unigene_braun_ever <- t(X_braun_ever[genes,])
	unigene_tcga <- t(X_tcga[genes,])

	## Save the minimal workspace for debugging below!
	save(file="generanking_workspace.RData", list=grep("dat_|X_|PFS_|OS_|RESP_|gsva_|unigene_", ls(), value=TRUE))
}




###
#
# oscar >= v0.7.x expected
#
###
library(oscar)
library(survival)

# Commented out due to runs done on June 24th
if(FALSE){

	# Current working directory ought to be inside ..\\PostChallenge\\ where the generanking_workspace.RData is located
	load("generanking_workspace.RData")

	set.seed(1)
	# 0.7.0 needs debugging - memory overflow in below calls:
	# try-function is not catching any errors; something goes wrong in Fortran

	# Prat et al.

	PFS_oscar_gsva_prat <- oscar::oscar(x = gsva_prat, y = PFS_prat, family = 'cox')
	PFS_oscar_unigene_prat <- oscar::oscar(x = unigene_prat, y = PFS_prat, family = 'cox')

	RESP_oscar_gsva_prat <- oscar::oscar(x = gsva_prat, y = RESP_prat, family = 'logistic')
	RESP_oscar_unigene_prat <- oscar::oscar(x = unigene_prat, y = RESP_prat, family = 'logistic')

	# Gide et al.

	PFS_oscar_gsva_gide <- oscar::oscar(x = gsva_gide, y = PFS_gide, family = 'cox')
	PFS_oscar_unigene_gide <- oscar::oscar(x = unigene_gide, y = PFS_gide, family = 'cox')

	OS_oscar_gsva_gide <- oscar::oscar(x = gsva_gide, y = OS_gide, family = 'cox')
	OS_oscar_unigene_gide <- oscar::oscar(x = unigene_gide, y = OS_gide, family = 'cox')

	RESP_oscar_gsva_gide <- oscar::oscar(x = gsva_gide, y = RESP_gide, family = 'logistic')
	RESP_oscar_unigene_gide <- oscar::oscar(x = unigene_gide, y = RESP_gide, family = 'logistic')

	# Hugo et al.

	OS_oscar_gsva_hugo <- oscar::oscar(x = gsva_hugo, y = OS_hugo, family = 'cox')
	OS_oscar_unigene_hugo <- oscar::oscar(x = unigene_hugo, y = OS_hugo, family = 'cox')

	RESP_oscar_gsva_hugo <- oscar::oscar(x = gsva_hugo, y = RESP_hugo, family = 'logistic')
	RESP_oscar_unigene_hugo <- oscar::oscar(x = unigene_hugo, y = RESP_hugo, family = 'logistic')

	# Lauss et al.

	PFS_oscar_gsva_lauss <- oscar::oscar(x = gsva_lauss, y = PFS_lauss, family = 'cox')
	PFS_oscar_unigene_lauss <- oscar::oscar(x = unigene_lauss, y = PFS_lauss, family = 'cox')

	OS_oscar_gsva_lauss <- oscar::oscar(x = gsva_lauss, y = OS_lauss, family = 'cox')
	OS_oscar_unigene_lauss <- oscar::oscar(x = unigene_lauss, y = OS_lauss, family = 'cox')

	RESP_oscar_gsva_lauss <- oscar::oscar(x = gsva_lauss, y = RESP_lauss, family = 'logistic')
	RESP_oscar_unigene_lauss <- oscar::oscar(x = unigene_lauss, y = RESP_lauss, family = 'logistic')

	# Kim et al.

	RESP_oscar_gsva_kim <- oscar::oscar(x = gsva_kim, y = RESP_kim, family = 'logistic')
	RESP_oscar_unigene_kim <- oscar::oscar(x = unigene_kim, y = RESP_kim, family = 'logistic')

	# Chen et al.

	RESP_oscar_gsva_chen <- oscar::oscar(x = gsva_chen, y = RESP_chen, family = 'logistic')
	RESP_oscar_unigene_chen <- oscar::oscar(x = unigene_chen, y = RESP_chen, family = 'logistic')

	# Westin et al.

	PFS_oscar_gsva_westin <- oscar::oscar(x = gsva_westin, y = PFS_westin, family = 'cox')
	PFS_oscar_unigene_westin <- oscar::oscar(x = unigene_westin, y = PFS_westin, family = 'cox')

	# Riaz et al.

	RESP_oscar_gsva_riaz <- oscar::oscar(x = gsva_riaz, y = RESP_riaz, family = 'logistic')
	RESP_oscar_unigene_riaz <- oscar::oscar(x = unigene_riaz, y = RESP_riaz, family = 'logistic')

	# Braun et al. (Nivo & Chemo)

	# Nivo-arm

	PFS_oscar_gsva_braun_nivo <- oscar::oscar(x = gsva_braun_nivo, y = PFS_braun_nivo, family = 'cox')
	PFS_oscar_unigene_braun_nivo <- oscar::oscar(x = unigene_braun_nivo, y = PFS_braun_nivo, family = 'cox')

	OS_oscar_gsva_braun_nivo <- oscar::oscar(x = gsva_braun_nivo, y = OS_braun_nivo, family = 'cox')
	OS_oscar_unigene_braun_nivo <- oscar::oscar(x = unigene_braun_nivo, y = OS_braun_nivo, family = 'cox')

	RESP_oscar_gsva_braun_nivo <- oscar::oscar(x = gsva_braun_nivo, y = RESP_braun_nivo, family = 'logistic')
	RESP_oscar_unigene_braun_nivo <- oscar::oscar(x = unigene_braun_nivo, y = RESP_braun_nivo, family = 'logistic')

	# Chemo-arm

	PFS_oscar_gsva_braun_ever <- oscar::oscar(x = gsva_braun_ever, y = PFS_braun_ever, family = 'cox')
	PFS_oscar_unigene_braun_ever <- oscar::oscar(x = unigene_braun_ever, y = PFS_braun_ever, family = 'cox')

	OS_oscar_gsva_braun_ever <- oscar::oscar(x = gsva_braun_ever, y = OS_braun_ever, family = 'cox')
	OS_oscar_unigene_braun_ever <- oscar::oscar(x = unigene_braun_ever, y = OS_braun_ever, family = 'cox')

	RESP_oscar_gsva_braun_ever <- oscar::oscar(x = gsva_braun_ever, y = RESP_braun_ever, family = 'logistic')
	RESP_oscar_unigene_braun_ever <- oscar::oscar(x = unigene_braun_ever, y = RESP_braun_ever, family = 'logistic')

	# TCGA

	PFS_oscar_gsva_tcga <- oscar::oscar(x = gsva_tcga, y = PFS_tcga, family = 'cox')
	PFS_oscar_unigene_tcga <- oscar::oscar(x = unigene_tcga, y = PFS_tcga, family = 'cox')

	OS_oscar_gsva_tcga <- oscar::oscar(x = gsva_tcga, y = OS_tcga, family = 'cox')
	OS_oscar_unigene_tcga <- oscar::oscar(x = unigene_tcga, y = OS_tcga, family = 'cox')

	RESP_oscar_gsva_tcga <- oscar::oscar(x = gsva_tcga, y = RESP_tcga, family = 'logistic')
	RESP_oscar_unigene_tcga <- oscar::oscar(x = unigene_tcga, y = RESP_tcga, family = 'logistic')
}

# Run and aggregate results from each dataset

load("DREAM_generanking_ash_240621.RData")

# Aggregate coefficient paths over PFS
## -> Weighting to be added with binarizing the beta-k matrices?

## GSVA

PFS_gsva_aggregate_nivo <- list(
	PFS_oscar_gsva_braun_nivo@bperk,
	PFS_oscar_gsva_gide@bperk,
	PFS_oscar_gsva_lauss@bperk,
	PFS_oscar_gsva_prat@bperk,
	PFS_oscar_gsva_westin@bperk
)

PFS_gsva_aggregate_chem <- list(
	PFS_oscar_gsva_braun_ever@bperk,
	PFS_oscar_gsva_tcga@bperk	
)


## unigene

PFS_unigene_aggregate_nivo <- list(
	PFS_oscar_unigene_braun_nivo@bperk,
	PFS_oscar_unigene_gide@bperk,
	PFS_oscar_unigene_lauss@bperk,
	PFS_oscar_unigene_prat@bperk,
	PFS_oscar_unigene_westin@bperk
)

PFS_unigene_aggregate_chem <- list(
	PFS_oscar_unigene_braun_ever@bperk,
	PFS_oscar_unigene_tcga@bperk	
)


# Aggregate coefficient paths over OS

## GSVA

OS_gsva_aggregate_nivo <- list(
	OS_oscar_gsva_braun_nivo@bperk,
	OS_oscar_gsva_gide@bperk,
	OS_oscar_gsva_hugo@bperk
)

OS_gsva_aggregate_chem <- list(
	OS_oscar_gsva_braun_ever@bperk,
	OS_oscar_gsva_tcga@bperk
)

## unigene

OS_unigene_aggregate_nivo <- list(
	OS_oscar_unigene_braun_nivo@bperk,
	OS_oscar_unigene_gide@bperk,
	OS_oscar_unigene_hugo@bperk
)

OS_unigene_aggregate_chem <- list(
	OS_oscar_unigene_braun_ever@bperk,
	OS_oscar_unigene_tcga@bperk
)


# Aggregate coefficient paths over RESP

## GSVA

RESP_gsva_aggregate_nivo <- list(
	RESP_oscar_gsva_braun_nivo@bperk,
	RESP_oscar_gsva_chen@bperk,
	RESP_oscar_gsva_gide@bperk,
	RESP_oscar_gsva_hugo@bperk,
	RESP_oscar_gsva_kim@bperk,
	RESP_oscar_gsva_lauss@bperk,
	RESP_oscar_gsva_prat@bperk,
	RESP_oscar_gsva_riaz@bperk
)
	
RESP_gsva_aggregate_chem <- list(
	RESP_oscar_gsva_braun_ever@bperk,
	RESP_oscar_gsva_tcga@bperk
)

## unigene

RESP_unigene_aggregate_nivo <- list(
	RESP_oscar_unigene_braun_nivo@bperk,
	RESP_oscar_unigene_chen@bperk,
	RESP_oscar_unigene_gide@bperk,
	RESP_oscar_unigene_hugo@bperk,
	RESP_oscar_unigene_kim@bperk,
	RESP_oscar_unigene_lauss@bperk,
	RESP_oscar_unigene_prat@bperk,
	RESP_oscar_unigene_riaz@bperk
)
	
RESP_unigene_aggregate_chem <- list(
	RESP_oscar_unigene_braun_ever@bperk,
	RESP_oscar_unigene_tcga@bperk
)

# Order of oscar feature selection paths:

# Nivo-arms

## gsva - ELIMINATION ORDER; inverse is importance

PFS_gsva_aggregate_nivo
# Braun nivo: TIGIT, CXCR6, PDCD1, CXCL9, CD274
## Braun et al. not really representative of NSCLC
# Gide: CXCR6, CD274, TIGIT, PDCD1, CXCL9
# Lauss: CXCR6, TIGIT, PDCD1, CD274, CXCL9
# Prat: CXCL9, CXCR6, TIGIT, PDCD1, CD274
# Westin: CXCL9, CXCR6, TIGIT, PDCD1, CD274
## Interestingly, gsva_all gets eliminated quite late (keeping C274 in?) - not reported here

OS_gsva_aggregate_nivo
# Braun nivo: TIGIT, CXCL9, CXCR6, CD274, PDCD1
## Braun et al. not really representative of NSCLC
# Gide: CD274, CXCR9, TIGIT, CXCL6, PDCD1
# Hugo: PDCD1, CXCL9, CD274, CXCR6/TIGIT

RESP_gsva_aggregate_nivo
# Braun nivo: Not converging
## Not very representative
# Chen: CXCR6, CD274, CXCL9/TIGIT, PDCD1
# Gide: CXCR6, CD274, PDCD1, CXCL9, TIGIT
# Hugo: PDCD1, TIGIT, CXCL9/CD274, CXCR6
# Kim: PDCD1, CXCL9/CXCR6, PDCD1/TIGIT, CD274
# Lauss: CXCR6, CXCL9/TIGIT, CD274, PDCD1
# Prat: PDCD1, CXCL9, TIGIT, CD274, CXCR6
# Riaz: PDCD1, CXCL9, CD274, CXCR6/TIGIT

## unigene - SELECTION ORDER; most important chosen first

PFS_unigene_aggregate_nivo
# Braun: CXCR6, CD274, CXCL9, TIGIT, PDCD1
## Not very representative
# Gide: CXCR6, CXCL9, CD274, PDCD1, TIGIT
# Lauss: CD274, TIGIT, CXCR6, PDCD1, CXCL9
# Prat: PDCD1, TIGIT, CXCL9, CXCR6, CD274
# Westin: CXCL9, CD274, TIGIT, CXCR6, PDCD1 

OS_unigene_aggregate_nivo
# Braun: PDCD1, CD274, CXCR6, TIGIT, CXCL9
## Not very representative
# Gide: CD274, CXCR6, PDCD1, CXCL9, TIGIT
# Hugo: CXCR6, TIGIT, CXCL9, PDCD1, CD274

RESP_unigene_aggregate_nivo
# Braun: Did not converge
# Chen: PDCD1, CD274, CXCR6, CXCL9
# Gide: CXCR6, CD274, TIGIT, CXCL9, PDCD1
# Hugo: CXCR6, TIGIT, PDCD1, CXCL9, CD274
# Kim: PDCD1, CD274, CXCL9, CXCR6, TIGIT
# Lauss: CD274, PDCD1, CXCR6, TIGIT, CXCL9
# Prat: CD274, PDCD1, TIGIT, CXCR6, CXCL9
# Riaz: TIGIT, CXCR9, CXCL9, CD274, PDCD1

# Chemo-arms

## gsva
PFS_gsva_aggregate_chem
# Braun ever: CXCL9, PDCD1, CD274, TIGIT, CXCR6
# TCGA: CXCR6, TIGIT, PDCD1, CXCL9, CD274

OS_gsva_aggregate_chem
# Braun ever: PDCD1, TIGIT, CD274/CXCL9, CXCR6
# TCGA: PDCD1, TIGIT, CD274/CXCL9, CXCR6

RESP_gsva_aggregate_chem
# Braun ever: Did not converge
# TCGA: Did not converge

## unigene
PFS_unigene_aggregate_chem
# Braun ever: TIGIT, CXCL9, CD274, CXCR6, PDCD1
# TCGA: CD274, CXCR6, CXCL9, PDCD1, TIGIT

OS_unigene_aggregate_chem
# Braun ever: CXCR6, CXCL9, PDCD1, TIGIT, CD274
# TCGA: CXCR6, CXCL9, PDCD1, TIGIT, CD274

RESP_unigene_aggregate_chem
# Braun ever: CD274, CXCL9, PDCD1, CXCR6, TIGIT
# TCGA: CD274, CXCL9, PDCD1, TIGIT, CXCR6

##
## Conflicting and very difficult to interpret
## Some consisting trends though... while considering also the chemo arms
##

## CXCR6 seems like a very strong candidate, for predictivity even better than CD274
## CXCL9 seems about at par with CD274, but is quite close with CXCR6
## PDCD1 and TIGIT clearly show the least signal, with a bit of dataset generalizability to order them as 4th or 5th

## With such a small gene set, priority is to be placed on the unigene feature selection rather than GSVA estimates (possibly misleading due to small gene set size etc)

## Best performing curation (OS) emphasis requested by organizers, in benchmarking by other teams in the challenge
## Importance curated from unigene results, weighting the most representative datasets, favoring nivo datasets over chemo arms in such a selective subset: 
## CXCR6 > CD274 >  CXCL9 > PDCD1 > TIGIT

ranks <- read.csv("gene_ranks_FICAN-OSCAR.csv")
#> sum(ranks$rank, na.rm=TRUE)
#[1] 0
ranks[match(c("CXCR6", "CD274", "CXCL9", "PDCD1", "TIGIT"), ranks$gene),"rank"] <- 1:5
#> sum(ranks$rank, na.rm=TRUE)
#[1] 15
#> sum(1:5)
#[1] 15
#> sum(is.na(ranks$rank))
#[1] 29182

## Write out the gene rankings in same format as the input file
write.csv(ranks, file="gene_ranks_FICAN-OSCAR.csv", quote=FALSE, row.names=FALSE)

