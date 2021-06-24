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

