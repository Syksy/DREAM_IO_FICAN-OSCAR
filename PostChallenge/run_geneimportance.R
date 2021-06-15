####
#
# Test gene importances by omitting them one-by-one from GSVA gene sets and then re-running predictions for them
# Output is ranking of genes, added to the CSV file provided by organizers
#
###

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
gsva_prat <- list(
	gsva_prat_all = t(GSVA::gsva(as.matrix(X_prat), gmt_all, verbose=FALSE))[,1],	
	gsva_prat_noCD274 = t(GSVA::gsva(as.matrix(X_prat), gmt_noCD274, verbose=FALSE))[,1],	
	gsva_prat_noCXCL9 = t(GSVA::gsva(as.matrix(X_prat), gmt_noCXCL9, verbose=FALSE))[,1],	
	gsva_prat_noCXCR6 = t(GSVA::gsva(as.matrix(X_prat), gmt_noCXCR6, verbose=FALSE))[,1],	
	gsva_prat_noPDCD1 = t(GSVA::gsva(as.matrix(X_prat), gmt_noPDCD1, verbose=FALSE))[,1],	
	gsva_prat_noTIGIT = t(GSVA::gsva(as.matrix(X_prat), gmt_noTIGIT, verbose=FALSE))[,1]	
)
gsva_gide <- list(
	gsva_gide_all = t(GSVA::gsva(as.matrix(X_gide), gmt_all, verbose=FALSE))[,1],	
	gsva_gide_noCD274 = t(GSVA::gsva(as.matrix(X_gide), gmt_noCD274, verbose=FALSE))[,1],
	gsva_gide_noCXCL9 = t(GSVA::gsva(as.matrix(X_gide), gmt_noCXCL9, verbose=FALSE))[,1],	
	gsva_gide_noCXCR6 = t(GSVA::gsva(as.matrix(X_gide), gmt_noCXCR6, verbose=FALSE))[,1],	
	gsva_gide_noPDCD1 = t(GSVA::gsva(as.matrix(X_gide), gmt_noPDCD1, verbose=FALSE))[,1],	
	gsva_gide_noTIGIT = t(GSVA::gsva(as.matrix(X_gide), gmt_noTIGIT, verbose=FALSE))[,1]	
)
gsva_hugo <- list(
	gsva_hugo_all = t(GSVA::gsva(as.matrix(X_hugo), gmt_all, verbose=FALSE))[,1],	
	gsva_hugo_noCD274 = t(GSVA::gsva(as.matrix(X_hugo), gmt_noCD274, verbose=FALSE))[,1],	
	gsva_hugo_noCXCL9 = t(GSVA::gsva(as.matrix(X_hugo), gmt_noCXCL9, verbose=FALSE))[,1],	
	gsva_hugo_noCXCR6 = t(GSVA::gsva(as.matrix(X_hugo), gmt_noCXCR6, verbose=FALSE))[,1],	
	gsva_hugo_noPDCD1 = t(GSVA::gsva(as.matrix(X_hugo), gmt_noPDCD1, verbose=FALSE))[,1],	
	gsva_hugo_noTIGIT = t(GSVA::gsva(as.matrix(X_hugo), gmt_noTIGIT, verbose=FALSE))[,1]	
)
gsva_lauss <- list(
	gsva_lauss_all = t(GSVA::gsva(as.matrix(X_lauss), gmt_all, verbose=FALSE))[,1],	
	gsva_lauss_noCD274 = t(GSVA::gsva(as.matrix(X_lauss), gmt_noCD274, verbose=FALSE))[,1],	
	gsva_lauss_noCXCL9 = t(GSVA::gsva(as.matrix(X_lauss), gmt_noCXCL9, verbose=FALSE))[,1],	
	gsva_lauss_noCXCR6 = t(GSVA::gsva(as.matrix(X_lauss), gmt_noCXCR6, verbose=FALSE))[,1],	
	gsva_lauss_noPDCD1 = t(GSVA::gsva(as.matrix(X_lauss), gmt_noPDCD1, verbose=FALSE))[,1],	
	gsva_lauss_noTIGIT = t(GSVA::gsva(as.matrix(X_lauss), gmt_noTIGIT, verbose=FALSE))[,1]	
)
gsva_kim <- list(
	gsva_kim_all = t(GSVA::gsva(as.matrix(X_kim), gmt_all, verbose=FALSE))[,1],	
	gsva_kim_noCD274 = t(GSVA::gsva(as.matrix(X_kim), gmt_noCD274, verbose=FALSE))[,1],	
	gsva_kim_noCXCL9 = t(GSVA::gsva(as.matrix(X_kim), gmt_noCXCL9, verbose=FALSE))[,1],	
	gsva_kim_noCXCR6 = t(GSVA::gsva(as.matrix(X_kim), gmt_noCXCR6, verbose=FALSE))[,1],	
	gsva_kim_noPDCD1 = t(GSVA::gsva(as.matrix(X_kim), gmt_noPDCD1, verbose=FALSE))[,1],	
	gsva_kim_noTIGIT = t(GSVA::gsva(as.matrix(X_kim), gmt_noTIGIT, verbose=FALSE))[,1]	
)
gsva_chen <- list(
	gsva_chen_all = t(GSVA::gsva(as.matrix(X_chen), gmt_all, verbose=FALSE))[,1],	
	gsva_chen_noCD274 = t(GSVA::gsva(as.matrix(X_chen), gmt_noCD274, verbose=FALSE))[,1],	
	gsva_chen_noCXCL9 = t(GSVA::gsva(as.matrix(X_chen), gmt_noCXCL9, verbose=FALSE))[,1],	
	gsva_chen_noCXCR6 = t(GSVA::gsva(as.matrix(X_chen), gmt_noCXCR6, verbose=FALSE))[,1],	
	gsva_chen_noPDCD1 = t(GSVA::gsva(as.matrix(X_chen), gmt_noPDCD1, verbose=FALSE))[,1],	
	gsva_chen_noTIGIT = t(GSVA::gsva(as.matrix(X_chen), gmt_noTIGIT, verbose=FALSE))[,1]	
)
gsva_westin <- list(
	gsva_westin_all = t(GSVA::gsva(as.matrix(X_westin), gmt_all, verbose=FALSE))[,1],	
	gsva_westin_noCD274 = t(GSVA::gsva(as.matrix(X_westin), gmt_noCD274, verbose=FALSE))[,1],	
	gsva_westin_noCXCL9 = t(GSVA::gsva(as.matrix(X_westin), gmt_noCXCL9, verbose=FALSE))[,1],	
	gsva_westin_noCXCR6 = t(GSVA::gsva(as.matrix(X_westin), gmt_noCXCR6, verbose=FALSE))[,1],	
	gsva_westin_noPDCD1 = t(GSVA::gsva(as.matrix(X_westin), gmt_noPDCD1, verbose=FALSE))[,1],	
	gsva_westin_noTIGIT = t(GSVA::gsva(as.matrix(X_westin), gmt_noTIGIT, verbose=FALSE))[,1]	
)
gsva_riaz <- list(
	gsva_riaz_all = t(GSVA::gsva(as.matrix(X_riaz), gmt_all, verbose=FALSE))[,1],	
	gsva_riaz_noCD274 = t(GSVA::gsva(as.matrix(X_riaz), gmt_noCD274, verbose=FALSE))[,1],	
	gsva_riaz_noCXCL9 = t(GSVA::gsva(as.matrix(X_riaz), gmt_noCXCL9, verbose=FALSE))[,1],	
	gsva_riaz_noCXCR6 = t(GSVA::gsva(as.matrix(X_riaz), gmt_noCXCR6, verbose=FALSE))[,1],	
	gsva_riaz_noPDCD1 = t(GSVA::gsva(as.matrix(X_riaz), gmt_noPDCD1, verbose=FALSE))[,1],	
	gsva_riaz_noTIGIT = t(GSVA::gsva(as.matrix(X_riaz), gmt_noTIGIT, verbose=FALSE))[,1]	
)

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
