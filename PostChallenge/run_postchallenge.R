###
#
# Post-challenge diagnostics etc
#
###

setwd("D:\\Gits\\DREAM_2020_IO\\")

curateX <- function(
	gex, # Gene expression matrix, with gene symbols for annotation and patient IDs for sample identification
	dat, # Combined clinical and non-GEX variables (e.g. TMB, IHC-stainings, age, ...)
	keygenes = unique(c(
		"CD274", "PDL1", # PD-L1
		"CD276", "B7-H3", "B7.H3",
		"PDCD1", "CD279",
		"EGFR",   # T790M mutation in particular, separate drugs used for squamous
		"CD246", "ALK",    # Mutation often seen in non-smoker, young, adenocarcinoma subtype
		"ROS1",   # Adenocarcinoma, often negative for ALK, KRAS and EGFR muts
		"BRAF",   # Mutation helps tumor to grow
		"RET",    # Mutation helps tumor to grow
		"MET",    # Mutation helps tumor to grow
		"NTRK",   # Mutation helps tumor to grow
		"HER2",   # NSLungCa PeerView podcast
		"KRAS",   # NSLungCa PeerView podcast
		"NRG1",    # NSLungCa PeerView podcast
		# Mikko's slides, anything related; https://www.frontiersin.org/articles/10.3389/fonc.2016.00112/full
		"VEGF",
		"HER1",
		"HER2",
		"HER3",
		"HER4",
		"PTEN",
		"PI3K",
		"RAS",
		"MEK",
		"ERK",
		"CDK4",
		"CDK6",
		"PDK-1",
		"AKT",
		# Mikko NGS syöpäpaneeli 1 lisäyksiä
		"PIK3CA",
		"KIT",
		"NRAS",
		"PDGFRA",
		# Other interesting? Cytokines, chemokines
		"CLCL10",
		"CSCL11",
		# https://www.cell.com/cancer-cell/pdfExtended/S1535-6108(19)30037-6 Table 2
		"TLR3", 
		"LAG3",
		"IDO1",
		"TIGIT", 
		"TNFAIP3",
		"ADORA2A",
		"ICOS",
		"TNFRSF9",
		"TNFRSF9",
		"CD52",		
		"BTLA",
		"TIGIT",
		"IDO1",
		"TLR8",
		# https://www.nature.com/articles/nm.3909
		"KLRB1",
		"CD161",
		"FOXM1",
		# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6487502/
		"SERPINB9",
		# Immuno marker, e.g. https://science.sciencemag.org/content/362/6411/eaar3593?rss=1#ref-15
		"CD8A",
		"CCL5",
		"CD27",
		"CMKLR1",
		"CXCL9",
		"CXCR6",
		"HLA-D1A1","HLA.D1A1",
		"HLA-DRB1","HLA.DRB1",
		"HLA-E", "HLA.E",
		"IDO1",
		"LAG3",
		"NKG7",
		"PDCD1LG2","PDL2",
		"PSMB10",
		"STAT1",
		"TIGIT",
		# https://www.cell.com/cancer-cell/pdfExtended/S1535-6108(19)30037-6
		"TLR3",
		"LAG3",
		"IDO1",
		"TIGIT",
		"TNFAIP3",
		"ADORA2A",
		"ICOS",
		"TNFRSF9",
		"CD52",
		"BTLA",
		"TLR8",
		# Call on Feb 22
		"CD80",
		"ATG5",
		"ATG7",
		"ARG1"
	)),
	normalize = FALSE, # Changed: by default, no normalization here
	normfunc = logz, # Function for normalizing gene values to be used as variables - could be e.g. z-score within sample? log-transform if normalized count data >0?
	...
){

	## Cherry-pick some relevant genes for modelling
	# Pick genes that are present, warn of non-present symbols, column bind into variables
	warning(paste("Following key gene symbols not found from gex: ", paste(keygenes[which(!keygenes %in% rownames(gex))], collapse=", ")), sep="")
	keygenes <- keygenes[which(keygenes %in% rownames(gex))]
	#X <- t(do.call("rbind", lapply(keygenes, FUN=function(z){
	#	print(paste(z, "..."))
	#	if(normalize){
	#		coln <- colnames(gex)
	#		tmp <- t(normfunc(gex[z,]))
	#		
	#	}else{
	#		gex[z,]
	#	}
	#})))
	# Rename
	#colnames(X) <- keygenes

	# Return X with newly derived variables
	t(as.matrix(gex[keygenes,]))
}
# Omit redundant columns that are either full of one single value or all NAs
omit.reducols <- function(mat){
	#mat <- mat[,-which(apply(mat, MARGIN=2, FUN=function(x) { all(x==unique(x)[1] | is.na(x)) }))]
	mat
}
# Omit any columns that include any amount of NAs
omit.nacols <- function(mat){
	#mat <- mat[,-which(apply(mat, MARGIN=2, FUN=function(x) { any(is.na(x)) }))]
	mat
}
# Do logz-transformation
logz <- function(x) { 
	tmp <- c(scale(log(x+1)))
	if(any(!is.finite(tmp))){
		tmp[!is.finite(tmp)] <- 0
	}
	tmp
}

# Load data and clinical variable matrices

load(".\\RData\\gex_tcga.RData")
load(".\\RData\\dat_tcga.RData")
X_tcga <- omit.reducols(curateX(gex=gex_tcga, dat=dat_tcga))
PFS_tcga <- survival::Surv(time = dat_tcga$PFS.time, event = dat_tcga$PFS.event)
OS_tcga <- survival::Surv(time = dat_tcga$OS.time, event = dat_tcga$OS.event)
RESP_tcga <- as.integer(dat_tcga$Responder)
# Hugo et al. (GEO)
load(".\\RData\\gex_hugo.RData")
load(".\\RData\\dat_hugo.RData")
X_hugo <- omit.reducols(curateX(gex=gex_hugo, dat=dat_hugo))
OS_hugo <- survival::Surv(time = dat_hugo$OS.time, event = dat_hugo$OS.event)
RESP_hugo <- as.integer(dat_hugo$Responder)
# Prat et al. (GEO)
load(".\\RData\\gex_prat.RData")
load(".\\RData\\dat_prat.RData")
X_prat <- omit.reducols(curateX(gex=gex_prat, dat=dat_prat, normalize=FALSE))
PFS_prat <- survival::Surv(time = dat_prat$PFS.time, event = dat_prat$PFS.event)
RESP_prat <- as.integer(dat_prat$Responder)
# Westin et al. (GEO)
load(".\\RData\\gex_westin.RData")
load(".\\RData\\dat_westin.RData")
X_westin <- omit.reducols(curateX(gex=gex_westin, dat=dat_westin))
PFS_westin <- survival::Surv(time = dat_westin$PFS.time, event = dat_westin$PFS.event)
# Riaz et al. (GEO)
load(".\\RData\\gex_riaz.RData")
load(".\\RData\\dat_riaz.RData")
X_riaz <- omit.reducols(curateX(gex=gex_riaz, dat=dat_riaz))
RESP_riaz <- dat_riaz[,"Responder"]
# Lauss et al. (TIDE)
load(".\\RData\\gex_lauss.RData")
load(".\\RData\\dat_lauss.RData")
X_lauss <- omit.reducols(curateX(gex=gex_lauss, dat=dat_lauss, normalize=FALSE))
PFS_lauss <- survival::Surv(time=dat_lauss[,"PFS.time"], event=dat_lauss[,"PFS.event"])
OS_lauss <- survival::Surv(time=dat_lauss[,"OS.time"], event=dat_lauss[,"OS.event"])
RESP_lauss <- dat_lauss[,"Responder"]
# Kim et al. (TIDE)
load(".\\RData\\gex_kim.RData")
load(".\\RData\\dat_kim.RData")
X_kim <- omit.reducols(curateX(gex=gex_kim, dat=dat_kim, normalize=FALSE))
RESP_kim <- dat_kim[,"Responder"]
# Chen et al. (TIDE)
load(".\\RData\\gex_chen.RData")
load(".\\RData\\dat_chen.RData")
X_chen <- omit.reducols(curateX(gex=gex_chen, dat=dat_chen, normalize=FALSE))
RESP_chen <- dat_chen[,"Responder"]
# Gide et al. (TIDE)
load(".\\RData\\gex_gide.RData")
load(".\\RData\\dat_gide.RData")
X_gide <- omit.reducols(curateX(gex=gex_gide, dat=dat_gide))
PFS_gide <- survival::Surv(time=dat_gide[,"PFS.time"], event=dat_gide[,"PFS.event"])
OS_gide <- survival::Surv(time=dat_gide[,"OS.time"], event=dat_gide[,"OS.event"])
RESP_gide <- dat_gide[,"Responder"]
# Braun et al. (Raw source, both Nivo and Chemo arms)
# NOTE! Possible differences in how TMB behaves:
# https://jitc.bmj.com/content/8/1/e000319
#
# "TMB scores ranged from 0.36 to 12.24 mutations/Mb (mean 2.83 mutations/Mb) with no significant difference between 
# the PD and DC groups (3.01 vs 2.63 mutations/Mb, respectively; p=0.7682)."
#
# "Overall, neither TMB nor PD-L1 correlated with ICI response and TMB was not significantly associated with 
# PD-L1 expression. The higher incidence of LOH-MHC in PD group suggests that loss of antigen presentation 
# may restrict response to ICIs. Separately, enrichment of HRR gene mutations in the DC group suggests potential 
# utility in predicting ICI response and a potential therapeutic target, warranting future studies."
load(".\\RData\\gex_braun_nivo.RData")
load(".\\RData\\dat_braun_nivo.RData")
load(".\\RData\\gex_braun_ever.RData")
load(".\\RData\\dat_braun_ever.RData")
X_braun_nivo <- omit.reducols(curateX(gex=gex_braun_nivo, dat=dat_braun_nivo))
X_braun_ever <- omit.reducols(curateX(gex=gex_braun_ever, dat=dat_braun_ever))
PFS_braun_nivo <- survival::Surv(time=dat_braun_nivo[,"PFS.time"], event=dat_braun_nivo[,"PFS.event"])
iPFS_braun_nivo <- survival::Surv(time=dat_braun_nivo[,"iPFS.time"], event=dat_braun_nivo[,"iPFS.event"])
OS_braun_nivo <- survival::Surv(time=dat_braun_nivo[,"OS.time"], event=dat_braun_nivo[,"OS.event"])
RESP_braun_nivo <- dat_braun_nivo[,"Responder"]
iRESP_braun_nivo <- dat_braun_nivo[,"iResponder"]
PFS_braun_ever <- survival::Surv(time=dat_braun_ever[,"PFS.time"], event=dat_braun_ever[,"PFS.event"])
iPFS_braun_ever <- survival::Surv(time=dat_braun_ever[,"iPFS.time"], event=dat_braun_ever[,"iPFS.event"])
OS_braun_ever <- survival::Surv(time=dat_braun_ever[,"OS.time"], event=dat_braun_ever[,"OS.event"])
RESP_braun_ever <- dat_braun_ever[,"Responder"]
iRESP_braun_ever <- dat_braun_ever[,"iResponder"]



## OSCAR main model fits



setwd("D:\\Gits\\DREAM_2020_IO\\PostChallenge")



# Modeling based on genes only - post-challenge diagnostics etc

# Anti-PD1 arms

OS_hugo_oscar <- oscar::oscar(x = X_hugo, y = OS_hugo, family="cox")
RESP_hugo_oscar <- oscar::oscar(x = X_hugo, y = RESP_hugo, family="logistic")

PFS_prat_oscar <- oscar::oscar(x = X_prat, y = PFS_prat, family="cox")
RESP_prat_oscar <- oscar::oscar(x = X_prat, y = RESP_prat, family="logistic")

PFS_westin_oscar <- oscar::oscar(x = X_westin, y = PFS_westin, family="cox")

RESP_riaz_oscar <- oscar::oscar(x = X_riaz, y = RESP_riaz, family="logistic")

PFS_lauss_oscar <- oscar::oscar(x = X_lauss, y = PFS_lauss, family="cox")
OS_lauss_oscar <- oscar::oscar(x = X_lauss, y = OS_lauss, family="cox")
RESP_lauss_oscar <- oscar::oscar(x = X_lauss, y = RESP_lauss, family="logistic")

PFS_gide_oscar <- oscar::oscar(x = X_gide, y = PFS_gide, family="cox")
OS_gide_oscar <- oscar::oscar(x = X_gide, y = OS_gide, family="cox")
RESP_gide_oscar <- oscar::oscar(x = X_gide, y = RESP_gide, family="logistic")

PFS_braun_nivo_oscar <- oscar::oscar(x = X_braun_nivo, y = PFS_braun_nivo, family="cox")
OS_braun_nivo_oscar <- oscar::oscar(x = X_braun_nivo, y = OS_braun_nivo, family="cox")
RESP_braun_nivo_oscar <- oscar::oscar(x = X_braun_nivo, y = RESP_braun_nivo, family="logistic")

RESP_kim_oscar <- oscar::oscar(x = X_kim, y = RESP_kim, family="logistic")

RESP_chen_oscar <- oscar::oscar(x = X_chen, y = RESP_chen, family="logistic")

save.image("fits.RData")

# Chemo control arms

PFS_tcga_oscar <- oscar::oscar(x = X_tcga, y = PFS_tcga, family="cox", verb=1)
OS_tcga_oscar <- oscar::oscar(x = X_tcga, y = OS_tcga, family="cox", verb=1)
RESP_tcga_oscar <- oscar::oscar(x = X_tcga, y = RESP_tcga, family="logistic", verb=1)

PFS_braun_ever_oscar <- oscar::oscar(x = X_braun_ever, y = PFS_braun_ever, family="cox")
OS_braun_ever_oscar <- oscar::oscar(x = X_braun_ever, y = OS_braun_ever, family="cox")
RESP_braun_ever_oscar <- oscar::oscar(x = X_braun_ever, y = RESP_braun_ever, family="logistic")

save.image("fits.RData")


load("fits.RData")

## OSCAR model fitting based on bootstrapped datasets

OS_hugo_oscar_bs <- oscar::bs.oscar(OS_hugo_oscar, bootstrap=20, seed=123)
RESP_hugo_oscar_bs <- oscar::bs.oscar(RESP_hugo_oscar, bootstrap=20, seed=123)

save.image("fits_bs_temp.RData")

PFS_prat_oscar_bs <- oscar::bs.oscar(PFS_prat_oscar, bootstrap=20, seed=123)
RESP_prat_oscar_bs <- oscar::bs.oscar(RESP_prat_oscar, bootstrap=20, seed=123)

save.image("fits_bs_temp.RData")

PFS_westin_oscar_bs <- oscar::bs.oscar(PFS_westin_oscar, bootstrap=20, seed=123)

save.image("fits_bs_temp.RData")

RESP_riaz_oscar_bs <- oscar::bs.oscar(RESP_riaz_oscar, bootstrap=20, seed=123)

save.image("fits_bs_temp.RData")

PFS_lauss_oscar_bs <- oscar::bs.oscar(PFS_lauss_oscar, bootstrap=20, seed=123)
OS_lauss_oscar_bs <- oscar::bs.oscar(OS_lauss_oscar, bootstrap=20, seed=123)
RESP_lauss_oscar_bs <- oscar::bs.oscar(RESP_lauss_oscar, bootstrap=20, seed=123)

save.image("fits_bs_temp.RData")

PFS_gide_oscar_bs <- oscar::bs.oscar(PFS_gide_oscar, bootstrap=20, seed=123)
OS_gide_oscar_bs <- oscar::bs.oscar(OS_gide_oscar, bootstrap=20, seed=123)
RESP_gide_oscar_bs <- oscar::bs.oscar(RESP_gide_oscar, bootstrap=20, seed=123)

save.image("fits_bs_temp.RData")

PFS_braun_nivo_oscar_bs <- oscar::bs.oscar(PFS_braun_nivo_oscar, bootstrap=20, seed=123)
OS_braun_nivo_oscar_bs <- oscar::bs.oscar(OS_braun_nivo_oscar, bootstrap=20, seed=123)
RESP_braun_nivo_oscar_bs <- oscar::bs.oscar(RESP_braun_nivo_oscar, bootstrap=20, seed=123)

save.image("fits_bs_temp.RData")

RESP_kim_oscar_bs <- oscar::bs.oscar(RESP_kim_oscar, bootstrap=20, seed=123)

save.image("fits_bs_temp.RData")

RESP_chen_oscar_bs <- oscar::bs.oscar(RESP_chen_oscar, bootstrap=20, seed=123)

save.image("fits_bs.RData")


# Chemo control arms

PFS_tcga_oscar_bs <- oscar::bs.oscar(PFS_tcga_oscar, bootstrap=20, seed=123)
OS_tcga_oscar_bs <- oscar::bs.oscar(OS_tcga_oscar, bootstrap=20, seed=123)
RESP_tcga_oscar_bs <- oscar::bs.oscar(RESP_tcga_oscar, bootstrap=20, seed=123)

save.image("fits_bs_temp.RData")

PFS_braun_ever_oscar_bs <- oscar::bs.oscar(PFS_braun_ever_oscar, bootstrap=20, seed=123)
OS_braun_ever_oscar_bs <- oscar::bs.oscar(OS_braun_ever_oscar, bootstrap=20, seed=123)
RESP_braun_ever_oscar_bs <- oscar::bs.oscar(RESP_braun_ever_oscar, bootstrap=20, seed=123)

save.image("fits_bs.RData")


### Cross-validations, 10-fold


load("fits_bs.RData")

## OSCAR model fitting based on bootstrapped datasets

OS_hugo_oscar_cv_seed123_fold10 <- oscar::cv.oscar(OS_hugo_oscar, fold=10, seed=123)
RESP_hugo_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_hugo_oscar, fold=10, seed=123)

save.image("fits_bs_cv_temp.RData")

PFS_prat_oscar_cv_seed123_fold10 <- oscar::cv.oscar(PFS_prat_oscar, fold=10, seed=123)
RESP_prat_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_prat_oscar, fold=10, seed=123)

save.image("fits_bs_cv_temp.RData")

PFS_westin_oscar_cv_seed123_fold10 <- oscar::cv.oscar(PFS_westin_oscar, fold=10, seed=123)

save.image("fits_bs_cv_temp.RData")

RESP_riaz_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_riaz_oscar, fold=10, seed=123)

save.image("fits_bs_cv_temp.RData")

PFS_lauss_oscar_cv_seed123_fold10 <- oscar::cv.oscar(PFS_lauss_oscar, fold=10, seed=123)
OS_lauss_oscar_cv_seed123_fold10 <- oscar::cv.oscar(OS_lauss_oscar, fold=10, seed=123)
RESP_lauss_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_lauss_oscar, fold=10, seed=123)

save.image("fits_bs_cv_temp.RData")

PFS_gide_oscar_cv_seed123_fold10 <- oscar::cv.oscar(PFS_gide_oscar, fold=10, seed=123)
OS_gide_oscar_cv_seed123_fold10 <- oscar::cv.oscar(OS_gide_oscar, fold=10, seed=123)
RESP_gide_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_gide_oscar, fold=10, seed=123)

save.image("fits_bs_cv_temp.RData")

PFS_braun_nivo_oscar_cv_seed123_fold10 <- oscar::cv.oscar(PFS_braun_nivo_oscar, fold=10, seed=123)
OS_braun_nivo_oscar_cv_seed123_fold10 <- oscar::cv.oscar(OS_braun_nivo_oscar, fold=10, seed=123)
RESP_braun_nivo_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_braun_nivo_oscar, fold=10, seed=123)

save.image("fits_bs_cv_temp.RData")

RESP_kim_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_kim_oscar, fold=10, seed=123)

save.image("fits_bs_cv_temp.RData")

RESP_chen_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_chen_oscar, fold=10, seed=123)

save.image("fits_bs_cv.RData")


# Chemo control arms

PFS_tcga_oscar_cv_seed123_fold10 <- oscar::cv.oscar(PFS_tcga_oscar, fold=10, seed=123)
OS_tcga_oscar_cv_seed123_fold10 <- oscar::cv.oscar(OS_tcga_oscar, fold=10, seed=123)
RESP_tcga_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_tcga_oscar, fold=10, seed=123)

save.image("fits_bs_cv_temp.RData")

PFS_braun_ever_oscar_cv_seed123_fold10 <- oscar::cv.oscar(PFS_braun_ever_oscar, fold=10, seed=123)
OS_braun_ever_oscar_cv_seed123_fold10 <- oscar::cv.oscar(OS_braun_ever_oscar, fold=10, seed=123)
RESP_braun_ever_oscar_cv_seed123_fold10 <- oscar::cv.oscar(RESP_braun_ever_oscar, fold=10, seed=123)

save.image("fits_bs_cv.RData")





### Cross-validations, 5-fold


load("fits_bs.RData")

## OSCAR model fitting based on bootstrapped datasets

OS_hugo_oscar_cv_seed234_fold5 <- oscar::cv.oscar(OS_hugo_oscar, fold=5, seed=234)
RESP_hugo_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_hugo_oscar, fold=5, seed=234)

save.image("fits_bs_cv_temp.RData")

PFS_prat_oscar_cv_seed234_fold5 <- oscar::cv.oscar(PFS_prat_oscar, fold=5, seed=234)
RESP_prat_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_prat_oscar, fold=5, seed=234)

save.image("fits_bs_cv_temp.RData")

PFS_westin_oscar_cv_seed234_fold5 <- oscar::cv.oscar(PFS_westin_oscar, fold=5, seed=234)

save.image("fits_bs_cv_temp.RData")

RESP_riaz_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_riaz_oscar, fold=5, seed=234)

save.image("fits_bs_cv_temp.RData")

PFS_lauss_oscar_cv_seed234_fold5 <- oscar::cv.oscar(PFS_lauss_oscar, fold=5, seed=234)
OS_lauss_oscar_cv_seed234_fold5 <- oscar::cv.oscar(OS_lauss_oscar, fold=5, seed=234)
RESP_lauss_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_lauss_oscar, fold=5, seed=234)

save.image("fits_bs_cv_temp.RData")

PFS_gide_oscar_cv_seed234_fold5 <- oscar::cv.oscar(PFS_gide_oscar, fold=5, seed=234)
OS_gide_oscar_cv_seed234_fold5 <- oscar::cv.oscar(OS_gide_oscar, fold=5, seed=234)
RESP_gide_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_gide_oscar, fold=5, seed=234)

save.image("fits_bs_cv_temp.RData")

PFS_braun_nivo_oscar_cv_seed234_fold5 <- oscar::cv.oscar(PFS_braun_nivo_oscar, fold=5, seed=234)
OS_braun_nivo_oscar_cv_seed234_fold5 <- oscar::cv.oscar(OS_braun_nivo_oscar, fold=5, seed=234)
RESP_braun_nivo_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_braun_nivo_oscar, fold=5, seed=234)

save.image("fits_bs_cv_temp.RData")

RESP_kim_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_kim_oscar, fold=5, seed=234)

save.image("fits_bs_cv_temp.RData")

RESP_chen_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_chen_oscar, fold=5, seed=234)

save.image("fits_bs_cv.RData")


# Chemo control arms

PFS_tcga_oscar_cv_seed234_fold5 <- oscar::cv.oscar(PFS_tcga_oscar, fold=5, seed=234)
OS_tcga_oscar_cv_seed234_fold5 <- oscar::cv.oscar(OS_tcga_oscar, fold=5, seed=234)
RESP_tcga_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_tcga_oscar, fold=5, seed=234)

save.image("fits_bs_cv_temp.RData")

PFS_braun_ever_oscar_cv_seed234_fold5 <- oscar::cv.oscar(PFS_braun_ever_oscar, fold=5, seed=234)
OS_braun_ever_oscar_cv_seed234_fold5 <- oscar::cv.oscar(OS_braun_ever_oscar, fold=5, seed=234)
RESP_braun_ever_oscar_cv_seed234_fold5 <- oscar::cv.oscar(RESP_braun_ever_oscar, fold=5, seed=234)

save.image("fits_bs_cv.RData")


