###
#
# Load readily curated datasets into workspace and continue working (i.e. modelling)
#
###

# Load temporary files (not stored in git due to size)
#load("temp.RData")

# For a given dataset, generate the finalized data matrix X where:
# Columns:
# - Selected genes and their expression (z-score normalized perhaps?)
# - Selected scoring indices (TIS, APM, IFN-gamma 1/2)
# - Immunedeconv immune cell compositions from xCell and MCP counter (maybe others still?)
# - Single sample GSEA's GSVA for at least HALLMARK gmt pathways
# Rows:
# - Unique sample (=typically patient) identifiers
curateX <- function(
	gex, # Gene expression matrix, with gene symbols for annotation and patient IDs for sample identification
	dat, # Combined clinical and non-GEX variables (e.g. TMB, IHC-stainings, age, ...)
	keygenes = unique(c(
		"CD274", "PDL1", # PD-L1
		"CD276", "B7-H3",
		"PDCD1", "CD279",
		"EGFR",   # T790M mutation in particular, separate drugs used for squamous
		"ALK",    # Mutation often seen in non-smoker, young, adenocarcinoma subtype
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
		"TIGIT"
	)),
	normfunc = function(input) { input }, # Function for normalizing gene values to be used as variables - could be e.g. z-score within sample? log-transform if normalized count data >0?
	gmts = c(1,2,4,5), # Hallmarks, oncology, custom self-made GMTs, filtered curated pathways from e.g. KEGG
	idcs = 1,
	clinvars = c("Age", "Smoking", "ECOG", "Squamous", "TMB", "SEX", "PDL1"),
	scores = c( # Included scoring metrics (preferably IO and lung cancer related)
		#"TIS", # Tumor inflammatory score
		#"MSI", # Microsatellite instability
		#"APM", # Antigen processing / presentation machinery
		#"IFN"  # Interferon gamma signalling
	)
){
	# Format data matrix X
	X <- matrix(NA, nrow=ncol(gex), ncol=0)
	
	## Cherry-pick some relevant genes for modelling
	# Pick genes that are present, warn of non-present symbols, column bind into variables
	warning(paste("Following key gene symbols not found from gex: ", paste(keygenes[which(!keygenes %in% rownames(gex))], collapse=", ")), sep="")
	keygenes <- keygenes[which(keygenes %in% rownames(gex))]
	X <- cbind(X, t(do.call("rbind", lapply(keygenes, FUN=function(z){
		print(paste(z, "..."))
		normfunc(gex[z,])
	}))))
	# Rename
	colnames(X) <- keygenes


	## Clinical variables
	# Age, integer
	if("Age" %in% clinvars){
		print("Age")
		X <- cbind(X, Age = as.integer(dat[,"AAGE"]))
		if(all(is.na(X[,"Age"]))){
			# Redundant information for age, imputation for zeroes
			X[,"Age"] <- 0
		}else{
			# Median imputation for patients without age reported
			X[is.na(X[,"Age"]),"Age"] <- median(X[,"Age"], na.rm=TRUE) 
		}
		
	}
	# Smoking status, binarized
	if("Smoking" %in% clinvars){
		print("Smoking")
		# Smoking, factor which will be turned into binary indicators
		X <- cbind(X, isSmokerCurrent = as.integer(dat[,"TOBACUSE"] == "CURRENT"))
		X <- cbind(X, isSmokerFormer = as.integer(dat[,"TOBACUSE"] == "FORMER"))
		X <- cbind(X, isSmokerNever = as.integer(dat[,"TOBACUSE"] == "NEVER"))
		# other case UNKNOWN omitted; basically NA
	}
	# ECOG status {0,1,2,3,4}, binarized
	if("ECOG" %in% clinvars){
		print("ECOG")
		X <- cbind(X, isECOG0 = as.integer(ifelse(is.na(dat[,"ECOGPS"]), 0, ifelse(dat[,"ECOGPS"] == 0, 1, 0))))
		X <- cbind(X, isECOG1 = as.integer(ifelse(is.na(dat[,"ECOGPS"]), 0, ifelse(dat[,"ECOGPS"] == 1, 1, 0))))
		X <- cbind(X, isECOG2 = as.integer(ifelse(is.na(dat[,"ECOGPS"]), 0, ifelse(dat[,"ECOGPS"] == 2, 1, 0))))
		X <- cbind(X, isECOG3 = as.integer(ifelse(is.na(dat[,"ECOGPS"]), 0, ifelse(dat[,"ECOGPS"] == 3, 1, 0))))
		X <- cbind(X, isECOG4 = as.integer(ifelse(is.na(dat[,"ECOGPS"]), 0, ifelse(dat[,"ECOGPS"] == 4, 1, 0))))
		X <- cbind(X, isECOG2orAbove = as.integer(ifelse(is.na(dat[,"ECOGPS"]), 0, ifelse(dat[,"ECOGPS"] >= 2, 1, 0))))
	}
	# If Squamous histologic subtype
	if("Squamous" %in% clinvars){
		print("Squamous")
		X <- cbind(X, isSquamous = as.integer(dat[,"CRFHIST"] == "SQUAMOUS"))
		# 0s, NAs and 1s -> code into -1, 0, +1
		if(NA %in% X[,"isSquamous"] & !all(is.na(X[,"isSquamous"]))){
			X[X[,"isSquamous"]==0,"isSquamous"] <- -1
			X[is.na(X[,"isSquamous"]),"isSquamous"] <- 0
		}
	}
	# TMB, coding into lower (-1), NA (0) or higher (+1)
	if("TMB" %in% clinvars){
		print("TMB")
		medianTMB <- median(dat[,"TMB"], na.rm=TRUE)
		# No finite TMB threshold could be determined
		if(!is.finite(medianTMB)){
			X <- cbind(X, TMB = 0)
		}else{
			X <- cbind(X, TMB = ifelse(is.na(dat[,"TMB"]), 0, ifelse(dat[,"TMB"]>=medianTMB, 1, -1)))
		}
	}
	# Sex, binary indicator for +1 that patient was male
	if("SEX" %in% clinvars){
		print("SEX")
		X <- cbind(X, isMale = as.integer(dat[,"SEX"] == "M"))
	}
	# PDL1, coding into lower (-1), NA (0) or higher (+1)
	if("PDL1" %in% clinvars){
		print("PDL1")
		medianPDL1 <- median(dat[,"PDL1"], na.rm=TRUE)
		# No finite TMB threshold could be determined
		if(!is.finite(medianPDL1)){
			X <- cbind(X, PDL1 = 0)
		}else{
			X <- cbind(X, PDL1 = ifelse(is.na(dat[,"PDL1"]), 0, ifelse(dat[,"PDL1"]>=medianPDL1, 1, -1)))
		}
	}
	
	# If PDL is >=5% expression; one threshold chosen from literature
	# Below <5% results in -1
	# NA values result in 0
	# Above >=5% results in +1
	if("PDL1" %in% clinvars){
		print("PDL1 5% threshold")
		if("PDL1" %in% colnames(dat)){
			if(any(is.finite(dat[,"PDL1"]))){
				X <- cbind(X, PDL1above5perc = ifelse(is.na(dat[,"PDL1"]), 0, ifelse(dat[,"PDL1"]>=5, 1, -1)))
			}
		}
	}

	# All scores prior to these are called "BASE" metrics; i.e. very basic patient characteristics, cherry-picked individual genes etc
	
	try({colnames(X) <- paste("BASE_", colnames(X), sep="") })
	

	## Custom gene-based scoring metrics
	# For starters, simply compute mean of expression as the "score" in all scoring methods (they were all very similar)
	# Later perhaps refine to rank-based, z-scored means or so?
	IFNGscore1 <- function(gex){
		genes <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")
		genes <- genes[which(genes %in% rownames(gex))]
		apply(gex[genes,], MARGIN=2, FUN=mean)
		
	}
	IFNGscore2 <- function(gex){
		genes <- c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB")
		genes <- genes[which(genes %in% rownames(gex))]
		apply(gex[genes,], MARGIN=2, FUN=mean)
	}
	APMscore <- function(gex){
		genes <- c("B2M", "CALR", "NLRC5", "PSMB9", "PSME1", "PSME3", "RFX5", "HSP90AB1")
		genes <- genes[which(genes %in% rownames(gex))]
		apply(gex[genes,], MARGIN=2, FUN=mean)
	}
	TISscore <- function(gex){
		# Taken from Fig 1 panel d in the publ.
		# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6829827/#MOESM2
		genes <- c("CD276", "HLA-DQA1", "CD274", "IDO1", "HLA-DRB1", "HLA-E", "CMKLR1", "PDCD1LG2", "PSMB10", "LAG3", "CXCL9", "STAT1", "CD8A", "CCL5", "NKG7", "TIGIT", "CD27", "CXCR6")
		genes <- genes[which(genes %in% rownames(gex))]
		apply(gex[genes,], MARGIN=2, FUN=mean)		
	}
	# MSI
	# ...
	
	
	
	# IFN gamma signaling
	if("IFN" %in% scores){
		print("IFN scores")
		X <- cbind(X, IFN1score = IFNGscore1(gex=gex))
		X <- cbind(X, IFN2score = IFNGscore2(gex=gex))
	}
	# Tumor inflammatory score
	if("TIS" %in% scores){
		print("TIS score")
		X <- cbind(X, TISscore = TISscore(gex=gex))
	}
	# Antigen processing / presentation machinery
	if("APM" %in% scores){
		print("APM score")
		X <- cbind(X, APMscore = APMscore(gex=gex))
	}
	# Microsallite instability
	# k-NN based prediction, 15 genes or so
	#if("MSI" %in% scores){
	#	# https://github.com/WangX-Lab/PreMSIm
	#	#if (!requireNamespace("devtools", quietly = TRUE))
	#	#    install.packages("devtools")
	#	#devtools::install_github("WangX-Lab/PreMSIm")
	#	library(PreMSIm)
	#	
	#}

	## GSVA
	library(GSVA)	
	# Hallmark gene sets
	gmt_h <- GSEABase::getGmt(".\\MSigDB\\h.all.v7.2.symbols.gmt")
	# Oncogenic
	gmt_c6 <- GSEABase::getGmt(".\\MSigDB\\c6.all.v7.2.symbols.gmt")
	# Immunogenic
	gmt_c7 <- GSEABase::getGmt(".\\MSigDB\\c7.all.v7.2.symbols.gmt")
	# Custom self-made
	gmt_custom <- GSEABase::getGmt(".\\MSigDB\\selfmade.gmt")
	# Curated GMTs, subset to interesting
	gmt_c2 <- GSEABase::getGmt(".\\MSigDB\\c2.all.v7.2.symbols.gmt")
	gmt_c2 <- gmt_c2[grep("NEUTROPH|LEUKOCYT|CD4|CD8|INTERLEUK|INFLAM|T.CELL|B.CELL|TCELL|BCELL|CHEMOK|CYTOK", unlist(lapply(gmt_c2, FUN=function(x) x@setName)))]

	# Create selected gmt-variables
	# Hallmarks
	if(1 %in% gmts){
		print("Hallmarks")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_h, verbose=FALSE))) # Hallmarks
			# Omit selected hallmarks based on a priori knowledge or educated guesses
			X <- X[,-grep("CHOLESTEROL|ESTROGEN|ANDROGEN|FATTY_ACID|OXIDATIVE|GLYCOLYSIS|REACTIVE_OXYGEN|UV_RESPONSE|BILE_ACID|ALLOGRAFT|SPERMATOGENESIS|PANCREAS", colnames(X))]
		})
	}
	# Oncogenic
	if(2 %in% gmts){
		print("Oncogenic")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_c6, verbose=FALSE))) # Oncogenic
		})
	}
	# Immunology
	if(3 %in% gmts){
		print("Immunology")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_c7, verbose=FALSE))) # Immunology
		})
	}
	# Custom self made GMTs
	if(4 %in% gmts){
		print("Selfmade GMTs")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE))) # Custom self made GMTs
		})
	}
	# Curated gene pathways
	if(5 %in% gmts){
		print("Curated GTMs from various databases")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_c2, verbose=FALSE))) # Curated GMTs
		})
	}

	## Immune deconvolution
	#install.packages("remotes")
	#remotes::install_github("icbi-lab/immunedeconv")
	library(immunedeconv)

	# xCell
	if(1 %in% idcs){
		print("xCell")
		try({
			tmp <- immunedeconv::deconvolute(gex, method="xcell")
			tmp <- as.matrix(tmp)
			rownames(tmp) <- paste("xce_", tmp[,1], sep="")
			tmp <- tmp[,-1]
			class(tmp) <- "numeric"
			X <- cbind(X, t(tmp))
		})
	}
	# MCP counter
	# NOTE: MCP counter is not functioning properly inside the cloud via Docker; it attempts to download
	#if(2 %in% idcs){
	#	print("MCP counter")
	#	try({
	#		tmp <- immunedeconv::deconvolute(gex, method="mcp_counter")
	#		X <- cbind(X, t(tmp))
	#	})
	#}

	# Return X with newly derived variables
	as.matrix(X)
}

# Omit redundant columns that are either full of one single value or all NAs
omit.reducols <- function(mat){
	mat <- mat[,-which(apply(mat, MARGIN=2, FUN=function(x) { all(x==unique(x)[1] | is.na(x)) }))]
}
# Omit any columns that include any amount of NAs
omit.nacols <- function(mat){
	mat <- mat[,-which(apply(mat, MARGIN=2, FUN=function(x) { any(is.na(x)) }))]
}

## LOAD DATA

###
#
# SYNTHETIC
#
###
gex_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
gex_synthetic <- gex_synthetic[order(rownames(gex_synthetic)),]
gex_synthetic <- as.matrix(gex_synthetic)
dat_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
X_synthetic <- curateX(gex=gex_synthetic, dat=dat_synthetic)

# Load premade GEX / DAT and generate X matrix
# Whole TCGA
load(".\\RData\\gex_tcga.RData")
load(".\\RData\\dat_tcga.RData")
X_tcga <- omit.reducols(curateX(gex=gex_tcga, dat=dat_tcga))
X_xce_tcga <- X_tcga[,grep("xce_", colnames(X_tcga))]
X_cus_tcga <- X_tcga[,grep("CUSTOM_", colnames(X_tcga))]
X_hal_tcga <- X_tcga[,grep("HALLMARK_", colnames(X_tcga))]
X_bas_tcga <- X_tcga[,grep("BASE_", colnames(X_tcga))]
PFS_tcga <- survival::Surv(time = dat_tcga$PFS.time, event = dat_tcga$PFS.event)
OS_tcga <- survival::Surv(time = dat_tcga$OS.time, event = dat_tcga$OS.event)
RESP_tcga <- as.integer(dat_tcga$Responder)
# Hugo et al. (GEO)
load(".\\RData\\gex_hugo.RData")
load(".\\RData\\dat_hugo.RData")
X_hugo <- omit.reducols(curateX(gex=gex_hugo, dat=dat_hugo))
X_xce_hugo <- X_hugo[,grep("xce_", colnames(X_hugo))]
X_cus_hugo <- X_hugo[,grep("CUSTOM_", colnames(X_hugo))]
X_hal_hugo <- X_hugo[,grep("HALLMARK_", colnames(X_hugo))]
X_bas_hugo <- X_hugo[,grep("BASE_", colnames(X_hugo))]
OS_hugo <- survival::Surv(time = dat_hugo$OS.time, event = dat_hugo$OS.event)
RESP_hugo <- as.integer(dat_hugo$Responder)
# Prat et al. (GEO)
load(".\\RData\\gex_prat.RData")
load(".\\RData\\dat_prat.RData")
X_prat <- omit.reducols(curateX(gex=gex_prat, dat=dat_prat))
X_xce_prat <- X_prat[,grep("xce_", colnames(X_prat))]
X_cus_prat <- X_prat[,grep("CUSTOM_", colnames(X_prat))]
X_hal_prat <- X_prat[,grep("HALLMARK_", colnames(X_prat))]
X_bas_prat <- X_prat[,grep("BASE_", colnames(X_prat))]
PFS_prat <- survival::Surv(time = dat_prat$PFS.time, event = dat_prat$PFS.event)
RESP_prat <- as.integer(dat_prat$Responder)
# Westin et al. (GEO)
load(".\\RData\\gex_westin.RData")
load(".\\RData\\dat_westin.RData")
X_westin <- omit.reducols(curateX(gex=gex_westin, dat=dat_westin))
X_xce_westin <- X_westin[,grep("xce_", colnames(X_westin))]
X_cus_westin <- X_westin[,grep("CUSTOM_", colnames(X_westin))]
X_hal_westin <- X_westin[,grep("HALLMARK_", colnames(X_westin))]
X_bas_westin <- X_westin[,grep("BASE_", colnames(X_westin))]
PFS_westin <- survival::Surv(time = dat_westin$PFS.time, event = dat_westin$PFS.event)
# Riaz et al. (GEO)
load(".\\RData\\gex_riaz.RData")
load(".\\RData\\dat_riaz.RData")
X_riaz <- omit.reducols(curateX(gex=gex_riaz, dat=dat_riaz))
X_xce_riaz <- X_riaz[,grep("xce_", colnames(X_riaz))]
X_cus_riaz <- X_riaz[,grep("CUSTOM_", colnames(X_riaz))]
X_hal_riaz <- X_riaz[,grep("HALLMARK_", colnames(X_riaz))]
X_bas_riaz <- X_riaz[,grep("BASE_", colnames(X_riaz))]
RESP_riaz <- dat_riaz[,"Responder"]
# Lauss et al. (TIDE)
load(".\\RData\\gex_lauss.RData")
load(".\\RData\\dat_lauss.RData")
X_lauss <- omit.reducols(curateX(gex=gex_lauss, dat=dat_lauss))
X_xce_lauss <- X_lauss[,grep("xce_", colnames(X_lauss))]
X_cus_lauss <- X_lauss[,grep("CUSTOM_", colnames(X_lauss))]
X_hal_lauss <- X_lauss[,grep("HALLMARK_", colnames(X_lauss))]
X_bas_lauss <- X_lauss[,grep("BASE_", colnames(X_lauss))]
PFS_lauss <- survival::Surv(time=dat_lauss[,"PFS.time"], event=dat_lauss[,"PFS.event"])
OS_lauss <- survival::Surv(time=dat_lauss[,"OS.time"], event=dat_lauss[,"OS.event"])
RESP_lauss <- dat_lauss[,"Responder"]
# Kim et al. (TIDE)
load(".\\RData\\gex_kim.RData")
load(".\\RData\\dat_kim.RData")
X_kim <- omit.reducols(curateX(gex=gex_kim, dat=dat_kim))
X_xce_kim <- X_kim[,grep("xce_", colnames(X_kim))]
X_cus_kim <- X_kim[,grep("CUSTOM_", colnames(X_kim))]
X_hal_kim <- X_kim[,grep("HALLMARK_", colnames(X_kim))]
X_bas_kim <- X_kim[,grep("BASE_", colnames(X_kim))]
RESP_kim <- dat_kim[,"Responder"]
# Chen et al. (TIDE)
load(".\\RData\\gex_chen.RData")
load(".\\RData\\dat_chen.RData")
X_chen <- omit.reducols(curateX(gex=gex_chen, dat=dat_chen))
X_xce_chen <- X_chen[,grep("xce_", colnames(X_chen))]
X_cus_chen <- X_chen[,grep("CUSTOM_", colnames(X_chen))]
X_hal_chen <- X_chen[,grep("HALLMARK_", colnames(X_chen))]
X_bas_chen <- X_chen[,grep("BASE_", colnames(X_chen))]
RESP_chen <- dat_chen[,"Responder"]
# Save image containing the GEXs, DATs, Xs, and various y-responses
#save.image("temp.RData")



## OSCAR + cross-validation
library(oscar)
# Seed for reproducibility
set.seed(1)

# TCGA Chemo control arm xCell
PFS_xce_tcga_oscar <- oscar::oscar(x = X_xce_tcga, y = PFS_tcga, family="cox", verb=1)
OS_xce_tcga_oscar <- oscar::oscar(x = X_xce_tcga, y = OS_tcga, family="cox", verb=1)
RESP_xce_tcga_oscar <- oscar::oscar(x = X_xce_tcga, y = RESP_tcga, family="logistic", verb=1)
# OSCAR CV
PFS_xce_cv_tcga <- oscar::cv.oscar(PFS_xce_tcga_oscar, fold=5)
OS_xce_cv_tcga <- oscar::cv.oscar(OS_xce_tcga_oscar, fold=5)
RESP_xce_cv_tcga <- oscar::cv.oscar(RESP_xce_tcga_oscar, fold=5)
# TCGA Chemo control arm CUSTOM GMTs
PFS_cus_tcga_oscar <- oscar::oscar(x = X_cus_tcga, y = PFS_tcga, family="cox", verb=1)
OS_cus_tcga_oscar <- oscar::oscar(x = X_cus_tcga, y = OS_tcga, family="cox", verb=1)
RESP_cus_tcga_oscar <- oscar::oscar(x = X_cus_tcga, y = RESP_tcga, family="logistic", verb=1)
# OSCAR CV
PFS_cus_cv_tcga <- oscar::cv.oscar(PFS_cus_tcga_oscar, fold=5)
OS_cus_cv_tcga <- oscar::cv.oscar(OS_cus_tcga_oscar, fold=5)
RESP_cus_cv_tcga <- oscar::cv.oscar(RESP_cus_tcga_oscar, fold=5)
# TCGA Chemo control arm HALLMARKS
PFS_hal_tcga_oscar <- oscar::oscar(x = X_hal_tcga, y = PFS_tcga, family="cox", verb=1)
OS_hal_tcga_oscar <- oscar::oscar(x = X_hal_tcga, y = OS_tcga, family="cox", verb=1)
RESP_hal_tcga_oscar <- oscar::oscar(x = X_hal_tcga, y = RESP_tcga, family="logistic", verb=1)
# OSCAR CV
PFS_hal_cv_tcga <- oscar::cv.oscar(PFS_hal_tcga_oscar, fold=5)
OS_hal_cv_tcga <- oscar::cv.oscar(OS_hal_tcga_oscar, fold=5)
RESP_hal_cv_tcga <- oscar::cv.oscar(RESP_hal_tcga_oscar, fold=5)
# TCGA Chemo control arm BASE MARKERS
PFS_bas_tcga_oscar <- oscar::oscar(x = X_bas_tcga, y = PFS_tcga, family="cox", verb=1)
OS_bas_tcga_oscar <- oscar::oscar(x = X_bas_tcga, y = OS_tcga, family="cox", verb=1)
RESP_bas_tcga_oscar <- oscar::oscar(x = X_bas_tcga, y = RESP_tcga, family="logistic", verb=1)
# OSCAR CV
PFS_bas_cv_tcga <- oscar::cv.oscar(PFS_bas_tcga_oscar, fold=5, seed=123)
OS_bas_cv_tcga <- oscar::cv.oscar(OS_bas_tcga_oscar, fold=5, seed=234)
RESP_bas_cv_tcga <- oscar::cv.oscar(RESP_bas_tcga_oscar, fold=5, seed=345)

# Visualize cross-validation and features as a function of k in TCGA chemo-arm
if(FALSE){
	cv.visu(PFS_xce_cv_tcga, main="PFS_xce_cv_tcga", ylab="C-index")
	cv.visu(OS_xce_cv_tcga, main="OS_xce_cv_tcga", ylab="C-index")
	cv.visu(RESP_xce_cv_tcga, main="RESP_xce_cv_tcga", ylab="Accuracy")
	
	cv.visu(PFS_cus_cv_tcga, main="PFS_cus_cv_tcga", ylab="C-index")
	cv.visu(OS_cus_cv_tcga, main="OS_cus_cv_tcga", ylab="C-index")
	cv.visu(RESP_cus_cv_tcga, main="RESP_cus_cv_tcga", ylab="Accuracy")

	cv.visu(PFS_hal_cv_tcga, main="PFS_hal_cv_tcga", ylab="C-index")
	cv.visu(OS_hal_cv_tcga, main="OS_hal_cv_tcga", ylab="C-index")
	cv.visu(RESP_hal_cv_tcga, main="RESP_hal_cv_tcga", ylab="Accuracy")

	cv.visu(PFS_bas_cv_tcga, main="PFS_bas_cv_tcga", ylab="C-index")
	cv.visu(OS_bas_cv_tcga, main="OS_bas_cv_tcga", ylab="C-index")
	cv.visu(RESP_bas_cv_tcga, main="RESP_bas_cv_tcga", ylab="Accuracy")
	
}

save.image("temprun_tcga.RData")

## GEO datasets

# Hugo et al. xCell
OS_xce_hugo_oscar <- oscar::oscar(x = X_xce_hugo, y = OS_hugo, family="cox")
RESP_xce_hugo_oscar <- oscar::oscar(x = X_xce_hugo, y = RESP_hugo, family="logistic")
# OSCAR CV
OS_xce_cv_hugo <- oscar::cv.oscar(OS_xce_hugo_oscar, fold=5, seed=1)
RESP_xce_cv_hugo <- oscar::cv.oscar(RESP_xce_hugo_oscar, fold=5, seed=2)
# Hugo et al. CUSTOM GMTs
OS_cus_hugo_oscar <- oscar::oscar(x = X_cus_hugo, y = OS_hugo, family="cox")
RESP_cus_hugo_oscar <- oscar::oscar(x = X_cus_hugo, y = RESP_hugo, family="logistic")
# OSCAR CV
OS_cus_cv_hugo <- oscar::cv.oscar(OS_cus_hugo_oscar, fold=5, seed=3)
RESP_cus_cv_hugo <- oscar::cv.oscar(RESP_cus_hugo_oscar, fold=5, seed=4)
# Hugo et al. HALLMARKS
OS_hal_hugo_oscar <- oscar::oscar(x = X_hal_hugo, y = OS_hugo, family="cox")
RESP_hal_hugo_oscar <- oscar::oscar(x = X_hal_hugo, y = RESP_hugo, family="logistic")
# OSCAR CV
OS_hal_cv_hugo <- oscar::cv.oscar(OS_hal_hugo_oscar, fold=5, seed=5)
RESP_hal_cv_hugo <- oscar::cv.oscar(RESP_hal_hugo_oscar, fold=5, seed=6)
# Hugo et al. BASE MARKERS
OS_bas_hugo_oscar <- oscar::oscar(x = X_bas_hugo, y = OS_hugo, family="cox")
RESP_bas_hugo_oscar <- oscar::oscar(x = X_bas_hugo, y = RESP_hugo, family="logistic")
# OSCAR CV
OS_bas_cv_hugo <- oscar::cv.oscar(OS_bas_hugo_oscar, fold=5, seed=7)
RESP_bas_cv_hugo <- oscar::cv.oscar(RESP_bas_hugo_oscar, fold=5, seed=8)

save.image("temprun.RData")

# Prat et al. xCell
PFS_xce_prat_oscar <- oscar::oscar(x = X_xce_prat, y = PFS_prat, family="cox")
RESP_xce_prat_oscar <- oscar::oscar(x = X_xce_prat, y = RESP_prat, family="logistic")
# OSCAR CV
PFS_xce_cv_prat <- oscar::cv.oscar(PFS_xce_prat_oscar, fold=5, seed=9)
RESP_xce_cv_prat <- oscar::cv.oscar(RESP_xce_prat_oscar, fold=5, seed=10)
# Prat et al. CUSTOM GMTs
PFS_cus_prat_oscar <- oscar::oscar(x = X_cus_prat, y = PFS_prat, family="cox")
RESP_cus_prat_oscar <- oscar::oscar(x = X_cus_prat, y = RESP_prat, family="logistic")
# OSCAR CV
PFS_cus_cv_prat <- oscar::cv.oscar(PFS_cus_prat_oscar, fold=5, seed=11)
RESP_cus_cv_prat <- oscar::cv.oscar(RESP_cus_prat_oscar, fold=5, seed=12)
# Prat et al. HALLMARKS
PFS_hal_prat_oscar <- oscar::oscar(x = X_hal_prat, y = PFS_prat, family="cox")
RESP_hal_prat_oscar <- oscar::oscar(x = X_hal_prat, y = RESP_prat, family="logistic")
# OSCAR CV
PFS_hal_cv_prat <- oscar::cv.oscar(PFS_hal_prat_oscar, fold=5, seed=13)
RESP_hal_cv_prat <- oscar::cv.oscar(RESP_hal_prat_oscar, fold=5, seed=14)
# Prat et al. BASE MARKERS
PFS_bas_prat_oscar <- oscar::oscar(x = X_bas_prat, y = PFS_prat, family="cox")
RESP_bas_prat_oscar <- oscar::oscar(x = X_bas_prat, y = RESP_prat, family="logistic")
# OSCAR CV
PFS_bas_cv_prat <- oscar::cv.oscar(PFS_bas_prat_oscar, fold=5, seed=15)
RESP_bas_cv_prat <- oscar::cv.oscar(RESP_bas_prat_oscar, fold=5, seed=16)

save.image("temprun.RData")

# Westin et al. xCell
PFS_xce_westin_oscar <- oscar::oscar(x = X_xce_westin, y = PFS_westin, family="cox")
# OSCAR CV
PFS_xce_cv_westin <- oscar::cv.oscar(PFS_xce_westin_oscar, fold=5, seed=17)
# Westin et al. CUSTOM GMTs
PFS_cus_westin_oscar <- oscar::oscar(x = X_cus_westin, y = PFS_westin, family="cox")
# OSCAR CV
PFS_cus_cv_westin <- oscar::cv.oscar(PFS_cus_westin_oscar, fold=5, seed=18)
# Westin et al. HALLMARKS
PFS_hal_westin_oscar <- oscar::oscar(x = X_hal_westin, y = PFS_westin, family="cox")
# OSCAR CV
PFS_hal_cv_westin <- oscar::cv.oscar(PFS_hal_westin_oscar, fold=5, seed=19)
# Westin et al. BASE MARKERS
PFS_bas_westin_oscar <- oscar::oscar(x = X_bas_westin, y = PFS_westin, family="cox")
# OSCAR CV
PFS_bas_cv_westin <- oscar::cv.oscar(PFS_bas_westin_oscar, fold=5, seed=20)

# Riaz et al. xCell
RESP_xce_riaz_oscar <- oscar::oscar(x = X_xce_riaz, y = RESP_riaz, family="logistic")
# OSCAR CV
RESP_xce_cv_riaz <- oscar::cv.oscar(RESP_xce_riaz_oscar, fold=5, seed=202)
# Riaz et al. CUSTOM GMTs
RESP_cus_riaz_oscar <- oscar::oscar(x = X_cus_riaz, y = RESP_riaz, family="logistic")
# OSCAR CV
RESP_cus_cv_riaz <- oscar::cv.oscar(RESP_cus_riaz_oscar, fold=5, seed=204)
# Riaz et al. HALLMARKS
RESP_hal_riaz_oscar <- oscar::oscar(x = X_hal_riaz, y = RESP_riaz, family="logistic")
# OSCAR CV
RESP_hal_cv_riaz <- oscar::cv.oscar(RESP_hal_riaz_oscar, fold=5, seed=206)
# Riaz et al. BASE MARKERS
RESP_bas_riaz_oscar <- oscar::oscar(x = X_bas_riaz, y = RESP_riaz, family="logistic")
# OSCAR CV
RESP_bas_cv_riaz <- oscar::cv.oscar(RESP_bas_riaz_oscar, fold=5, seed=208)

save.image("temprun.RData")

### TIDE datasets

# Lauss et al. xCell
PFS_xce_lauss_oscar <- oscar::oscar(x = X_xce_lauss, y = PFS_lauss, family="cox")
OS_xce_lauss_oscar <- oscar::oscar(x = X_xce_lauss, y = OS_lauss, family="cox")
RESP_xce_lauss_oscar <- oscar::oscar(x = X_xce_lauss, y = RESP_lauss, family="logistic")
# OSCAR CV
PFS_xce_cv_lauss <- oscar::cv.oscar(PFS_xce_lauss_oscar, fold=5, seed=21)
OS_xce_cv_lauss <- oscar::cv.oscar(OS_xce_lauss_oscar, fold=5, seed=22)
RESP_xce_cv_lauss <- oscar::cv.oscar(RESP_xce_lauss_oscar, fold=5, seed=23)
# Lauss et al. CUSTOM GMTs
PFS_cus_lauss_oscar <- oscar::oscar(x = X_cus_lauss, y = PFS_lauss, family="cox")
OS_cus_lauss_oscar <- oscar::oscar(x = X_cus_lauss, y = OS_lauss, family="cox")
RESP_cus_lauss_oscar <- oscar::oscar(x = X_cus_lauss, y = RESP_lauss, family="logistic")
# OSCAR CV
PFS_cus_cv_lauss <- oscar::cv.oscar(PFS_cus_lauss_oscar, fold=5, seed=24)
OS_cus_cv_lauss <- oscar::cv.oscar(OS_cus_lauss_oscar, fold=5, seed=25)
RESP_cus_cv_lauss <- oscar::cv.oscar(RESP_cus_lauss_oscar, fold=5, seed=26)
# Lauss et al. HALLMARKS
PFS_hal_lauss_oscar <- oscar::oscar(x = X_hal_lauss, y = PFS_lauss, family="cox")
OS_hal_lauss_oscar <- oscar::oscar(x = X_hal_lauss, y = OS_lauss, family="cox")
RESP_hal_lauss_oscar <- oscar::oscar(x = X_hal_lauss, y = RESP_lauss, family="logistic")
# OSCAR CV
PFS_hal_cv_lauss <- oscar::cv.oscar(PFS_hal_lauss_oscar, fold=5, seed=27)
OS_hal_cv_lauss <- oscar::cv.oscar(OS_hal_lauss_oscar, fold=5, seed=28)
RESP_hal_cv_lauss <- oscar::cv.oscar(RESP_hal_lauss_oscar, fold=5, seed=29)
# Lauss et al. BASE MARKERS
PFS_bas_lauss_oscar <- oscar::oscar(x = X_bas_lauss, y = PFS_lauss, family="cox")
OS_bas_lauss_oscar <- oscar::oscar(x = X_bas_lauss, y = OS_lauss, family="cox")
RESP_bas_lauss_oscar <- oscar::oscar(x = X_bas_lauss, y = RESP_lauss, family="logistic")
# OSCAR CV
PFS_bas_cv_lauss <- oscar::cv.oscar(PFS_bas_lauss_oscar, fold=5, seed=30)
OS_bas_cv_lauss <- oscar::cv.oscar(OS_bas_lauss_oscar, fold=5, seed=31)
RESP_bas_cv_lauss <- oscar::cv.oscar(RESP_bas_lauss_oscar, fold=5, seed=32)

save.image("temprun.RData")

# Kim et al. xCell
RESP_xce_kim_oscar <- oscar::oscar(x = X_xce_kim, y = RESP_kim, family="logistic")
# OSCAR CV
RESP_xce_cv_kim <- oscar::cv.oscar(RESP_xce_kim_oscar, fold=5, seed=33)
# Kim et al. CUSTOM GMTs
RESP_cus_kim_oscar <- oscar::oscar(x = X_cus_kim, y = RESP_kim, family="logistic")
# OSCAR CV
RESP_cus_cv_kim <- oscar::cv.oscar(RESP_cus_kim_oscar, fold=5, seed=34)
# Kim et al. HALLMARKS
RESP_hal_kim_oscar <- oscar::oscar(x = X_hal_kim, y = RESP_kim, family="logistic")
# OSCAR CV
RESP_hal_cv_kim <- oscar::cv.oscar(RESP_hal_kim_oscar, fold=5, seed=35)
# Kim et al. BASE MARKERS
RESP_bas_kim_oscar <- oscar::oscar(x = X_bas_kim, y = RESP_kim, family="logistic")
# OSCAR CV
RESP_bas_cv_kim <- oscar::cv.oscar(RESP_bas_kim_oscar, fold=5, seed=36)

save.image("temprun.RData")

# Chen et al. xCell
RESP_xce_chen_oscar <- oscar::oscar(x = X_xce_chen, y = RESP_chen, family="logistic")
# OSCAR CV
RESP_xce_cv_chen <- oscar::cv.oscar(RESP_xce_chen_oscar, fold=5, seed=37)
# Chen et al. CUSTOM GMTs
RESP_cus_chen_oscar <- oscar::oscar(x = X_cus_chen, y = RESP_chen, family="logistic")
# OSCAR CV
RESP_cus_cv_chen <- oscar::cv.oscar(RESP_cus_chen_oscar, fold=5, seed=38)
# Chen et al. HALLMARKS
RESP_hal_chen_oscar <- oscar::oscar(x = X_hal_chen, y = RESP_chen, family="logistic")
# OSCAR CV
RESP_hal_cv_chen <- oscar::cv.oscar(RESP_hal_chen_oscar, fold=5, seed=39)
# Chen et al. BASE MARKERS
RESP_bas_chen_oscar <- oscar::oscar(x = X_bas_chen, y = RESP_chen, family="logistic")
# OSCAR CV
RESP_bas_cv_chen <- oscar::cv.oscar(RESP_bas_chen_oscar, fold=5, seed=40)

save.image("temprun.RData")





if(FALSE){
	# Visualizing AIC and/or cross-validation
	visu(OS_bas_hugo_oscar, y="AIC", main="AIC OS, Basic variables, Hugo et al.")
	cv.visu(OS_bas_cv_hugo, main="Cross-validated OS, Basic variables, Hugo et al.", ylab="C-index")
	visu(RESP_bas_hugo_oscar, y="AIC", main="AIC RESP, Basic variables, Hugo et al.")
	cv.visu(RESP_bas_cv_hugo, main="Cross-validated RESP, Basic variables, Hugo et al.", ylab="C-index")



	visu(OS_cus_hugo_oscar, y="AIC", main="AIC OS, Custom GMTs, Hugo et al.")
	cv.visu(OS_cus_cv_hugo, main="Cross-validated OS, Custom GMTs, Hugo et al.", ylab="C-index")
	visu(RESP_cus_hugo_oscar, y="AIC", main="AIC RESP, Custom GMTs, Hugo et al.")
	cv.visu(RESP_cus_cv_hugo, main="Cross-validated RESP, Custom GMTs, Hugo et al.", ylab="Accuracy")


	visu(RESP_cus_riaz_oscar, y="AIC", main="AIC RESP, Custom GMTs, Riaz et al.")
	cv.visu(RESP_cus_cv_riaz, main="Cross-validated RESP, Custom GMTs, Riaz et al.", ylab="Accuracy")

	# cat for printing out features at various ks
	# cat(paste(names(feat(OS_xce_hugo_oscar, k=3)), ", "), "\n")

	# Fails:
	#cv.visu(OS_xce_cv_hugo, main="Cross-validated OS, xCell, Hugo et al.")
	visu(OS_xce_hugo_oscar, y="AIC", main="AIC OS, xCell, Hugo et al.")
	visu(RESP_xce_hugo_oscar, y="AIC", main="AIC RESP, xCell, Hugo et al.")
	cv.visu(RESP_xce_cv_hugo, main="Cross-validated RESP, xCell, Hugo et al.", ylab="Accuracy")
}

if(FALSE){
	visu(OS_xle_prat_oscar, y="AIC", main="AIC OS, xCell, Prat et al.")
	cv.visu(OS_xle_cv_prat, main="Cross-validated OS, xCell, Prat et al.", ylab="C-index")
	visu(RESP_xle_prat_oscar, y="AIC", main="AIC RESP, xCell, Prat et al.")
	cv.visu(RESP_xle_cv_prat, main="Cross-validated RESP, xCell, Prat et al.", ylab="C-index")

	
	cv.visu(OS_bas_cv_hugo, main="Cross-validated OS, Basic variables, Hugo et al.", ylab="C-index")
	cv.visu(OS_bas_cv_lauss, main="Cross-validated OS, Basic variables, Lauss et al.", ylab="C-index")
	cv.visu(RESP_bas_cv_hugo, main="Cross-validated RESP, Basic variables, Hugo et al.", ylab="Accuracy")
	
	cv.visu(RESP_bas_cv_kim, main="Cross-validated RESP, Basic variables, Kim et al.", ylab="Accuracy")
	
	cv.visu(RESP_bas_cv_chen, main="Cross-validated RESP, Basic variables, Chen et al.", ylab="Accuracy")
	
	cv.visu(PFS_bas_cv_westin, main="Cross-validated PFS, Basic variables, Westin et al.", ylab="C-index")
	#cv.visu(PFS_xce_cv_westin, main="Cross-validated PFS, xCell, Westin et al.", ylab="C-index")
	
	cv.visu(PFS_hal_cv_westin, main="Cross-validated PFS, Hallmark GMTs, Westin et al.", ylab="C-index")
	
	cv.visu(PFS_hal_cv_prat, main="Cross-validated PFS, Hallmark GMTs, Prat et al.", ylab="C-index")
	cv.visu(RESP_hal_cv_prat, main="Cross-validated RESP, Hallmark GMTs, Prat et al.", ylab="Accuracy")
	cv.visu(RESP_bas_cv_prat, main="Cross-validated RESP, Basic variables, Prat et al.", ylab="Accuracy")
	
	cv.visu(PFS_cus_cv_lauss, main="Cross-validated PFS, Custom GMTs, Lauss et al.", ylab="C-index")
	cv.visu(OS_cus_cv_lauss, main="Cross-validated OS, Custom GMTs, Lauss et al.", ylab="C-index")
	cv.visu(RESP_cus_cv_lauss, main="Cross-validated RESP, Custom GMTs, Lauss et al.", ylab="Accuracy")
	
	cv.visu(RESP_cus_cv_chen, main="Cross-validated RESP, Custom GMTs, Chen et al.", ylab="Accuracy")
	cv.visu(RESP_cus_cv_kim, main="Cross-validated RESP, Custom GMTs, Kim et al.", ylab="Accuracy")
	
	cv.visu(RESP_hal_cv_chen, main="Cross-validated RESP, Hallmark GMTs, Chen et al.", ylab="Accuracy")
	cv.visu(RESP_hal_cv_kim, main="Cross-validated RESP, Hallmark GMTs, Kim et al.", ylab="Accuracy")
	
	cv.visu(RESP_xce_cv_lauss, main="Cross-validated RESP, xCell, Lauss et al.", ylab="Accuracy")
	cv.visu(RESP_xce_cv_kim, main="Cross-validated RESP, xCell, Kim et al.", ylab="Accuracy")
	
	visu(PFS_xle_lauss_oscar, y="AIC", main="AIC PFS, xCell, Lauss et al.")
	cv.visu(PFS_xle_cv_lauss, main="Cross-validated PFS, xCell, Lauss et al.", ylab="C-index")
	visu(OS_xle_lauss_oscar, y="AIC", main="AIC OS, xCell, Lauss et al.")
	cv.visu(OS_xle_cv_lauss, main="Cross-validated OS, xCell, Lauss et al.", ylab="C-index")
	visu(RESP_xle_lauss_oscar, y="AIC", main="AIC RESP, xCell, Lauss et al.")
	cv.visu(RESP_xle_cv_lauss, main="Cross-validated RESP, xCell, Lauss et al.", ylab="C-index")


	cv.visu(OS_hal_cv_lauss, main="Cross-validated OS, Custom GMTs, Lauss et al.", ylab="C-index")
}




# Plot CVs
#cv.visu(PFS_hal_cv_lauss)
#cv.visu(OS_hal_cv_lauss)
#cv.visu(RESP_hal_cv_lauss)
#cv.visu(RESP_hal_cv_kim)
#cv.visu(RESP_hal_cv_chen)


# Potential LASSO-based benchmarking
if(FALSE){

	# LASSO benchmarking
	library(glmnet)
	# TCGA
	PFS_lasso_tcga <- glmnet(x=X_tcga[!is.na(PFS_tcga),], y=PFS_tcga[!is.na(PFS_tcga)], family="cox")
	PFS_cv_lasso_tcga <- cv.glmnet(x=X_tcga[!is.na(PFS_tcga),], y=PFS_tcga[!is.na(PFS_tcga)], family="cox", nfolds=3)
	colnames(X_tcga)[predict(PFS_lasso_tcga, type="nonzero", s=PFS_cv_lasso_tcga$lambda.min)[,1]]

	OS_lasso_tcga <- glmnet(x=X_tcga[!is.na(OS_tcga),], y=OS_tcga[!is.na(OS_tcga)], family="cox")
	OS_cv_lasso_tcga <- cv.glmnet(x=X_tcga[!is.na(OS_tcga),], y=OS_tcga[!is.na(OS_tcga)], family="cox", nfolds=3)
	colnames(X_tcga)[predict(OS_lasso_tcga, type="nonzero", s=OS_cv_lasso_tcga$lambda.min)[,1]]

	RESP_lasso_tcga <- glmnet(x=X_tcga[!is.na(RESP_tcga),], y=RESP_tcga[!is.na(RESP_tcga)], family="binomial")
	RESP_cv_lasso_tcga <- cv.glmnet(x=X_tcga[!is.na(RESP_tcga),], y=RESP_tcga[!is.na(RESP_tcga)], family="binomial", nfolds=3)
	colnames(X_tcga)[predict(RESP_lasso_tcga, type="nonzero", s=RESP_cv_lasso_tcga$lambda.min)[,1]]

	# Hugo
	OS_lasso_hugo <- glmnet(x=X_hugo[!is.na(OS_hugo),], y=OS_hugo[!is.na(OS_hugo)], family="cox")
	OS_cv_lasso_hugo <- cv.glmnet(x=X_hugo[!is.na(OS_hugo),], y=OS_hugo[!is.na(OS_hugo)], family="cox", nfolds=3)
	colnames(X_hugo)[predict(OS_lasso_hugo, type="nonzero", s=OS_cv_lasso_hugo$lambda.min)[,1]]

	RESP_lasso_hugo <- glmnet(x=X_hugo[!is.na(RESP_hugo),], y=RESP_hugo[!is.na(RESP_hugo)], family="binomial")
	RESP_cv_lasso_hugo <- cv.glmnet(x=X_hugo[!is.na(RESP_hugo),], y=RESP_hugo[!is.na(RESP_hugo)], family="binomial", nfolds=3)
	colnames(X_hugo)[predict(RESP_lasso_hugo, type="nonzero", s=RESP_cv_lasso_hugo$lambda.min)[,1]]

	# Prat
	PFS_lasso_prat <- glmnet(x=X_prat[!is.na(PFS_prat) & as.matrix(PFS_prat)[,1]>0,], y=PFS_prat[!is.na(PFS_prat) & as.matrix(PFS_prat)[,1]>0], family="cox")
	PFS_cv_lasso_prat <- cv.glmnet(x=X_prat[!is.na(PFS_prat) & as.matrix(PFS_prat)[,1]>0,], y=PFS_prat[!is.na(PFS_prat) & as.matrix(PFS_prat)[,1]>0], family="cox", nfolds=3)
	colnames(X_prat)[predict(PFS_lasso_prat, type="nonzero", s=PFS_cv_lasso_prat$lambda.min)[,1]]

	RESP_lasso_prat <- glmnet(x=X_prat[!is.na(RESP_prat),], y=RESP_prat[!is.na(RESP_prat)], family="binomial")
	RESP_cv_lasso_prat <- cv.glmnet(x=X_prat[!is.na(RESP_prat),], y=RESP_prat[!is.na(RESP_prat)], family="binomial", nfolds=3)
	colnames(X_prat)[predict(RESP_lasso_prat, type="nonzero", s=RESP_cv_lasso_prat$lambda.min)[,1]]

	# Westin
	PFS_lasso_westin <- glmnet(x=X_westin[!is.na(PFS_westin) & as.matrix(PFS_westin)[,1]>0,], y=PFS_westin[!is.na(PFS_westin) & as.matrix(PFS_westin)[,1]>0], family="cox")
	PFS_cv_lasso_westin <- cv.glmnet(x=X_westin[!is.na(PFS_westin) & as.matrix(PFS_westin)[,1]>0,], y=PFS_westin[!is.na(PFS_westin) & as.matrix(PFS_westin)[,1]>0], family="cox", nfolds=3)
	colnames(X_westin)[predict(PFS_lasso_westin, type="nonzero", s=PFS_cv_lasso_westin$lambda.min)[,1]]

	# Lauss et al.
	PFS_lasso_lauss <- glmnet(x=X_lauss[!is.na(PFS_lauss),], y=PFS_lauss[!is.na(PFS_lauss)], family="cox")
	PFS_cv_lasso_lauss <- cv.glmnet(x=X_lauss[!is.na(PFS_lauss),], y=PFS_lauss[!is.na(PFS_lauss)], family="cox", nfolds=3)
	colnames(X_lauss)[predict(PFS_lasso_lauss, type="nonzero", s=PFS_cv_lasso_lauss$lambda.min)[,1]]

	OS_lasso_lauss <- glmnet(x=X_lauss[!is.na(OS_lauss),], y=OS_lauss[!is.na(OS_lauss)], family="cox")
	OS_cv_lasso_lauss <- cv.glmnet(x=X_lauss[!is.na(OS_lauss),], y=OS_lauss[!is.na(OS_lauss)], family="cox", nfolds=3)
	colnames(X_lauss)[predict(OS_lasso_lauss, type="nonzero", s=OS_cv_lasso_lauss$lambda.min)[,1]]

	RESP_lasso_lauss <- glmnet(x=X_lauss[!is.na(RESP_lauss),], y=RESP_lauss[!is.na(RESP_lauss)], family="binomial")
	RESP_cv_lasso_lauss <- cv.glmnet(x=X_lauss[!is.na(RESP_lauss),], y=RESP_lauss[!is.na(RESP_lauss)], family="binomial", nfolds=3)
	colnames(X_lauss)[predict(RESP_lasso_lauss, type="nonzero", s=RESP_cv_lasso_lauss$lambda.min)[,1]]

}





# Focusing on type-specific subsets, FALSE-commented out for time being
if(FALSE){


	## MODEL DATA

	# TCGA
	# OSCAR
	PFS_tcga_oscar <- oscar::oscar(x = X_tcga, y = PFS_tcga, family="cox", start=1, verb=1)
	OS_tcga_oscar <- oscar::oscar(x = X_tcga, y = OS_tcga, family="cox", start=1, verb=1)
	RESP_tcga_oscar <- oscar::oscar(x = X_tcga, y = RESP_tcga, family="logistic", start=1, verb=1)
	# OSCAR CV
	PFS_cv_tcga <- oscar::cv.oscar(PFS_tcga_oscar, fold=5, seed=1)
	OS_cv_tcga <- oscar::cv.oscar(OS_tcga_oscar, fold=5, seed=2)
	RESP_cv_tcga <- oscar::cv.oscar(RESP_tcga_oscar, fold=5, seed=3)


	# Lauss et al. (GEO)
	# OSCAR
	PFS_lauss_oscar <- oscar::oscar(x = X_lauss, y = PFS_lauss, family="cox")
	OS_lauss_oscar <- oscar::oscar(x = X_lauss, y = OS_lauss, family="cox")
	RESP_lauss_oscar <- oscar::oscar(x = X_lauss, y = RESP_lauss, family="logistic")
	# OSCAR CV
	PFS_cv_lauss <- oscar::cv.oscar(PFS_lauss_oscar, fold=5, seed=1)
	OS_cv_lauss <- oscar::cv.oscar(OS_lauss_oscar, fold=5, seed=2)
	RESP_cv_lauss <- oscar::cv.oscar(RESP_lauss_oscar, fold=5, seed=3)

	# Westin et al. (GEO)
	# OSCAR
	PFS_westin_oscar <- oscar::oscar(x = X_westin, y = PFS_westin, family="cox", start=1)



	library(survival)
	library(oscar)


	###
	#
	# TCGA (luad and lusc), chemo arm
	#
	###
	# Load premade data 
	load(".\\RData\\gex_tcga.RData")
	load(".\\RData\\dat_tcga.RData")
	# Create X for TCGA
	X_tcga <- curateX(gex=gex_tcga, dat=dat_tcga)
	# Remove redundant columns; should be added to oscar as debugging
	#X_tcga <- X_tcga[,-which(apply(X_tcga, MARGIN=2, FUN=function(x) { all(x==unique(x)[1]) }))]
	X_tcga <- omit.reducols(X_tcga)

	# 3 columns get omitted
	#> dim(X_tcga)
	#[1] 314 117

	#> dim(X_tcga)
	#[1] 314 475

	# Model lusc and luad separately
	X_lusc_tcga <- X_tcga[which(X_tcga[,"isSquamous"] == 1),]
	X_lusc_tcga <- omit.reducols(X_lusc_tcga)
	X_luad_tcga <- X_tcga[which(X_tcga[,"isSquamous"] == 0),]
	X_luad_tcga <- omit.reducols(X_luad_tcga)

	# Whole cohort responses
	PFS_tcga <- survival::Surv(time = dat_tcga$PFS.time, event = dat_tcga$PFS.event)
	OS_tcga <- survival::Surv(time = dat_tcga$OS.time, event = dat_tcga$OS.event)
	RESP_tcga <- as.integer(dat_tcga$Responder)

	# LUSC
	PFS_lusc_tcga <- survival::Surv(time = dat_tcga$PFS.time, event = dat_tcga$PFS.event)[which(X_tcga[,"isSquamous"] == 1)]
	OS_lusc_tcga <- survival::Surv(time = dat_tcga$OS.time, event = dat_tcga$OS.event)[which(X_tcga[,"isSquamous"] == 1)]
	RESP_lusc_tcga <- as.integer(dat_tcga$Responder)[which(X_tcga[,"isSquamous"] == 1)]

	# LUAD
	PFS_luad_tcga <- survival::Surv(time = dat_tcga$PFS.time, event = dat_tcga$PFS.event)[which(X_tcga[,"isSquamous"] == 0)]
	OS_luad_tcga <- survival::Surv(time = dat_tcga$OS.time, event = dat_tcga$OS.event)[which(X_tcga[,"isSquamous"] == 0)]
	RESP_luad_tcga <- as.integer(dat_tcga$Responder)[which(X_tcga[,"isSquamous"] == 0)]

	# Whole cohort
	set.seed(1)
	#PFS_tcga_oscar <- oscar::oscar(x = X_tcga[!is.na(PFS_tcga),], y = PFS_tcga[!is.na(PFS_tcga)], family = "cox", kmax=20, verb=1)
	PFS_tcga_oscar <- oscar::oscar(x = X_tcga[!is.na(PFS_tcga),], y = PFS_tcga[!is.na(PFS_tcga)], family = "cox", verb=1)
	#OS_tcga_oscar <- oscar::oscar(x = X_tcga[!is.na(OS_tcga),], y = OS_tcga[!is.na(OS_tcga)], family = "cox", kmax=20, verb=1)
	OS_tcga_oscar <- oscar::oscar(x = X_tcga[!is.na(OS_tcga),], y = OS_tcga[!is.na(OS_tcga)], family = "cox", verb=1)
	#RESP_tcga_oscar <- oscar::oscar(x = X_tcga[!is.na(RESP_tcga),], y = RESP_tcga[!is.na(RESP_tcga)], family = "logistic", kmax=20, verb=1)
	RESP_tcga_oscar <- oscar::oscar(x = X_tcga[!is.na(RESP_tcga),], y = RESP_tcga[!is.na(RESP_tcga)], family = "logistic", verb=1)

	# Manual bugfix from earlier runs for AIC
	#> OS_tcga_oscar@AIC <- unlist(lapply(OS_tcga_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] }))
	#> PFS_tcga_oscar@AIC <- unlist(lapply(PFS_tcga_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] }))
	#> RESP_tcga_oscar@AIC <- unlist(lapply(RESP_tcga_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] }))

	save(PFS_tcga_oscar, file=".\\RData\\PFS_tcga_oscar.RData")
	save(OS_tcga_oscar, file=".\\RData\\OS_tcga_oscar.RData")
	save(RESP_tcga_oscar, file=".\\RData\\RESP_tcga_oscar.RData")

	par(mfrow=c(1,2))
	plot(OS_tcga_oscar)
	plot(OS_tcga_oscar@AIC, type="l", xlab="Cardinality 'k'", ylab="AIC", main="Overall survival, TCGA")

	par(mfrow=c(1,2))
	plot(PFS_tcga_oscar)
	plot(PFS_tcga_oscar@AIC, type="l", xlab="Cardinality 'k'", ylab="AIC", main="Progression free survival, TCGA")

	par(mfrow=c(1,2))
	plot(RESP_tcga_oscar)
	plot(RESP_tcga_oscar@AIC, type="l", xlab="Cardinality 'k'", ylab="AIC", main="Responder, TCGA")



	# LUSC
	set.seed(1)
	#PFS_lusc_tcga_oscar <- oscar::oscar(x = X_lusc_tcga[!is.na(PFS_lusc_tcga),], y = PFS_lusc_tcga[!is.na(PFS_lusc_tcga)], family = "cox", kmax=20, verb=1)
	PFS_lusc_tcga_oscar <- oscar::oscar(x = X_lusc_tcga[!is.na(PFS_lusc_tcga),], y = PFS_lusc_tcga[!is.na(PFS_lusc_tcga)], family = "cox", verb=1)
	#OS_lusc_tcga_oscar <- oscar::oscar(x = X_lusc_tcga[!is.na(OS_lusc_tcga),], y = OS_lusc_tcga[!is.na(OS_lusc_tcga)], family = "cox", kmax=20, verb=1)
	OS_lusc_tcga_oscar <- oscar::oscar(x = X_lusc_tcga[!is.na(OS_lusc_tcga),], y = OS_lusc_tcga[!is.na(OS_lusc_tcga)], family = "cox", verb=1)
	#RESP_lusc_tcga_oscar <- oscar::oscar(x = X_lusc_tcga[!is.na(RESP_lusc_tcga),], y = RESP_lusc_tcga[!is.na(RESP_lusc_tcga)], family = "logistic", kmax=20, verb=1)
	RESP_lusc_tcga_oscar <- oscar::oscar(x = X_lusc_tcga[!is.na(RESP_lusc_tcga),], y = RESP_lusc_tcga[!is.na(RESP_lusc_tcga)], family = "logistic", verb=1)

	save(PFS_lusc_tcga_oscar, file=".\\RData\\PFS_lusc_tcga_oscar.RData")
	save(OS_lusc_tcga_oscar, file=".\\RData\\OS_lusc_tcga_oscar.RData")
	save(RESP_lusc_tcga_oscar, file=".\\RData\\RESP_lusc_tcga_oscar.RData")

	# LUAD 
	set.seed(1)
	#PFS_luad_tcga_oscar <- oscar::oscar(x = X_luad_tcga[!is.na(PFS_luad_tcga),], y = PFS_luad_tcga[!is.na(PFS_luad_tcga)], family = "cox", kmax=20, verb=1)
	PFS_luad_tcga_oscar <- oscar::oscar(x = X_luad_tcga[!is.na(PFS_luad_tcga),], y = PFS_luad_tcga[!is.na(PFS_luad_tcga)], family = "cox", verb=1)
	#OS_luad_tcga_oscar <- oscar::oscar(x = X_luad_tcga[!is.na(OS_luad_tcga),], y = OS_luad_tcga[!is.na(OS_luad_tcga)], family = "cox", kmax=20, verb=1)
	OS_luad_tcga_oscar <- oscar::oscar(x = X_luad_tcga[!is.na(OS_luad_tcga),], y = OS_luad_tcga[!is.na(OS_luad_tcga)], family = "cox", verb=1)
	#RESP_luad_tcga_oscar <- oscar::oscar(x = X_luad_tcga[!is.na(RESP_luad_tcga),], y = RESP_luad_tcga[!is.na(RESP_luad_tcga)], family = "logistic", kmax=20, verb=1)
	RESP_luad_tcga_oscar <- oscar::oscar(x = X_luad_tcga[!is.na(RESP_luad_tcga),], y = RESP_luad_tcga[!is.na(RESP_luad_tcga)], family = "logistic", verb=1)

	save(PFS_luad_tcga_oscar, file=".\\RData\\PFS_luad_tcga_oscar.RData")
	save(OS_luad_tcga_oscar, file=".\\RData\\OS_luad_tcga_oscar.RData")
	save(RESP_luad_tcga_oscar, file=".\\RData\\RESP_luad_tcga_oscar.RData")


	#par(mfrow=c(1,3))
	#plot(unlist(lapply(PFS_tcga_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] })), type="l", xlab="k", ylab="AIC", main="PFS OSCAR")
	#plot(unlist(lapply(OS_tcga_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] })), type="l", xlab="k", ylab="AIC", main="OS OSCAR")
	#plot(unlist(lapply(RESP_tcga_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] })), type="l", xlab="k", ylab="AIC", main="RESP OSCAR")

	#PFS_tcga_cv_oscar <- oscar::cv.oscar(fit = PFS_tcga_oscar, fold=5, seed=1, verb=0)
	#OS_tcga_cv_oscar <- oscar::cv.oscar(fit = OS_tcga_oscar, fold=5, seed=2, verb=0)
	#RESP_tcga_cv_oscar <- oscar::cv.oscar(fit = RESP_tcga_oscar, fold=5, seed=3, verb=0)

	PFS_luad_tcga_cv_oscar <- oscar::cv.oscar(fit = PFS_luad_tcga_oscar, fold=5, seed=1, verb=0)
	OS_luad_tcga_cv_oscar <- oscar::cv.oscar(fit = OS_luad_tcga_oscar, fold=5, seed=1, verb=0)
	RESP_luad_tcga_cv_oscar <- oscar::cv.oscar(fit = RESP_luad_tcga_oscar, fold=5, seed=1, verb=0)

	PFS_lusc_tcga_cv_oscar <- oscar::cv.oscar(fit = PFS_lusc_tcga_oscar, fold=5, seed=1, verb=0)
	OS_lusc_tcga_cv_oscar <- oscar::cv.oscar(fit = OS_lusc_tcga_oscar, fold=5, seed=1, verb=0)
	RESP_lusc_tcga_cv_oscar <- oscar::cv.oscar(fit = RESP_lusc_tcga_oscar, fold=5, seed=1, verb=0)







	###
	#
	# Hugo et al., melanoma metas
	#
	###
	library(survival); library(oscar)
	# Load premade data 
	load(".\\RData\\gex_hugo.RData")
	load(".\\RData\\dat_hugo.RData")
	# Create X for Hugo et al.
	X_hugo <- curateX(gex=gex_hugo, dat=dat_hugo)
	# Remove redundant columns; should be added to oscar as debugging
	#X_hugo <- X_hugo[,-which(apply(X_hugo, MARGIN=2, FUN=function(x) { all(x==unique(x)[1] | is.na(x)) }))]
	X_hugo <- omit.reducols(X_hugo)
	# from 120, 11 variables are omitted
	#> dim(X_hugo)
	#[1]  27 109

	OS_hugo <- survival::Surv(time = dat_hugo$OS.time, event = dat_hugo$OS.event)
	RESP_hugo <- as.integer(dat_hugo$Responder)

	OS_hugo_oscar <- oscar::oscar(x = X_hugo[!is.na(OS_hugo),], y = OS_hugo[!is.na(OS_hugo)], family = "cox", verb=1)
	RESP_hugo_oscar <- oscar::oscar(x = X_hugo[!is.na(RESP_hugo),], y = RESP_hugo[!is.na(RESP_hugo)], family = "logistic", verb=1)

	# Manual bugfix for the older version of oscar AIC
	#> OS_hugo_oscar@AIC <- unlist(lapply(OS_hugo_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] }))
	#> RESP_hugo_oscar@AIC <- unlist(lapply(RESP_hugo_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] }))

	par(mfrow=c(1,2))
	plot(OS_hugo_oscar, main="Overall survival, Hugo et al.")
	plot(OS_hugo_oscar@AIC, type="l", xlab="Cardinality 'k'", ylab="AIC", main="Overall survival, Hugo et al.")

	par(mfrow=c(1,2))
	plot(RESP_hugo_oscar, main="Responder, Hugo et al.")
	plot(RESP_hugo_oscar@AIC, type="l", xlab="Cardinality 'k'", ylab="AIC", main="Responder, Hugo et al.")

	save(OS_hugo_oscar, file=".\\RData\\OS_hugo_oscar.RData")
	save(RESP_hugo_oscar, file=".\\RData\\RESP_hugo_oscar.RData")

	OS_hugo_cv_oscar <- oscar::cv.oscar(fit = OS_hugo_oscar, fold=5, seed=1, verb=0)
	RESP_hugo_cv_oscar <- oscar::cv.oscar(fit = RESP_hugo_oscar, fold=5, seed=2, verb=0)

	save(OS_hugo_cv_oscar, file=".\\RData\\OS_hugo_cv_oscar.RData")
	save(RESP_hugo_cv_oscar, file=".\\RData\\RESP_hugo_cv_oscar.RData")



	###
	#
	# Prat et al., small amount of measured genes but relevant phenotypes
	#
	###
	library(survival); library(oscar)
	# Load premade data 
	load(".\\RData\\gex_prat.RData")
	load(".\\RData\\dat_prat.RData")
	# Create X for Hugo et al.
	X_prat <- curateX(gex=gex_prat, dat=dat_prat)
	# Remove redundant columns; should be added to oscar as debugging
	X_prat <- omit.reducols(X_prat)
	# A lot of variables that could not be estimated
	#> dim(X_prat)
	#[1] 65 62
	X_prat <- omit.nacols(X_prat)
	#> dim(X_prat)
	#[1] 65 61

	PFS_prat <- survival::Surv(time = dat_prat$PFS.time, event = dat_prat$PFS.event)
	RESP_prat <- as.integer(dat_prat$Responder)

	PFS_prat_oscar <- oscar::oscar(x = X_prat[!is.na(PFS_prat),], y = PFS_prat[!is.na(PFS_prat)], family = "cox", verb=1)
	RESP_prat_oscar <- oscar::oscar(x = X_prat[!is.na(RESP_prat),], y = RESP_prat[!is.na(RESP_prat)], family = "logistic", verb=1)

	PFS_prat_oscar@AIC <- unlist(lapply(PFS_prat_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] }))
	RESP_prat_oscar@AIC <- unlist(lapply(RESP_prat_oscar@fits, FUN=function(z) { stats::extractAIC(z)[2] }))

	save(PFS_prat_oscar, file=".\\RData\\PFS_prat_oscar.RData")
	save(RESP_prat_oscar, file=".\\RData\\RESP_prat_oscar.RData")


	par(mfrow=c(1,2))
	plot(PFS_prat_oscar, main="Overall survival, Hugo et al.")
	plot(PFS_prat_oscar@AIC, type="l", xlab="Cardinality 'k'", ylab="AIC", main="Progression free survival, Prat et al.")

	par(mfrow=c(1,2))
	plot(RESP_prat_oscar, main="Responder, Hugo et al.")
	plot(RESP_prat_oscar@AIC, type="l", xlab="Cardinality 'k'", ylab="AIC", main="Responder, Prat et al.")



	PFS_prat_cv_oscar <- oscar::cv.oscar(fit = PFS_prat_oscar, fold=5, seed=1, verb=0)
	RESP_prat_cv_oscar <- oscar::cv.oscar(fit = RESP_prat_oscar, fold=5, seed=2, verb=0)

	save(PFS_prat_cv_oscar, file=".\\RData\\PFS_prat_cv_oscar.RData")
	save(RESP_prat_cv_oscar, file=".\\RData\\RESP_prat_cv_oscar.RData")

}




