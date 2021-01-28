###
#
# Load readily curated datasets into workspace and continue working (i.e. modelling)
#
###

# Road pre-curated datasets for use
load(".\\RData\\gex_tcga.RData")
load(".\\RData\\dat_tcga.RData")

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
	keygenes = c(
		"CD274", "PDL1", # PD-L1
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
		"mTOR",
		"PDK-1",
		"AKT",
		# Mikko NGS syöpäpaneeli 1 lisäyksiä
		"PIK3CA",
		"KIT",
		"NRAS",
		"PDGFRA",
		# Other interesting? Cytokines, chemokines
		"CLCL10",
		"CSCL11"
	),
	normfunc = function(input) { input }, # Function for normalizing gene values to be used as variables - could be e.g. z-score within sample? log-transform if normalized count data >0?
	gmts = 1,
	idcs = 1:2,
	clinvars = c("Age", "Smoking", "ECOG", "Squamous", "TMB", "SEX", "PDL1"),
	scores = c( # Included scoring metrics (preferably IO and lung cancer related)
		"TIS", # Tumor inflammatory score
		"MSI", # Microsatellite instability
		"APM", # Antigen processing / presentation machinery
		"IFN"  # Interferon gamma signalling
	)
){
	# Format data matrix X
	X <- matrix(NA, nrow=ncol(gex), ncol=0)
	
	## Cherry-pick some relevant genes for modelling
	# Pick genes that are present, warn of non-present symbols, column bind into variables
	warning(paste("Following key gene symbols not found from gex: ", paste(keygenes[which(!keygenes %in% rownames(gex))], collapse=", ")), sep="")
	keygenes <- keygenes[which(keygenes %in% rownames(gex))]
	X <- cbind(X, do.call("cbind", lapply(keygenes, FUN=function(z){
		print(paste(z, "..."))
		normfunc(gex[z,])
	})))
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
	}
	# If Squamous histologic subtype
	if("Squamous" %in% clinvars){
		print("Squamous")
		X <- cbind(X, isSquamous = as.integer(dat[,"CRFHIST"] == "SQUAMOUS"))
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

	# Create selected gmt-variables
	# Hallmarks
	if(1 %in% gmts){
		print("Hallmarks")
		X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_h))) # Hallmarks
		# Omit selected hallmarks based on a priori knowledge or educated guesses
		X <- X[,-grep("CHOLESTEROL|ESTROGEN|ANDROGEN|FATTY_ACID|OXIDATIVE|GLYCOLYSIS|REACTIVE_OXYGEN|UV_RESPONSE|BILE_ACID|ALLOGRAFT|SPERMATOGENESIS|PANCREAS", colnames(X))]
	}
	# Oncogenic
	if(2 %in% gmts){
		print("Oncogenic")
		X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_c6))) # Oncogenic
	}
	# Immunology
	if(3 %in% gmts){
		print("Immunology")
		X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_c7))) # Immunology
	}

	## Immune deconvolution
	#install.packages("remotes")
	#remotes::install_github("icbi-lab/immunedeconv")
	library(immunedeconv)
	rename_idc <- function(x, method="idc_"){
		tmp <- gsub(" ", "_", gsub("\\-", ".", c(x[,1])[[1]]))
		x <- as.matrix(x[,-1])
		rownames(x) <- paste(method, "_", tmp, sep="")
		x
	}

	# xCell
	if(1 %in% idcs){
		print("xCell")
		tmp <- rename_idc(immunedeconv::deconvolute(gex, method="xcell"), method="xce")
		X <- cbind(X, t(tmp))
	}
	# MCP counter
	if(2 %in% idcs){
		print("MCP counter")
		tmp <- rename_idc(immunedeconv::deconvolute(gex, method="mcp_counter"), method="mcp")
		X <- cbind(X, t(tmp))
	}
	
	
	# Return X with newly derived variables
	as.matrix(X)
}

# Create X for TCGA
X_tcga <- curateX(gex=gex_tcga, dat=dat_tcga)
# Remove redundant columns; should be added to oscar as debugging
X_tcga <- X_tcga[,-which(apply(X_tcga, MARGIN=2, FUN=function(x) { all(x==unique(x)[1]) }))]
# 3 columns get omitted

#> dim(X_tcga)
#[1] 314 120
library(survival)
library(oscar)

PFS_tcga <- survival::Surv(time = dat_tcga$PFS.time, event = dat_tcga$PFS.event)
OS_tcga <- survival::Surv(time = dat_tcga$OS.time, event = dat_tcga$OS.event)
RESP_tcga <- as.integer(dat_tcga$Responder)

set.seed(1)
PFS_tcga_oscar <- oscar::oscar(x = X_tcga[!is.na(PFS_tcga),], y = PFS_tcga[!is.na(PFS_tcga)], family = "cox", kmax=20, verb=1)
OS_tcga_oscar <- oscar::oscar(x = X_tcga[!is.na(OS_tcga),], y = OS_tcga[!is.na(OS_tcga)], family = "cox", kmax=20, verb=1)
RESP_tcga_oscar <- oscar::oscar(x = X_tcga[!is.na(RESP_tcga),], y = RESP_tcga[!is.na(RESP_tcga)], family = "logistic", kmax=20, verb=1)

PFS_tcga_cv_oscar <- oscar::cv.oscar(fit = PFS_tcga_oscar, fold=10, seed=1, verb=0)
OS_tcga_cv_oscar <- oscar::cv.oscar(fit = OS_tcga_oscar, fold=10, seed=2, verb=0)
RESP_tcga_cv_oscar <- oscar::cv.oscar(fit = RESP_tcga_oscar, fold=10, seed=3, verb=0)




