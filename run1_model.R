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
		"FOXM1"
		
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
	rename_idc <- function(x, method="idc_"){
		tmp <- gsub(" ", "_", gsub("\\-", "minus", c(x[,1])[[1]]))
		x <- as.matrix(x[,-1])
		rownames(x) <- paste(method, "_", tmp, sep="")
		x
	}

	# xCell
	if(1 %in% idcs){
		print("xCell")
		try({
			tmp <- rename_idc(immunedeconv::deconvolute(gex, method="xcell"), method="xce")
			X <- cbind(X, t(tmp))
		})
	}
	# MCP counter
	# NOTE: MCP counter is not functioning properly inside the cloud via Docker; it attempts to download latest 
	if(2 %in% idcs){
		print("MCP counter")
		try({
			tmp <- rename_idc(immunedeconv::deconvolute(gex, method="mcp_counter"), method="mcp")
			X <- cbind(X, t(tmp))
		})
	}
	
	# Sanitize '+' symbol
	colnames(X) <- gsub("\\+", "plus", colnames(X))
	
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

# Load premade GEX / DAT and generate X matrix
# Whole TCGA
load(".\\RData\\gex_tcga.RData")
load(".\\RData\\dat_tcga.RData")
X_tcga <- curateX(gex=gex_tcga, dat=dat_tcga)
X_tcga <- omit.reducols(X_tcga)
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
X_prat <- omit.reducols(curateX(gex=gex_prat, dat=dat_prat))
PFS_prat <- survival::Surv(time = dat_prat$PFS.time, event = dat_prat$PFS.event)
RESP_prat <- as.integer(dat_prat$Responder)
# Westin et al. (GEO)
load(".\\RData\\gex_westin.RData")
load(".\\RData\\dat_westin.RData")
X_westin <- omit.reducols(curateX(gex=gex_westin, dat=dat_westin))
PFS_westin <- survival::Surv(time = dat_westin$PFS.time, event = dat_westin$PFS.event)
# Lauss et al. (TIDE)
load(".\\RData\\gex_lauss.RData")
load(".\\RData\\dat_lauss.RData")
X_lauss <- omit.reducols(curateX(gex=gex_lauss, dat=dat_lauss))
PFS_lauss <- survival::Surv(time=dat_lauss[,"PFS.time"], event=dat_lauss[,"PFS.event"])
OS_lauss <- survival::Surv(time=dat_lauss[,"OS.time"], event=dat_lauss[,"OS.event"])
RESP_lauss <- dat_lauss[,"Responder"]
# Kim et al. (TIDE)
load(".\\RData\\gex_kim.RData")
load(".\\RData\\dat_kim.RData")
X_kim <- omit.reducols(curateX(gex=gex_kim, dat=dat_kim))
RESP_kim <- dat_kim[,"Responder"]
# Chen et al. (TIDE)
load(".\\RData\\gex_chen.RData")
load(".\\RData\\dat_chen.RData")
X_chen <- omit.reducols(curateX(gex=gex_chen, dat=dat_chen))
RESP_chen <- dat_chen[,"Responder"]
# Save image containing the GEXs, DATs, Xs, and various y-responses
save.image("temp.RData")

## MODEL DATA

# TCGA
# OSCAR
PFS_tcga_oscar <- oscar::oscar(x = X_tcga, y = PFS_tcga, family="cox")
OS_tcga_oscar <- oscar::oscar(x = X_tcga, y = OS_tcga, family="cox")
RESP_tcga_oscar <- oscar::oscar(x = X_tcga, y = RESP_tcga, family="logistic")
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
PFS_westin_oscar <- oscar::oscar(x = X_westin, y = PFS_westin, family="cox")



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



