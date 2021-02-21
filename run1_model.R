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
		"CD276", "B7-H3", "B7.H3",
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
		"TLR8"
	)),
	normfunc = function(input) { input }, # Function for normalizing gene values to be used as variables - could be e.g. z-score within sample? log-transform if normalized count data >0?
	gmts = c(1,2,4,5), # Hallmarks, oncology, custom self-made GMTs, filtered curated pathways from e.g. KEGG
	idcs = 1,
	clinvars = c("Age", "Smoking", "ECOG", "Squamous", "TMB", "SEX", "PDL1")
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
	# If Squamous histologic subtype
	if("NonSquamous" %in% clinvars){
		print("NonSquamous")
		X <- cbind(X, isSquamous = as.integer(dat[,"CRFHIST"] == "NON-SQUAMOUS"))
		# 0s, NAs and 1s -> code into -1, 0, +1
		if(NA %in% X[,"isNonSquamous"] & !all(is.na(X[,"isNonSquamous"]))){
			X[X[,"isNonSquamous"]==0,"isNonSquamous"] <- -1
			X[is.na(X[,"isNonSquamous"]),"isNonSquamous"] <- 0
		}
	}
	# TMB median, coding into lower (-1), NA (0) or higher (+1), and raw TMB values
	if("TMB" %in% clinvars){
		print("TMB")
		medianTMB <- median(dat[,"TMB"], na.rm=TRUE)
		# No finite TMB threshold could be determined
		if(!is.finite(medianTMB)){
			X <- cbind(X, TMBmedian = 0)
			X <- cbind(X, TMB = 0)
			X <- cbind(X, log10TMB = 0)
		}else{
			X <- cbind(X, TMBmedian = ifelse(is.na(dat[,"TMB"]), 0, ifelse(dat[,"TMB"]>=medianTMB, 1, -1)))
			X <- cbind(X, TMB = dat[,"TMB"])
			X <- cbind(X, log10TMB = log10(dat[,"TMB"]))
			#X[is.na(X[,"TMB"]),"TMB"] <- medianTMB
		}		
	}
	# Sex, binary indicator for +1 that patient was male
	if("SEX" %in% clinvars){
		print("SEX")
		X <- cbind(X, isMale = as.integer(dat[,"SEX"] %in% c("M", "Male")))
	}
	# PDL1, coding into lower (-1), NA (0) or higher (+1)
	if("PDL1" %in% clinvars){
		print("PDL1")
		medianPDL1 <- median(dat[,"PDL1"], na.rm=TRUE)
		# No finite TMB threshold could be determined
		if(!is.finite(medianPDL1)){
			X <- cbind(X, PDL1median = 0)
			X <- cbind(X, PDL1 = 0)
		}else{
			X <- cbind(X, PDL1median = ifelse(is.na(dat[,"PDL1"]), 0, ifelse(dat[,"PDL1"]>=medianPDL1, 1, -1)))
			X <- cbind(X, PDL1 = ifelse(is.na(dat[,"PDL1"]), medianPDL1, dat[,"PDL1"]))
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
	colnames(X) <- paste("BASE_", colnames(X), sep="")

	## GSVA
	library(GSVA)	
	# Hallmark gene sets
	gmt_h <- GSEABase::getGmt("h.all.v7.2.symbols.gmt")
	# Oncogenic
	gmt_c6 <- GSEABase::getGmt("c6.all.v7.2.symbols.gmt")
	# Immunogenic
	gmt_c7 <- GSEABase::getGmt("c7.all.v7.2.symbols.gmt")
	# Custom self-made
	gmt_custom <- GSEABase::getGmt("selfmade.gmt")
	# Curated GMTs, subset to interesting
	gmt_c2 <- GSEABase::getGmt("c2.all.v7.2.symbols.gmt")
	gmt_c2 <- gmt_c2[grep("NEUTROPH|LEUKOCYT|CD4|CD8|INTERLEUK|INFLAM|T.CELL|B.CELL|TCELL|BCELL|CHEMOK|CYTOK", unlist(lapply(gmt_c2, FUN=function(x) x@setName)))]

	# Create selected gmt-variables
	# Add parameter mx.diff = FALSE because it's unsure whether the test statistic is uni- or bivariate (genes expressed in both extremes within gene set)
	# Hallmarks
	if(1 %in% gmts){
		print("Hallmarks")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_h, verbose=FALSE, mx.diff=FALSE)))
			# Omit selected hallmarks based on a priori knowledge or educated guesses
			X <- X[,-grep("CHOLESTEROL|ESTROGEN|ANDROGEN|FATTY_ACID|OXIDATIVE|GLYCOLYSIS|REACTIVE_OXYGEN|UV_RESPONSE|BILE_ACID|ALLOGRAFT|SPERMATOGENESIS|PANCREAS", colnames(X))]
		})
	}
	# Oncogenic
	if(2 %in% gmts){
		print("Oncogenic")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_c6, verbose=FALSE, mx.diff=FALSE))) # Oncogenic
		})
	}
	# Immunology
	if(3 %in% gmts){
		print("Immunology")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_c7, verbose=FALSE, mx.diff=FALSE))) # Immunology
		})
	}
	# Custom self made GMTs
	if(4 %in% gmts){
		print("Selfmade GMTs")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE, mx.diff=FALSE))) # Custom self made GMTs
		})
	}
	# Curated gene pathways
	if(5 %in% gmts){
		print("Curated GTMs from various databases")
		try({
			X <- cbind(X, t(GSVA::gsva(as.matrix(gex), gmt_c2, verbose=FALSE, mx.diff=FALSE))) # Curated GMTs
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

# Do logz-transformation for the synthetic (or target) data
logz <- function(x) { 
	tmp <- scale(log(x+1)) 
	if(any(!is.finite(tmp))){
		tmp[!is.finite(tmp)] <- 0
	}
	colnames(tmp) <- rownames(x)
	tmp
}

###
#
# SYNTHETIC
#
###
gex_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
gex_synthetic <- gex_synthetic[order(rownames(gex_synthetic)),]
gex_synthetic <- as.matrix(gex_synthetic)
dat_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
#tmp <- rownames(gex_synthetic)
#gex_synthetic <- t(apply(gex_synthetic, MARGIN=1, FUN=logz))
#rownames(gex_synthetic) <- tmp
#rm(tmp)
X_synthetic <- cbind(omit.reducols(curateX(gex=gex_synthetic, dat=dat_synthetic)), isIO = 1, isChemo = 0)

# Normalize using quantile normalization
if(FALSE){
	geq_hugo <- preprocessCore::normalize.quantiles.use.target(gex_hugo, target=c(gex_synthetic))
	rownames(geq_hugo) <- rownames(gex_hugo)
	colnames(geq_hugo) <- colnames(gex_hugo)
	geq_prat <- preprocessCore::normalize.quantiles.use.target(gex_prat, target=c(gex_synthetic))
	rownames(geq_prat) <- rownames(gex_prat)
	colnames(geq_prat) <- colnames(gex_prat)
	geq_westin <- preprocessCore::normalize.quantiles.use.target(gex_westin, target=c(gex_synthetic))
	rownames(geq_westin) <- rownames(gex_westin)
	colnames(geq_westin) <- colnames(gex_westin)
	geq_lauss <- preprocessCore::normalize.quantiles.use.target(gex_lauss, target=c(gex_synthetic))
	rownames(geq_lauss) <- rownames(gex_lauss)
	colnames(geq_lauss) <- colnames(gex_lauss)
	geq_gide <- preprocessCore::normalize.quantiles.use.target(gex_gide, target=c(gex_synthetic))
	rownames(geq_gide) <- rownames(gex_gide)
	colnames(geq_gide) <- colnames(gex_gide)
	geq_chen <- preprocessCore::normalize.quantiles.use.target(gex_chen, target=c(gex_synthetic))
	rownames(geq_chen) <- rownames(gex_chen)
	colnames(geq_chen) <- colnames(gex_chen)
	geq_kim <- preprocessCore::normalize.quantiles.use.target(gex_kim, target=c(gex_synthetic))
	rownames(geq_kim) <- rownames(gex_kim)
	colnames(geq_kim) <- colnames(gex_kim)
	geq_riaz <- preprocessCore::normalize.quantiles.use.target(gex_riaz, target=c(gex_synthetic))
	rownames(geq_riaz) <- rownames(gex_riaz)
	colnames(geq_riaz) <- colnames(gex_riaz)
	geq_braun_nivo <- preprocessCore::normalize.quantiles.use.target(gex_braun_nivo, target=c(gex_synthetic))
	rownames(geq_braun_nivo) <- rownames(gex_braun_nivo)
	colnames(geq_braun_nivo) <- colnames(gex_braun_nivo)
	geq_braun_ever <- preprocessCore::normalize.quantiles.use.target(gex_braun_ever, target=c(gex_synthetic))
	rownames(geq_braun_ever) <- rownames(gex_braun_ever)
	colnames(geq_braun_ever) <- colnames(gex_braun_ever)
	geq_tcga <- preprocessCore::normalize.quantiles.use.target(gex_tcga, target=c(gex_synthetic))
	rownames(geq_tcga) <- rownames(gex_tcga)
	colnames(geq_tcga) <- colnames(gex_tcga)
})

# Load premade GEX / DAT and generate X matrix
# Whole TCGA
load(".\\RData\\gex_tcga.RData")
load(".\\RData\\dat_tcga.RData")
gex_tcga <- logz(gex_tcga)
X_tcga <- cbind(omit.reducols(curateX(gex=gex_tcga, dat=dat_tcga)), isIO = 0, isChemo = 1)
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
gex_hugo <- logz(gex_hugo)
X_hugo <- cbind(omit.reducols(curateX(gex=gex_hugo, dat=dat_hugo)), isIO = 1, isChemo = 0)
X_xce_hugo <- X_hugo[,grep("xce_", colnames(X_hugo))]
X_cus_hugo <- X_hugo[,grep("CUSTOM_", colnames(X_hugo))]
X_hal_hugo <- X_hugo[,grep("HALLMARK_", colnames(X_hugo))]
X_bas_hugo <- X_hugo[,grep("BASE_", colnames(X_hugo))]
OS_hugo <- survival::Surv(time = dat_hugo$OS.time, event = dat_hugo$OS.event)
RESP_hugo <- as.integer(dat_hugo$Responder)
# Prat et al. (GEO)
load(".\\RData\\gex_prat.RData")
load(".\\RData\\dat_prat.RData")
gex_hugo <- logz(gex_hugo)
X_prat <- cbind(omit.reducols(curateX(gex=gex_prat, dat=dat_prat)), isIO = 1, isChemo = 0)
X_xce_prat <- X_prat[,grep("xce_", colnames(X_prat))]
X_cus_prat <- X_prat[,grep("CUSTOM_", colnames(X_prat))]
X_hal_prat <- X_prat[,grep("HALLMARK_", colnames(X_prat))]
X_bas_prat <- X_prat[,grep("BASE_", colnames(X_prat))]
PFS_prat <- survival::Surv(time = dat_prat$PFS.time, event = dat_prat$PFS.event)
RESP_prat <- as.integer(dat_prat$Responder)
# Westin et al. (GEO)
load(".\\RData\\gex_westin.RData")
load(".\\RData\\dat_westin.RData")
gex_hugo <- logz(gex_hugo)
X_westin <- cbind(omit.reducols(curateX(gex=gex_westin, dat=dat_westin)), isIO = 1, isChemo = 0)
X_xce_westin <- X_westin[,grep("xce_", colnames(X_westin))]
X_cus_westin <- X_westin[,grep("CUSTOM_", colnames(X_westin))]
X_hal_westin <- X_westin[,grep("HALLMARK_", colnames(X_westin))]
X_bas_westin <- X_westin[,grep("BASE_", colnames(X_westin))]
PFS_westin <- survival::Surv(time = dat_westin$PFS.time, event = dat_westin$PFS.event)
# Riaz et al. (GEO)
load(".\\RData\\gex_riaz.RData")
load(".\\RData\\dat_riaz.RData")
gex_hugo <- logz(gex_hugo)
X_riaz <- cbind(omit.reducols(curateX(gex=gex_riaz, dat=dat_riaz)), isIO = 1, isChemo = 0)
X_xce_riaz <- X_riaz[,grep("xce_", colnames(X_riaz))]
X_cus_riaz <- X_riaz[,grep("CUSTOM_", colnames(X_riaz))]
X_hal_riaz <- X_riaz[,grep("HALLMARK_", colnames(X_riaz))]
X_bas_riaz <- X_riaz[,grep("BASE_", colnames(X_riaz))]
RESP_riaz <- dat_riaz[,"Responder"]
# Lauss et al. (TIDE)
load(".\\RData\\gex_lauss.RData")
load(".\\RData\\dat_lauss.RData")
gex_lauss <- logz(gex_lauss)
X_lauss <- cbind(omit.reducols(curateX(gex=gex_lauss, dat=dat_lauss)), isIO = 1, isChemo = 0)
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
gex_kim <- logz(gex_kim)
X_kim <- cbind(omit.reducols(curateX(gex=gex_kim, dat=dat_kim)), isIO = 1, isChemo = 0)
X_xce_kim <- X_kim[,grep("xce_", colnames(X_kim))]
X_cus_kim <- X_kim[,grep("CUSTOM_", colnames(X_kim))]
X_hal_kim <- X_kim[,grep("HALLMARK_", colnames(X_kim))]
X_bas_kim <- X_kim[,grep("BASE_", colnames(X_kim))]
RESP_kim <- dat_kim[,"Responder"]
# Chen et al. (TIDE)
load(".\\RData\\gex_chen.RData")
load(".\\RData\\dat_chen.RData")
gex_chen <- logz(gex_chen)
X_chen <- cbind(omit.reducols(curateX(gex=gex_chen, dat=dat_chen)), isIO = 1, isChemo = 0)
X_xce_chen <- X_chen[,grep("xce_", colnames(X_chen))]
X_cus_chen <- X_chen[,grep("CUSTOM_", colnames(X_chen))]
X_hal_chen <- X_chen[,grep("HALLMARK_", colnames(X_chen))]
X_bas_chen <- X_chen[,grep("BASE_", colnames(X_chen))]
RESP_chen <- dat_chen[,"Responder"]
# Gide et al. (TIDE)
load(".\\RData\\gex_gide.RData")
load(".\\RData\\dat_gide.RData")
gex_gide <- logz(gex_gide)
X_gide <- cbind(omit.reducols(curateX(gex=gex_gide, dat=dat_gide)), isIO = 1, isChemo = 0)
X_xce_gide <- X_gide[,grep("xce_", colnames(X_gide))]
X_cus_gide <- X_gide[,grep("CUSTOM_", colnames(X_gide))]
X_hal_gide <- X_gide[,grep("HALLMARK_", colnames(X_gide))]
X_bas_gide <- X_gide[,grep("BASE_", colnames(X_gide))]
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
gex_braun_nivo <- logz(gex_braun_nivo)
gex_braun_ever <- logz(gex_braun_ever)
X_braun_nivo <- cbind(omit.reducols(curateX(gex=gex_braun_nivo, dat=dat_braun_nivo)), isIO = 1, isChemo = 0)
X_braun_ever <- cbind(omit.reducols(curateX(gex=gex_braun_ever, dat=dat_braun_ever)), isIO = 0, isChemo = 1)
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

# Save image containing the GEXs, DATs, Xs, and various y-responses
save.image("datas.RData")



## OSCAR + cross-validation
library(oscar)


### OSCAR FITS

# TCGA Chemo control arm 
PFS_xce_tcga_oscar <- oscar::oscar(x = X_xce_tcga, y = PFS_tcga, family="cox", verb=1)
OS_xce_tcga_oscar <- oscar::oscar(x = X_xce_tcga, y = OS_tcga, family="cox", verb=1)
RESP_xce_tcga_oscar <- oscar::oscar(x = X_xce_tcga, y = RESP_tcga, family="logistic", verb=1)
PFS_cus_tcga_oscar <- oscar::oscar(x = X_cus_tcga, y = PFS_tcga, family="cox", verb=1)
OS_cus_tcga_oscar <- oscar::oscar(x = X_cus_tcga, y = OS_tcga, family="cox", verb=1)
RESP_cus_tcga_oscar <- oscar::oscar(x = X_cus_tcga, y = RESP_tcga, family="logistic", verb=1)
PFS_hal_tcga_oscar <- oscar::oscar(x = X_hal_tcga, y = PFS_tcga, family="cox", verb=1)
OS_hal_tcga_oscar <- oscar::oscar(x = X_hal_tcga, y = OS_tcga, family="cox", verb=1)
RESP_hal_tcga_oscar <- oscar::oscar(x = X_hal_tcga, y = RESP_tcga, family="logistic", verb=1)
PFS_bas_tcga_oscar <- oscar::oscar(x = X_bas_tcga, y = PFS_tcga, family="cox", verb=1)
OS_bas_tcga_oscar <- oscar::oscar(x = X_bas_tcga, y = OS_tcga, family="cox", verb=1)
RESP_bas_tcga_oscar <- oscar::oscar(x = X_bas_tcga, y = RESP_tcga, family="logistic", verb=1)

save.image("temprun_oscar.RData")

# Hugo et al. xCell
OS_xce_hugo_oscar <- oscar::oscar(x = X_xce_hugo, y = OS_hugo, family="cox")
RESP_xce_hugo_oscar <- oscar::oscar(x = X_xce_hugo, y = RESP_hugo, family="logistic")
OS_cus_hugo_oscar <- oscar::oscar(x = X_cus_hugo, y = OS_hugo, family="cox")
RESP_cus_hugo_oscar <- oscar::oscar(x = X_cus_hugo, y = RESP_hugo, family="logistic")
OS_hal_hugo_oscar <- oscar::oscar(x = X_hal_hugo, y = OS_hugo, family="cox")
RESP_hal_hugo_oscar <- oscar::oscar(x = X_hal_hugo, y = RESP_hugo, family="logistic")
OS_bas_hugo_oscar <- oscar::oscar(x = X_bas_hugo, y = OS_hugo, family="cox")
RESP_bas_hugo_oscar <- oscar::oscar(x = X_bas_hugo, y = RESP_hugo, family="logistic")

save.image("temprun_oscar.RData")

# Prat et al. 
PFS_xce_prat_oscar <- oscar::oscar(x = X_xce_prat, y = PFS_prat, family="cox")
RESP_xce_prat_oscar <- oscar::oscar(x = X_xce_prat, y = RESP_prat, family="logistic")
PFS_cus_prat_oscar <- oscar::oscar(x = X_cus_prat, y = PFS_prat, family="cox")
RESP_cus_prat_oscar <- oscar::oscar(x = X_cus_prat, y = RESP_prat, family="logistic")
PFS_hal_prat_oscar <- oscar::oscar(x = X_hal_prat, y = PFS_prat, family="cox")
RESP_hal_prat_oscar <- oscar::oscar(x = X_hal_prat, y = RESP_prat, family="logistic")
PFS_bas_prat_oscar <- oscar::oscar(x = X_bas_prat, y = PFS_prat, family="cox")
RESP_bas_prat_oscar <- oscar::oscar(x = X_bas_prat, y = RESP_prat, family="logistic")

save.image("temprun_oscar.RData")

# Westin et al. 
PFS_xce_westin_oscar <- oscar::oscar(x = X_xce_westin, y = PFS_westin, family="cox")
PFS_cus_westin_oscar <- oscar::oscar(x = X_cus_westin, y = PFS_westin, family="cox")
PFS_hal_westin_oscar <- oscar::oscar(x = X_hal_westin, y = PFS_westin, family="cox")
PFS_bas_westin_oscar <- oscar::oscar(x = X_bas_westin, y = PFS_westin, family="cox")

save.image("temprun_oscar.RData")

# Riaz et al. 
RESP_xce_riaz_oscar <- oscar::oscar(x = X_xce_riaz, y = RESP_riaz, family="logistic")
RESP_cus_riaz_oscar <- oscar::oscar(x = X_cus_riaz, y = RESP_riaz, family="logistic")
RESP_hal_riaz_oscar <- oscar::oscar(x = X_hal_riaz, y = RESP_riaz, family="logistic")
RESP_bas_riaz_oscar <- oscar::oscar(x = X_bas_riaz, y = RESP_riaz, family="logistic")

save.image("temprun_oscar.RData")

# Lauss et al. 
PFS_xce_lauss_oscar <- oscar::oscar(x = X_xce_lauss, y = PFS_lauss, family="cox")
OS_xce_lauss_oscar <- oscar::oscar(x = X_xce_lauss, y = OS_lauss, family="cox")
RESP_xce_lauss_oscar <- oscar::oscar(x = X_xce_lauss, y = RESP_lauss, family="logistic")
PFS_cus_lauss_oscar <- oscar::oscar(x = X_cus_lauss, y = PFS_lauss, family="cox")
OS_cus_lauss_oscar <- oscar::oscar(x = X_cus_lauss, y = OS_lauss, family="cox")
RESP_cus_lauss_oscar <- oscar::oscar(x = X_cus_lauss, y = RESP_lauss, family="logistic")
PFS_hal_lauss_oscar <- oscar::oscar(x = X_hal_lauss, y = PFS_lauss, family="cox")
OS_hal_lauss_oscar <- oscar::oscar(x = X_hal_lauss, y = OS_lauss, family="cox")
RESP_hal_lauss_oscar <- oscar::oscar(x = X_hal_lauss, y = RESP_lauss, family="logistic")
PFS_bas_lauss_oscar <- oscar::oscar(x = X_bas_lauss, y = PFS_lauss, family="cox")
OS_bas_lauss_oscar <- oscar::oscar(x = X_bas_lauss, y = OS_lauss, family="cox")
RESP_bas_lauss_oscar <- oscar::oscar(x = X_bas_lauss, y = RESP_lauss, family="logistic")

save.image("temprun_oscar.RData")

# Gide et al. 
PFS_xce_gide_oscar <- oscar::oscar(x = X_xce_gide, y = PFS_gide, family="cox")
OS_xce_gide_oscar <- oscar::oscar(x = X_xce_gide, y = OS_gide, family="cox")
RESP_xce_gide_oscar <- oscar::oscar(x = X_xce_gide, y = RESP_gide, family="logistic")
PFS_cus_gide_oscar <- oscar::oscar(x = X_cus_gide, y = PFS_gide, family="cox")
OS_cus_gide_oscar <- oscar::oscar(x = X_cus_gide, y = OS_gide, family="cox")
RESP_cus_gide_oscar <- oscar::oscar(x = X_cus_gide, y = RESP_gide, family="logistic")
PFS_hal_gide_oscar <- oscar::oscar(x = X_hal_gide, y = PFS_gide, family="cox")
OS_hal_gide_oscar <- oscar::oscar(x = X_hal_gide, y = OS_gide, family="cox")
RESP_hal_gide_oscar <- oscar::oscar(x = X_hal_gide, y = RESP_gide, family="logistic")
PFS_bas_gide_oscar <- oscar::oscar(x = X_bas_gide, y = PFS_gide, family="cox")
OS_bas_gide_oscar <- oscar::oscar(x = X_bas_gide, y = OS_gide, family="cox")
RESP_bas_gide_oscar <- oscar::oscar(x = X_bas_gide, y = RESP_gide, family="logistic")

save.image("temprun_oscar.RData")

# Braun et al. Nivo  
PFS_xce_braun_nivo_oscar <- oscar::oscar(x = X_xce_braun_nivo, y = PFS_braun_nivo, family="cox")
OS_xce_braun_nivo_oscar <- oscar::oscar(x = X_xce_braun_nivo, y = OS_braun_nivo, family="cox")
RESP_xce_braun_nivo_oscar <- oscar::oscar(x = X_xce_braun_nivo, y = RESP_braun_nivo, family="logistic")
PFS_cus_braun_nivo_oscar <- oscar::oscar(x = X_cus_braun_nivo, y = PFS_braun_nivo, family="cox")
OS_cus_braun_nivo_oscar <- oscar::oscar(x = X_cus_braun_nivo, y = OS_braun_nivo, family="cox")
RESP_cus_braun_nivo_oscar <- oscar::oscar(x = X_cus_braun_nivo, y = RESP_braun_nivo, family="logistic")
PFS_hal_braun_nivo_oscar <- oscar::oscar(x = X_hal_braun_nivo, y = PFS_braun_nivo, family="cox")
OS_hal_braun_nivo_oscar <- oscar::oscar(x = X_hal_braun_nivo, y = OS_braun_nivo, family="cox")
RESP_hal_braun_nivo_oscar <- oscar::oscar(x = X_hal_braun_nivo, y = RESP_braun_nivo, family="logistic")
PFS_bas_braun_nivo_oscar <- oscar::oscar(x = X_bas_braun_nivo, y = PFS_braun_nivo, family="cox")
OS_bas_braun_nivo_oscar <- oscar::oscar(x = X_bas_braun_nivo, y = OS_braun_nivo, family="cox")
RESP_bas_braun_nivo_oscar <- oscar::oscar(x = X_bas_braun_nivo, y = RESP_braun_nivo, family="logistic")

save.image("temprun_oscar.RData")

# Braun et al. Ever  xCell
PFS_xce_braun_ever_oscar <- oscar::oscar(x = X_xce_braun_ever, y = PFS_braun_ever, family="cox")
OS_xce_braun_ever_oscar <- oscar::oscar(x = X_xce_braun_ever, y = OS_braun_ever, family="cox")
RESP_xce_braun_ever_oscar <- oscar::oscar(x = X_xce_braun_ever, y = RESP_braun_ever, family="logistic")
PFS_cus_braun_ever_oscar <- oscar::oscar(x = X_cus_braun_ever, y = PFS_braun_ever, family="cox")
OS_cus_braun_ever_oscar <- oscar::oscar(x = X_cus_braun_ever, y = OS_braun_ever, family="cox")
RESP_cus_braun_ever_oscar <- oscar::oscar(x = X_cus_braun_ever, y = RESP_braun_ever, family="logistic")
PFS_hal_braun_ever_oscar <- oscar::oscar(x = X_hal_braun_ever, y = PFS_braun_ever, family="cox")
OS_hal_braun_ever_oscar <- oscar::oscar(x = X_hal_braun_ever, y = OS_braun_ever, family="cox")
RESP_hal_braun_ever_oscar <- oscar::oscar(x = X_hal_braun_ever, y = RESP_braun_ever, family="logistic")
PFS_bas_braun_ever_oscar <- oscar::oscar(x = X_bas_braun_ever, y = PFS_braun_ever, family="cox")
OS_bas_braun_ever_oscar <- oscar::oscar(x = X_bas_braun_ever, y = OS_braun_ever, family="cox")
RESP_bas_braun_ever_oscar <- oscar::oscar(x = X_bas_braun_ever, y = RESP_braun_ever, family="logistic")

save.image("temprun_oscar.RData")

# Kim et al.
RESP_xce_kim_oscar <- oscar::oscar(x = X_xce_kim, y = RESP_kim, family="logistic")
RESP_cus_kim_oscar <- oscar::oscar(x = X_cus_kim, y = RESP_kim, family="logistic")
RESP_hal_kim_oscar <- oscar::oscar(x = X_hal_kim, y = RESP_kim, family="logistic")
RESP_bas_kim_oscar <- oscar::oscar(x = X_bas_kim, y = RESP_kim, family="logistic")

save.image("temprun_oscar.RData")

# Chen et al. 
RESP_xce_chen_oscar <- oscar::oscar(x = X_xce_chen, y = RESP_chen, family="logistic")
RESP_cus_chen_oscar <- oscar::oscar(x = X_cus_chen, y = RESP_chen, family="logistic")
RESP_hal_chen_oscar <- oscar::oscar(x = X_hal_chen, y = RESP_chen, family="logistic")
RESP_bas_chen_oscar <- oscar::oscar(x = X_bas_chen, y = RESP_chen, family="logistic")

save.image("temprun_oscar.RData")


### OSCAR CVS, seed 1

load("temprun_oscar.RData")

# Hugo et al. CV
OS_xce_cv_hugo_seed1 <- oscar::cv.oscar(OS_xce_hugo_oscar, fold=5, seed=1)
RESP_xce_cv_hugo_seed1 <- oscar::cv.oscar(RESP_xce_hugo_oscar, fold=5, seed=1)
OS_cus_cv_hugo_seed1 <- oscar::cv.oscar(OS_cus_hugo_oscar, fold=5, seed=1)
RESP_cus_cv_hugo_seed1 <- oscar::cv.oscar(RESP_cus_hugo_oscar, fold=5, seed=1)
OS_hal_cv_hugo_seed1 <- oscar::cv.oscar(OS_hal_hugo_oscar, fold=5, seed=1)
RESP_hal_cv_hugo_seed1 <- oscar::cv.oscar(RESP_hal_hugo_oscar, fold=5, seed=1)
OS_bas_cv_hugo_seed1 <- oscar::cv.oscar(OS_bas_hugo_oscar, fold=5, seed=1)
RESP_bas_cv_hugo_seed1 <- oscar::cv.oscar(RESP_bas_hugo_oscar, fold=5, seed=1)

# Prat et al. CV
PFS_xce_cv_prat_seed1 <- oscar::cv.oscar(PFS_xce_prat_oscar, fold=5, seed=1)
RESP_xce_cv_prat_seed1 <- oscar::cv.oscar(RESP_xce_prat_oscar, fold=5, seed=1)
PFS_cus_cv_prat_seed1 <- oscar::cv.oscar(PFS_cus_prat_oscar, fold=5, seed=1)
RESP_cus_cv_prat_seed1 <- oscar::cv.oscar(RESP_cus_prat_oscar, fold=5, seed=1)
PFS_hal_cv_prat_seed1 <- oscar::cv.oscar(PFS_hal_prat_oscar, fold=5, seed=1)
RESP_hal_cv_prat_seed1 <- oscar::cv.oscar(RESP_hal_prat_oscar, fold=5, seed=1)
PFS_bas_cv_prat_seed1 <- oscar::cv.oscar(PFS_bas_prat_oscar, fold=5, seed=1)
RESP_bas_cv_prat_seed1 <- oscar::cv.oscar(RESP_bas_prat_oscar, fold=5, seed=1)

# Westin et al. CV
PFS_xce_cv_westin_seed1 <- oscar::cv.oscar(PFS_xce_westin_oscar, fold=5, seed=1)
PFS_cus_cv_westin_seed1 <- oscar::cv.oscar(PFS_cus_westin_oscar, fold=5, seed=1)
PFS_hal_cv_westin_seed1 <- oscar::cv.oscar(PFS_hal_westin_oscar, fold=5, seed=1)
PFS_bas_cv_westin_seed1 <- oscar::cv.oscar(PFS_bas_westin_oscar, fold=5, seed=1)

# Riaz et al. CV
RESP_xce_cv_riaz_seed1 <- oscar::cv.oscar(RESP_xce_riaz_oscar, fold=5, seed=1)
RESP_cus_cv_riaz_seed1 <- oscar::cv.oscar(RESP_cus_riaz_oscar, fold=5, seed=1)
RESP_hal_cv_riaz_seed1 <- oscar::cv.oscar(RESP_hal_riaz_oscar, fold=5, seed=1)
RESP_bas_cv_riaz_seed1 <- oscar::cv.oscar(RESP_bas_riaz_oscar, fold=5, seed=1)

# Lauss et al CV
PFS_xce_cv_lauss_seed1 <- oscar::cv.oscar(PFS_xce_lauss_oscar, fold=5, seed=1)
OS_xce_cv_lauss_seed1 <- oscar::cv.oscar(OS_xce_lauss_oscar, fold=5, seed=1)
RESP_xce_cv_lauss_seed1 <- oscar::cv.oscar(RESP_xce_lauss_oscar, fold=5, seed=1)
PFS_cus_cv_lauss_seed1 <- oscar::cv.oscar(PFS_cus_lauss_oscar, fold=5, seed=1)
OS_cus_cv_lauss_seed1 <- oscar::cv.oscar(OS_cus_lauss_oscar, fold=5, seed=1)
RESP_cus_cv_lauss_seed1 <- oscar::cv.oscar(RESP_cus_lauss_oscar, fold=5, seed=1)
PFS_hal_cv_lauss_seed1 <- oscar::cv.oscar(PFS_hal_lauss_oscar, fold=5, seed=1)
OS_hal_cv_lauss_seed1 <- oscar::cv.oscar(OS_hal_lauss_oscar, fold=5, seed=1)
RESP_hal_cv_lauss_seed1 <- oscar::cv.oscar(RESP_hal_lauss_oscar, fold=5, seed=1)
PFS_bas_cv_lauss_seed1 <- oscar::cv.oscar(PFS_bas_lauss_oscar, fold=5, seed=1)
OS_bas_cv_lauss_seed1 <- oscar::cv.oscar(OS_bas_lauss_oscar, fold=5, seed=1)
RESP_bas_cv_lauss_seed1 <- oscar::cv.oscar(RESP_bas_lauss_oscar, fold=5, seed=1)

# Gide et al. CV
PFS_xce_cv_gide_seed1 <- oscar::cv.oscar(PFS_xce_gide_oscar, fold=5, seed=1)
OS_xce_cv_gide_seed1 <- oscar::cv.oscar(OS_xce_gide_oscar, fold=5, seed=1)
RESP_xce_cv_gide_seed1 <- oscar::cv.oscar(RESP_xce_gide_oscar, fold=5, seed=1)
PFS_cus_cv_gide_seed1 <- oscar::cv.oscar(PFS_cus_gide_oscar, fold=5, seed=1)
OS_cus_cv_gide_seed1 <- oscar::cv.oscar(OS_cus_gide_oscar, fold=5, seed=1)
RESP_cus_cv_gide_seed1 <- oscar::cv.oscar(RESP_cus_gide_oscar, fold=5, seed=1)
PFS_hal_cv_gide_seed1 <- oscar::cv.oscar(PFS_hal_gide_oscar, fold=5, seed=1)
OS_hal_cv_gide_seed1 <- oscar::cv.oscar(OS_hal_gide_oscar, fold=5, seed=1)
RESP_hal_cv_gide_seed1 <- oscar::cv.oscar(RESP_hal_gide_oscar, fold=5, seed=1)
PFS_bas_cv_gide_seed1 <- oscar::cv.oscar(PFS_bas_gide_oscar, fold=5, seed=1)
OS_bas_cv_gide_seed1 <- oscar::cv.oscar(OS_bas_gide_oscar, fold=5, seed=1)
RESP_bas_cv_gide_seed1 <- oscar::cv.oscar(RESP_bas_gide_oscar, fold=5, seed=1)

# Kim et al. CV
RESP_xce_cv_kim_seed1 <- oscar::cv.oscar(RESP_xce_kim_oscar, fold=5, seed=1)
RESP_cus_cv_kim_seed1 <- oscar::cv.oscar(RESP_cus_kim_oscar, fold=5, seed=1)
RESP_hal_cv_kim_seed1 <- oscar::cv.oscar(RESP_hal_kim_oscar, fold=5, seed=1)
RESP_bas_cv_kim_seed1 <- oscar::cv.oscar(RESP_bas_kim_oscar, fold=5, seed=1)

# Chen et al. CV
RESP_xce_cv_chen_seed1 <- oscar::cv.oscar(RESP_xce_chen_oscar, fold=5, seed=1)
RESP_cus_cv_chen_seed1 <- oscar::cv.oscar(RESP_cus_chen_oscar, fold=5, seed=1)
RESP_hal_cv_chen_seed1 <- oscar::cv.oscar(RESP_hal_chen_oscar, fold=5, seed=1)
RESP_bas_cv_chen_seed1 <- oscar::cv.oscar(RESP_bas_chen_oscar, fold=5, seed=1)

# Braun et al. Nivo CV
PFS_xce_cv_braun_nivo_seed1 <- oscar::cv.oscar(PFS_xce_braun_nivo_oscar, fold=5, seed=1)
OS_xce_cv_braun_nivo_seed1 <- oscar::cv.oscar(OS_xce_braun_nivo_oscar, fold=5, seed=1)
RESP_xce_cv_braun_nivo_seed1 <- oscar::cv.oscar(RESP_xce_braun_nivo_oscar, fold=5, seed=1)
PFS_cus_cv_braun_nivo_seed1 <- oscar::cv.oscar(PFS_cus_braun_nivo_oscar, fold=5, seed=1)
OS_cus_cv_braun_nivo_seed1 <- oscar::cv.oscar(OS_cus_braun_nivo_oscar, fold=5, seed=1)
RESP_cus_cv_braun_nivo_seed1 <- oscar::cv.oscar(RESP_cus_braun_nivo_oscar, fold=5, seed=1)
PFS_hal_cv_braun_nivo_seed1 <- oscar::cv.oscar(PFS_hal_braun_nivo_oscar, fold=5, seed=1)
OS_hal_cv_braun_nivo_seed1 <- oscar::cv.oscar(OS_hal_braun_nivo_oscar, fold=5, seed=1)
RESP_hal_cv_braun_nivo_seed1 <- oscar::cv.oscar(RESP_hal_braun_nivo_oscar, fold=5, seed=1)
PFS_bas_cv_braun_nivo_seed1 <- oscar::cv.oscar(PFS_bas_braun_nivo_oscar, fold=5, seed=1)
OS_bas_cv_braun_nivo_seed1 <- oscar::cv.oscar(OS_bas_braun_nivo_oscar, fold=5, seed=1)
RESP_bas_cv_braun_nivo_seed1 <- oscar::cv.oscar(RESP_bas_braun_nivo_oscar, fold=5, seed=1)

# Braun et al. Chemo CV
PFS_xce_cv_braun_ever_seed1 <- oscar::cv.oscar(PFS_xce_braun_ever_oscar, fold=5, seed=1)
OS_xce_cv_braun_ever_seed1 <- oscar::cv.oscar(OS_xce_braun_ever_oscar, fold=5, seed=1)
RESP_xce_cv_braun_ever_seed1 <- oscar::cv.oscar(RESP_xce_braun_ever_oscar, fold=5, seed=1)
PFS_cus_cv_braun_ever_seed1 <- oscar::cv.oscar(PFS_cus_braun_ever_oscar, fold=5, seed=1)
OS_cus_cv_braun_ever_seed1 <- oscar::cv.oscar(OS_cus_braun_ever_oscar, fold=5, seed=1)
RESP_cus_cv_braun_ever_seed1 <- oscar::cv.oscar(RESP_cus_braun_ever_oscar, fold=5, seed=1)
PFS_hal_cv_braun_ever_seed1 <- oscar::cv.oscar(PFS_hal_braun_ever_oscar, fold=5, seed=1)
OS_hal_cv_braun_ever_seed1 <- oscar::cv.oscar(OS_hal_braun_ever_oscar, fold=5, seed=1)
RESP_hal_cv_braun_ever_seed1 <- oscar::cv.oscar(RESP_hal_braun_ever_oscar, fold=5, seed=1)
PFS_bas_cv_braun_ever_seed1 <- oscar::cv.oscar(PFS_bas_braun_ever_oscar, fold=5, seed=1)
OS_bas_cv_braun_ever_seed1 <- oscar::cv.oscar(OS_bas_braun_ever_oscar, fold=5, seed=1)
RESP_bas_cv_braun_ever_seed1 <- oscar::cv.oscar(RESP_bas_braun_ever_oscar, fold=5, seed=1)

# TCGA CV
PFS_xce_cv_tcga_seed1 <- oscar::cv.oscar(PFS_xce_tcga_oscar, fold=5, seed=1)
OS_xce_cv_tcga_seed1 <- oscar::cv.oscar(OS_xce_tcga_oscar, fold=5, seed=1)
RESP_xce_cv_tcga_seed1 <- oscar::cv.oscar(RESP_xce_tcga_oscar, fold=5, seed=1)
PFS_cus_cv_tcga_seed1 <- oscar::cv.oscar(PFS_cus_tcga_oscar, fold=5, seed=1)
OS_cus_cv_tcga_seed1 <- oscar::cv.oscar(OS_cus_tcga_oscar, fold=5, seed=1)
RESP_cus_cv_tcga_seed1 <- oscar::cv.oscar(RESP_cus_tcga_oscar, fold=5, seed=1)
PFS_hal_cv_tcga_seed1 <- oscar::cv.oscar(PFS_hal_tcga_oscar, fold=5, seed=1)
OS_hal_cv_tcga_seed1 <- oscar::cv.oscar(OS_hal_tcga_oscar, fold=5, seed=1)
RESP_hal_cv_tcga_seed1 <- oscar::cv.oscar(RESP_hal_tcga_oscar, fold=5, seed=1)
PFS_bas_cv_tcga_seed1 <- oscar::cv.oscar(PFS_bas_tcga_oscar, fold=5, seed=1)
OS_bas_cv_tcga_seed1 <- oscar::cv.oscar(OS_bas_tcga_oscar, fold=5, seed=1)
RESP_bas_cv_tcga_seed1 <- oscar::cv.oscar(RESP_bas_tcga_oscar, fold=5, seed=1)



save.image("temprun_cv_seed1.RData")

### OSCAR CVS, seed 2

load("temprun_oscar.RData")

# Hugo et al. CV
OS_xce_cv_hugo_seed2 <- oscar::cv.oscar(OS_xce_hugo_oscar, fold=5, seed=2)
RESP_xce_cv_hugo_seed2 <- oscar::cv.oscar(RESP_xce_hugo_oscar, fold=5, seed=2)
OS_cus_cv_hugo_seed2 <- oscar::cv.oscar(OS_cus_hugo_oscar, fold=5, seed=2)
RESP_cus_cv_hugo_seed2 <- oscar::cv.oscar(RESP_cus_hugo_oscar, fold=5, seed=2)
OS_hal_cv_hugo_seed2 <- oscar::cv.oscar(OS_hal_hugo_oscar, fold=5, seed=2)
RESP_hal_cv_hugo_seed2 <- oscar::cv.oscar(RESP_hal_hugo_oscar, fold=5, seed=2)
OS_bas_cv_hugo_seed2 <- oscar::cv.oscar(OS_bas_hugo_oscar, fold=5, seed=2)
RESP_bas_cv_hugo_seed2 <- oscar::cv.oscar(RESP_bas_hugo_oscar, fold=5, seed=2)

# Prat et al. CV
PFS_xce_cv_prat_seed2 <- oscar::cv.oscar(PFS_xce_prat_oscar, fold=5, seed=2)
RESP_xce_cv_prat_seed2 <- oscar::cv.oscar(RESP_xce_prat_oscar, fold=5, seed=2)
PFS_cus_cv_prat_seed2 <- oscar::cv.oscar(PFS_cus_prat_oscar, fold=5, seed=2)
RESP_cus_cv_prat_seed2 <- oscar::cv.oscar(RESP_cus_prat_oscar, fold=5, seed=2)
PFS_hal_cv_prat_seed2 <- oscar::cv.oscar(PFS_hal_prat_oscar, fold=5, seed=2)
RESP_hal_cv_prat_seed2 <- oscar::cv.oscar(RESP_hal_prat_oscar, fold=5, seed=2)
PFS_bas_cv_prat_seed2 <- oscar::cv.oscar(PFS_bas_prat_oscar, fold=5, seed=2)
RESP_bas_cv_prat_seed2 <- oscar::cv.oscar(RESP_bas_prat_oscar, fold=5, seed=2)

# Westin et al. CV
PFS_xce_cv_westin_seed2 <- oscar::cv.oscar(PFS_xce_westin_oscar, fold=5, seed=2)
PFS_cus_cv_westin_seed2 <- oscar::cv.oscar(PFS_cus_westin_oscar, fold=5, seed=2)
PFS_hal_cv_westin_seed2 <- oscar::cv.oscar(PFS_hal_westin_oscar, fold=5, seed=2)
PFS_bas_cv_westin_seed2 <- oscar::cv.oscar(PFS_bas_westin_oscar, fold=5, seed=2)

# Riaz et al. CV
RESP_xce_cv_riaz_seed2 <- oscar::cv.oscar(RESP_xce_riaz_oscar, fold=5, seed=2)
RESP_cus_cv_riaz_seed2 <- oscar::cv.oscar(RESP_cus_riaz_oscar, fold=5, seed=2)
RESP_hal_cv_riaz_seed2 <- oscar::cv.oscar(RESP_hal_riaz_oscar, fold=5, seed=2)
RESP_bas_cv_riaz_seed2 <- oscar::cv.oscar(RESP_bas_riaz_oscar, fold=5, seed=2)

# Lauss et al CV
PFS_xce_cv_lauss_seed2 <- oscar::cv.oscar(PFS_xce_lauss_oscar, fold=5, seed=2)
OS_xce_cv_lauss_seed2 <- oscar::cv.oscar(OS_xce_lauss_oscar, fold=5, seed=2)
RESP_xce_cv_lauss_seed2 <- oscar::cv.oscar(RESP_xce_lauss_oscar, fold=5, seed=2)
PFS_cus_cv_lauss_seed2 <- oscar::cv.oscar(PFS_cus_lauss_oscar, fold=5, seed=2)
OS_cus_cv_lauss_seed2 <- oscar::cv.oscar(OS_cus_lauss_oscar, fold=5, seed=2)
RESP_cus_cv_lauss_seed2 <- oscar::cv.oscar(RESP_cus_lauss_oscar, fold=5, seed=2)
PFS_hal_cv_lauss_seed2 <- oscar::cv.oscar(PFS_hal_lauss_oscar, fold=5, seed=2)
OS_hal_cv_lauss_seed2 <- oscar::cv.oscar(OS_hal_lauss_oscar, fold=5, seed=2)
RESP_hal_cv_lauss_seed2 <- oscar::cv.oscar(RESP_hal_lauss_oscar, fold=5, seed=2)
PFS_bas_cv_lauss_seed2 <- oscar::cv.oscar(PFS_bas_lauss_oscar, fold=5, seed=2)
OS_bas_cv_lauss_seed2 <- oscar::cv.oscar(OS_bas_lauss_oscar, fold=5, seed=2)
RESP_bas_cv_lauss_seed2 <- oscar::cv.oscar(RESP_bas_lauss_oscar, fold=5, seed=2)

# Gide et al. CV
PFS_xce_cv_gide_seed2 <- oscar::cv.oscar(PFS_xce_gide_oscar, fold=5, seed=2)
OS_xce_cv_gide_seed2 <- oscar::cv.oscar(OS_xce_gide_oscar, fold=5, seed=2)
RESP_xce_cv_gide_seed2 <- oscar::cv.oscar(RESP_xce_gide_oscar, fold=5, seed=2)
PFS_cus_cv_gide_seed2 <- oscar::cv.oscar(PFS_cus_gide_oscar, fold=5, seed=2)
OS_cus_cv_gide_seed2 <- oscar::cv.oscar(OS_cus_gide_oscar, fold=5, seed=2)
RESP_cus_cv_gide_seed2 <- oscar::cv.oscar(RESP_cus_gide_oscar, fold=5, seed=2)
PFS_hal_cv_gide_seed2 <- oscar::cv.oscar(PFS_hal_gide_oscar, fold=5, seed=2)
OS_hal_cv_gide_seed2 <- oscar::cv.oscar(OS_hal_gide_oscar, fold=5, seed=2)
RESP_hal_cv_gide_seed2 <- oscar::cv.oscar(RESP_hal_gide_oscar, fold=5, seed=2)
PFS_bas_cv_gide_seed2 <- oscar::cv.oscar(PFS_bas_gide_oscar, fold=5, seed=2)
OS_bas_cv_gide_seed2 <- oscar::cv.oscar(OS_bas_gide_oscar, fold=5, seed=2)
RESP_bas_cv_gide_seed2 <- oscar::cv.oscar(RESP_bas_gide_oscar, fold=5, seed=2)

# Kim et al. CV
RESP_xce_cv_kim_seed2 <- oscar::cv.oscar(RESP_xce_kim_oscar, fold=5, seed=2)
RESP_cus_cv_kim_seed2 <- oscar::cv.oscar(RESP_cus_kim_oscar, fold=5, seed=2)
RESP_hal_cv_kim_seed2 <- oscar::cv.oscar(RESP_hal_kim_oscar, fold=5, seed=2)
RESP_bas_cv_kim_seed2 <- oscar::cv.oscar(RESP_bas_kim_oscar, fold=5, seed=2)

# Chen et al. CV
RESP_xce_cv_chen_seed2 <- oscar::cv.oscar(RESP_xce_chen_oscar, fold=5, seed=2)
RESP_cus_cv_chen_seed2 <- oscar::cv.oscar(RESP_cus_chen_oscar, fold=5, seed=2)
RESP_hal_cv_chen_seed2 <- oscar::cv.oscar(RESP_hal_chen_oscar, fold=5, seed=2)
RESP_bas_cv_chen_seed2 <- oscar::cv.oscar(RESP_bas_chen_oscar, fold=5, seed=2)

# Braun et al. Nivo CV
PFS_xce_cv_braun_nivo_seed2 <- oscar::cv.oscar(PFS_xce_braun_nivo_oscar, fold=5, seed=2)
OS_xce_cv_braun_nivo_seed2 <- oscar::cv.oscar(OS_xce_braun_nivo_oscar, fold=5, seed=2)
RESP_xce_cv_braun_nivo_seed2 <- oscar::cv.oscar(RESP_xce_braun_nivo_oscar, fold=5, seed=2)
PFS_cus_cv_braun_nivo_seed2 <- oscar::cv.oscar(PFS_cus_braun_nivo_oscar, fold=5, seed=2)
OS_cus_cv_braun_nivo_seed2 <- oscar::cv.oscar(OS_cus_braun_nivo_oscar, fold=5, seed=2)
RESP_cus_cv_braun_nivo_seed2 <- oscar::cv.oscar(RESP_cus_braun_nivo_oscar, fold=5, seed=2)
PFS_hal_cv_braun_nivo_seed2 <- oscar::cv.oscar(PFS_hal_braun_nivo_oscar, fold=5, seed=2)
OS_hal_cv_braun_nivo_seed2 <- oscar::cv.oscar(OS_hal_braun_nivo_oscar, fold=5, seed=2)
RESP_hal_cv_braun_nivo_seed2 <- oscar::cv.oscar(RESP_hal_braun_nivo_oscar, fold=5, seed=2)
PFS_bas_cv_braun_nivo_seed2 <- oscar::cv.oscar(PFS_bas_braun_nivo_oscar, fold=5, seed=2)
OS_bas_cv_braun_nivo_seed2 <- oscar::cv.oscar(OS_bas_braun_nivo_oscar, fold=5, seed=2)
RESP_bas_cv_braun_nivo_seed2 <- oscar::cv.oscar(RESP_bas_braun_nivo_oscar, fold=5, seed=2)

# Braun et al. Chemo CV
PFS_xce_cv_braun_ever_seed2 <- oscar::cv.oscar(PFS_xce_braun_ever_oscar, fold=5, seed=2)
OS_xce_cv_braun_ever_seed2 <- oscar::cv.oscar(OS_xce_braun_ever_oscar, fold=5, seed=2)
RESP_xce_cv_braun_ever_seed2 <- oscar::cv.oscar(RESP_xce_braun_ever_oscar, fold=5, seed=2)
PFS_cus_cv_braun_ever_seed2 <- oscar::cv.oscar(PFS_cus_braun_ever_oscar, fold=5, seed=2)
OS_cus_cv_braun_ever_seed2 <- oscar::cv.oscar(OS_cus_braun_ever_oscar, fold=5, seed=2)
RESP_cus_cv_braun_ever_seed2 <- oscar::cv.oscar(RESP_cus_braun_ever_oscar, fold=5, seed=2)
PFS_hal_cv_braun_ever_seed2 <- oscar::cv.oscar(PFS_hal_braun_ever_oscar, fold=5, seed=2)
OS_hal_cv_braun_ever_seed2 <- oscar::cv.oscar(OS_hal_braun_ever_oscar, fold=5, seed=2)
RESP_hal_cv_braun_ever_seed2 <- oscar::cv.oscar(RESP_hal_braun_ever_oscar, fold=5, seed=2)
PFS_bas_cv_braun_ever_seed2 <- oscar::cv.oscar(PFS_bas_braun_ever_oscar, fold=5, seed=2)
OS_bas_cv_braun_ever_seed2 <- oscar::cv.oscar(OS_bas_braun_ever_oscar, fold=5, seed=2)
RESP_bas_cv_braun_ever_seed2 <- oscar::cv.oscar(RESP_bas_braun_ever_oscar, fold=5, seed=2)

# TCGA CV
PFS_xce_cv_tcga_seed2 <- oscar::cv.oscar(PFS_xce_tcga_oscar, fold=5, seed=2)
OS_xce_cv_tcga_seed2 <- oscar::cv.oscar(OS_xce_tcga_oscar, fold=5, seed=2)
RESP_xce_cv_tcga_seed2 <- oscar::cv.oscar(RESP_xce_tcga_oscar, fold=5, seed=2)
PFS_cus_cv_tcga_seed2 <- oscar::cv.oscar(PFS_cus_tcga_oscar, fold=5, seed=2)
OS_cus_cv_tcga_seed2 <- oscar::cv.oscar(OS_cus_tcga_oscar, fold=5, seed=2)
RESP_cus_cv_tcga_seed2 <- oscar::cv.oscar(RESP_cus_tcga_oscar, fold=5, seed=2)
PFS_hal_cv_tcga_seed2 <- oscar::cv.oscar(PFS_hal_tcga_oscar, fold=5, seed=2)
OS_hal_cv_tcga_seed2 <- oscar::cv.oscar(OS_hal_tcga_oscar, fold=5, seed=2)
RESP_hal_cv_tcga_seed2 <- oscar::cv.oscar(RESP_hal_tcga_oscar, fold=5, seed=2)
PFS_bas_cv_tcga_seed2 <- oscar::cv.oscar(PFS_bas_tcga_oscar, fold=5, seed=2)
OS_bas_cv_tcga_seed2 <- oscar::cv.oscar(OS_bas_tcga_oscar, fold=5, seed=2)
RESP_bas_cv_tcga_seed2 <- oscar::cv.oscar(RESP_bas_tcga_oscar, fold=5, seed=2)



save.image("temprun_cv_seed2.RData")

### OSCAR CVS, seed 3

load("temprun_oscar.RData")

# Hugo et al. CV
OS_xce_cv_hugo_seed3 <- oscar::cv.oscar(OS_xce_hugo_oscar, fold=5, seed=3)
RESP_xce_cv_hugo_seed3 <- oscar::cv.oscar(RESP_xce_hugo_oscar, fold=5, seed=3)
OS_cus_cv_hugo_seed3 <- oscar::cv.oscar(OS_cus_hugo_oscar, fold=5, seed=3)
RESP_cus_cv_hugo_seed3 <- oscar::cv.oscar(RESP_cus_hugo_oscar, fold=5, seed=3)
OS_hal_cv_hugo_seed3 <- oscar::cv.oscar(OS_hal_hugo_oscar, fold=5, seed=3)
RESP_hal_cv_hugo_seed3 <- oscar::cv.oscar(RESP_hal_hugo_oscar, fold=5, seed=3)
OS_bas_cv_hugo_seed3 <- oscar::cv.oscar(OS_bas_hugo_oscar, fold=5, seed=3)
RESP_bas_cv_hugo_seed3 <- oscar::cv.oscar(RESP_bas_hugo_oscar, fold=5, seed=3)

# Prat et al. CV
PFS_xce_cv_prat_seed3 <- oscar::cv.oscar(PFS_xce_prat_oscar, fold=5, seed=3)
RESP_xce_cv_prat_seed3 <- oscar::cv.oscar(RESP_xce_prat_oscar, fold=5, seed=3)
PFS_cus_cv_prat_seed3 <- oscar::cv.oscar(PFS_cus_prat_oscar, fold=5, seed=3)
RESP_cus_cv_prat_seed3 <- oscar::cv.oscar(RESP_cus_prat_oscar, fold=5, seed=3)
PFS_hal_cv_prat_seed3 <- oscar::cv.oscar(PFS_hal_prat_oscar, fold=5, seed=3)
RESP_hal_cv_prat_seed3 <- oscar::cv.oscar(RESP_hal_prat_oscar, fold=5, seed=3)
PFS_bas_cv_prat_seed3 <- oscar::cv.oscar(PFS_bas_prat_oscar, fold=5, seed=3)
RESP_bas_cv_prat_seed3 <- oscar::cv.oscar(RESP_bas_prat_oscar, fold=5, seed=3)

# Westin et al. CV
PFS_xce_cv_westin_seed3 <- oscar::cv.oscar(PFS_xce_westin_oscar, fold=5, seed=3)
PFS_cus_cv_westin_seed3 <- oscar::cv.oscar(PFS_cus_westin_oscar, fold=5, seed=3)
PFS_hal_cv_westin_seed3 <- oscar::cv.oscar(PFS_hal_westin_oscar, fold=5, seed=3)
PFS_bas_cv_westin_seed3 <- oscar::cv.oscar(PFS_bas_westin_oscar, fold=5, seed=3)

# Riaz et al. CV
RESP_xce_cv_riaz_seed3 <- oscar::cv.oscar(RESP_xce_riaz_oscar, fold=5, seed=3)
RESP_cus_cv_riaz_seed3 <- oscar::cv.oscar(RESP_cus_riaz_oscar, fold=5, seed=3)
RESP_hal_cv_riaz_seed3 <- oscar::cv.oscar(RESP_hal_riaz_oscar, fold=5, seed=3)
RESP_bas_cv_riaz_seed3 <- oscar::cv.oscar(RESP_bas_riaz_oscar, fold=5, seed=3)

# Lauss et al CV
PFS_xce_cv_lauss_seed3 <- oscar::cv.oscar(PFS_xce_lauss_oscar, fold=5, seed=3)
OS_xce_cv_lauss_seed3 <- oscar::cv.oscar(OS_xce_lauss_oscar, fold=5, seed=3)
RESP_xce_cv_lauss_seed3 <- oscar::cv.oscar(RESP_xce_lauss_oscar, fold=5, seed=3)
PFS_cus_cv_lauss_seed3 <- oscar::cv.oscar(PFS_cus_lauss_oscar, fold=5, seed=3)
OS_cus_cv_lauss_seed3 <- oscar::cv.oscar(OS_cus_lauss_oscar, fold=5, seed=3)
RESP_cus_cv_lauss_seed3 <- oscar::cv.oscar(RESP_cus_lauss_oscar, fold=5, seed=3)
PFS_hal_cv_lauss_seed3 <- oscar::cv.oscar(PFS_hal_lauss_oscar, fold=5, seed=3)
OS_hal_cv_lauss_seed3 <- oscar::cv.oscar(OS_hal_lauss_oscar, fold=5, seed=3)
RESP_hal_cv_lauss_seed3 <- oscar::cv.oscar(RESP_hal_lauss_oscar, fold=5, seed=3)
PFS_bas_cv_lauss_seed3 <- oscar::cv.oscar(PFS_bas_lauss_oscar, fold=5, seed=3)
OS_bas_cv_lauss_seed3 <- oscar::cv.oscar(OS_bas_lauss_oscar, fold=5, seed=3)
RESP_bas_cv_lauss_seed3 <- oscar::cv.oscar(RESP_bas_lauss_oscar, fold=5, seed=3)

# Gide et al. CV
PFS_xce_cv_gide_seed3 <- oscar::cv.oscar(PFS_xce_gide_oscar, fold=5, seed=3)
OS_xce_cv_gide_seed3 <- oscar::cv.oscar(OS_xce_gide_oscar, fold=5, seed=3)
RESP_xce_cv_gide_seed3 <- oscar::cv.oscar(RESP_xce_gide_oscar, fold=5, seed=3)
PFS_cus_cv_gide_seed3 <- oscar::cv.oscar(PFS_cus_gide_oscar, fold=5, seed=3)
OS_cus_cv_gide_seed3 <- oscar::cv.oscar(OS_cus_gide_oscar, fold=5, seed=3)
RESP_cus_cv_gide_seed3 <- oscar::cv.oscar(RESP_cus_gide_oscar, fold=5, seed=3)
PFS_hal_cv_gide_seed3 <- oscar::cv.oscar(PFS_hal_gide_oscar, fold=5, seed=3)
OS_hal_cv_gide_seed3 <- oscar::cv.oscar(OS_hal_gide_oscar, fold=5, seed=3)
RESP_hal_cv_gide_seed3 <- oscar::cv.oscar(RESP_hal_gide_oscar, fold=5, seed=3)
PFS_bas_cv_gide_seed3 <- oscar::cv.oscar(PFS_bas_gide_oscar, fold=5, seed=3)
OS_bas_cv_gide_seed3 <- oscar::cv.oscar(OS_bas_gide_oscar, fold=5, seed=3)
RESP_bas_cv_gide_seed3 <- oscar::cv.oscar(RESP_bas_gide_oscar, fold=5, seed=3)

# Kim et al. CV
RESP_xce_cv_kim_seed3 <- oscar::cv.oscar(RESP_xce_kim_oscar, fold=5, seed=3)
RESP_cus_cv_kim_seed3 <- oscar::cv.oscar(RESP_cus_kim_oscar, fold=5, seed=3)
RESP_hal_cv_kim_seed3 <- oscar::cv.oscar(RESP_hal_kim_oscar, fold=5, seed=3)
RESP_bas_cv_kim_seed3 <- oscar::cv.oscar(RESP_bas_kim_oscar, fold=5, seed=3)

# Chen et al. CV
RESP_xce_cv_chen_seed3 <- oscar::cv.oscar(RESP_xce_chen_oscar, fold=5, seed=3)
RESP_cus_cv_chen_seed3 <- oscar::cv.oscar(RESP_cus_chen_oscar, fold=5, seed=3)
RESP_hal_cv_chen_seed3 <- oscar::cv.oscar(RESP_hal_chen_oscar, fold=5, seed=3)
RESP_bas_cv_chen_seed3 <- oscar::cv.oscar(RESP_bas_chen_oscar, fold=5, seed=3)

# Braun et al. Nivo CV
PFS_xce_cv_braun_nivo_seed3 <- oscar::cv.oscar(PFS_xce_braun_nivo_oscar, fold=5, seed=3)
OS_xce_cv_braun_nivo_seed3 <- oscar::cv.oscar(OS_xce_braun_nivo_oscar, fold=5, seed=3)
RESP_xce_cv_braun_nivo_seed3 <- oscar::cv.oscar(RESP_xce_braun_nivo_oscar, fold=5, seed=3)
PFS_cus_cv_braun_nivo_seed3 <- oscar::cv.oscar(PFS_cus_braun_nivo_oscar, fold=5, seed=3)
OS_cus_cv_braun_nivo_seed3 <- oscar::cv.oscar(OS_cus_braun_nivo_oscar, fold=5, seed=3)
RESP_cus_cv_braun_nivo_seed3 <- oscar::cv.oscar(RESP_cus_braun_nivo_oscar, fold=5, seed=3)
PFS_hal_cv_braun_nivo_seed3 <- oscar::cv.oscar(PFS_hal_braun_nivo_oscar, fold=5, seed=3)
OS_hal_cv_braun_nivo_seed3 <- oscar::cv.oscar(OS_hal_braun_nivo_oscar, fold=5, seed=3)
RESP_hal_cv_braun_nivo_seed3 <- oscar::cv.oscar(RESP_hal_braun_nivo_oscar, fold=5, seed=3)
PFS_bas_cv_braun_nivo_seed3 <- oscar::cv.oscar(PFS_bas_braun_nivo_oscar, fold=5, seed=3)
OS_bas_cv_braun_nivo_seed3 <- oscar::cv.oscar(OS_bas_braun_nivo_oscar, fold=5, seed=3)
RESP_bas_cv_braun_nivo_seed3 <- oscar::cv.oscar(RESP_bas_braun_nivo_oscar, fold=5, seed=3)

# Braun et al. Chemo CV
PFS_xce_cv_braun_ever_seed3 <- oscar::cv.oscar(PFS_xce_braun_ever_oscar, fold=5, seed=3)
OS_xce_cv_braun_ever_seed3 <- oscar::cv.oscar(OS_xce_braun_ever_oscar, fold=5, seed=3)
RESP_xce_cv_braun_ever_seed3 <- oscar::cv.oscar(RESP_xce_braun_ever_oscar, fold=5, seed=3)
PFS_cus_cv_braun_ever_seed3 <- oscar::cv.oscar(PFS_cus_braun_ever_oscar, fold=5, seed=3)
OS_cus_cv_braun_ever_seed3 <- oscar::cv.oscar(OS_cus_braun_ever_oscar, fold=5, seed=3)
RESP_cus_cv_braun_ever_seed3 <- oscar::cv.oscar(RESP_cus_braun_ever_oscar, fold=5, seed=3)
PFS_hal_cv_braun_ever_seed3 <- oscar::cv.oscar(PFS_hal_braun_ever_oscar, fold=5, seed=3)
OS_hal_cv_braun_ever_seed3 <- oscar::cv.oscar(OS_hal_braun_ever_oscar, fold=5, seed=3)
RESP_hal_cv_braun_ever_seed3 <- oscar::cv.oscar(RESP_hal_braun_ever_oscar, fold=5, seed=3)
PFS_bas_cv_braun_ever_seed3 <- oscar::cv.oscar(PFS_bas_braun_ever_oscar, fold=5, seed=3)
OS_bas_cv_braun_ever_seed3 <- oscar::cv.oscar(OS_bas_braun_ever_oscar, fold=5, seed=3)
RESP_bas_cv_braun_ever_seed3 <- oscar::cv.oscar(RESP_bas_braun_ever_oscar, fold=5, seed=3)

# TCGA CV
PFS_xce_cv_tcga_seed3 <- oscar::cv.oscar(PFS_xce_tcga_oscar, fold=5, seed=3)
OS_xce_cv_tcga_seed3 <- oscar::cv.oscar(OS_xce_tcga_oscar, fold=5, seed=3)
RESP_xce_cv_tcga_seed3 <- oscar::cv.oscar(RESP_xce_tcga_oscar, fold=5, seed=3)
PFS_cus_cv_tcga_seed3 <- oscar::cv.oscar(PFS_cus_tcga_oscar, fold=5, seed=3)
OS_cus_cv_tcga_seed3 <- oscar::cv.oscar(OS_cus_tcga_oscar, fold=5, seed=3)
RESP_cus_cv_tcga_seed3 <- oscar::cv.oscar(RESP_cus_tcga_oscar, fold=5, seed=3)
PFS_hal_cv_tcga_seed3 <- oscar::cv.oscar(PFS_hal_tcga_oscar, fold=5, seed=3)
OS_hal_cv_tcga_seed3 <- oscar::cv.oscar(OS_hal_tcga_oscar, fold=5, seed=3)
RESP_hal_cv_tcga_seed3 <- oscar::cv.oscar(RESP_hal_tcga_oscar, fold=5, seed=3)
PFS_bas_cv_tcga_seed3 <- oscar::cv.oscar(PFS_bas_tcga_oscar, fold=5, seed=3)
OS_bas_cv_tcga_seed3 <- oscar::cv.oscar(OS_bas_tcga_oscar, fold=5, seed=3)
RESP_bas_cv_tcga_seed3 <- oscar::cv.oscar(RESP_bas_tcga_oscar, fold=5, seed=3)

save.image("temprun_cv_seed3.RData")
