# Estimate tumor mutational burden coefficients and/or X matrix

aggregateX <- function(
	gex,
	dat,
	keyGenes = unique(c(
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
	))
){
	library(GSVA)
	library(immunedeconv)
	
	# Hard-coded Carbone et al. thresholds
	# 242 mutations threshold for 'high' mutatation; impute NAs as if being not highly mutated
	X <- data.frame(TMB = as.integer(dat[,"TMB"]))
	rownames(X) <- rownames(dat)
	# Impute zero indicators if there are NA values
	if(any(is.na(X[,"TMB"]))) X[is.na(X[,"TMB"]),"TMB"] <- median(X[,"TMB"], na.rm=TRUE)
	# Arbitrary threshold, based on tertiary from Carbone et al
	X <- cbind(X, TMBhigh = X[,"TMB"]>242)
	
	# There was no significant association between tumor-mutation burden and 
	# PD-L1 expression level (Pearson’s correlation coefficient = 0.059).	
	# --> Complementary information
	
	# However, in the nivolumab group, patients with both a high tumor-mutation burden and a 
	# PD-L1 expression level of 50% or more had a higher response rate (75%) than those with 
	# only one of these factors (32% among patients with a high tumor-mutation burden only and 
	# 34% among those with a PD-L1 expression level of >50% only) or neither factor (16%).
	
	# 75% -> 3/4
	# 32%, 34 % ~ 1/3
	# 16% ~ 1/7 - 1/8
	# looks like a non-linear trend?
	
	# Thirds for quantiles
	#TMBq3s <- quantile(dat[,"TMB"], probs=c(0, 1/3, 2/3, 1), na.rm=TRUE)
	#dat[is.na(dat[,"TMB"]),"TMB"] <- median(dat[,"TMB"], na.rm=TRUE)
	#PDL1q3s <- quantile(dat[,"PDL1"], probs=c(0, 1/3, 2/3, 1), na.rm=TRUE)
	#X <- cbind(X, TMBq3s = findInterval(x=dat[,"TMB"], vec=TMBq3s, rightmost.closed = TRUE))
	#X <- cbind(X, PDL1q3s = findInterval(x=dat[,"PDL1"], vec=PDL1q3s, rightmost.closed = TRUE))
	#X <- cbind(X, TMB.PDL1 = X[,"TMBq3s"] * X[,"PDL1q3s"])
	# CD274 expression level modelled as a surrogate for PD-L1 IHC
	# Normalized expressions between various platforms and their respective distributional characteristics
	X <- cbind(X, CD274 = gex["CD274",])
	
	# GSVA for CUSTOM_IFNG3
	gmt_custom <- GSEABase::getGmt(".\\selfmade.gmt")
	res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE)) # Custom GMTs
	#X <- cbind(X, IFNG = res_gsva[,grep("CUSTOM_IFNG3", colnames(res_gsva))])
	X <- cbind(X, res_gsva)
	
	# Interaction with T-cell inflammatory signal Gene Expression Profile (GEP) and high tumor mutational burden
	#GEPq2s <- median(res_gsva[,"CUSTOM_GEP"])
	#X <- cbind(X, GEPq2s = as.integer(res_gsva[,"CUSTOM_GEP"] > GEPq2s))
	# TMB x GEP interaction, upper tertile of TMB vs. upper median of GEP signature
	#X <- cbind(X, TMB.GEP = as.integer(X[,"TMBq3s"] == 2) * X[,"GEPq2s"])
	
	# GSVA for Hallmarks
	gmt_hallmarks <- GSEABase::getGmt(".\\h.all.v7.2.symbols.gmt")
	res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_hallmarks, verbose=FALSE)) # Custom GMTs
	#X <- cbind(X, IFNG = res_gsva[,grep("CUSTOM_IFNG3", colnames(res_gsva))])
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

setwd("D:\\Gits\\DREAM_2020_IO")
#load("temp.RData")
#load("temprun_tcga.RData")
#load("tempspace_1.RData")
#load("tempspace_2.RData")
#load("tempspace_3.RData")
#load("tempspace_4.RData")
#load("tempspace_5.RData")
#load("tempspace_6.RData")
#load("tempspace_7.RData")
load("datas.RData")

# Plot various gene expression distributions to see if Cox PH coefficients would generalize between studies

png("GEXs_notransforms.png", width=800, height=800)
par(mfrow=c(3,4))
gexplot <- function(x, main=""){
	plot.new()
	par(mar=c(4,4,4,1))
	# Ranges based on trimmed quantiles
	xr <- quantile(x, probs=seq(0,1,.1))[c(2,10)]
	plot.window(xlim=xr, ylim=c(0,1))
	box(); axis(1); axis(2)
	title(xlab="GEX values (trimmed, first 100)", ylab="Kernel density", main=main)
	lapply(1:min(100, nrow(x)), FUN=function(z){
		points(density(x[z,]), type="l", col=z)
	})
}
gexplot(gex_synthetic, main="GEX Synthetic (orig.)")
gexplot(gex_tcga, main="GEX TCGA (orig.)")
gexplot(gex_chen, main="GEX Chen et al. (orig.)")
gexplot(gex_hugo, main="GEX Hugo et al. (orig.)")
gexplot(gex_kim, main="GEX Kim et al. (orig.)")
gexplot(gex_lauss, main="GEX Lauss et al. (orig.)")
gexplot(gex_hugo, main="GEX Prat et al. (orig.)")
gexplot(gex_westin, main="GEX Westin et al. (orig.)")
gexplot(gex_riaz, main="GEX Riaz et al. (orig.)")
gexplot(gex_gide, main="GEX Gide et al. (orig.)")
gexplot(gex_braun_nivo, main="GEX Braun et al. Nivo (orig.)")
gexplot(gex_braun_ever, main="GEX Braun et al. Chemo (orig.)")
dev.off()

logz <- function(x) { 
	tmp <- scale(log(x+1)) 
	if(any(!is.finite(tmp))){
		tmp[!is.finite(tmp)] <- 0
	}
	tmp
}

if(FALSE){
	png("GEXs_logzs.png", width=800, height=800)
	par(mfrow=c(3,4))
	gexplot <- function(x, main=""){
		plot.new()
		par(mar=c(4,4,4,1))
		# Ranges based on trimmed quantiles
		xr <- quantile(x, probs=seq(0,1,.1))[c(2,10)]
		plot.window(xlim=xr, ylim=c(0,1))
		box(); axis(1); axis(2)
		title(xlab="GEX values (trimmed, first 100)", ylab="Kernel density", main=main)
		lapply(1:min(100, nrow(x)), FUN=function(z){
			points(density(x[z,]), type="l", col=z)
		})
	}
	gexplot(apply(gex_synthetic, MARGIN=1, FUN=logz), main="GEX Synthetic (log-z)")
	gexplot(apply(gex_tcga, MARGIN=1, FUN=logz), main="GEX TCGA (log-z)")
	gexplot(gex_chen, main="GEX Chen et al. (orig.)")
	gexplot(apply(gex_hugo, MARGIN=1, FUN=logz), main="GEX Hugo et al. (log-z)")
	gexplot(gex_kim, main="GEX Kim et al. (orig.)")
	gexplot(gex_lauss, main="GEX Lauss et al. (orig.)")
	gexplot(gex_prat, main="GEX Prat et al. (orig.)")
	gexplot(apply(gex_westin, MARGIN=1, FUN=logz), main="GEX Westin et al. (log-z)")
	gexplot(apply(gex_riaz, MARGIN=1, FUN=logz), main="GEX Riaz et al. (log-z)")
	gexplot(gex_gide, main="GEX Gide et al. (orig.)")
	gexplot(apply(gex_braun_nivo, MARGIN=1, FUN=logz), main="GEX Braun et al. Nivo (log-z)")
	gexplot(apply(gex_braun_ever, MARGIN=1, FUN=logz), main="GEX Braun et al. Chemo (log-z)")
	dev.off()
}

# Same the log-z standardized gene expression patterns
#gxz_synthetic <- t(apply(gex_synthetic, MARGIN=1, FUN=logz))
# Normalize using quantile normalization
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

png("GEXs_quantnorm.png", width=800, height=800)
par(mfrow=c(3,4))
gexplot <- function(x, main=""){
	plot.new()
	par(mar=c(4,4,4,1))
	# Ranges based on trimmed quantiles
	xr <- quantile(x, probs=seq(0,1,.1))[c(2,10)]
	plot.window(xlim=xr, ylim=c(0,1))
	box(); axis(1); axis(2)
	title(xlab="GEX values (trimmed, first 100)", ylab="Kernel density", main=main)
	lapply(1:min(100, nrow(x)), FUN=function(z){
		points(density(x[z,]), type="l", col=z)
	})
}
gexplot(gex_synthetic, main="GEX Synthetic (orig.)")
gexplot(geq_tcga, main="GEX TCGA (quant.)")
gexplot(geq_chen, main="GEX Chen et al. (quant.)")
gexplot(geq_hugo, main="GEX Hugo et al. (quant.)")
gexplot(geq_kim, main="GEX Kim et al. (quant.)")
gexplot(geq_lauss, main="GEX Lauss et al. (quant.)")
gexplot(geq_hugo, main="GEX Prat et al. (quant.)")
gexplot(geq_westin, main="GEX Westin et al. (quant.)")
gexplot(geq_riaz, main="GEX Riaz et al. (quant.)")
gexplot(geq_gide, main="GEX Gide et al. (quant.)")
gexplot(geq_braun_nivo, main="GEX Braun et al. Nivo (quant.)")
gexplot(geq_braun_ever, main="GEX Braun et al. Chemo (quant.)")
dev.off()

if(FALSE){
	rownames(gxz_synthetic) <- rownames(gex_synthetic)
	gxz_tcga <- t(apply(gex_tcga, MARGIN=1, FUN=logz))
	rownames(gxz_tcga) <- rownames(gex_tcga)
	gxz_hugo <- t(apply(gex_hugo, MARGIN=1, FUN=logz))
	rownames(gxz_hugo) <- rownames(gex_hugo)
	gxz_westin <- t(apply(gex_westin, MARGIN=1, FUN=logz))
	rownames(gxz_westin) <- rownames(gex_westin)
	gxz_riaz <- t(apply(gex_riaz, MARGIN=1, FUN=logz))
	rownames(gxz_riaz) <- rownames(gex_riaz)
	gxz_braun_nivo <- t(apply(gex_braun_nivo, MARGIN=1, FUN=logz))
	rownames(gxz_braun_nivo) <- rownames(gex_braun_nivo)
	gxz_braun_ever <- t(apply(gex_braun_ever, MARGIN=1, FUN=logz))
	rownames(gxz_braun_ever) <- rownames(gex_braun_ever)
}

# Compute synthetic data 
Xs_synthetic <- aggregateX(gex_synthetic, dat=dat_synthetic)
# Test if logz has an effect
Xz_synthetic <- aggregateX(gxz_synthetic, dat=dat_synthetic)

# Some interesting fields for OS, interactions

OScols <- c("TMB.PDL1", "TMB.GEP", "TMBq3s", "PDL1q3s", "GEPq2s")
Xs_synthetic[,OScols]
Xz_synthetic[,OScols]

# Aggregate X matrices
# PFS:	(Prat et al. &) Lauss et al. & Westin et al.
#	+ Gide et al. & Braun et al (Nivo)
# OS:	Hugo et al. & Lauss et al.
#	+ Gide et al. & Braun et al (Nivo)
# RESP:	Hugo et a. & Prat et al. & Lauss et al. & Kim et al. & Chen et al. & Riaz et al.
#	+ Gide et al. & Braun et al (Nivo)

# Chemo-arm
Xz_tcga <- aggregateX(gex=gxz_tcga, dat=dat_tcga)
# IOs
# Xz_prat <- aggregateX(gex=gxz_prat, dat=dat_prat)
Xz_hugo <- aggregateX(gex=gxz_hugo, dat=dat_hugo)
Xz_lauss <- aggregateX(gex=gex_lauss, dat=dat_lauss)
Xz_westin <- aggregateX(gex=gxz_westin, dat=dat_westin)
Xz_kim <- aggregateX(gex=gex_kim, dat=dat_kim)
#Xz_chen <- aggregateX(gex=gex_chen, dat=dat_chen)
Xz_riaz <- aggregateX(gex=gex_riaz, dat=dat_riaz)
Xz_gide <- aggregateX(gex=gex_gide, dat=dat_gide)
# Nivo arm
Xz_braun_nivo <- aggregateX(gex=gex_braun_nivo, dat=dat_braun_nivo)
# Chemo arm
Xz_braun_ever <- aggregateX(gex=gex_braun_ever, dat=dat_braun_ever)


# Combine studies to generate aggregate HR estimates
combcols <- c("CD274", "CUSTOM_IFNG3", "xce_Endothelial_cell")

# xCell fails for Prat et al. due to too narrow gene panel
library(survival)

## NOTE!
# Responses not necessarily in same unit; e.g. Hugo et al. OS in days, Lauss in months
PFS2_westin <- PFS_westin
PFS2_westin[,1] <- round(PFS2_westin[,1]/30.5,0)

PFS2_gide <- PFS_gide
PFS2_gide[,1] <- round(PFS2_gide[,1]/30.5,0)

OS2_gide <- OS_gide
OS2_gide[,1] <- round(OS2_gide[,1]/30.5,0)

OS2_hugo <- OS_hugo
OS2_hugo[,1] <- round(OS2_hugo[,1]/30.5,0)

# Coefficient for linear combination for Carbone et al.
# for tumor high burden (based on absolute tumor mutation count)
# HR: 0.62
# -> 
#> log(0.62)
#[1] -0.4780358

######## COMBINED PROGRESSION FREE SURVIVAL ######

Xz_PFS <- rbind(Xz_lauss[,combcols], Xz_westin[,combcols], Xz_gide[,combcols], Xz_braun_nivo[,combcols])
Yz_PFS <- c(PFS_lauss, PFS2_westin, PFS2_gide, PFS_braun_nivo)
survival::coxph(Yz_PFS ~ ., data = as.data.frame(Xz_PFS))
#> survival::coxph(Yz_PFS ~ ., data = as.data.frame(Xz_PFS))
#Call:
#survival::coxph(formula = Yz_PFS ~ ., data = as.data.frame(Xz_PFS))
#
#                        coef exp(coef) se(coef)      z      p
#CD274                -0.5326    0.5871   0.2388 -2.231 0.0257
#IFNG                  0.1753    1.1917   0.5129  0.342 0.7325
#xce_Endothelial_cell  0.8812    2.4138   0.8040  1.096 0.2731
#
#Likelihood ratio test=6.25  on 3 df, p=0.1
#n= 43, number of events= 31

######## COMBINED OVERALL SURVIVAL #######

Xz_OS <- rbind(Xz_hugo[,combcols], Xz_lauss[,combcols])
Yz_OS <- c(OS2_hugo, OS_lauss)
survival::coxph(Yz_OS ~ ., data = as.data.frame(Xz_OS))
#> survival::coxph(Yz_OS ~ ., data = as.data.frame(Xz_OS))
#Call:
#survival::coxph(formula = Yz_OS ~ ., data = as.data.frame(Xz_OS))
#
#                        coef exp(coef) se(coef)      z      p
#CD274                -0.7039    0.4947   0.3685 -1.910 0.0561
#IFNG                  0.4233    1.5269   0.5263  0.804 0.4213
#xce_Endothelial_cell  1.7892    5.9849   0.7229  2.475 0.0133
#
#Likelihood ratio test=7.51  on 3 df, p=0.05726
#n= 51, number of events= 28 
#   (1 observation deleted due to missingness)

# Anni coefs
anniOS <- c("xce_Myeloid_dendritic_cell_activated", "xce_Common_lymphoid_progenitor", "xce_T_cell_regulatory_(Tregs)")

Xz_OS <- rbind(Xz_hugo[,anniOS], Xz_lauss[,anniOS])
Yz_OS <- c(OS2_hugo, OS_lauss)
survival::coxph(Yz_OS ~ ., data = as.data.frame(Xz_OS))

   
######### COMBINED RESPONSE ########

#RESPcols <- c("CUSTOM_MCP_ENDOTHELIAL.CELLS", "CUSTOM_IFNG3", "xce_Myeloid_dendritic_cell", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "CUSTOM_MCP_MYELOID.DENDRITIC.CELLS")
RESPcols <- c("CUSTOM_MCP_ENDOTHELIAL.CELLS", "CUSTOM_IFNG3", "CD274", "CUSTOM_GEP", "HALLMARK_INTERFERON_GAMMA_RESPONSE")

#Xz_RESP <- rbind(Xz_hugo[,combcols], Xz_lauss[,combcols], Xz_kim[,combcols], Xz_riaz[,combcols])
Xz_RESP <- rbind(Xz_hugo[,RESPcols], Xz_lauss[,RESPcols], Xz_kim[,RESPcols], Xz_riaz[,RESPcols], Xz_gide[,RESPcols], Xz_braun_nivo[,RESPcols])
Yz_RESP <- c(RESP_hugo, RESP_lauss, RESP_kim, RESP_riaz, RESP_gide, RESP_braun_nivo)
summary(stats::glm(Yz_RESP ~ ., data = as.data.frame(Xz_RESP), family="binomial"))
### Updated with more variable and Riaz et al.



### OLD MODEL
#> summary(stats::glm(Yz_RESP ~ ., data = as.data.frame(Xz_RESP), family="binomial"))
#
#Call:
#stats::glm(formula = Yz_RESP ~ ., family = "binomial", data = as.data.frame(Xz_RESP))
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-2.0400  -0.8642  -0.4856   0.9134   2.1618  
#
#Coefficients:
#                     Estimate Std. Error z value Pr(>|z|)   
#(Intercept)            0.1932     0.3285   0.588  0.55647   
#CD274                  0.7570     0.2996   2.527  0.01151 * 
#IFNG                   0.8468     0.5784   1.464  0.14322   
#xce_Endothelial_cell  -3.1055     1.0381  -2.992  0.00277 **
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 128.96  on 96  degrees of freedom
#Residual deviance: 102.03  on 93  degrees of freedom
#AIC: 110.03
#
#Number of Fisher Scoring iterations: 5

####### INDIVIDUAL PROGRESSION FREE SURVIVAL ########

# Braun et al Nivo vs. Chemo
survival::coxph(PFS_braun_nivo ~ ., data = as.data.frame(X_braun_nivo[,c("BASE_TMB", "BASE_CD274", "CUSTOM_GEP")]))
survival::coxph(OS_braun_nivo ~ ., data = as.data.frame(X_braun_nivo[,c("BASE_TMB", "BASE_CD274", "CUSTOM_GEP")]))

survival::coxph(PFS_braun_ever ~ ., data = as.data.frame(X_braun_ever[,c("BASE_TMB", "BASE_CD274")]))
survival::coxph(OS_braun_ever ~ ., data = as.data.frame(X_braun_ever[,c("BASE_TMB", "BASE_CD274")]))




# Individual models   
survival::coxph(PFS_lauss ~ ., data = as.data.frame(Xz_lauss))
#> survival::coxph(PFS_lauss ~ ., data = as.data.frame(Xz_lauss))
#Call:
#survival::coxph(formula = PFS_lauss ~ ., data = as.data.frame(Xz_lauss))
#
#                         coef exp(coef) se(coef)      z     p
#isTMBhigh                  NA        NA  0.00000     NA    NA
#CD274                -0.65527   0.51930  0.56371 -1.162 0.245
#IFNG                  0.03933   1.04011  0.77718  0.051 0.960
#xce_Endothelial_cell  0.51908   1.68049  0.95249  0.545 0.586
#
#Likelihood ratio test=3.4  on 3 df, p=0.3344
#n= 25, number of events= 20
survival::coxph(PFS_westin ~ ., data = as.data.frame(Xz_westin))
#> survival::coxph(PFS_westin ~ ., data = as.data.frame(Xz_westin))
#Call:
#survival::coxph(formula = PFS_westin ~ ., data = as.data.frame(Xz_westin))
#
#                        coef exp(coef) se(coef)      z      p
#isTMBhigh                 NA        NA   0.0000     NA     NA
#CD274                -1.0275    0.3579   0.4424 -2.323 0.0202
#IFNG                  1.4974    4.4700   1.0696  1.400 0.1615
#xce_Endothelial_cell -0.7019    0.4956   3.4139 -0.206 0.8371
#
#Likelihood ratio test=7.77  on 3 df, p=0.05098
#n= 18, number of events= 11

####### INDIVIDUAL OVERALL SURVIVAL ########


survival::coxph(OS_hugo ~ ., data = as.data.frame(Xz_hugo))
#> survival::coxph(OS_hugo ~ ., data = as.data.frame(Xz_hugo))
#Call:
#survival::coxph(formula = OS_hugo ~ ., data = as.data.frame(Xz_hugo))
#
#                        coef exp(coef) se(coef)      z      p
#isTMBhigh                 NA        NA   0.0000     NA     NA
#CD274                -0.3827    0.6820   0.5657 -0.677 0.4987
#IFNG                  0.7265    2.0678   0.8581  0.847 0.3972
#xce_Endothelial_cell  3.0756   21.6633   1.1843  2.597 0.0094
#
#Likelihood ratio test=6.45  on 3 df, p=0.09156
#n= 26, number of events= 11 
#   (1 observation deleted due to missingness)
survival::coxph(OS_lauss ~ ., data = as.data.frame(Xz_lauss[,c("CD274")]))
#> survival::coxph(OS_lauss ~ ., data = as.data.frame(Xz_lauss))
#Call:
#survival::coxph(formula = OS_lauss ~ ., data = as.data.frame(Xz_lauss))
#
#                        coef exp(coef) se(coef)      z       p
#isTMBhigh                 NA        NA   0.0000     NA      NA
#CD274                -1.7585    0.1723   0.6015 -2.924 0.00346
#IFNG                  0.9853    2.6785   0.7117  1.384 0.16623
#xce_Endothelial_cell  2.0898    8.0830   0.9973  2.095 0.03614
#
#Likelihood ratio test=9.76  on 3 df, p=0.02073
#n= 25, number of events= 17

summary(stats::glm(RESP_hugo ~ ., data = as.data.frame(Xz_hugo), family="binomial"))
#> summary(stats::glm(RESP_hugo ~ ., data = as.data.frame(Xz_hugo), family="binomial"))
#
#Call:
#stats::glm(formula = RESP_hugo ~ ., family = "binomial", data = as.data.frame(Xz_hugo))
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-2.0615  -1.0030   0.6423   0.8309   1.6578  
#
#Coefficients: (1 not defined because of singularities)
#                     Estimate Std. Error z value Pr(>|z|)  
#(Intercept)            1.3906     0.6803   2.044   0.0409 *
#isTMBhigh                  NA         NA      NA       NA  
#CD274                  0.4800     0.4898   0.980   0.3271  
#IFNG                  -0.0469     0.9750  -0.048   0.9616  
#xce_Endothelial_cell  -4.7414     2.2492  -2.108   0.0350 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 37.096  on 26  degrees of freedom
#Residual deviance: 30.370  on 23  degrees of freedom
#AIC: 38.37
#
#Number of Fisher Scoring iterations: 4
summary(stats::glm(RESP_lauss ~ ., data = as.data.frame(Xz_lauss), family="binomial"))
#> summary(stats::glm(RESP_lauss ~ ., data = as.data.frame(Xz_lauss), family="binomial"))
#
#Call:
#stats::glm(formula = RESP_lauss ~ ., family = "binomial", data = as.data.frame(Xz_lauss))
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.7534  -0.8164  -0.5206   0.9739   1.7780  
#
#Coefficients: (1 not defined because of singularities)
#                     Estimate Std. Error z value Pr(>|z|)
#(Intercept)           -0.3112     0.6678  -0.466    0.641
#isTMBhigh                  NA         NA      NA       NA
#CD274                  0.4949     0.9427   0.525    0.600
#IFNG                   1.6109     1.4151   1.138    0.255
#xce_Endothelial_cell  -1.0358     1.8955  -0.546    0.585
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 33.651  on 24  degrees of freedom
#Residual deviance: 27.951  on 21  degrees of freedom
#AIC: 35.951
#
#Number of Fisher Scoring iterations: 4
summary(stats::glm(RESP_kim ~ ., data = as.data.frame(Xz_kim), family="binomial"))
#> summary(stats::glm(RESP_kim ~ ., data = as.data.frame(Xz_kim), family="binomial"))
#
#Call:
#stats::glm(formula = RESP_kim ~ ., family = "binomial", data = as.data.frame(Xz_kim))
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-1.3225  -0.5395  -0.2091   0.1691   2.9275  
#
#Coefficients: (1 not defined because of singularities)
#                     Estimate Std. Error z value Pr(>|z|)  
#(Intercept)           -0.3752     0.6111  -0.614   0.5392  
#isTMBhigh                  NA         NA      NA       NA  
#CD274                  1.2036     0.5559   2.165   0.0304 *
#IFNG                   1.2933     1.1547   1.120   0.2627  
#xce_Endothelial_cell  -5.0684     2.4167  -2.097   0.0360 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 52.192  on 44  degrees of freedom
#Residual deviance: 29.530  on 41  degrees of freedom
#AIC: 37.53
#
#Number of Fisher Scoring iterations: 6

## Raw data version for Hugo et al.
survival::coxph(OS_hugo ~ ., data = as.data.frame(X_hugo[,c("BASE_CD274", "xce_Endothelial cell", "CUSTOM_IFNG3")]))
#> survival::coxph(OS_hugo ~ ., data = as.data.frame(X_hugo[,c("BASE_CD274", "xce_Endothelial cell", "CUSTOM_IFNG3")]))
#Call:
#survival::coxph(formula = OS_hugo ~ ., data = as.data.frame(X_hugo[, 
#    c("BASE_CD274", "xce_Endothelial cell", "CUSTOM_IFNG3")]))
#
#                            coef exp(coef)  se(coef)     z       p
#BASE_CD274             4.940e-02 1.051e+00 1.705e-01 0.290 0.77201
#`xce_Endothelial cell` 2.268e+01 7.102e+09 7.845e+00 2.892 0.00383
#CUSTOM_IFNG3           1.845e-01 1.203e+00 8.186e-01 0.225 0.82167
#
#Likelihood ratio test=7.67  on 3 df, p=0.05334
#n= 26, number of events= 11 
#   (1 observation deleted due to missingness)
   


#### TCGA

# Chemo-arm
survival::coxph(PFS_tcga ~ ., data = as.data.frame(Xz_tcga))
#> survival::coxph(PFS_tcga ~ ., data = as.data.frame(Xz_tcga))
#Call:
#survival::coxph(formula = PFS_tcga ~ ., data = as.data.frame(Xz_tcga))
#
#                        coef exp(coef) se(coef)      z     p
#isTMBhigh                 NA        NA   0.0000     NA    NA
#CD274                -0.1082    0.8974   0.1078 -1.004 0.315
#IFNG                  0.1260    1.1343   0.2234  0.564 0.573
#xce_Endothelial_cell  0.1018    1.1072   0.2737  0.372 0.710
#
#Likelihood ratio test=1.08  on 3 df, p=0.7822
#n= 278, number of events= 118 
#   (36 observations deleted due to missingness)
survival::coxph(OS_tcga ~ ., data = as.data.frame(Xz_tcga))
#> survival::coxph(OS_tcga ~ ., data = as.data.frame(Xz_tcga))
#Call:
#survival::coxph(formula = OS_tcga ~ ., data = as.data.frame(Xz_tcga))
#
#                        coef exp(coef) se(coef)      z        p
#isTMBhigh                 NA        NA   0.0000     NA       NA
#CD274                -0.0180    0.9822   0.2011 -0.089 0.928685
#IFNG                  0.5324    1.7030   0.3790  1.405 0.160118
#xce_Endothelial_cell  1.5172    4.5597   0.3981  3.811 0.000138
#
#Likelihood ratio test=17.03  on 3 df, p=0.0006984
#n= 311, number of events= 41 
#   (3 observations deleted due to missingness)
summary(stats::glm(RESP_tcga ~ ., data = as.data.frame(Xz_tcga), family="binomial"))
### NEW VARS
summary(stats::glm(RESP_tcga ~ ., data = as.data.frame(Xz_tcga[,RESPcols]), family="binomial"))
#> summary(stats::glm(RESP_tcga ~ ., data = as.data.frame(Xz_tcga[,RESPcols]), family="binomial"))
#
#Call:
#stats::glm(formula = RESP_tcga ~ ., family = "binomial", data = as.data.frame(Xz_tcga[, 
#    RESPcols]))
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-2.1884   0.4653   0.6066   0.7067   1.0939  
#
#Coefficients:
#                                   Estimate Std. Error z value Pr(>|z|)  
#(Intercept)                          0.9070     0.4573   1.983   0.0473 *
#CUSTOM_MCP_ENDOTHELIAL.CELLS        -0.6622     0.4895  -1.353   0.1761  
#CUSTOM_IFNG3                         0.2560     0.6666   0.384   0.7009  
#xce_Myeloid_dendritic_cell           2.4199     2.5203   0.960   0.3370  
#HALLMARK_INTERFERON_GAMMA_RESPONSE  -1.5045     1.0232  -1.470   0.1415  
#CUSTOM_MCP_MYELOID.DENDRITIC.CELLS  -0.3224     0.7705  -0.418   0.6756  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 165.52  on 159  degrees of freedom
#Residual deviance: 159.19  on 154  degrees of freedom
#  (154 observations deleted due to missingness)
#AIC: 171.19
#
#Number of Fisher Scoring iterations: 4

### OLD MODEL
#> summary(stats::glm(RESP_tcga ~ ., data = as.data.frame(Xz_tcga), family="binomial"))
#
#Call:
#stats::glm(formula = RESP_tcga ~ ., family = "binomial", data = as.data.frame(Xz_tcga))
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-2.0275   0.4291   0.6080   0.7325   1.0940  
#
#Coefficients: (1 not defined because of singularities)
#                     Estimate Std. Error z value Pr(>|z|)    
#(Intercept)            1.7805     0.3060   5.818 5.95e-09 ***
#isTMBhigh                  NA         NA      NA       NA    
#CD274                  0.3688     0.2477   1.489   0.1365    
#IFNG                  -0.7461     0.5096  -1.464   0.1431    
#xce_Endothelial_cell  -1.2613     0.6178  -2.042   0.0412 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#    Null deviance: 165.52  on 159  degrees of freedom
#Residual deviance: 158.88  on 156  degrees of freedom
#  (154 observations deleted due to missingness)
#AIC: 166.88
#
#Number of Fisher Scoring iterations: 4

load("datas.RData")

# Days to months
dat_tcga[,"PFS.time"] <- dat_tcga[,"PFS.time"] / 30.5
dat_tcga[,"OS.time"] <- dat_tcga[,"OS.time"] / 30.5

# Impute <0 gex in gex_prat with min(prat>0)
gex_prat[gex_prat<0] <- min(gex_prat[gex_prat>0])

# Bind Prat et al. vs. TCGA
gex_prat_tcga <- cbind(
	gex_prat[intersect(rownames(gex_prat), rownames(gex_tcga)),],
	gex_tcga[intersect(rownames(gex_prat), rownames(gex_tcga)),]
)
dat_prat_tcga <- rbind(
	cbind(dat_prat[,intersect(colnames(dat_prat), colnames(dat_tcga))], isIO = 1, isChemo = 0),
	cbind(dat_tcga[,intersect(colnames(dat_prat), colnames(dat_tcga))], isIO = 0, isChemo = 1)
)
PFS_prat_tcga <- c(Surv(time=dat_prat[,"PFS.time"], event=dat_prat[,"PFS.event"]), 
	Surv(time=dat_tcga[,"PFS.time"], event=dat_tcga[,"PFS.event"]))
# TODO: Update TCGA RESP to comply with Responder ~ !PD
RESP_prat_tcga <- c(RESP_prat, RESP_tcga)

# Bind Braun Nivo vs. Ever
gex_brauns <- cbind(
	gex_braun_nivo,
	gex_braun_ever
)
dat_brauns <- rbind(
	cbind(dat_braun_nivo, isIO = 1, isChemo = 0),
	cbind(dat_braun_ever, isIO = 0, isChemo = 0)
)
PFS_brauns <- c(Surv(time=dat_braun_nivo[,"PFS.time"], event=dat_braun_nivo[,"PFS.event"]), 
	Surv(time=dat_braun_ever[,"PFS.time"], event=dat_braun_ever[,"PFS.event"]))


# logz the combined datasets
logz <- function(x) { 
	tmp <- scale(log(x+1)) 
	if(any(!is.finite(tmp))){
		tmp[!is.finite(tmp)] <- 0
	}
	tmp
}
translogz <- function(gex){
	rown <- colnames(gex)
	gex <- apply(gex, MARGIN=1, FUN=logz)
	rownames(gex) <- rown
	gex
}
scaleonly <- function(gex){
	rown <- colnames(gex)
	gex <- apply(gex, MARGIN=1, FUN=scale)
	rownames(gex) <- rown
	gex
}
# Combined
gex_prat_tcga <- translogz(gex_prat_tcga)
gex_brauns <- translogz(gex_brauns)
# Separate
#gex_prat <- translogz(gex_prat) 
gex_prat <- scaleonly(gex_prat)
gex_tcga <- translogz(gex_tcga)
gex_braun_nivo <- translogz(gex_braun_nivo)
gex_braun_ever <- translogz(gex_braun_ever)

# Omit melanoma from Prat et al.
gex_prat <- gex_prat[-which(dat_prat[,"Type"]=="MELANOMA"),]
PFS_prat <- PFS_prat[-which(dat_prat[,"Type"]=="MELANOMA"),]
RESP_prat <- RESP_prat[-which(dat_prat[,"Type"]=="MELANOMA")]
dat_prat <- dat_prat[-which(dat_prat[,"Type"]=="MELANOMA"),]

## Datasets look roughly equal for distributions
#
#> quantile(gex_prat_tcga, probs=seq(0,1,.1))
#         0%         10%         20%         30%         40%         50%         60%         70%         80%         90%        100% 
#-4.15229061 -1.37559702 -0.63504017 -0.29961864 -0.00954894  0.18894548  0.34669507  0.50259304  0.69657637  1.01222756  9.75493651 
#> quantile(gex_brauns, probs=seq(0,1,.1))
#          0%          10%          20%          30%          40%          50%          60%          70%          80%          90%         100% 
#-13.46884098  -0.86833959  -0.51815395  -0.33348711  -0.21501742  -0.11617387   0.00222266   0.14766584   0.55049372   1.03722375  17.55783962

#> dim(gex_prat_tcga)
#[1] 354 716
#> dim(gex_brauns)
#[1]   311 43864

coxphs <- lapply(1:ncol(gex_prat), FUN=function(z){
	try({
		cox_prat <- survival::coxph(PFS ~ gex_prat, data = data.frame(PFS = PFS_prat, gex_prat = gex_prat[,colnames(gex_prat)[z]]))
		cox_tcga <- survival::coxph(PFS ~ gex_tcga, data = data.frame(PFS = PFS_tcga, gex_tcga = gex_tcga[,colnames(gex_prat)[z]]))
		rbind(coef(summary(cox_prat)), coef(summary(cox_tcga)))
	})
})
names(coxphs) <- colnames(gex_prat)
coxphs[[1]]

biom_prat <- unlist(lapply(coxphs, FUN=function(z){
	if(!class(z)[1]=="try-error"){
		# Extract if Nivo-arm is significantly associated while Chemo-armis not
		z[1,5]<0.1 & z[2,5]>=0.2
	}else{
		FALSE
	}
}))
names(biom_prat) <- colnames(gex_prat)
names(biom_prat)[unlist(biom_prat)]
length(names(biom_prat)[unlist(biom_prat)])

# Identified as having connection to NSCLC or IO in literature, checking for overlap
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
		"IFNG",	# Interferon gamma
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
		"CD8A", "CD8B",
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
	)
)

candidates <- intersect(names(biom_prat)[unlist(biom_prat)], keygenes)

#> intersect(names(biom_prat)[unlist(biom_prat)], keygenes)
# [1] "ADORA2A" "BTLA"    "CCL5"    "CD27"    "CD274"   "CD8A"    "CXCL9"   "CXCR6"   "HLA-E"   "ICOS"    "LAG3"    "PDCD1"   "TIGIT"   "TLR8"    "TNFAIP3" "TNFRSF9"

#> all(candidates %in% rownames(gex_synthetic))
#[1] TRUE

library(oscar)
pfs_prat_oscar_candidates <- oscar(x = gex_prat[,candidates], y = PFS_prat, family = "cox")
cv_pfs_prat_oscar_candidates_seed1 <- cv.oscar(pfs_prat_oscar_candidates, fold=5, seed=1)
cv_pfs_prat_oscar_candidates_seed2 <- cv.oscar(pfs_prat_oscar_candidates, fold=5, seed=2)
cv_pfs_prat_oscar_candidates_seed3 <- cv.oscar(pfs_prat_oscar_candidates, fold=5, seed=3)
cv_pfs_prat_oscar_candidates_seed4 <- cv.oscar(pfs_prat_oscar_candidates, fold=5, seed=4)
cv_pfs_prat_oscar_candidates_seed5 <- cv.oscar(pfs_prat_oscar_candidates, fold=5, seed=5)
cv_pfs_prat_oscar_candidates_seed6 <- cv.oscar(pfs_prat_oscar_candidates, fold=5, seed=6)

# Earlier candidates before omitting melanomas
# ICOS: https://cancerres.aacrjournals.org/content/80/14/3023
# TIGIT: https://jitc.bmj.com/content/8/2/e000957

#> predict(pfs_prat_oscar_candidates, k=2, type="nonzero")
#   ADORA2A    TNFAIP3 
#-0.3083733 -0.3441822 
#> predict(pfs_prat_oscar_candidates, k=3, type="nonzero")
#   ADORA2A      CD274    TNFAIP3 
#-0.2894373 -0.1897381 -0.2903873 
#> predict(pfs_prat_oscar_candidates, k=4, type="nonzero")
#   ADORA2A      CD274    TNFAIP3    TNFRSF9 
#-0.3326448 -0.2189064 -0.3226417  0.1225185


coefs <- predict(pfs_prat_oscar_candidates, k=4, type="nonzero")
names(coefs) <- gsub("minus", "-", names(coefs))

x11()
par(mfrow=c(3,2))
cv.visu(cv_pfs_prat_oscar_candidates_seed1)
cv.visu(cv_pfs_prat_oscar_candidates_seed2)
cv.visu(cv_pfs_prat_oscar_candidates_seed3)
cv.visu(cv_pfs_prat_oscar_candidates_seed4)
cv.visu(cv_pfs_prat_oscar_candidates_seed5)
cv.visu(cv_pfs_prat_oscar_candidates_seed6)

# k = 1-4 seem optimal

unlist(lapply(1:5, FUN=function(k){
	pred_prat <- predict(pfs_prat_oscar_candidates, k=k, type="response", newdata = gex_prat[,match(candidates, colnames(gex_prat))])
	survival::coxph(PFS_prat ~ pred_prat[,1])$concordance["concordance"]
}))
#concordance concordance concordance concordance concordance 
#  0.6723549   0.6979522   0.7013652   0.7081911   0.7235495
  

unlist(lapply(1:5, FUN=function(k){
	pred_tcga <- predict(pfs_prat_oscar_candidates, k=k, type="response", newdata = gex_tcga[,match(candidates, colnames(gex_tcga))])
	survival::coxph(PFS_tcga ~ pred_tcga[,1])$concordance["concordance"]
}))
#concordance concordance concordance concordance concordance 
#  0.5028199   0.4945319   0.5205238   0.5173851   0.4912461
  
#> survival::coxph(PFS_tcga ~ pred_tcga[,1])$concordance["concordance"]
#concordance 
#  0.4971311
#
# Nice!
  
#pROC::auc(response = PFS_tcga, predictor = pred_tcga[,1])

fill.na0 <- function(x){
	x[is.na(x)] <- 0
	x
}

# No real response in Braun et al., perhaps due to cancer differences

unlist(lapply(1:5, FUN=function(k){
	pred_braun_nivo <- predict(pfs_prat_oscar_candidates, k=k, type="response", newdata = fill.na0(gex_braun_nivo[,match(candidates, colnames(gex_braun_nivo))]))
	survival::coxph(PFS_braun_nivo ~ pred_braun_nivo[,1])$concordance["concordance"]
}))

unlist(lapply(1:5, FUN=function(k){
	pred_braun_ever <- predict(pfs_prat_oscar_candidates, k=k, type="response", newdata = fill.na0(gex_braun_ever[,match(candidates, colnames(gex_braun_ever))]))
	survival::coxph(PFS_braun_ever ~ pred_braun_ever[,1])$concordance["concordance"]
}))

cor(cbind(gex_tcga, TMB = dat_tcga[,"TMB"])[,c("CD274", "ADORA2A", "TNFAIP3", "TNFRSF9", "TMB")], use="pairwise.complete.obs")
#> cor(cbind(gex_tcga, TMB = dat_tcga[,"TMB"])[,c("CD274", "ADORA2A", "TNFAIP3", "TNFRSF9", "TMB")], use="pairwise.complete.obs")
#            CD274     ADORA2A     TNFAIP3   TNFRSF9         TMB
#CD274   1.0000000  0.29780668  0.42624139 0.4099023  0.11554263
#ADORA2A 0.2978067  1.00000000  0.44378796 0.3436178 -0.04046758
#TNFAIP3 0.4262414  0.44378796  1.00000000 0.4419534 -0.07185906
#TNFRSF9 0.4099023  0.34361782  0.44195341 1.0000000  0.18156038
#TMB     0.1155426 -0.04046758 -0.07185906 0.1815604  1.00000000

## Similar results with spearman

## TMB independent of all the genetic markers, ICOS ~ TIGIT correlate, HLA-E somewhat correlates with the rest




# Gender and Squamous as potential markers

summary(coxph(PFS_prat ~ dat_prat[,"SEX"] == "M"))

#> coef(summary(coxph(PFS_prat ~ dat_prat[,"SEX"] == "M")))
#                                   coef exp(coef)  se(coef)        z   Pr(>|z|)
#dat_prat[, "SEX"] == "M"TRUE -0.9438411 0.3891303 0.4311981 -2.18888 0.02860553

## -> Being male significantly improves PFS


coef(summary(coxph(PFS_prat ~ dat_prat[,"CRFHIST"] == "SQUAMOUS")))

#> coef(summary(coxph(PFS_prat ~ dat_prat[,"CRFHIST"] == "SQUAMOUS")))
#                                               coef exp(coef)  se(coef)          z  Pr(>|z|)
#dat_prat[, "CRFHIST"] == "SQUAMOUS"TRUE -0.06328496 0.9386759 0.4100476 -0.1543356 0.8773451

## -> No direct link with being squamous,although small coefficient possibly exists based on literature

# Limiting search to only key genes (no TCGA sanity checking for univariate here)
library(oscar)
pfs_prat_oscar_keygenes <- oscar(x = gex_prat[,intersect(colnames(gex_prat), keygenes)], y = PFS_prat, family = "cox")

library(oscar)
resp_prat_oscar_keygenes <- oscar(x = gex_prat[,intersect(colnames(gex_prat), keygenes)], y = RESP_prat, family = "logistic")

cv_resp_prat_oscar_keygenes_seed1 <- cv.oscar(resp_prat_oscar_keygenes, fold=5, seed=1)
cv_resp_prat_oscar_keygenes_seed2 <- cv.oscar(resp_prat_oscar_keygenes, fold=5, seed=2)
cv_resp_prat_oscar_keygenes_seed3 <- cv.oscar(resp_prat_oscar_keygenes, fold=5, seed=3)
cv_resp_prat_oscar_keygenes_seed4 <- cv.oscar(resp_prat_oscar_keygenes, fold=5, seed=4)
cv_resp_prat_oscar_keygenes_seed5 <- cv.oscar(resp_prat_oscar_keygenes, fold=5, seed=5)
cv_resp_prat_oscar_keygenes_seed6 <- cv.oscar(resp_prat_oscar_keygenes, fold=5, seed=6)

x11()
par(mfrow=c(3,2))
cv.visu(cv_resp_prat_oscar_keygenes_seed1)
cv.visu(cv_resp_prat_oscar_keygenes_seed2)
cv.visu(cv_resp_prat_oscar_keygenes_seed3)
cv.visu(cv_resp_prat_oscar_keygenes_seed4)
cv.visu(cv_resp_prat_oscar_keygenes_seed5)
cv.visu(cv_resp_prat_oscar_keygenes_seed6)

# Close to perfect ROC-AUC at k=4
lapply(1:4, FUN=function(k) { predict(resp_prat_oscar_keygenes, k=k, type="nonzero") })
#> lapply(1:4, FUN=function(k) { predict(resp_prat_oscar_keygenes, k=k, type="nonzero") })
#[[1]]
#(Intercept)     TNFAIP3 
#  0.4358063   1.5371501 
#
#[[2]]
#(Intercept)       STAT1     TNFAIP3 
#  0.4187963  -0.8576741   2.0492448 
#
#[[3]]
#(Intercept)       PDCD1       STAT1     TNFAIP3 
#  0.5536882   2.1052032  -1.8249023   1.8744194 
#
#[[4]]
#(Intercept)        ICOS       PDCD1       STAT1     TNFAIP3 
#  0.6710472  -1.5739444   3.5773745  -2.0436992   1.9443965

# Very similar as for PFS

# PCDC1 = CD279 = Programmed Cell Death 1
# "Inhibitory receptor on antigen activated T-cells that plays a critical role in induction and 
# maintenance of immune tolerance to self (PubMed:21276005). Delivers inhibitory signals upon 
# binding to ligands CD274/PDCD1L1 and CD273/PDCD1LG2 (PubMed:21276005). Following T-cell receptor 
# (TCR) engagement, PDCD1 associates with CD3-TCR in the immunological synapse and directly inhibits 
# T-cell activation (By similarity)."
# (Gene Card)

# TNFAIP3: 
# Inhibition of TNFAIP3 increases invasiveness etc of NSCLC
# https://onlinelibrary.wiley.com/doi/abs/10.1002/jcb.29323

