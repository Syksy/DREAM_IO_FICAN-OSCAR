# Various tests combined from OSCAR and literature
load("datas.RData")


# IO arms
coef(summary(stats::glm(RESP_hugo ~ ., data = as.data.frame(X_hugo[,c("BASE_CD274", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))
coef(summary(stats::glm(RESP_prat ~ ., data = as.data.frame(X_prat[,c("BASE_CD274", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))
coef(summary(stats::glm(RESP_kim ~ ., data = as.data.frame(X_kim[,c("BASE_CD274", "BASE_log10TMB", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))
coef(summary(stats::glm(RESP_kim ~ ., data = as.data.frame(X_kim[,c("BASE_CD274", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))
coef(summary(stats::glm(RESP_kim ~ ., data = as.data.frame(X_kim[,c("BASE_CD274", "BASE_log10TMB")]), family="binomial")))
coef(summary(stats::glm(RESP_gide ~ ., data = as.data.frame(X_gide[,c("BASE_CD274",  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))
coef(summary(stats::glm(RESP_braun_nivo ~ ., data = as.data.frame(X_braun_nivo[,c("BASE_CD274", "BASE_log10TMB", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))
coef(summary(stats::glm(RESP_braun_nivo ~ ., data = as.data.frame(X_braun_nivo[,c("BASE_CD274", "BASE_log10TMB")]), family="binomial")))
coef(summary(stats::glm(RESP_braun_nivo ~ ., data = as.data.frame(X_braun_nivo[,c("BASE_CD274", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))
coef(summary(stats::glm(RESP_lauss ~ ., data = as.data.frame(X_lauss[,c("BASE_CD274", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))
# Chemo arms
coef(summary(stats::glm(RESP_tcga ~ ., data = as.data.frame(X_tcga[,c("BASE_CD274", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))
coef(summary(stats::glm(RESP_braun_ever ~ ., data = as.data.frame(X_braun_ever[,c("BASE_CD274", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]), family="binomial")))

#tempcols <- c("BASE_CD274", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "BASE_CD8A")
tempcols <- c("BASE_CD274", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

X_RESP_comb <- rbind(
	X_hugo[,tempcols],
	X_prat[,tempcols],
	X_kim[,tempcols],
	X_gide[,tempcols],
	X_braun_nivo[,tempcols],
	X_lauss[,tempcols],
	X_riaz[,tempcols]
)
Y_RESP_comb <- c(RESP_hugo, RESP_prat, RESP_kim, RESP_gide, RESP_braun_nivo, RESP_lauss, RESP_riaz)

coef(summary(stats::glm(Y_RESP_comb ~ ., data = as.data.frame(X_RESP_comb[,"BASE_CD274"]), family="binomial")))
#> coef(summary(stats::glm(Y_RESP_comb ~ ., data = as.data.frame(X_RESP_comb[,"BASE_CD274"]), family="binomial")))
#                                Estimate Std. Error  z value    Pr(>|z|)
#(Intercept)                   0.16884019 0.11615423 1.453586 0.146061082
#`X_RESP_comb[, "BASE_CD274"]` 0.03635064 0.01386794 2.621200 0.008762083
coef(summary(stats::glm(Y_RESP_comb ~ ., data = as.data.frame(X_RESP_comb), family="binomial")))
#> coef(summary(stats::glm(Y_RESP_comb ~ ., data = as.data.frame(X_RESP_comb), family="binomial")))
#                                              Estimate Std. Error     z value    Pr(>|z|)
#(Intercept)                                 0.16845068 0.11628082  1.44865401 0.147434229
#BASE_CD274                                  0.03643122 0.01391818  2.61752701 0.008856948
#HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION -0.02435734 0.31816680 -0.07655526 0.938977354

# Create subset consisting of relevant cancer subtypes
# Eligible (IO):
# - PFS & RESP: NSCLC and Head&Neck from Prat et al.
# - RESP: Gastric (Kim et al)
# - PFS, OS & RESP: CCRCC from Braun et al.
# Eligible (Chemo):
# - TCGA (NSCLC, both Sq and non-Sq)
# - Braun et al CCRCC Everlumitab chemo arm

setwd("D:\\Gits\\DREAM_2020_IO\\")
load("datas.RData")
dat_prat <- dat_prat[!dat_prat[,"Type"]=="MELANOMA",]
gex_prat <- gex_prat[,rownames(dat_prat)]
#> dim(dat_prat)
#[1] 40 20
#> dim(gex_prat)
#[1] 725  40

#> table(dat_prat$Responder)
#
# 0  1 
#15 25

# Damotte et al. J Transl Med 2019, Fig1d genes
TIL_genes <- c("CD276","HLA-DQA1","CD274","IDO1","HLA-DRB1","HLA-E","CMKLR1","PDCD1LG2","PSMB10","LAG3","CXCL9","STAT1","CD8A","CCL5","NKG7","TIGIT","CD27","CXCR6")
TIL_genes[TIL_genes %in% rownames(gex_prat)]
TIL_genes[!TIL_genes %in% rownames(gex_prat)]
#> TIL_genes[TIL_genes %in% rownames(gex_prat)]
# [1] "CD276"    "HLA-DQA1" "CD274"    "IDO1"     "HLA-E"    "CMKLR1"   "PDCD1LG2" "PSMB10"   "LAG3"     "CXCL9"    "STAT1"    "CD8A"     "CCL5"     "TIGIT"   
#[15] "CD27"     "CXCR6"   
#> TIL_genes[!TIL_genes %in% rownames(gex_prat)]
#[1] "HLA-DRB1" "NKG7"
TIL_genes <- TIL_genes[TIL_genes %in% rownames(gex_prat)]
# Aliases for NKG7: GMP-17, GIG1
# Aliases for HLA-DRB1: HLA-DR1B, DLA-DRB, DRB1, SS1
# Aliases not found in Prat et al.

library(oscar)
library(survival)
RESP_Prat_TIL_OSCAR <- oscar(x = t(gex_prat[TIL_genes,]), y = dat_prat[,"Responder"], family = "logistic")
RESP_Prat_TIL_OSCAR_cv1 <- cv.oscar(RESP_Prat_TIL_OSCAR, fold=5, seed=1)
RESP_Prat_TIL_OSCAR_cv2 <- cv.oscar(RESP_Prat_TIL_OSCAR, fold=5, seed=2)
RESP_Prat_TIL_OSCAR_cv3 <- cv.oscar(RESP_Prat_TIL_OSCAR, fold=5, seed=3)

PFS_Prat_TIL_OSCAR <- oscar(x = t(gex_prat[TIL_genes,]), y = Surv(time = dat_prat[,"PFS.time"], event = as.integer(dat_prat[,"PFS.event"])), family = "cox")
PFS_Prat_TIL_OSCAR_cv1 <- cv.oscar(PFS_Prat_TIL_OSCAR, fold=5, seed=1)
PFS_Prat_TIL_OSCAR_cv2 <- cv.oscar(PFS_Prat_TIL_OSCAR, fold=5, seed=2)
PFS_Prat_TIL_OSCAR_cv3 <- cv.oscar(PFS_Prat_TIL_OSCAR, fold=5, seed=3)

par(mfrow=c(2,3))
cv.visu(RESP_Prat_TIL_OSCAR_cv1, main="CV RESP, Prat, seed=1")
cv.visu(RESP_Prat_TIL_OSCAR_cv2, main="CV RESP, Prat, seed=2")
cv.visu(RESP_Prat_TIL_OSCAR_cv3, main="CV RESP, Prat, seed=3")
cv.visu(PFS_Prat_TIL_OSCAR_cv1, main="CV PFS, Prat, seed=1")
cv.visu(PFS_Prat_TIL_OSCAR_cv2, main="CV PFS, Prat, seed=2")
cv.visu(PFS_Prat_TIL_OSCAR_cv3, main="CV PFS, Prat, seed=3")

#> predict(RESP_Prat_TIL_OSCAR, k=5, type="nonzero")
#(Intercept)       CD274   HLAminusE      PSMB10       TIGIT        CD27 
# -27.686713    1.493569    3.345066   -2.366040   -2.785168    2.401625 
#> predict(PFS_Prat_TIL_OSCAR, k=5, type="nonzero")
#     CD276      CD274  HLAminusE     PSMB10      STAT1 
# 0.8745616 -0.5877714 -1.3864352  0.9427290  0.3036566

RESP_markers <- c("CD274", "HLA-E", "PSMB10", "TIGIT", "CD27")
PFS_markers <- c("CD276", "CD274", "HLA-E", "PSMB10", "STAT1")

coef(summary(stats::glm(dat_prat[,"Responder"] ~ ., data = as.data.frame(t(gex_prat)[,c("CD274", "HLA-E", "PSMB10", "TIGIT", "CD27")]), family="binomial")))
#> coef(summary(stats::glm(dat_prat[,"Responder"] ~ ., data = as.data.frame(t(gex_prat)[,c("CD274", "HLA-E", "PSMB10", "TIGIT", "CD27")]), family="binomial")))
#              Estimate Std. Error   z value   Pr(>|z|)
#(Intercept) -27.691991 13.1325348 -2.108655 0.03497433
#CD274         1.493976  0.7057675  2.116810 0.03427598
#`HLA-E`       3.345868  1.4882162  2.248241 0.02456085
#PSMB10       -2.366827  1.3716342 -1.725553 0.08442786
#TIGIT        -2.786026  1.2341974 -2.257359 0.02398567
#CD27          2.402406  1.0905251  2.202981 0.02759606

# Sanity check matches with the OSCAR solution; all except one coef significant, and even that is borderline

RESP_kim_pred <- c(predict(RESP_Prat_TIL_OSCAR, k=5, newdata=t(gex_kim[TIL_genes,])))
RESP_braun_nivo_pred <- c(predict(RESP_Prat_TIL_OSCAR, k=5, newdata=t(gex_braun_nivo[TIL_genes,])))
RESP_braun_ever_pred <- c(predict(RESP_Prat_TIL_OSCAR, k=5, newdata=t(gex_braun_ever[TIL_genes,])))
RESP_tcga_pred <- c(predict(RESP_Prat_TIL_OSCAR, k=5, newdata=t(gex_tcga[TIL_genes,])))

RESP_kim_pred <- c(predict(RESP_Prat_TIL_OSCAR, k=1, newdata=t(gex_kim[TIL_genes,])))
RESP_braun_nivo_pred <- c(predict(RESP_Prat_TIL_OSCAR, k=1, newdata=t(gex_braun_nivo[TIL_genes,])))
RESP_braun_ever_pred <- c(predict(RESP_Prat_TIL_OSCAR, k=1, newdata=t(gex_braun_ever[TIL_genes,])))
RESP_tcga_pred <- c(predict(RESP_Prat_TIL_OSCAR, k=1, newdata=t(gex_tcga[TIL_genes,])))

pROC::auc(response=RESP_kim, predictor=RESP_kim_pred)
#> pROC::auc(response=RESP_kim, predictor=RESP_kim_pred)
#Setting levels: control = 0, case = 1
#Setting direction: controls > cases
#Area under the curve: 0.5905
pROC::auc(response=RESP_tcga, predictor=RESP_tcga_pred)

#> pROC::auc(response=RESP_tcga, predictor=RESP_tcga_pred)
#Setting levels: control = 0, case = 1
#Setting direction: controls > cases
#Area under the curve: 0.6361





## Top 3 genes in Prat for RESP and PFS: CD274, HLA-E, and STAT1
#
#> predict(RESP_Prat_TIL_OSCAR, k=3, type="nonzero")
#(Intercept)       CD274   HLAminusE       STAT1 
#-12.5281776   0.9504497   1.7603145  -1.3460683 
#> predict(PFS_Prat_TIL_OSCAR, k=3, type="nonzero")
#     CD274  HLAminusE      STAT1 
#-0.3070573 -0.7347333  0.4828231

#> pROC::auc(response=RESP_kim, predictor=RESP_kim_pred)
#Setting levels: control = 0, case = 1
#Setting direction: controls > cases
#Area under the curve: 0.5206
#> pROC::auc(response=RESP_tcga, predictor=RESP_tcga_pred)
#Setting levels: control = 0, case = 1
#Setting direction: controls > cases
#Area under the curve: 0.637
#> predict(RESP_Prat_TIL_OSCAR, k=1, type="nonzero")
#(Intercept)   HLAminusE 
# -17.634693    1.492741

# -> HLA-E predicts in Chemo arm, eliminate it from the model



#### Creating submission 2 for PFS

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


coef(summary(coxph(PFS_prat ~ .,
	data = data.frame(
		CD274 = gex_prat[,"CD274"],
		CD8A = gex_prat[,"CD8A"],
		isMale = as.integer(dat_prat[,"SEX"]=="M"),
		isEversmoker = as.integer(dat_prat[,"TOBACUSE"] %in% c("FORMER", "CURRENT"))
		#isSquamous = as.integer(dat_prat[,"CRFHIST"]=="M"),		
	)
)))
#                   coef exp(coef)  se(coef)         z    Pr(>|z|)
#CD274        -0.3579477 0.6991096 0.2804826 -1.276185 0.201890099
#CD8A         -0.3910739 0.6763301 0.2902005 -1.347599 0.177787328
#isMale       -1.6070974 0.2004687 0.5860794 -2.742115 0.006104487
#isEversmoker  1.1134446 3.0448285 0.8415577  1.323076 0.185810211

# Slight trend in aggregating non-responders to low values in both
png("CD274_CD8A_vs_RESP.png", width=800, height=800)
plot(gex_prat[,"CD274"], gex_prat[,"CD8A"], pch=16, col=1+RESP_prat, xlab="CD274", ylab="CD8A", main="Prat et al.,\nn-SQ NSCLC n=22, SQ NSCLC n=13; H&NCa n=5")
legend("topleft", pch=16, col=1:2, legend=c("Non-responder (PD)", "Responder (non-PD)"))
dev.off()


## isEversmoker value is in opposite direction from what would've been expected from literature

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

