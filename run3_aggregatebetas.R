# Estimate tumor mutational burden coefficients and/or X matrix

aggregateX <- function(
	gex,
	dat
){
	# Carbone et al.
	# For efficacy analyses, patients were grouped in thirds according to tumor-mutation burden. 
	# The boundaries for these three
	# groups were a tumor-mutation burden of 0 to less than 100 mutations (low burden), 100 to
	# 242 mutations (medium burden), and 243 or more mutations (high burden)
	# ...
	# Among the patients with a high tumor-mutation burden, the response rate was higher in the
	# nivolumab group than in the chemotherapy group (47% vs. 28%), and progression-free
	# survival was longer (median, 9.7 vs. 5.8 months; hazard ratio for disease progression or
	# death, 0.62; 95% CI, 0.38 to 1.00) (Fig. 2C).
	# ...
	# However, in the nivolumab group, patients with both a high
	# tumor-mutation burden and a PD-L1 expression level of 50% or more had a higher response
	# rate (75%) than those with only one of these factors
	# ...
	# Patients with both a high tumor-mutation
	# burden and a PD-L1 expression level of 50% or more may have a greater likelihood of re
	# sponse to nivolumab than those with only one or neither of these factors. Overall, the current
	# findings are consistent with the hypothesis that immunotherapy may have enhanced activity
	# in patients with a high tumor-mutation burden.14 However, because this was an exploratory
	#
	# ---> Fig 2
	#TMBclass <- c("Low", "Medium", "High")[findInterval(dat[,"TMB"], c(-Inf, 100, 242, Inf)]
	
	library(GSVA)
	library(immunedeconv)
	
	# 242 mutations threshold for 'high' mutatation; impute NAs as if being not highly mutated
	X <- data.frame(isTMBhigh = as.integer(dat[,"TMB"] > 242))
	rownames(X) <- rownames(dat)
	# Impute zero indicators if there are NA values
	if(any(is.na(X[,"isTMBhigh"]))) X[is.na(X[,"isTMBhigh"]),"isTMBhigh"] <- 0
	
	# CD274 expression level modelled as a surrogate for PD-L1 IHC
	# Normalized expressions between various platforms and their respective distributional characteristics
	X <- cbind(X, CD274 = gex["CD274",])
	
	# GSVA for CUSTOM_IFNG3
	gmt_custom <- GSEABase::getGmt(".\\selfmade.gmt")
	res_gsva <- t(GSVA::gsva(as.matrix(gex), gmt_custom, verbose=FALSE)) # Custom GMTs
	X <- cbind(X, IFNG = res_gsva[,grep("CUSTOM_IFNG3", colnames(res_gsva))])
	
	# Epithelial cell expression as reported by xCell
	tmp <- immunedeconv::deconvolute(gex, method="xcell")
	tmp <- as.matrix(tmp)
	rownames(tmp) <- paste("xce_", gsub(" ", "_", tmp[,1]), sep="")
	tmp <- tmp[,-1]
	class(tmp) <- "numeric"
	X <- cbind(X, xce_Endothelial_cell = t(tmp)[,"xce_Endothelial_cell"])	
	
	# Return X as a matrix
	as.matrix(X)
}

setwd("D:\\Gits\\DREAM_2020_IO")
load("temp.RData")
load("temprun_tcga.RData")
load("tempspace_1.RData")
load("tempspace_2.RData")
load("tempspace_3.RData")
load("tempspace_4.RData")
load("tempspace_5.RData")
load("tempspace_6.RData")

# Plot various gene expression distributions to see if Cox PH coefficients would generalize between studies

png("GEXs_notransforms.png", width=800, height=800)
par(mfrow=c(3,3))
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
dev.off()

logz <- function(x) { 
	tmp <- scale(log(x+1)) 
	if(any(!is.finite(tmp))){
		tmp[!is.finite(tmp)] <- 0
	}
	tmp
}

png("GEXs_logzs.png", width=800, height=800)
par(mfrow=c(3,3))
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
dev.off()

# Same the log-z standardized gene expression patterns
gxz_synthetic <- t(apply(gex_synthetic, MARGIN=1, FUN=logz))
rownames(gxz_synthetic) <- rownames(gex_synthetic)
gxz_tcga <- t(apply(gex_tcga, MARGIN=1, FUN=logz))
rownames(gxz_tcga) <- rownames(gex_tcga)
gxz_hugo <- t(apply(gex_hugo, MARGIN=1, FUN=logz))
rownames(gxz_hugo) <- rownames(gex_hugo)
gxz_westin <- t(apply(gex_westin, MARGIN=1, FUN=logz))
rownames(gxz_westin) <- rownames(gex_westin)

# Aggregate X matrices
# PFS:	(Prat et al. &) Lauss et al. & Westin et al.
# OS:	Hugo et al. & Lauss et al.
# RESP	Hugo et a. & Prat et al. & Lauss et al. & Kim et al. & Chen et al.

# Chemo-arm
Xz_tcga <- aggregateX(gex=gxz_tcga, dat=dat_tcga)
# IOs
# Xz_prat <- aggregateX(gex=gxz_prat, dat=dat_prat)
Xz_hugo <- aggregateX(gex=gxz_hugo, dat=dat_hugo)
Xz_lauss <- aggregateX(gex=gex_lauss, dat=dat_lauss)
Xz_westin <- aggregateX(gex=gxz_westin, dat=dat_westin)
Xz_kim <- aggregateX(gex=gex_kim, dat=dat_kim)
#Xz_chen <- aggregateX(gex=gex_chen, dat=dat_chen)


# Combine studies to generate aggregate HR estimates

combcols <- c("CD274", "IFNG", "xce_Endothelial_cell")

# xCell fails for Prat et al. due to too narrow gene panel
library(survival)

## NOTE!
# Responses not necessarily in same unit; e.g. Hugo et al. OS in days, Lauss in months
PFS2_westin <- PFS_westin
PFS2_westin[,1] <- PFS2_westin[,1]/30.5

OS2_hugo <- OS_hugo
OS2_hugo[,1] <- OS2_hugo[,1]/30.5

# Coefficient for linear combination for Carbone et al.
# for tumor high burden (based on absolute tumor mutation count)
# HR: 0.62
# -> 
#> log(0.62)
#[1] -0.4780358

######## COMBINED PROGRESSION FREE SURVIVAL ######

Xz_PFS <- rbind(Xz_lauss[,combcols], Xz_westin[,combcols])
Yz_PFS <- c(PFS_lauss, PFS2_westin)
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
   
######### COMBINED RESPONSE ########

Xz_RESP <- rbind(Xz_hugo[,combcols], Xz_lauss[,combcols], Xz_kim[,combcols])
Yz_RESP <- c(RESP_hugo, RESP_lauss, RESP_kim)
summary(stats::glm(Yz_RESP ~ ., data = as.data.frame(Xz_RESP), family="binomial"))
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
survival::coxph(OS_lauss ~ ., data = as.data.frame(Xz_lauss))
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







