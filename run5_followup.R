####
#
# Model testing, prediction, visualizations etc for Wiki and post-community phase in DREAM IO
#
###

# Local working directory
setwd("D:/Gits/DREAM_2020_IO/")

# Raw data, as well as normalized training X matrices with PFS, OS, and RESP y-variables
load("datas.RData")
# Large scale OSCAR runs spanning hundreds of candidate features
load("temprun_oscar.RData")
# Much more selective OSCAR runs (including curated FO-panels as candidate features) with ~50 features
load("temprun_oscar_selective.RData")
# Cross validation for the selective runs
load("temprun_oscar_selective_cv.RData")

library(ComplexHeatmap)
library(circlize)

####
#### HEATMAP FIGURE
####

## Following transformation was used primarily to bring count-normalized data closer to e.g. normal distribution
# Transformation for gene expression data
logz <- function(x) { 
	tmp <- scale(log(x+1)) 
	if(any(!is.finite(tmp))){
		tmp[!is.finite(tmp)] <- 0
	}
	tmp
}
# Count and represent gene expression as tertiles {low, medium, high} mean over the whole panel
tertiles <- function(gex, genes, summaryfunc = function(x) { mean(x, na.rm=TRUE)} ){
	# Omit genes that are not present in primary data
	genes <- genes[which(genes %in% rownames(gex))]
	# Apply summary statistic over genes (by default mean expression)
	gem <- apply(gex, MARGIN=2, FUN=summaryfunc)
	# Split into tertiles and return a character representation
	qs <- quantile(gem, probs=c(0, 1/3, 2/3, 1), na.rm=TRUE)
	factor(c("Low", "Medium", "High")[findInterval(gem, qs, rightmost.closed = TRUE)], levels=c("Low", "Medium", "High"))
}
# Handle Surv objects to summarize at fixed time point (PFS and OS statuses at 1 year for individuals)
pfs_summary <- function(pfs, time=12){ # Assume time in months by default
	ifelse(pfs[,1]<time,
		ifelse(pfs[,2]==1, "Progressor", "Censored"),
		"Non-progressor"
	)
}
os_summary <- function(os, time=12){ # Assume time in months by default
	ifelse(os[,1]<time,
		ifelse(os[,2]==1, "Dead", "Censored"),
		"Alive"
	)
}

# Tertiles of gene panel averages; one for each: base, extended, and extensive
FO1 <- c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6")
FO2 <- c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CCL5", "BASE_CD8A")
FO3 <- c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CCL5", "BASE_CD8A", "BASE_ALK", "BASE_IDO1")

# Prat et al.
Xz_prat <- t(apply(t(X_prat[,c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_CCL5", "BASE_IDO1")]), MARGIN=1, FUN=scale))
Tz_prat_FO1 <- tertiles(Xz_prat, genes=FO1)
#Tz_prat_FO2 <- tertiles(Xz_prat, genes=FO2)
#Tz_prat_FO3 <- tertiles(Xz_prat, genes=FO3)

ht_prat_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_prat[,"Responder"]],
	"Progression at 1y" = pfs_summary(PFS_prat),
	Type = dat_prat[,"Type"],
	Hist = dat_prat[,"CRFHIST"],
	Drug = dat_prat[,"Drug"],
	Sex  = dat_prat[,"SEX"],
	"Mean-tertiles base-panel" = anno_simple(Tz_prat_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk"),
		"Progression at 1y" = c("Progressor" = "black", "Censored" = "grey", "Non-progressor" = "cornsilk"),
		Type = c("HEADNECK" = rainbow(4)[1], "LUNG NON-SQUAMOUS CANCER" = rainbow(4)[2], "MELANOMA" = rainbow(4)[3], "SQUAMOUS LUNG CANCER" = rainbow(4)[4]),
		Drug = c("NIVOLUMAB" = "yellow", "PEMBROLIZUMAB" = "orange"),
		Hist = c("NON-SQUAMOUS" = rainbow(4)[2], "SQUAMOUS" = rainbow(4)[4]),
		Sex = c("F" = "pink", "M" = "brown")
	)
)
ht_prat_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=5), rep("Extended FO-panel", times=2), rep("Extensive FO-panel", times=1)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)
ht_prat_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_prat <- ComplexHeatmap::Heatmap(
	Xz_prat, 
	top_annotation = ht_prat_col, 
	right_annotation=ht_prat_row, 
	name="GEX (z-score)", 
	column_title = "Prat et al. (Mixed types, primarily NSCLC)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_prat[,"Responder"],
	row_dend_reorder=1:nrow(Xz_prat),
	row_labels = gsub("BASE_", "", rownames(Xz_prat))
)

draw(ht_prat, annotation_legend_list = packLegend(ht_prat_leg0))


# Kim et al.

# Scale
Xz_kim <- t(apply(t(X_kim[,c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_CCL5", "BASE_IDO1", "BASE_ALK")]), MARGIN=1, FUN=scale))
Tz_kim_FO1 <- tertiles(Xz_kim, genes=FO1)

ht_kim_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_kim[,"Responder"]],
	TMB  = anno_simple(dat_kim[,"TMB"], col=circlize::colorRamp2(quantile(dat_kim[,"TMB"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	PDL1 = anno_simple(dat_kim[,"PDL1"], col=circlize::colorRamp2(quantile(dat_kim[,"PDL1"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	"Mean-tertiles base-panel" = anno_simple(Tz_kim_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk")
	)
)
ht_kim_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=5), rep("Extended FO-panel", times=2), rep("Extensive FO-panel", times=2)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)
ht_kim_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_kim_leg1 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_kim[,"TMB"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "TMB")
ht_kim_leg2 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_kim[,"PDL1"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "PDL1")
ht_kim <- ComplexHeatmap::Heatmap(
	Xz_kim,  
	top_annotation = ht_kim_col, name="GEX (z-score)", 
	right_annotation=ht_kim_row, 
	column_title = "Kim et al. (Gastric cancer, metastatic)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_kim[,"Responder"],
	row_dend_reorder=1:nrow(Xz_kim),
	row_labels = gsub("BASE_", "", rownames(Xz_kim))
)

draw(ht_kim, annotation_legend_list = packLegend(ht_kim_leg0, ht_kim_leg1, ht_kim_leg2))

# Gide et al.
Xz_gide <- t(apply(t(X_gide[,c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_CCL5", "BASE_IDO1", "BASE_ALK")]), MARGIN=1, FUN=scale))
Tz_gide_FO1 <- tertiles(Xz_gide, genes=FO1)

ht_gide_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_gide[,"Responder"]],
	"Progression at 1y" = pfs_summary(PFS_gide, time=365),
	"Survival at 1y" = os_summary(OS_gide, time=365),
	Drug = dat_gide[,"Treatment"],
	Age  = anno_simple(dat_gide[,"AAGE"], col=circlize::colorRamp2(quantile(dat_gide[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	Sex  = dat_gide[,"SEX"],
	Site = dat_gide[,"Site"],
	"Mean-tertiles base-panel" = anno_simple(Tz_gide_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk"),
		"Progression at 1y" = c("Progressor" = "black", "Censored" = "grey", "Non-progressor" = "cornsilk"),
		"Survival at 1y" = c("Dead" = "black", "Censored" = "grey", "Alive" = "cornsilk"),
		Drug = c("Nivolumab" = "yellow", "Pembrolizumab" = "orange"),
		Sex = c("Female" = "pink", "Male" = "brown"),
		Site = c("LN" = rainbow(3)[1], "Mucosa" = rainbow(3)[2], "SQ" = rainbow(3)[3])
	)
)
ht_gide_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=5), rep("Extended FO-panel", times=2), rep("Extensive FO-panel", times=2)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)
ht_gide_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_gide_leg1 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_gide[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "Age")
ht_gide <- ComplexHeatmap::Heatmap(
	Xz_gide, 
	top_annotation = ht_gide_col, 
	right_annotation=ht_gide_row, 
	name="GEX (z-score)", 
	column_title = "Gide et al. (Melanoma)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_gide[,"Responder"],
	row_dend_reorder=1:nrow(Xz_gide),
	row_labels = gsub("BASE_", "", rownames(Xz_gide))
)

draw(ht_gide, annotation_legend_list = packLegend(ht_gide_leg0, ht_gide_leg1))


# Hugo et al.
Xz_hugo <- t(apply(t(X_hugo[,c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_CCL5", "BASE_IDO1", "BASE_ALK")]), MARGIN=1, FUN=scale))
Tz_hugo_FO1 <- tertiles(Xz_hugo, genes=FO1)

ht_hugo_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_hugo[,"Responder"]],
	"Survival at 1y" = os_summary(OS_hugo, time=365),
	Age  = anno_simple(dat_hugo[,"AAGE"], col=circlize::colorRamp2(quantile(dat_hugo[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	Sex  = dat_hugo[,"SEX"],
	"Mean-tertiles base-panel" = anno_simple(Tz_hugo_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk"),
		"Survival at 1y" = c("Dead" = "black", "Censored" = "grey", "Alive" = "cornsilk"),
		Sex = c("F" = "pink", "M" = "brown")
	)
)
ht_hugo_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=5), rep("Extended FO-panel", times=2), rep("Extensive FO-panel", times=2)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)
ht_hugo_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_hugo_leg1 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_hugo[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "Age")
ht_hugo <- ComplexHeatmap::Heatmap(
	Xz_hugo, 
	top_annotation = ht_hugo_col, 
	right_annotation=ht_hugo_row, 
	name="GEX (z-score)", 
	column_title = "Hugo et al. (Melanoma)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_hugo[,"Responder"],
	row_dend_reorder=1:nrow(Xz_hugo),
	row_labels = gsub("BASE_", "", rownames(Xz_hugo))
)

draw(ht_hugo, annotation_legend_list = packLegend(ht_hugo_leg0, ht_hugo_leg1))


# Lauss et al.
Xz_lauss <- t(apply(t(X_lauss[,c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_CCL5", "BASE_IDO1", "BASE_ALK")]), MARGIN=1, FUN=scale))
Tz_lauss_FO1 <- tertiles(Xz_lauss, genes=FO1)

ht_lauss_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_lauss[,"Responder"]],
	"Progression at 1y" = pfs_summary(PFS_lauss),
	"Survival at 1y" = os_summary(OS_lauss),
	"Mean-tertiles base-panel" = anno_simple(Tz_lauss_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk"),
		"Progression at 1y" = c("Progressor" = "black", "Censored" = "grey", "Non-progressor" = "cornsilk"),
		"Survival at 1y" = c("Dead" = "black", "Censored" = "grey", "Alive" = "cornsilk")
	)
)
ht_lauss_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=5), rep("Extended FO-panel", times=2), rep("Extensive FO-panel", times=2)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)
ht_lauss_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_lauss <- ComplexHeatmap::Heatmap(
	Xz_lauss, 
	top_annotation = ht_lauss_col, 
	right_annotation=ht_lauss_row, 
	name="GEX (z-score)", 
	column_title = "Lauss et al. (Melanoma)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_lauss[,"Responder"],
	row_dend_reorder=1:nrow(Xz_lauss),
	row_labels = gsub("BASE_", "", rownames(Xz_lauss))
)

draw(ht_lauss, annotation_legend_list = packLegend(ht_riaz_leg0))


# Riaz et al.
Xz_riaz <- t(apply(t(X_riaz[,c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_CCL5", "BASE_IDO1", "BASE_ALK")]), MARGIN=1, FUN=scale))
Tz_riaz_FO1 <- tertiles(Xz_riaz, genes=FO1)

ht_riaz_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_riaz[,"Responder"]],
	"Mean-tertiles base-panel" = anno_simple(Tz_riaz_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk")
	)
)
ht_riaz_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=5), rep("Extended FO-panel", times=2), rep("Extensive FO-panel", times=2)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)
ht_riaz_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_riaz <- ComplexHeatmap::Heatmap(
	Xz_riaz, 
	top_annotation = ht_riaz_col, 
	right_annotation=ht_riaz_row, 
	name="GEX (z-score)", 
	column_title = "Riaz et al. (Advanced melanoma)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_riaz[,"Responder"],
	row_dend_reorder=1:nrow(Xz_riaz),
	row_labels = gsub("BASE_", "", rownames(Xz_riaz))
)

draw(ht_riaz, annotation_legend_list = packLegend(ht_riaz_leg0))


# Chen et al.
Xz_chen <- t(apply(t(X_chen[,c("BASE_CD274", "BASE_PDCD1", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_CCL5", "BASE_IDO1")]), MARGIN=1, FUN=scale))
Tz_chen_FO1 <- tertiles(Xz_chen, genes=FO1)

ht_chen_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_chen[,"Responder"]],
	"Mean-tertiles base-panel" = anno_simple(Tz_chen_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk")
	)
)
ht_chen_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=4), rep("Extended FO-panel", times=2), rep("Extensive FO-panel", times=1)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)
ht_chen_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_chen <- ComplexHeatmap::Heatmap(
	Xz_chen, 
	top_annotation = ht_chen_col, 
	right_annotation=ht_chen_row, 
	name="GEX (z-score)", 
	column_title = "Chen et al. (Melanoma, metastatic)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_chen[,"Responder"],
	row_dend_reorder=1:nrow(Xz_chen),
	row_labels = gsub("BASE_", "", rownames(Xz_chen))
)

draw(ht_chen, annotation_legend_list = packLegend(ht_chen_leg0))


# Westin et al.
Xz_westin <- t(apply(t(X_westin[,c("BASE_CD274", "BASE_PDCD1", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_CCL5", "BASE_IDO1")]), MARGIN=1, FUN=scale))
Tz_westin_FO1 <- tertiles(Xz_westin, genes=FO1)

ht_westin_col <- ComplexHeatmap::columnAnnotation(
	"Progression at 1y" = pfs_summary(PFS_westin, time=365),
	"Mean-tertiles base-panel" = anno_simple(Tz_westin_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		"Progression at 1y" = c("Progressor" = "black", "Censored" = "grey", "Non-progressor" = "cornsilk")
	)
)
ht_westin_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=4), rep("Extended FO-panel", times=2), rep("Extensive FO-panel", times=1)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)
ht_westin_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")

ht_westin <- ComplexHeatmap::Heatmap(
	Xz_westin, 
	top_annotation = ht_westin_col, 
	right_annotation=ht_westin_row, 
	name="GEX (z-score)", 
	column_title = "Westin et al. (Relapsed follicular lymphoma)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_westin[,"Responder"],
	row_dend_reorder=1:nrow(Xz_westin),
	row_labels = gsub("BASE_", "", rownames(Xz_westin))
)

draw(ht_westin, annotation_legend_list = packLegend(ht_westin_leg0))


# Braun et al. (Nivolumab)

# Scale
Xz_braun_nivo <- t(apply(t(X_braun_nivo[,c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_IDO1", "BASE_ALK")]), MARGIN=1, FUN=scale))
Tz_braun_nivo_FO1 <- tertiles(Xz_braun_nivo, genes=FO1)

ht_braun_nivo_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_braun_nivo[,"Responder"]],
	"Progression at 1y" = pfs_summary(PFS_braun_nivo),
	"Survival at 1y" = os_summary(OS_braun_nivo),
	Age  = anno_simple(dat_braun_nivo[,"AAGE"], col=circlize::colorRamp2(quantile(dat_braun_nivo[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	SEX  = dat_braun_nivo[,"SEX"],
	TMB  = anno_simple(dat_braun_nivo[,"TMB"], col=circlize::colorRamp2(quantile(dat_braun_nivo[,"TMB"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	"Mean-tertiles base-panel" = anno_simple(Tz_braun_nivo_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk"),
		"Progression at 1y" = c("Progressor" = "black", "Censored" = "grey", "Non-progressor" = "cornsilk"),
		"Survival at 1y" = c("Dead" = "black", "Censored" = "grey", "Alive" = "cornsilk"),
		SEX = c("Female" = "pink", "Male" = "brown"),
		"Primary vs. Metas" = dat_braun_nivo[,"PrimaryOrMetastasis"],
		"Mean-tertiles base-panel" = c("Low" = "blue", "Medium" = "black", "High" = "red")
	)
)
ht_braun_nivo_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=5), rep("Extended FO-panel", times=1), rep("Extensive FO-panel", times=2)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)

ht_braun_nivo_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_braun_nivo_leg1 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_braun_nivo[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "Age")
ht_braun_nivo_leg2 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_braun_nivo[,"TMB"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "TMB")

ht_braun_nivo <- ComplexHeatmap::Heatmap(
	Xz_braun_nivo,  
	top_annotation = ht_braun_nivo_col, name="GEX (z-score)", 
	right_annotation=ht_braun_nivo_row, 
	column_title = "Braun et al., Nivolumab chemo-arm (Clear cell renal cell carcinoma)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_braun_nivo[,"Responder"],
	row_dend_reorder=1:nrow(Xz_braun_nivo),
	row_labels = gsub("BASE_", "", rownames(Xz_braun_nivo))
)

draw(ht_braun_nivo, annotation_legend_list = packLegend(ht_braun_nivo_leg0, ht_braun_nivo_leg1, ht_braun_nivo_leg2))




# Braun et al. (Everolimus)

# Scale
Xz_braun_ever <- t(apply(t(X_braun_ever[,c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_IDO1", "BASE_ALK")]), MARGIN=1, FUN=scale))
Tz_braun_ever_FO1 <- tertiles(Xz_braun_ever, genes=FO1)

ht_braun_ever_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_braun_ever[,"Responder"]],
	"Progression at 1y" = pfs_summary(PFS_braun_ever),
	"Survival at 1y" = os_summary(OS_braun_ever),
	Age  = anno_simple(dat_braun_ever[,"AAGE"], col=circlize::colorRamp2(quantile(dat_braun_ever[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	SEX  = dat_braun_ever[,"SEX"],
	TMB  = anno_simple(dat_braun_ever[,"TMB"], col=circlize::colorRamp2(quantile(dat_braun_ever[,"TMB"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	"Mean-tertiles base-panel" = anno_simple(Tz_braun_ever_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk"),
		"Progression at 1y" = c("Progressor" = "black", "Censored" = "grey", "Non-progressor" = "cornsilk"),
		"Survival at 1y" = c("Dead" = "black", "Censored" = "grey", "Alive" = "cornsilk"),
		SEX = c("Female" = "pink", "Male" = "brown"),
		"Primary vs. Metas" = dat_braun_ever[,"PrimaryOrMetastasis"],
		"Mean-tertiles base-panel" = c("Low" = "blue", "Medium" = "black", "High" = "red")
	)
)
ht_braun_ever_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=5), rep("Extended FO-panel", times=1), rep("Extensive FO-panel", times=2)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)

ht_braun_ever_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_braun_ever_leg1 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_braun_ever[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "Age")
ht_braun_ever_leg2 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_braun_ever[,"TMB"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "TMB")

ht_braun_ever <- ComplexHeatmap::Heatmap(
	Xz_braun_ever,  
	top_annotation = ht_braun_ever_col, name="GEX (z-score)", 
	right_annotation=ht_braun_ever_row, 
	column_title = "Braun et al., Everolimus chemo-arm (Clear cell renal cell carcinoma)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_braun_ever[,"Responder"],
	row_dend_reorder=1:nrow(Xz_braun_ever),
	row_labels = gsub("BASE_", "", rownames(Xz_braun_ever))
)

draw(ht_braun_ever, annotation_legend_list = packLegend(ht_braun_ever_leg0, ht_braun_ever_leg1, ht_braun_ever_leg2))




# TCGA (Chemo arm)

# Scale
Xz_tcga <- t(apply(t(X_tcga[,c("BASE_CD274", "BASE_PDCD1", "BASE_TIGIT", "BASE_CXCL9","BASE_CXCR6", "BASE_CD8A", "BASE_CCL5", "BASE_IDO1", "BASE_ALK")]), MARGIN=1, FUN=scale))
Tz_tcga_FO1 <- tertiles(Xz_tcga, genes=FO1)

ht_tcga_col <- ComplexHeatmap::columnAnnotation(
	Responder = c("PD", "CR, PR, SD")[1+dat_tcga[,"Responder"]],
	"Progression at 1y" = pfs_summary(PFS_tcga, time=365),
	"Survival at 1y" = os_summary(OS_tcga, time=365),
	Age  = anno_simple(dat_tcga[,"AAGE"], col=circlize::colorRamp2(quantile(dat_tcga[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	Sex  = dat_tcga[,"SEX"],
	TMB  = anno_simple(dat_tcga[,"TMB"], col=circlize::colorRamp2(quantile(dat_tcga[,"TMB"], na.rm=TRUE), c("black", "black", "orange", "orange", "red"))),
	"Tobacco use" = dat_tcga[,"TOBACUSE"],
	Hist = dat_tcga[,"CRFHIST"],
	"Mean-tertiles base-panel" = anno_simple(Tz_tcga_FO1, col=c("Low" = "darkolivegreen1", "Medium" = "darkgoldenrod1", "High" = "darkorange2")),
	col = list(
		Responder = c("PD" = "black", "CR, PR, SD" = "cornsilk", "Missing" = "grey"),
		"Progression at 1y" = c("Progressor" = "black", "Censored" = "grey", "Non-progressor" = "cornsilk"),
		"Survival at 1y" = c("Dead" = "black", "Censored" = "grey", "Alive" = "cornsilk"),
		Sex = c("F" = "pink", "M" = "brown"),
		"Tobacco use" = c("NEVER" = "cornsilk", "FORMER" = "tan", "CURRENT" = "brown", "UNKNOWN" = "grey"),
		Hist = c("NON-SQUAMOUS" = rainbow(4)[2], "SQUAMOUS" = rainbow(4)[4])		
	)
)

ht_tcga_row <- ComplexHeatmap::rowAnnotation(
	# Adjusted according to what genes were available in the particular data
	PANEL = c(rep("Base FO-panel", times=5), rep("Extended FO-panel", times=2), rep("Extensive FO-panel", times=2)),
	col = list(
		PANEL = c("Base FO-panel" = "cadetblue2", "Extended FO-panel" = "cornflowerblue", "Extensive FO-panel" = "darkgreen")
	)
)

ht_tcga_leg0 = ComplexHeatmap::Legend(labels = c("Low", "Medium", "High"), legend_gp = gpar(fill=c("darkolivegreen1", "darkgoldenrod1", "darkorange2")), title="Base-panel\nexpression\ntertiles")
ht_tcga_leg1 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_tcga[,"AAGE"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "Age")
ht_tcga_leg2 = ComplexHeatmap::Legend(col_fun = circlize::colorRamp2(quantile(dat_tcga[,"TMB"], na.rm=TRUE), c("black", "black", "orange", "orange", "red")), title = "TMB")

ht_tcga <- ComplexHeatmap::Heatmap(
	Xz_tcga,  
	top_annotation = ht_tcga_col, name="GEX (z-score)", 
	right_annotation=ht_tcga_row, 
	column_title = "TCGA, chemo-arm (LUAD & LUSC, NSCLC)", 
	#rect_gp = gpar(col = "cornsilk", lwd = 0), 
	col=colorRampPalette(c("cyan","blue","black","red","orange"))(1000), 
	column_dend_reorder=1+dat_tcga[,"Responder"],
	row_dend_reorder=1:nrow(Xz_tcga),
	row_labels = gsub("BASE_", "", rownames(Xz_tcga))
)

draw(ht_tcga, annotation_legend_list = packLegend(ht_tcga_leg0, ht_tcga_leg1, ht_tcga_leg2))




####
#### KAPLAN-MEIERS / ROC-CURVES
####


library(survminer)






####
#### Example of LASSO penalization coefficients as a function of lambda vs OSCAR unbiased estimates as a function of cardinality in Prat et al.
####

# Oscar fits and CV-runs
load("temprun_oscar_selective.RData")

#par(mfrow=c(2,2))
par(mfrow=c(2,2))
set.seed(0)
library(glmnet)
RESP_prat_lasso <- glmnet::glmnet(y = RESP_prat, x = X_prat[,intersect(selective, colnames(X_prat))], family = "binomial")
RESP_prat_lasso_cv <- glmnet::cv.glmnet(y = RESP_prat, x = X_prat[,intersect(selective, colnames(X_prat))], family = "binomial", nfolds=5, type.measure="auc")
plot(RESP_prat_lasso, xvar="lambda")
abline(v=log(RESP_prat_lasso_cv$lambda.min), col="grey", lwd=2)
plot.new(); plot.window(xlim=c(1,15), ylim=c(-50,50)); box(); axis(1); axis(2)
plot(RESP_prat_oscar, add=TRUE) # Stored object in the .RData
abline(v=which.max(apply(RESP_prat_oscar_cv_f5s1, MARGIN=2, FUN=mean)), col="grey", lwd=2)
title(main="\nOSCAR coefficients per cardinality")
plot(RESP_prat_lasso_cv)
oscar::cv.visu(RESP_prat_oscar_cv_f5s1, xlab="Cardinality 'k'", ylab="ROC-AUC") # 5-fold, seed=1 suffix





