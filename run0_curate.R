###
#
# DREAM 2020; Anti-PD1 Response Prediction DREAM Challenge
#
###

# First setwd to correct project root!
# setwd("...")
setwd("D:\\Gits\\DREAM_2020_IO\\")

# Just in case ...
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.12")
#update.packages(ask=FALSE, checkBuilt=TRUE)

# Generate biomaRt gene mapping data.frame for convenience
#BiocManager::install("biomaRt")
library("biomaRt")
# Generate gene names
#> grep("entrez", listAttributes(biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"))[,1], value=TRUE)
#[1] "entrezgene_trans_name"  "entrezgene_description" "entrezgene_accession"   "entrezgene_id"
# Fetch gene names for various annotations
genes <- biomaRt::getBM(
	attributes = 
		c(
			# ENSEMBL
			'ensembl_gene_id', 'ensembl_transcript_id',
			# entrez
			'entrezgene_id',
			# Hugo
			'hgnc_symbol',
			# RefSeq
			'refseq_mrna',
			# Chromosomal information
			'chromosome_name','start_position','end_position',
			# Description
			'description'
		),
	mart = biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
)
# Sort using chromosomes and then bp locations
genes <- genes[order(genes$chromosome_name, genes$start_position, genes$end_position),]
# Omit row names (wrong order indices)
rownames(genes) <- NULL

### Synthetic dataset from DREAM IO 2020 (no signal, must right composition for data)

# Read in synthetic data and examine fields
cli_synthetic <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\clinical_data.csv", row.names=1)
gex_synthetic_ensembl75_genes_count <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_ensembl75_genes_count.csv", row.names=1)
gex_synthetic_ensembl75_genes_tpm <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_ensembl75_genes_tpm.csv", row.names=1)
gex_synthetic_ensembl75_isoforms_count <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_ensembl75_isoforms_count.csv", row.names=1)
gex_synthetic_ensembl75_isoforms_tpm <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_ensembl75_isoforms_tpm.csv", row.names=1)
gex_synthetic_refseq105_genes_count <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_count.csv", row.names=1)
gex_synthetic_refseq105_genes_tpm <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_genes_tpm.csv", row.names=1)
gex_synthetic_refseq105_isoforms_count <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_isoforms_count.csv", row.names=1)
gex_synthetic_refseq105_isoforms_tpm <- read.csv(".\\CM_026_formatted_synthetic_data_subset\\GRCh37ERCC_refseq105_isoforms_tpm.csv", row.names=1)

## Ideal format of a clinical information table (strings as factors): 
#> head(cli_synthetic)
#     SEX     AAGE      CRFHIST TOBACUSE ECOGPS PDL1 TMB TCR_Shannon TCR_Richness TCR_Evenness BCR_Shannon BCR_Richness BCR_Evenness
#p227   M 49.10398 NON-SQUAMOUS  CURRENT      1   55  76    2.806938           21   0.12289143    6.721297          903  0.046129329
#p359   M 76.11855 NON-SQUAMOUS    NEVER      1   69  NA    1.829465           11   0.20272961    6.048394         3083  0.027932308
#p533   M 55.32944 NON-SQUAMOUS    NEVER      1   11 482    3.576665           40   0.04420034    5.655366          978  0.009695415
#p149   F 64.66889 NON-SQUAMOUS    NEVER      1   13  16    1.623662           27   0.07185817    4.025997         4092  0.015561859
#p160   M 50.30217 NON-SQUAMOUS   FORMER      1   93 171    3.239417           36   0.15662204    7.523416          240  0.001109798
#p2     F 75.05628 NON-SQUAMOUS   FORMER      1   18  71    3.111728           24   0.11196773    5.702087         5213  0.005632520

## Ideal format (especially row names) for gene expression (gex) data:
#> gex_synthetic_ensembl75_genes_count[1:3,1:3]
#                   p227    p359   p533
#ENSG00000203782    0.00    0.00   0.00
#ENSG00000163431 2189.24 3642.99 981.26
#ENSG00000239696    0.30    0.22   0.05
#> gex_synthetic_ensembl75_genes_tpm[1:3,1:3]
#                p227  p359 p533
#ENSG00000203782 0.00  0.01 0.00
#ENSG00000163431 0.67 10.80 5.05
#ENSG00000239696 0.01  0.00 0.02
#> gex_synthetic_ensembl75_isoforms_count[1:3,1:3]
#                  p227   p359     p533
#ENST00000415270 898.94 332.74    302.3
#ENST00000509684  44.70  38.92 131366.3
#ENST00000555703   0.01   0.00      0.0
#> gex_synthetic_ensembl75_isoforms_tpm[1:3,1:3]
#                p227  p359 p533
#ENST00000415270 3.14 10.68 1.55
#ENST00000509684 0.16  0.06 0.33
#ENST00000555703 0.00  0.01 0.03
#> gex_synthetic_refseq105_genes_count[1:3,1:3]
#           p227    p359   p533
#PHF20L1 4055.00 6892.00 3911.0
#OGDHL     56.00   34.00    5.0
#ZCCHC18   41.75   88.78   73.7
#> gex_synthetic_refseq105_genes_tpm[1:3,1:3]
#         p227  p359 p533
#PHF20L1 27.12 12.76 9.64
#OGDHL    0.00  0.20 0.30
#ZCCHC18  1.18  0.39 0.84
#> gex_synthetic_refseq105_isoforms_count[1:3,1:3]
#                 p227   p359    p533
#NR_049841.1      0.00   0.00    0.00
#NM_001039584.1 256.59 463.02 1005.26
#NM_001350039.2   0.00   0.00    0.00
#> gex_synthetic_refseq105_isoforms_tpm[1:3,1:3]
#               p227 p359 p533
#NR_049841.1    0.00 0.00 0.00
#NM_001039584.1 0.42 5.45 1.44
#NM_001350039.2 0.00 0.00 0.00

## Note especially G -> T in the ENSx###-codes when 'genes' -> 'isoforms'
## Intuitive explanation for G -> T mapping (i.e. one gene ENSGx can have multiple corresponding transcript IDs ENSTx):
## https://www.biostars.org/p/199073/
## Refseq uses gene names when in 'genes' not 'isoforms'


# Sanity check of dimensions
#> dim(gex_synthetic_ensembl75_genes_count)
#[1] 57997    56
#> dim(gex_synthetic_refseq105_genes_count)
#[1] 29182    56
#> dim(gex_synthetic_ensembl75_genes_tpm)
#[1] 57997    56
#> dim(gex_synthetic_refseq105_genes_tpm)
#[1] 29182    56
#> dim(gex_synthetic_ensembl75_isoforms_count)
#[1] 196593     56
#> dim(gex_synthetic_refseq105_isoforms_count)
#[1] 74375    56

## Non-comforming dimensions, except when normalized count -> tpm

# '-' symbols may require some sanitizing along the way (?), e.g. manual replacement with '.'?
#> grep("HLA", rownames(gex_synthetic_refseq105_genes_tpm), value=TRUE)
# [1] "SCHLAP1"      "HHLA3"        "HLA-V"        "HLA-DPA1"     "HLA-B"        "HLA-G"        "HLA-DMA"     
# [8] "HLA-L"        "HLA-DRB1"     "HLA-DQB1"     "HLA-E"        "HLA-A"        "HLA-DOA"      "HLA-DMB"     
#[15] "HLA-DQA2"     "HHLA2"        "HLA-DPB1"     "HLA-DQB1-AS1" "HLA-DRB6"     "HLA-DRB5"     "HLA-DQB2"    
#[22] "HLA-J"        "HLA-DPB2"     "HLA-F"        "HLA-F-AS1"    "HLA-C"        "HLA-DOB"      "HHLA1"       
#[29] "HLA-DRA"      "HLA-H"        "HLA-DQA1"


### Corresponding datasets from TCGA

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("RTCGA")
#library("RTCGA")
#RTCGA::checkTCGA(what="Dates", c("LUAD", "LUSC"))
#RTCGA::checkTCGA(what="DataSets", "LUAD")
#RTCGA::checkTCGA(what="DataSets", "LUSC")
# Latest "2016-01-28"
# Clinical annotations
#RTCGA::downloadTCGA(c("LUAD", "LUSC"), "Merge_Clinical.Level_1", "TCGA", date="2016-01-28", untarFile=TRUE, removeTar=TRUE)
#> list.files("TCGA")
#[1] "gdac.broadinstitute.org_LUAD.Merge_Clinical.Level_1.2016012800.0.0"
#[2] "gdac.broadinstitute.org_LUSC.Merge_Clinical.Level_1.2016012800.0.0"
#RTCGA::downloadTCGA(c("LUAD", "LUSC"), "mRNAseq_Preprocess.Level_3", "TCGA", date="2016-01-28", untarFile=TRUE, removeTar=TRUE)

# curatedTCGAdata from BioconductoR
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("curatedTCGAData")

tmp_luad <- curatedTCGAData::curatedTCGAData("LUAD", c("RNASeq2GeneNorm"), FALSE)
# Prompts possibly a question to create a temporary directory
tmp_lusc <- curatedTCGAData::curatedTCGAData("LUSC", c("RNASeq2GeneNorm"), FALSE)
# Prompts possibly a question to create a temporary directory

## use curatedTCGAData instead:
gex_luad <- tmp_luad[[1]]@assays$data@listData[[1]]
gex_lusc <- tmp_lusc[[1]]@assays$data@listData[[1]]
# Generate easier unique colnames
colnames(gex_luad) <- substr(colnames(gex_luad), 1, 12)
colnames(gex_lusc) <- substr(colnames(gex_lusc), 1, 12)


## Generation of clinical info
## Specific fields

#> colData(tmp_luad)[1:2,grep("age_at|tobac|performance|cell|purity|targeted|therapy", colnames(colData(tmp_luad)), value=TRUE)]
#DataFrame with 2 rows and 130 columns
#             radiation_therapy karnofsky_performance_score year_of_tobacco_smoking_onset patient.age_at_initial_pathologic_diagnosis patient.drugs.drug.10.days_to_drug_therapy_end patient.drugs.drug.10.days_to_drug_therapy_start
#                   <character>                   <integer>                     <integer>                                   <integer>                                      <integer>                                        <integer>
#TCGA-05-4249                no                          NA                            NA                                          67                                             NA                                               NA
#TCGA-05-4382               yes                          NA                            NA                                          68                                             NA                                               NA
#             patient.drugs.drug.10.therapy_ongoing patient.drugs.drug.10.therapy_types.therapy_type patient.drugs.drug.2.days_to_drug_therapy_end patient.drugs.drug.2.days_to_drug_therapy_start patient.drugs.drug.2.therapy_ongoing
#                                       <character>                                      <character>                                     <integer>                                       <integer>                          <character>
#TCGA-05-4249                                    NA                                               NA                                            NA                                              NA                                   NA
#TCGA-05-4382                                    NA                                               NA                                            NA                                              NA                                   NA
#             patient.drugs.drug.2.therapy_types.therapy_type patient.drugs.drug.2.therapy_types.therapy_type_notes patient.drugs.drug.3.days_to_drug_therapy_end patient.drugs.drug.3.days_to_drug_therapy_start
#                                                 <character>                                           <character>                                     <integer>                                       <integer>
#TCGA-05-4249                                              NA                                                    NA                                            NA                                              NA
#TCGA-05-4382                                              NA                                                    NA                                            NA                                              NA
#             patient.drugs.drug.3.therapy_ongoing patient.drugs.drug.3.therapy_types.therapy_type patient.drugs.drug.4.days_to_drug_therapy_end patient.drugs.drug.4.days_to_drug_therapy_start patient.drugs.drug.4.therapy_ongoing
#                                      <character>                                     <character>                                     <integer>                                       <integer>                          <character>
#TCGA-05-4249                                   NA                                              NA                                            NA                                              NA                                   NA
#TCGA-05-4382                                   NA                                              NA                                            NA                                              NA                                   NA
#             patient.drugs.drug.4.therapy_types.therapy_type patient.drugs.drug.4.therapy_types.therapy_type_notes patient.drugs.drug.5.days_to_drug_therapy_end patient.drugs.drug.5.days_to_drug_therapy_start
#                                                 <character>                                           <character>                                     <integer>                                       <integer>
#TCGA-05-4249                                              NA                                                    NA                                            NA                                              NA
#TCGA-05-4382                                              NA                                                    NA                                            NA                                              NA
#             patient.drugs.drug.5.therapy_ongoing patient.drugs.drug.5.therapy_types.therapy_type patient.drugs.drug.6.days_to_drug_therapy_end patient.drugs.drug.6.days_to_drug_therapy_start patient.drugs.drug.6.therapy_ongoing
#                                      <character>                                     <character>                                     <integer>                                       <integer>                          <character>
#TCGA-05-4249                                   NA                                              NA                                            NA                                              NA                                   NA
#TCGA-05-4382                                   NA                                              NA                                            NA                                              NA                                   NA
#             patient.drugs.drug.6.therapy_types.therapy_type patient.drugs.drug.7.days_to_drug_therapy_end patient.drugs.drug.7.days_to_drug_therapy_start patient.drugs.drug.7.therapy_ongoing
#                                                 <character>                                     <integer>                                       <integer>                          <character>
#TCGA-05-4249                                              NA                                            NA                                              NA                                   NA
#TCGA-05-4382                                              NA                                            NA                                              NA                                   NA
#             patient.drugs.drug.7.therapy_types.therapy_type patient.drugs.drug.8.days_to_drug_therapy_end patient.drugs.drug.8.days_to_drug_therapy_start patient.drugs.drug.8.therapy_ongoing
#                                                 <character>                                     <integer>                                       <integer>                          <character>
#TCGA-05-4249                                              NA                                            NA                                              NA                                   NA
#TCGA-05-4382                                              NA                                            NA                                              NA                                   NA
#             patient.drugs.drug.8.therapy_types.therapy_type patient.drugs.drug.9.days_to_drug_therapy_end patient.drugs.drug.9.days_to_drug_therapy_start patient.drugs.drug.9.therapy_ongoing
#                                                 <character>                                     <integer>                                       <integer>                          <character>
#TCGA-05-4249                                              NA                                            NA                                              NA                                   NA
#TCGA-05-4382                                              NA                                            NA                                              NA                                   NA
#             patient.drugs.drug.9.therapy_types.therapy_type patient.drugs.drug.days_to_drug_therapy_end patient.drugs.drug.days_to_drug_therapy_start patient.drugs.drug.therapy_ongoing patient.drugs.drug.therapy_types.therapy_type
#                                                 <character>                                   <integer>                                     <integer>                        <character>                                   <character>
#TCGA-05-4249                                              NA                                          NA                                            NA                                 NA                                            NA
#TCGA-05-4382                                              NA                                          NA                                            NA                                 NA                                            NA
#             patient.drugs.drug.therapy_types.therapy_type_notes patient.follow_ups.follow_up.2.additional_pharmaceutical_therapy patient.follow_ups.follow_up.2.additional_radiation_therapy
#                                                     <character>                                                      <character>                                                 <character>
#TCGA-05-4249                                                  NA                                                               NA                                                          NA
#TCGA-05-4382                                                  NA                                                               NA                                                          NA
#             patient.follow_ups.follow_up.2.karnofsky_performance_score patient.follow_ups.follow_up.2.performance_status_scale_timing patient.follow_ups.follow_up.2.primary_therapy_outcome_success
#                                                              <integer>                                                    <character>                                                    <character>
#TCGA-05-4249                                                         NA                                                             NA                                                             NA
#TCGA-05-4382                                                         NA                                                             NA                                                             NA
#             patient.follow_ups.follow_up.2.radiation_therapy patient.follow_ups.follow_up.2.targeted_molecular_therapy patient.follow_ups.follow_up.3.additional_pharmaceutical_therapy
#                                                  <character>                                               <character>                                                      <character>
#TCGA-05-4249                                               NA                                                        NA                                                               NA
#TCGA-05-4382                                               NA                                                        NA                                                               NA
#             patient.follow_ups.follow_up.3.additional_radiation_therapy patient.follow_ups.follow_up.3.karnofsky_performance_score patient.follow_ups.follow_up.3.performance_status_scale_timing
#                                                             <character>                                                  <integer>                                                    <character>
#TCGA-05-4249                                                          NA                                                         NA                                                             NA
#TCGA-05-4382                                                          NA                                                         NA                                                             NA
#             patient.follow_ups.follow_up.3.primary_therapy_outcome_success patient.follow_ups.follow_up.3.radiation_therapy patient.follow_ups.follow_up.3.targeted_molecular_therapy
#                                                                <character>                                      <character>                                               <character>
#TCGA-05-4249                                                             NA                                               NA                                                        NA
#TCGA-05-4382                                                             NA                                               NA                                                        NA
#             patient.follow_ups.follow_up.4.additional_pharmaceutical_therapy patient.follow_ups.follow_up.4.additional_radiation_therapy patient.follow_ups.follow_up.4.performance_status_scale_timing
#                                                                  <character>                                                 <character>                                                    <character>
#TCGA-05-4249                                                               NA                                                          NA                                                             NA
#TCGA-05-4382                                                               NA                                                          NA                                                             NA
#             patient.follow_ups.follow_up.4.primary_therapy_outcome_success patient.follow_ups.follow_up.4.radiation_therapy patient.follow_ups.follow_up.4.targeted_molecular_therapy
#                                                                <character>                                      <character>                                               <character>
#TCGA-05-4249                                                             NA                                               NA                                                        NA
#TCGA-05-4382                                                             NA                                               NA                                                        NA
#             patient.follow_ups.follow_up.additional_pharmaceutical_therapy patient.follow_ups.follow_up.additional_radiation_therapy patient.follow_ups.follow_up.karnofsky_performance_score
#                                                                <character>                                               <character>                                                <integer>
#TCGA-05-4249                                                             NA                                                        NA                                                       NA
#TCGA-05-4382                                                             NA                                                       yes                                                       NA
#             patient.follow_ups.follow_up.performance_status_scale_timing patient.follow_ups.follow_up.primary_therapy_outcome_success patient.follow_ups.follow_up.radiation_therapy
#                                                              <character>                                                  <character>                                    <character>
#TCGA-05-4249                                                           NA                                                           NA                                             no
#TCGA-05-4382                                                           NA                                       complete remission/r..                                            yes
#             patient.follow_ups.follow_up.targeted_molecular_therapy patient.karnofsky_performance_score patient.new_tumor_events.new_tumor_event.additional_pharmaceutical_therapy
#                                                         <character>                           <integer>                                                                <character>
#TCGA-05-4249                                                      no                                  NA                                                                         NA
#TCGA-05-4382                                                      no                                  NA                                                                         NA
#             patient.new_tumor_events.new_tumor_event.additional_radiation_therapy patient.performance_status_scale_timing patient.primary_therapy_outcome_success patient.radiation_therapy
#                                                                       <character>                             <character>                             <character>               <character>
#TCGA-05-4249                                                                    NA                                      NA                                      NA                        NA
#TCGA-05-4382                                                                    NA                                      NA                                      NA                        NA
#             patient.radiations.radiation.2.days_to_radiation_therapy_end patient.radiations.radiation.2.days_to_radiation_therapy_start patient.radiations.radiation.3.days_to_radiation_therapy_end
#                                                                <integer>                                                      <integer>                                                    <integer>
#TCGA-05-4249                                                           NA                                                             NA                                                           NA
#TCGA-05-4382                                                           NA                                                             NA                                                           NA
#             patient.radiations.radiation.3.days_to_radiation_therapy_start patient.radiations.radiation.4.days_to_radiation_therapy_end patient.radiations.radiation.4.days_to_radiation_therapy_start
#                                                                  <integer>                                                    <integer>                                                      <integer>
#TCGA-05-4249                                                             NA                                                           NA                                                             NA
#TCGA-05-4382                                                             NA                                                           NA                                                             NA
#             patient.radiations.radiation.days_to_radiation_therapy_end patient.radiations.radiation.days_to_radiation_therapy_start patient.targeted_molecular_therapy patient.tobacco_smoking_history
#                                                              <integer>                                                    <integer>                        <character>                     <character>
#TCGA-05-4249                                                         NA                                                           NA                                 NA          current reformed smo..
#TCGA-05-4382                                                        393                                                          365                                 NA          current reformed smo..
#             patient.year_of_tobacco_smoking_onset patient.samples.sample.2.portions.portion.2.slides.slide.2.percent_normal_cells patient.samples.sample.2.portions.portion.2.slides.slide.2.percent_stromal_cells
#                                         <integer>                                                                       <integer>                                                                        <integer>
#TCGA-05-4249                                    NA                                                                              NA                                                                               NA
#TCGA-05-4382                                    NA                                                                              NA                                                                               NA
#             patient.samples.sample.2.portions.portion.2.slides.slide.2.percent_tumor_cells patient.samples.sample.2.portions.portion.2.slides.slide.percent_normal_cells
#                                                                                  <integer>                                                                     <integer>
#TCGA-05-4249                                                                             NA                                                                            NA
#TCGA-05-4382                                                                             NA                                                                            NA
#             patient.samples.sample.2.portions.portion.2.slides.slide.percent_stromal_cells patient.samples.sample.2.portions.portion.2.slides.slide.percent_tumor_cells
#                                                                                  <integer>                                                                    <integer>
#TCGA-05-4249                                                                             NA                                                                           NA
#TCGA-05-4382                                                                             NA                                                                           NA
#             patient.samples.sample.2.portions.portion.slides.slide.2.percent_normal_cells patient.samples.sample.2.portions.portion.slides.slide.2.percent_stromal_cells
#                                                                                 <integer>                                                                      <integer>
#TCGA-05-4249                                                                            NA                                                                             NA
#TCGA-05-4382                                                                            NA                                                                             NA
#             patient.samples.sample.2.portions.portion.slides.slide.2.percent_tumor_cells patient.samples.sample.2.portions.portion.slides.slide.percent_normal_cells
#                                                                                <integer>                                                                   <integer>
#TCGA-05-4249                                                                           NA                                                                          NA
#TCGA-05-4382                                                                           NA                                                                          NA
#             patient.samples.sample.2.portions.portion.slides.slide.percent_stromal_cells patient.samples.sample.2.portions.portion.slides.slide.percent_tumor_cells
#                                                                                <integer>                                                                  <integer>
#TCGA-05-4249                                                                           NA                                                                         NA
#TCGA-05-4382                                                                           NA                                                                         NA
#             patient.samples.sample.3.portions.portion.slides.slide.percent_normal_cells patient.samples.sample.3.portions.portion.slides.slide.percent_stromal_cells
#                                                                               <integer>                                                                    <integer>
#TCGA-05-4249                                                                          NA                                                                           NA
#TCGA-05-4382                                                                          NA                                                                           NA
#             patient.samples.sample.3.portions.portion.slides.slide.percent_tumor_cells patient.samples.sample.4.portions.portion.2.slides.slide.percent_normal_cells
#                                                                              <integer>                                                                     <integer>
#TCGA-05-4249                                                                         NA                                                                            NA
#TCGA-05-4382                                                                         NA                                                                            NA
#             patient.samples.sample.4.portions.portion.2.slides.slide.percent_stromal_cells patient.samples.sample.4.portions.portion.2.slides.slide.percent_tumor_cells
#                                                                                  <integer>                                                                    <integer>
#TCGA-05-4249                                                                             NA                                                                           NA
#TCGA-05-4382                                                                             NA                                                                           NA
#             patient.samples.sample.4.portions.portion.3.slides.slide.percent_normal_cells patient.samples.sample.4.portions.portion.3.slides.slide.percent_stromal_cells
#                                                                                 <integer>                                                                      <integer>
#TCGA-05-4249                                                                            NA                                                                             NA
#TCGA-05-4382                                                                            NA                                                                             NA
#             patient.samples.sample.4.portions.portion.3.slides.slide.percent_tumor_cells patient.samples.sample.4.portions.portion.slides.slide.percent_normal_cells
#                                                                                <integer>                                                                   <integer>
#TCGA-05-4249                                                                           NA                                                                          NA
#TCGA-05-4382                                                                           NA                                                                          NA
#             patient.samples.sample.4.portions.portion.slides.slide.percent_stromal_cells patient.samples.sample.4.portions.portion.slides.slide.percent_tumor_cells
#                                                                                <integer>                                                                  <integer>
#TCGA-05-4249                                                                           NA                                                                         NA
#TCGA-05-4382                                                                           NA                                                                         NA
#             patient.samples.sample.portions.portion.2.slides.slide.2.percent_normal_cells patient.samples.sample.portions.portion.2.slides.slide.2.percent_stromal_cells
#                                                                                 <integer>                                                                      <integer>
#TCGA-05-4249                                                                            NA                                                                             NA
#TCGA-05-4382                                                                            NA                                                                             NA
#             patient.samples.sample.portions.portion.2.slides.slide.2.percent_tumor_cells patient.samples.sample.portions.portion.2.slides.slide.percent_normal_cells
#                                                                                <integer>                                                                   <integer>
#TCGA-05-4249                                                                           NA                                                                          NA
#TCGA-05-4382                                                                           NA                                                                          NA
#             patient.samples.sample.portions.portion.2.slides.slide.percent_stromal_cells patient.samples.sample.portions.portion.2.slides.slide.percent_tumor_cells
#                                                                                <integer>                                                                  <integer>
#TCGA-05-4249                                                                           NA                                                                         NA
#TCGA-05-4382                                                                           NA                                                                         NA
#             patient.samples.sample.portions.portion.3.slides.slide.percent_normal_cells patient.samples.sample.portions.portion.3.slides.slide.percent_stromal_cells
#                                                                               <integer>                                                                    <integer>
#TCGA-05-4249                                                                          NA                                                                           NA
#TCGA-05-4382                                                                          NA                                                                           NA
#             patient.samples.sample.portions.portion.3.slides.slide.percent_tumor_cells patient.samples.sample.portions.portion.slides.slide.2.percent_normal_cells
#                                                                              <integer>                                                                   <integer>
#TCGA-05-4249                                                                         NA                                                                           0
#TCGA-05-4382                                                                         NA                                                                           0
#             patient.samples.sample.portions.portion.slides.slide.2.percent_stromal_cells patient.samples.sample.portions.portion.slides.slide.2.percent_tumor_cells
#                                                                                <integer>                                                                  <integer>
#TCGA-05-4249                                                                           20                                                                         80
#TCGA-05-4382                                                                           10                                                                         NA
#             patient.samples.sample.portions.portion.slides.slide.percent_normal_cells patient.samples.sample.portions.portion.slides.slide.percent_stromal_cells patient.samples.sample.portions.portion.slides.slide.percent_tumor_cells
#                                                                             <integer>                                                                  <integer>                                                                <integer>
#TCGA-05-4249                                                                         0                                                                         20                                                                       80
#TCGA-05-4382                                                                         0                                                                         10                                                                       NA

table(colData(tmp_luad)[,"radiation_therapy"])
#> table(colData(tmp_luad)[,"radiation_therapy"])
#
# no yes 
#409  61

## Clinical info
cli_luad <- as.data.frame(colData(tmp_luad))
cli_lusc <- as.data.frame(colData(tmp_lusc))

# Key fields...

# Smoking
#> table(cli_luad$Smoking.Status)
#
#Current reformed smoker for < or = 15 years      Current reformed smoker for > 15 years                              Current smoker 
#                                         73                                          69                                          45 
#                        Lifelong Non-smoker 
#                                         32 
#> table(cli_lusc$Smoking.Status)
#
#Current reformed smoker for < or = 15 years      Current reformed smoker for > 15 years                              Current smoker 
#                                         87                                          50                                          28 
#                        Lifelong Non-smoker                                         N/A 
#                                          7                                           6
                                          
# ECOG-like performance
#> colnames(cli_lusc)[grep("performance", colnames(cli_lusc))]
#[1] "karnofsky_performance_score"                                    "patient.follow_ups.follow_up.2.karnofsky_performance_score"    
#[3] "patient.follow_ups.follow_up.2.performance_status_scale_timing" "patient.follow_ups.follow_up.3.karnofsky_performance_score"    
#[5] "patient.follow_ups.follow_up.karnofsky_performance_score"       "patient.follow_ups.follow_up.performance_status_scale_timing"  
#[7] "patient.karnofsky_performance_score"                            "patient.performance_status_scale_timing"

# -> mapping https://oncologypro.esmo.org/oncology-in-practice/practice-tools/performance-scales
# As published in Am J Clin. Oncol.: Oken, M.M., Creech, R.H., Tormey, D.C., Horton, J., Davis, T.E., McFadden, E.T., Carbone, P.P.: Toxicity And Response Criteria Of The Eastern Cooperative Oncology Group. Am J Clin Oncol 5:649-655, 1982.The Eastern Cooperative Oncology Group, Robert Comis M.D., Group Chair.
#
# Karnofsky 100 >= x >= 90 -> ECOG = 0
# Karnofsky  80 >= x >= 70 -> ECOG = 1
# Karnofsky  60 >= x >= 50 -> ECOG = 2
# Karnofsky  40 >= x >= 30 -> ECOG = 3
# Karnosfky  20 >= x >= 10 -> ECOG = 4
# Karnofsky  0, dead       -> ECOG = 5
           
#> table(cli_lusc$patient.performance_status_scale_timing, cli_lusc$patient.karnofsky_performance_score, useNA="ifany")
#                       
#                          0  20  40  50  70  80  90 100 <NA>
#  other                   3   0   0   0   0   1   0   0   11
#  post-adjuvant therapy   0   0   0   2   0   1   2   0    7
#  pre-adjuvant therapy    0   0   0   0   0   0   2   0   21
#  pre-operative           0   0   0   0   5  10  27  10   90
#  <NA>                   30   1   1   0   3   4   1   0  269
#> table(cli_luad$patient.performance_status_scale_timing, cli_luad$patient.karnofsky_performance_score, useNA="ifany")
#                       
#                          0  40  60  70  80  90 100 <NA>
#  other                   2   1   0   0   1   3   0   16
#  post-adjuvant therapy   0   0   0   0   1   1   0    9
#  pre-adjuvant therapy    0   0   0   0   0   0   0   20
#  pre-operative           0   1   1   5  16  27  32   94
#  <NA>                    2   0   1   0   5   1   0  277

#ifelse(cli_luad$patient.performance_status_scale_timing == "pre-operative", 7-findInterval(cli_luad$patient.karnofsky_performance_score, c(-1, 9, 29, 49, 69, 89, 100)), NA)
#7-findInterval(cli_luad$patient.karnofsky_performance_score, c(-1, 9, 20, 40, 60, 79, 100))                           
table(ECOG = ifelse(cli_luad$patient.performance_status_scale_timing == "pre-operative", 6-findInterval(cli_luad$patient.karnofsky_performance_score, c(-1, 9, 29, 49, 69, 89, 101)), NA),
	Karnofsky = cli_luad$patient.karnofsky_performance_score,
	timing = cli_luad$patient.performance_status_scale_timing
)
# ...
#
#, , timing = pre-operative
#
#    Karnofsky
#ECOG  0 40 60 70 80 90 100
#   0  0  0  0  0  0 27  32
#   1  0  0  0  5 16  0   0
#   2  0  0  1  0  0  0   0
#   3  0  1  0  0  0  0   0   

  
# Gender:
#> table(cli_lusc$patient.gender, useNA="ifany")
#
#female   male 
#   130    371

#> table(cli_luad$patient.person_neoplasm_cancer_status)
#
#tumor free with tumor 
#       310        113

#> table(cli_luad$patient.new_tumor_events.new_tumor_event_after_initial_treatment, useNA="ifany")
#
#  no  yes <NA> 
# 155   30  331

#> table(cli_luad$patient.follow_ups.follow_up.new_tumor_event_after_initial_treatment, useNA="ifany")
#
#  no  yes <NA> 
# 276  144   96

#> table(cli_luad$patient.follow_ups.follow_up.new_tumor_event_after_initial_treatment, cli_luad$patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment, useNA="ifany")
#      
#        15  18  21  29  35  42  54  60  96 127 132 139 150 153 158 162 164 177 179 182 183 184 195 199 209 213 214 216 221 224 231 232 238 245 246 251
#  no     0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#  yes    1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2   2   1   1   1   1   1   1   1   1   1   2   1   1   2   1   1
#  <NA>   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
#
#       925 1018 1083 1090 1202 1255 1433 1500 1568 2045 2218 3521 4812 <NA>
#  no     0    0    0    0    0    0    0    0    0    0    0    0    0  276
#  yes    1    1    1    1    1    1    1    1    1    1    1    1    1    9
#  <NA>   0    0    0    0    0    0    0    0    0    0    0    0    0   96
#
## ... -> PFS


# For subchallenge 3, in https://www.synapse.org/#!Synapse:syn18404605/wiki/607473
#
# - For Sub-challenge 3, do we just predict the patient will be PD or not PD? Or PD, Complete Response and Partial Response?
#
# > Here the response variable is binary. PD (progressive disease) will be treated as non-responder, 
# NE (not evaluable) will be set as NA, and so, not used in the scoring calculation. All others are considered responders.

#> table(Responder = c(0, 1, 1, 1, NA)[match(cli_luad[,"patient.drugs.drug.measure_of_response"], c("clinical progressive disease", "complete response", "partial response", "stable disease", NA))], raw = cli_luad[,"patient.drugs.drug.measure_of_response"], useNA="ifany")
#         raw
#Responder clinical progressive disease complete response partial response stable disease <NA>
#     0                              22                 0                0              0    0
#     1                               0                56                7             10    0
#     <NA>                            0                 0                0              0  421

#> table(Responder = c(0, 1, 1, 1, NA)[match(cli_lusc[,"patient.drugs.drug.measure_of_response"], c("clinical progressive disease", "complete response", "partial response", "stable disease", NA))], raw = cli_lusc[,"patient.drugs.drug.measure_of_response"], useNA="ifany")
#         raw
#Responder clinical progressive disease complete response partial response stable disease <NA>
#     0                              13                 0                0              0    0
#     1                               0                48                2              3    0
#     <NA>                            0                 0                0              0  435

# Clean up clinical information from cli_luad and cli_lusc
# No obvious immunohistochemistry or PD* available for IO
#> grep("IHC", colnames(cli_luad))
#integer(0)
#> grep("PD", colnames(cli_luad))
#integer(0)
#> grep("IHC", colnames(cli_lusc))
#integer(0)
#> grep("PD", colnames(cli_lusc))
#[1] 2196 2212 2241
#> colnames(cli_lusc)[grep("PD", colnames(cli_lusc))]
#[1] "Mutation_PDYN"   "Mutation_PDGFRA" "CNA_PDGFRA"


# Summarize given treatments as additional fields
#> table(unlist(cli_luad[,grep("therapy_type", colnames(cli_luad), value=TRUE)]))
#
#                                 ancillary aurora kinase inhibitor - protocol therapy                               chemotherapy                     gsk mage vaccine study                              immunotherapy 
#                                         6                                          1                                        429                                          1                                          6 
#                            other, specify            phase ii clinical trial krw2170                 targeted molecular therapy                                    vaccine 
#                                         3                                          1                                         11                                          1


# Follow-up times scattered along multiple columns
# Event: cli_lusc$patient.clinical_cqcf.consent_or_death_status
# Time somewhere possibly in: cbind(DTD = cli_lusc$patient.clinical_cqcf.days_to_death, FU.1 = cli_lusc$patient.follow_ups.follow_up.days_to_death, FU.2 = cli_lusc$patient.follow_ups.follow_up.2.days_to_death, FU.3 = cli_lusc$patient.follow_ups.follow_up.3.days_to_death, FU.4 = cli_lusc$patient.follow_ups.follow_up.4.days_to_death)
# or rather in any of the "days_to"-fields, with varying specificity
omit.infinite <- function(x) { ifelse(is.finite(x), x, NA) }

dat_luad <- data.frame(
	patientID = as.character(paste(rownames(cli_luad), ".01", sep="")),
	SEX = as.factor(ifelse(cli_luad$gender == "female", "F", "M")),
	AAGE = as.numeric(cli_luad$patient.age_at_initial_pathologic_diagnosis),
	CRFHIST = factor("NON-SQUAMOUS", levels=c("NON-SQUAMOUS", "SQUAMOUS")),
	TOBACUSE = factor(c("CURRENT", "FORMER", "FORMER", "FORMER", "NEVER", "UNKNOWN")[match(cli_luad$patient.tobacco_smoking_history, c("current smoker", "current reformed smoker for < or = 15 years", "current reformed smoker for > 15 years", "current reformed smoker, duration not specified", "lifelong non-smoker", NA))]),
	ECOGPS = as.numeric(ifelse(cli_luad$patient.performance_status_scale_timing == "pre-operative", 6-findInterval(cli_luad$patient.karnofsky_performance_score, c(-1, 9, 29, 49, 69, 89, 101)), NA)),
	PDL1 = as.numeric(NA),
	TMB = as.numeric(cli_luad[,"Nonsilent.Mutations.per.Mb"]),
	TCR_Shannon = as.numeric(NA),
	TCR_Richness = as.numeric(NA),
	TCR_Evenness = as.numeric(NA),
	BCR_Shannon = as.numeric(NA),
	BCR_Richness = as.numeric(NA),
	BCR_Evenness = as.numeric(NA),
	PFS.time = as.numeric(cli_luad$patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment),
	PFS.event = as.numeric(ifelse(cli_luad$patient.follow_ups.follow_up.new_tumor_event_after_initial_treatment == "yes", 1, 0)),
	#OS = as.numeric(cli_luad$days_to_last_followup),
	OS.time = omit.infinite(as.numeric(apply(cli_luad[,grep("days_to", colnames(cli_luad), value=TRUE)[1:3]], MARGIN=1, FUN=function(x) { max(x, na.rm=TRUE) }))),
	OS.event = ifelse(cli_luad$patient.clinical_cqcf.consent_or_death_status == "deceased", 1, 0),
	Responder = c(0, 1, 1, 1, NA)[match(cli_luad[,"patient.drugs.drug.measure_of_response"], c("clinical progressive disease", "complete response", "partial response", "stable disease", NA))],
	# Reported chemotherapy therapy
	tChemo = as.integer(apply(cli_luad[,grep("therapy_type", colnames(cli_luad), value=TRUE)], MARGIN=1, FUN=function(z) { any(z == "chemotherapy") })),
	# Reported radiation therapy
	#tRadia = as.integer(cli_luad$radiation_therapy == "yes"),
	tRadia = as.integer(apply(cli_luad[,grep("radiation_therapy|radiation_treatment_ongoing", colnames(cli_luad), value=TRUE)], MARGIN=1, FUN=function(x) { any(x == "yes") })),	
	# Reported other therapy
	#tOther = as.integer(apply(cli_luad[,grep("therapy_type", colnames(cli_luad), value=TRUE)], MARGIN=1, FUN=function(z) { any(!z == "chemotherapy") }))
	tOther = apply(cli_luad[,grep("therapy_type|targeted_molecular_therapy", colnames(cli_luad), value=TRUE)], MARGIN=1, FUN=function(z) { as.integer(!all(z %in% c("chemotherapy", "no", NA))) })
)
# PFS times missing for patients that did not have a reported progression, filling with the longest known survival in those cases
dat_luad[,"PFS.time"] <- ifelse(is.na(dat_luad$PFS.event), NA, omit.infinite(apply(dat_luad[,c("PFS.time", "OS.time")], MARGIN=1, FUN=function(x) { min(x, na.rm=TRUE) })))
# Fill in NAs in treatments as non-observed
dat_luad[is.na(dat_luad$tChemo),"tChemo"] <- 0
dat_luad[is.na(dat_luad$tRadia),"tRadia"] <- 0
dat_luad[is.na(dat_luad$tOther),"tOther"] <- 0
# Shorten IDs
rownames(dat_luad) <- dat_luad$patientID <- substr(dat_luad$patientID, 1, 12)

#> head(dat_luad)
#        patientID SEX AAGE      CRFHIST TOBACUSE ECOGPS PDL1   TMB TCR_Shannon TCR_Richness TCR_Evenness BCR_Shannon BCR_Richness BCR_Evenness PFS.time PFS.event OS.time OS.event Responder
#1 TCGA-05-4249.01   M   67 NON-SQUAMOUS   FORMER     NA   NA  8.47          NA           NA           NA          NA           NA           NA     1523         0    1523        0        NA
#2 TCGA-05-4382.01   M   68 NON-SQUAMOUS   FORMER     NA   NA 37.78          NA           NA           NA          NA           NA           NA      334         1     607        0        NA
#3 TCGA-05-4384.01   M   66 NON-SQUAMOUS   FORMER     NA   NA  3.50          NA           NA           NA          NA           NA           NA      183         1     426        0         0
#4 TCGA-05-4389.01   M   70 NON-SQUAMOUS   FORMER     NA   NA  6.43          NA           NA           NA          NA           NA           NA     1369         0    1369        0        NA
#5 TCGA-05-4390.01   F   58 NON-SQUAMOUS   FORMER     NA   NA 14.41          NA           NA           NA          NA           NA           NA      395         1    1126        0         1
#6 TCGA-05-4395.01   M   76 NON-SQUAMOUS   FORMER     NA   NA  6.34          NA           NA           NA          NA           NA           NA       NA        NA       0        0        NA

dat_lusc <- data.frame(
	patientID = as.character(paste(rownames(cli_lusc), ".01", sep="")),
	SEX = as.factor(ifelse(cli_lusc$gender == "female", "F", "M")),
	AAGE = as.numeric(cli_lusc$patient.age_at_initial_pathologic_diagnosis),
	CRFHIST = factor("SQUAMOUS", levels=c("NON-SQUAMOUS", "SQUAMOUS")),
	TOBACUSE = factor(c("CURRENT", "FORMER", "FORMER", "FORMER", "NEVER", "UNKNOWN")[match(cli_lusc$patient.tobacco_smoking_history, c("current smoker", "current reformed smoker for < or = 15 years", "current reformed smoker for > 15 years", "current reformed smoker, duration not specified", "lifelong non-smoker", NA))]),
	ECOGPS = as.numeric(ifelse(cli_lusc$patient.performance_status_scale_timing == "pre-operative", 6-findInterval(cli_lusc$patient.karnofsky_performance_score, c(-1, 9, 29, 49, 69, 89, 101)), NA)),
	PDL1 = as.numeric(NA),
	TMB = as.numeric(cli_lusc[,"Nonsilent.Mutatios.per.Mb"]), # Typo [sic]
	TCR_Shannon = as.numeric(NA),
	TCR_Richness = as.numeric(NA),
	TCR_Evenness = as.numeric(NA),
	BCR_Shannon = as.numeric(NA),
	BCR_Richness = as.numeric(NA),
	BCR_Evenness = as.numeric(NA),
	PFS.time = as.numeric(cli_lusc$patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment),
	PFS.event = as.numeric(ifelse(cli_lusc$patient.follow_ups.follow_up.new_tumor_event_after_initial_treatment == "yes", 1, 0)),
	#OS = as.numeric(cli_lusc$days_to_last_followup),
	OS.time = omit.infinite(as.numeric(apply(cli_lusc[,grep("days_to", colnames(cli_lusc), value=TRUE)[1:3]], MARGIN=1, FUN=function(x) { max(x, na.rm=TRUE) }))),
	OS.event = ifelse(cli_lusc$patient.clinical_cqcf.consent_or_death_status == "deceased", 1, 0),
	Responder = c(0, 1, 1, 1, NA)[match(cli_lusc[,"patient.drugs.drug.measure_of_response"], c("clinical progressive disease", "complete response", "partial response", "stable disease", NA))],
	# Reported chemotherapy therapy
	tChemo = as.integer(apply(cli_lusc[,grep("therapy_type", colnames(cli_lusc), value=TRUE)], MARGIN=1, FUN=function(z) { any(z == "chemotherapy") })),
	# Reported radiation therapy
	#tRadia = as.integer(cli_lusc$radiation_therapy == "yes"),
	tRadia = as.integer(apply(cli_lusc[,grep("radiation_therapy|radiation_treatment_ongoing", colnames(cli_lusc), value=TRUE)], MARGIN=1, FUN=function(x) { any(x == "yes") })),
	# Reported other therapy
	#tOther = as.integer(apply(cli_lusc[,grep("therapy_type", colnames(cli_lusc), value=TRUE)], MARGIN=1, FUN=function(z) { any(!z == "chemotherapy") }))
	tOther = apply(cli_lusc[,grep("therapy_type|targeted_molecular_therapy", colnames(cli_lusc), value=TRUE)], MARGIN=1, FUN=function(z) { as.integer(!all(z %in% c("chemotherapy", "no", NA))) })
	
)	
# PFS times missing for patients that did not have a reported progression, filling with the longest known survival in those cases
dat_lusc[,"PFS.time"] <- ifelse(is.na(dat_lusc$PFS.event), NA, omit.infinite(apply(dat_lusc[,c("PFS.time", "OS.time")], MARGIN=1, FUN=function(x) { min(x, na.rm=TRUE) })))
# Fill in NAs in treatments as non-observed
dat_lusc[is.na(dat_lusc$tChemo),"tChemo"] <- 0
dat_lusc[is.na(dat_lusc$tRadia),"tRadia"] <- 0
dat_lusc[is.na(dat_lusc$tOther),"tOther"] <- 0
# Shorten IDs
rownames(dat_lusc) <- dat_lusc$patientID <- substr(dat_lusc$patientID, 1, 12)

#> head(dat_lusc)
#        patientID SEX AAGE  CRFHIST TOBACUSE ECOGPS PDL1       TMB TCR_Shannon TCR_Richness TCR_Evenness BCR_Shannon BCR_Richness BCR_Evenness PFS.time PFS.event OS.time OS.event Responder
#1 TCGA-18-3406.01   M   67 SQUAMOUS   FORMER     NA   NA  8.068481          NA           NA           NA          NA           NA           NA      357         1     371        1        NA
#2 TCGA-18-3407.01   M   72 SQUAMOUS   FORMER     NA   NA  6.766357          NA           NA           NA          NA           NA           NA      136         0     136        1        NA
#3 TCGA-18-3408.01   F   77 SQUAMOUS   FORMER     NA   NA  2.965731          NA           NA           NA          NA           NA           NA     1793         1    2304        0        NA
#4 TCGA-18-3409.01   M   74 SQUAMOUS   FORMER     NA   NA 69.597030          NA           NA           NA          NA           NA           NA     2291         1    3747        0        NA
#5 TCGA-18-3410.01   M   81 SQUAMOUS   FORMER     NA   NA 10.451920          NA           NA           NA          NA           NA           NA      146         0     146        1        NA
#6 TCGA-18-3411.01   F   63 SQUAMOUS  CURRENT     NA   NA 10.561520          NA           NA           NA          NA           NA           NA       NA        NA    3576        0        NA

# Combine datasets to produce combined SQUAMOUS & NON-SQUAMOUS TCGA
# row-bind
dat_tcga <- rbind(dat_luad, dat_lusc)
# column-bind and transpose gene expression
gex_tcga <- cbind(gex_luad, gex_lusc)
# Matching rownames, everything's looking good!
#> sum(rownames(gex_tcga) %in% rownames(dat_tcga))
#[1] 1128
#> sum(!rownames(gex_tcga) %in% rownames(dat_tcga))
#[1] 0
#> sum(!rownames(dat_tcga) %in% rownames(gex_tcga))
#[1] 0

## Updated: subset TCGA cohort to only patients that have been reportedly treated with chemo:
#
#> table(Chemo = dat_tcga$tChemo, Squamous = dat_tcga$CRFHIST, useNA="ifany")
#     Squamous
#Chemo NON-SQUAMOUS SQUAMOUS
#    0          342      361
#    1          174      140
#
## Relatively well represented in both cohorts
## Subsetting based on binary chemo indicator
gex_tcga <- gex_tcga[,rownames(dat_tcga[which(dat_tcga[,"tChemo"]==1),])]
#> dim(gex_tcga)
#[1] 20501  1128
#> gex_tcga <- gex_tcga[,rownames(dat_tcga[which(dat_tcga[,"tChemo"]==1),])]
#> dim(gex_tcga)
#[1] 20501   314
## N down to 314
dat_tcga <- dat_tcga[which(dat_tcga[,"tChemo"]==1),]
#> dim(dat_tcga)
#[1] 314  22

# Save TCGA temporary datasets
setwd("./RData")
save(dat_tcga, file="dat_tcga.RData")
save(gex_tcga, file="gex_tcga.RData")
setwd("..")


## Datasets from GEO
library(GEOquery)
library(DESeq2)
# Follow guidelines e.g. in http://genomicsclass.github.io/book/pages/GEOquery.html
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Subdirectory for GEO stuff
setwd("GEO")

# Regarding gene length normalization in DESeq2: https://www.biostars.org/p/140090/

# Extract GPLs https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/

# GSE115821 - Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115821
gse_auslander <- GEOquery::getGEO("GSE115821", GSEMatrix = TRUE)
sup_auslander <- GEOquery::getGEOSuppFiles("GSE115821")
GEOquery::gunzip(rownames(sup_auslander)[1])
sup_auslander <- read.csv("GSE115821/GSE115821_MGH_counts.csv")
pData(gse_auslander[[1]])
## -> Only one responder, rest are non-responders

# GSE121810 - Neoadjuvant anti-PD-1 immunotherapy promotes a survival benefit with intratumoral and systemic immune responses in recurrent glioblastoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121810
gse_cloughesy <- GEOquery::getGEO("GSE121810", GSEMatrix = TRUE)
sup_cloughesy <- GEOquery::getGEOSuppFiles("GSE121810")
sup_cloughesy <- readxl::read_excel(rownames(sup_cloughesy)[1])
pData(gse_cloughesy[[1]])
## -> No endpoints available??

# GSE78220 - mRNA expressions in pre-treatment melanomas undergoing anti-PD-1 checkpoint inhibition therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220
gse_hugo <-  GEOquery::getGEO("GSE78220", GSEMatrix = TRUE)
sup_hugo <- GEOquery::getGEOSuppFiles("GSE78220")
sup_hugo <- readxl::read_excel(rownames(sup_hugo)[1])
# sup_hugo -> tibble table, cast to matrix preferably
gex_hugo <- as.matrix(sup_hugo)
rownames(gex_hugo) <- gex_hugo[,"Gene"]
gex_hugo <- gex_hugo[,-1]
class(gex_hugo) <- "numeric"
# Clinical data
cli_hugo <- pData(gse_hugo[[1]])
# Omit patient who had already been on treatment
cli_hugo <- cli_hugo[-grep("on-treatment", cli_hugo[,"biopsy time:ch1"]),]
gex_hugo <- gex_hugo[,-grep("OnTx", colnames(gex_hugo))]
colnames(gex_hugo) <- gsub(".baseline", "", colnames(gex_hugo))
# Create a dat info frame
dat_hugo <- data.frame(
	patientID = cli_hugo[,"title"],
	SEX = as.factor(cli_hugo[,"gender:ch1"]),
	AAGE = as.integer(cli_hugo[,"age (yrs):ch1"]),
	CRFHIST = factor(NA, levels=c("NON-SQUAMOUS", "SQUAMOUS")),
	TOBACUSE = factor("UNKNOWN", levels = c("CURRENT", "FORMER", "NEVER", "UNKNOWN")),
	ECOGPS = as.integer(NA),
	PDL1 = as.integer(NA),
	TMB = as.numeric(NA),
	TCR_Shannon = as.numeric(NA),
	TCR_Richness = as.numeric(NA),
	TCR_Evenness = as.numeric(NA),
	BCR_Shannon = as.numeric(NA),
	BCR_Richness = as.numeric(NA),
	BCR_Evenness = as.numeric(NA),	
	PFS.time = NA,
	PFS.event = NA,
	OS.time = as.integer(cli_hugo[,"overall survival (days):ch1"]),
	OS.event = as.integer(cli_hugo[,"vital status:ch1"]=="Dead"),
	#### NOTE! Only progressive disease --> 0, others -> 1
	#Responder = as.integer(cli_hugo[,"anti-pd-1 response:ch1"] %in% c("Complete Response", "Partial Response")),
	#Responder = as.integer(cli_hugo[,"anti-pd-1 response:ch1"] == "Complete Response"),
	Responder = 1 - as.integer(cli_hugo[,"anti-pd-1 response:ch1"] == "Progressive Disease"),
	MetasLocation = cli_hugo[,"anatomical location:ch1"]
)
rownames(dat_hugo) <- dat_hugo$patientID
# > all(colnames(gex_hugo) == rownames(dat_hugo))
# [1] TRUE
# Save RDatas
save(gex_hugo, file="./RData/gex_hugo.RData")
save(dat_hugo, file="./RData/dat_hugo.RData")


# GSE52562 - Gene expression profiling of tumor biopsies before and after pidilizumab therapy in patients with relapsed follicular lymphoma grade 1 or grade 2.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52562
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558
gse_westin <-  GEOquery::getGEO("GSE52562", GSEMatrix = TRUE, getGPL = TRUE)
# Mapping/collapsing:
# gse_westin[[1]]@featureData@data[1:2,]
# column "ID" <-> "ILMN_Gene"
cli_westin <- pData(gse_westin[[1]])
# Quantile normalized signals
gex_westin <- as.matrix(exprs(gse_westin[[1]]))
gpl_westin <- read.table("./Annotations/GPL10558-50081.txt", sep="\t", skip=30, quote="§", comment="§", header=TRUE)
gpl_westin <- gpl_westin[,c("ID", "Transcript", "ILMN_Gene", "RefSeq_ID", "Unigene_ID")]
#> all(rownames(gex_westin) %in% gpl_westin[,"ID"])
#[1] TRUE
## Rather old bead array; issues in detecting with probes, mapping to gene symbols
gpl_westin <- gpl_westin[match(rownames(gex_westin), gpl_westin[,"ID"]),]
gex_westin <- gex_westin[!gpl_westin[,"ILMN_Gene"]=="",]
gpl_westin <- gpl_westin[!gpl_westin[,"ILMN_Gene"]=="",]
rownames(gex_westin) <- gpl_westin[,"ILMN_Gene"]
gex_westin <- gex_westin[order(rownames(gex_westin)),]
gex_westin <- gex_westin[-grep("-Mar|-Dec", rownames(gex_westin)),]
# Collapse same symbols using mean expression for probes
gex_westin <- do.call("rbind", by(gex_westin, INDICES=rownames(gex_westin), FUN=function(z){ apply(z, MARGIN=2, FUN=mean) }))
# Create dat
dat_westin <- data.frame(
	patientID = cli_westin[,"title"],
	SEX = as.factor(cli_westin[,"gender:ch1"]),
	AAGE = as.integer(cli_westin[,"age:ch1"]),
	CRFHIST = factor(NA, levels=c("NON-SQUAMOUS", "SQUAMOUS")),
	TOBACUSE = factor("UNKNOWN", levels = c("CURRENT", "FORMER", "NEVER", "UNKNOWN")),
	ECOGPS = as.integer(NA),
	PDL1 = as.integer(NA),
	TMB = as.numeric(NA),
	TCR_Shannon = as.numeric(NA),
	TCR_Richness = as.numeric(NA),
	TCR_Evenness = as.numeric(NA),
	BCR_Shannon = as.numeric(NA),
	BCR_Richness = as.numeric(NA),
	BCR_Evenness = as.numeric(NA),	
	PFS.time = as.integer(cli_westin[,"pfs.days:ch1"]),
	PFS.event = as.integer(cli_westin[,"pfs.status.censorship:ch1"]),
	OS.time = as.integer(NA),
	OS.event = as.integer(NA),
	Responder = as.integer(NA)
)
rownames(dat_westin) <- rownames(cli_westin)
# Only pick pre-treatment samples
gex_westin <- gex_westin[,grep("pre", dat_westin[,"patientID"])]
dat_westin <- dat_westin[grep("pre", dat_westin[,"patientID"]),]
#> all(colnames(gex_westin) == rownames(dat_westin))
#[1] TRUE
# Save RDatas
save(gex_westin, file="./RData/gex_westin.RData")
save(dat_westin, file="./RData/dat_westin.RData")

# GSE79691 - Transcriptional mechanisms of resistance to anti-PD-1 therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67501
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL18281
gse_ascierto <- GEOquery::getGEO("GSE67501", GSEMatrix = TRUE, getGPL = TRUE)
# Mapping/collapsing:
# gse_ascierto[[1]]@featureData@data[1:2,]
# column "ID" <-> "ILMN_Gene"

# GSE79691 - Transcriptional mechanisms of resistance to anti-PD-1 therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79691
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL14951
gse_ascierto2 <- GEOquery::getGEO("GSE79691", GSEMatrix = TRUE, getGPL = TRUE)
# Mapping/collapsing:
# gse_ascierto2[[1]]@featureData@data[1:2,]
# column "ID" <-> "ILMN_Gene"

# GSE91061 - Molecular portraits of tumor mutational and micro-environmental sculpting by immune checkpoint blockade therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
# Note: 'Gene IDs in the processed data files are NCBI Entrez Gene IDs.'
gse_riaz <- GEOquery::getGEO("GSE91061", GSEMatrix = TRUE, getGPL = TRUE)
sup_riaz <- GEOquery::getGEOSuppFiles("GSE91061")
lapply(rownames(sup_riaz), FUN=GEOquery::gunzip)
sup_riaz1 <- read.csv(gsub(".gz", "", rownames(sup_riaz)[1]))
sup_riaz2 <- read.csv(gsub(".gz", "", rownames(sup_riaz)[2]))
sup_riaz3 <- read.csv(gsub(".gz", "", rownames(sup_riaz)[3]))
sup_riaz4 <- read.csv(gsub(".gz", "", rownames(sup_riaz)[4]))
# sup_riaz2 are raw counts, using those with DESeq2
gpl_riaz <- getGEO('GPL9052', destdir=".")
##
## Riaz et al.
##
#rownames(sup_riaz2) <- genes[match(sup_riaz2[,1], genes[,"entrezgene_id"]),"hgnc_symbol"]
#> table(table(genes[match(sup_riaz2[,1], genes[,"entrezgene_id"]),"hgnc_symbol"]))
#
#    1     2     3   217 
#20584    25     1     1
#
# 20584 had unique hugo symbol mapping, filtering the rest ("" occurred 217 times)
sup_riaz2[,1] <- genes[match(sup_riaz2[,1], genes[,"entrezgene_id"]),"hgnc_symbol"]
sup_riaz2 <- sup_riaz2[-which(sup_riaz2[,1] %in% names(table(sup_riaz2[,1])>1)[which(table(sup_riaz2[,1])>1)]),]
sup_riaz2 <- sup_riaz2[-which(is.na(sup_riaz2[,1])),]
rownames(sup_riaz2) <- sup_riaz2[,1]
sup_riaz2 <- sup_riaz2[,-1]
# Substrata
str_riaz <- pData(gse_riaz[[1]])$`characteristics_ch1`
# DESeq2
gex_riaz <- DESeq2::DESeqDataSetFromMatrix(countData=sup_riaz2,
	colData = data.frame("str_riaz" = str_riaz),
	design = ~ str_riaz)
gex_riaz <- estimateSizeFactors(gex_riaz)
gex_riaz <- counts(gex_riaz, normalized=TRUE)
gex_riaz <- gex_riaz[order(rownames(gex_riaz)),]
## Clinical info
cli_riaz <- pData(gse_riaz[[1]])
rownames(cli_riaz) <- gsub("-", ".", cli_riaz[,1])
cli_riaz <- cli_riaz[match(colnames(gex_riaz), rownames(cli_riaz)),]
cli_riaz <- cli_riaz[grep("Pre", colnames(gex_riaz)),]
gex_riaz <- gex_riaz[,grep("Pre", colnames(gex_riaz))]
# Omit patients for whom responder status was unknown
gex_riaz <- gex_riaz[,-which(cli_riaz[,"response:ch1"]=="UNK")]
cli_riaz <- cli_riaz[-which(cli_riaz[,"response:ch1"]=="UNK"),]
# Genereta patient data frame
dat_riaz <- data.frame(
	patientID = rownames(cli_riaz),
	SEX = NA,
	AAGE = NA,
	CRFHIST = NA,
	TOBACUSE = NA,
	ECOGPS = NA,
	PDL1 = NA,
	TMB = NA,
	TCR_Shannon = NA,
	TCR_Richness = NA,
	TCR_Evenness = NA,
	BCR_Shannon = NA,
	BCR_Richness = NA,
	BCR_Evenness = NA,
	#Responder = as.integer(cli_riaz[,"response:ch1"] == "PRCR") # Partial/complete response, other options PD (progressive dis) or SD (stable dis)
	# Update! Responders should be defined as any non-progressive disease
	Responder = 1 - as.integer(cli_riaz[,"response:ch1"] == "PD")
)
rownames(dat_riaz) <- dat_riaz[,1]
#table(dat_riaz[,"Responder"])
setwd("..")
save(gex_riaz, file="./RData/gex_riaz.RData")
save(dat_riaz, file="./RData/dat_riaz.RData")
setwd("GEO")

# GSE93157 - Programmed death 1 receptor blockade and immune-related gene expression profiling in non-small cell lung carcinoma, head and neck squamous cell carcinoma and melanoma
# Very small gene panel ("A minimum of 100 ng of total RNA was used to measure the expression of 105 breast cancer-related genes and 5 house-keeping genes (ACTB, MRPL19, PSMC4, RLP0 and SF3A1)")
gse_prat <- GEOquery::getGEO("GSE93157", GSEMatrix = TRUE)
#sup_prat <- GEOquery::getGEOSuppFiles("GSE93157")
#GEOquery::gunzip(rownames(sup_prat)[1])
#gex_prat <- read.table("GSE93157\GSE93157_raw_data_values.txt")
gex_prat <- as.matrix(exprs(gse_prat[[1]]))
gex_prat <- gex_prat[-which(apply(gex_prat, MARGIN=1, FUN=function(z) { all(is.na(z)) })),]
#> dim(gex_prat)
#[1] 725  65
cli_prat <- pData(gse_prat[[1]])
# Study:
# Twenty-three immune-related genes or signatures were linked to response. The intra- and inter-biopsy variability of PD1, PDL1, CD8A and CD4 mRNA.
dat_prat <- data.frame(
	#patientID = cli_prat[,"title"],
	patientID = rownames(cli_prat),
	#SEX = factor(ifelse(cli_prat[,"characteristics_ch1.9"] == "Sex: M", "M", "F")),
	SEX = factor(cli_prat[,"Sex:ch1"]),
	AAGE = as.integer(cli_prat[,"age:ch1"]),
	CRFHIST = factor(ifelse(cli_prat[,"source_name_ch1"]=="SQUAMOUS LUNG CANCER", "SQUAMOUS", ifelse(cli_prat[,"source_name_ch1"]=="LUNG NON-SQUAMOUS CANCER", "NON-SQUAMOUS", NA)), levels=c("NON-SQUAMOUS", "SQUAMOUS")),
	TOBACUSE = c("NEVER", "FORMER", "CURRENT")[match(cli_prat[,"characteristics_ch1.7"], c("smoking: NS", "smoking: FS", "smoking: CS"))],
	ECOGPS = c(0,1,2)[match(cli_prat[,"characteristics_ch1.8"], c("ecog: 0", "ecog: 1", "ecog: 2"))],
	PDL1 = as.integer(NA),
	TMB = as.numeric(NA),
	TCR_Shannon = as.numeric(NA),
	TCR_Richness = as.numeric(NA),
	TCR_Evenness = as.numeric(NA),
	BCR_Shannon = as.numeric(NA),
	BCR_Richness = as.numeric(NA),
	BCR_Evenness = as.numeric(NA),	
	#PFS.event = as.integer(ifelse(cli_prat[,"characteristics_ch1.13"]=="pfse: 1")),
	PFS.event = as.integer(cli_prat[,"pfse:ch1"]),
	#PFS.time = as.integer(gsub("pfs: ", "", cli_prat[,"characteristics_ch1.14"])),
	PFS.time = as.numeric(cli_prat[,"pfs:ch1"]),
	Responder = as.integer(cli_prat[,"characteristics_ch1.11"] %in% c("best.resp: CR", "best.resp: PR")),
	Type = cli_prat[,"source_name_ch1"],
	#Drug = factor(gsub("drug: ", "", cli_prat[,"characteristics_ch1.10"]))
	Drug = cli_prat[,"drug:ch1"],
	LungCa = ifelse(cli_prat[,"source_name_ch1"] %in% c("LUNG NON-SQUAMOUS CANCER", "SQUAMOUS LUNG CANCER"), 1, 0)
)
rownames(dat_prat) <- dat_prat$patientID
#> all(rownames(dat_prat) == colnames(gex_prat))
#[1] TRUE
save(gex_prat, file="./RData/gex_prat.RData")
save(dat_prat, file="./RData/dat_prat.RData")

# Response abbreviations:
#
# From e.g. Hugo et al.
# PD: Progressive disease
# CR: Complete response
# PR: Partial response


# Create gene expression matrices

##
## Auslander et al.
##
exprs(gse_auslander[[1]])
# Botched
# 0 features, 23 + 14 samples? Raw counts in supplementary files

## First batch
colnames(sup_auslander) <- gsub("X|.bam", "", colnames(sup_auslander))
# Substrata
str_auslander1 <- pData(gse_auslander[[1]])$`treatment state:ch1`
gex_auslander1 <- sup_auslander[,gsub("X|.bam", "", gsub("-", ".", pData(gse_auslander[[1]])$title))]
# DESeq2
gex_auslander1 <- DESeq2::DESeqDataSetFromMatrix(countData=gex_auslander1, 
	colData = data.frame("str_auslander1" = str_auslander1), 
	design = ~ str_auslander1)
gex_auslander1 <- estimateSizeFactors(gex_auslander1)
gex_auslander1 <- counts(gex_auslander1, normalized=TRUE)
rownames(gex_auslander1) <- make.unique(sup_auslander[,1])

## Second batch

## Substrata
str_auslander2 <- pData(gse_auslander[[2]])$`treatment state:ch1`
# DESeq2
gex_auslander2 <- sup_auslander[,gsub(".bam", "", gsub("-", ".", pData(gse_auslander[[2]])$title))]
gex_auslander2 <- DESeq2::DESeqDataSetFromMatrix(countData=gex_auslander2, 
	colData = data.frame("str_auslander2" = pData(gse_auslander[[2]])$`treatment state:ch1`), 
	design = ~ str_auslander2)
gex_auslander2 <- estimateSizeFactors(gex_auslander2)
gex_auslander2 <- counts(gex_auslander2, normalized=TRUE)
rownames(gex_auslander2) <- make.unique(sup_auslander[,1])
# Combine the two batches
gex_auslander <- cbind(gex_auslander1, gex_auslander2)
# Spurious gene names such as '1-Dec', '1-Mar', '1-Sep', '10-Mar', '10-Sep', ... ??
#> sort(table(rownames(gex_auslander)))[1:10]
#
#  1-Dec   1-Mar 1-Mar.1   1-Sep  10-Mar  10-Sep  11-Mar  11-Sep  12-Sep  14-Sep 
#      1       1       1       1       1       1       1       1       1       1
# Gene names shouldn't be repeated
# Collapse over unique gene IDs, use mean
gex_auslander <- do.call("rbind", by(gex_auslander, INDICES=sup_auslander[,1], FUN=function(z) { apply(z, MARGIN=2, FUN=mean) }))
# Include only ones included in genes from biomaRt
gex_auslander <- gex_auslander[which(rownames(gex_auslander) %in% genes$hgnc_symbol),]
#> dim(gex_auslander)
#[1] 20245    37
# Much more reasonable row count - previously over 70k?!
# contains e.g. NONHSAG000315 etc (noncoding genes)
## Clinical data
cli_auslander <- rbind(pData(gse_auslander[[1]]), pData(gse_auslander[[2]]))
## TODO: Filter down to PRE samples only!
grep("PRE", cli_auslander$`treatment state:ch1`)
#> grep("PRE", cli_auslander$`treatment state:ch1`)
# [1]  1  4  5  7  8 10 17 18 24 26 28 30 32 35
## -> 14 samples usable
## Should also filter down to only antibody:ch1 being "anti-PD-1" or is "anti-CTLA-4" or combination representative as well? 
 
 
##
## Westin et al.
##
#str(gse_westin[[1]]@featureData)
#head(gse_westin[[1]]@featureData@data)


##
## Ascierto et al.
##

##
## Ascierto et al. 2
##


#setwd("..")








## Datasets downloaded from TIDE

# Lauss et al., Nat Commun 2017
# PFS, OS and 0/1 Response available
gex_lauss <- read.table(".\\TIDE\\Lauss2017_ACT_Melanoma_RNASeq\\ICB.Lauss2017_ACT_Melanoma.self_subtract", header=TRUE)
cli_lauss <- read.table(".\\TIDE\\Lauss2017_ACT_Melanoma_RNASeq\\ICB.Lauss2017_ACT_Melanoma.clinical", header=TRUE)
# Map to ensembl/hugo/etc
gex_lauss[,1] <- genes[match(gex_lauss[,1], genes$entrezgene_id),"hgnc_symbol"]
gex_lauss <- do.call("rbind", by(gex_lauss, INDICES=gex_lauss[,1], FUN=function(z){ apply(z[,-1], MARGIN=2, FUN=mean) }))
# Take intersection of patients present in both
pat_lauss <- intersect(rownames(cli_lauss), colnames(gex_lauss))
gex_lauss <- gex_lauss[,pat_lauss]
cli_lauss <- cli_lauss[pat_lauss,]
# Create dat_lauss
dat_lauss <- data.frame(
	patientID = rownames(cli_lauss),
	SEX = NA,
	AAGE = NA,
	CRFHIST = NA,
	TOBACUSE = NA,
	ECOGPS = as.integer(NA),
	PDL1 = as.integer(NA),
	TMB = as.numeric(NA),
	TCR_Shannon = as.numeric(NA),
	TCR_Richness = as.numeric(NA),
	TCR_Evenness = as.numeric(NA),
	BCR_Shannon = as.numeric(NA),
	BCR_Richness = as.numeric(NA),
	BCR_Evenness = as.numeric(NA),	
	PFS.time = as.integer(cli_lauss[,"PFS"]),
	PFS.event = as.integer(cli_lauss[,"PFS.Event"]),
	OS.time = as.integer(cli_lauss[,"OS"]),
	OS.event = as.integer(cli_lauss[,"OS.Event"]),
	Responder = as.integer(cli_lauss[,"Response"])
)
rownames(dat_lauss) <- dat_lauss$patientID
#> all(colnames(gex_lauss) == rownames(dat_lauss))
#[1] TRUE
save(gex_lauss, file="./RData/gex_lauss.RData")
save(dat_lauss, file="./RData/dat_lauss.RData")

# Kim et al., Nat Medicine 2018
# (Pembrolizumab)
# Response 0/1 available

gex_kim <- read.table(".\\TIDE\\Kim2018_PD1_Gastric_RNASeq\\ICB.Kim2018_Pembrolizumab_Gastric.self_subtract", header=TRUE)
cli_kim <- read.table(".\\TIDE\\Kim2018_PD1_Gastric_RNASeq\\ICB.Kim2018_Pembrolizumab_Gastric.clinical", header=TRUE, nrows=57)
# Map to ensembl/hugo/etc
gex_kim[,1] <- genes[match(gex_kim[,1], genes$entrezgene_id),"hgnc_symbol"]
gex_kim <- do.call("rbind", by(gex_kim, INDICES=gex_kim[,1], FUN=function(z){ apply(z[,-1], MARGIN=2, FUN=mean) }))
# cli_kim missing some responses in the end
# '-' symbols transform to '.' in column names
cli_kim$patient <- gsub('-', '.', cli_kim$patient)
rownames(cli_kim) <- cli_kim[,1]
cli_kim <- cli_kim[,-1,drop=FALSE]
# Take intersection of patients present in both
pat_kim <- intersect(rownames(cli_kim), colnames(gex_kim))
gex_kim <- gex_kim[,pat_kim]
cli_kim <- cli_kim[pat_kim,,drop=FALSE]
# Create dat_kim
dat_kim <- data.frame(
	patientID = rownames(cli_kim),
	SEX = NA,
	AAGE = NA,
	CRFHIST = NA,
	TOBACUSE = NA,
	ECOGPS = as.integer(NA),
	PDL1 = as.integer(NA),
	TMB = as.numeric(NA),
	TCR_Shannon = as.numeric(NA),
	TCR_Richness = as.numeric(NA),
	TCR_Evenness = as.numeric(NA),
	BCR_Shannon = as.numeric(NA),
	BCR_Richness = as.numeric(NA),
	BCR_Evenness = as.numeric(NA),	
	PFS.time = as.integer(NA),
	PFS.event = as.integer(NA),
	OS.time = as.integer(NA),
	OS.event = as.integer(NA),
	Responder = as.integer(cli_kim[,"Response"])
)
rownames(dat_kim) <- dat_kim$patientID
#> all(colnames(gex_kim) == rownames(dat_kim))
#[1] TRUE
save(gex_kim, file="./RData/gex_kim.RData")
save(dat_kim, file="./RData/dat_kim.RData")


# Chen et al., Cancer Discov 2016
# Response 0/1 available
gex_chen <- read.table(".\\TIDE\\Chen2016_PD1_Melanoma_Nanostring_Ipi.Prog\\ICB.Chen2016_PD1_Melanoma_Ipi.Prog.self_subtract", header=TRUE)
cli_chen <- read.table(".\\TIDE\\Chen2016_PD1_Melanoma_Nanostring_Ipi.Prog\\ICB.Chen2016_PD1_Melanoma_Ipi.Prog.clinical", header=TRUE)
# P6 duplicate omitted
gex_chen <- gex_chen[,-6]
cli_chen <- cli_chen[-5,]
# Map to ensembl/hugo/etc
gex_chen[,1] <- genes[match(gex_chen[,1], genes$entrezgene_id),"hgnc_symbol"]
gex_chen <- do.call("rbind", by(gex_chen, INDICES=gex_chen[,1], FUN=function(z){ apply(z[,-1], MARGIN=2, FUN=mean) }))
rownames(cli_chen) <- cli_chen[,1]
cli_chen <- cli_chen[,-1,drop=FALSE]
# Take intersection of patients present in both
pat_chen <- intersect(rownames(cli_chen), colnames(gex_chen))
gex_chen <- gex_chen[,pat_chen]
cli_chen <- cli_chen[pat_chen,,drop=FALSE]
# Create dat_chen
dat_chen <- data.frame(
	patientID = rownames(cli_chen),
	SEX = NA,
	AAGE = NA,
	CRFHIST = NA,
	TOBACUSE = NA,
	ECOGPS = as.integer(NA),
	PDL1 = as.integer(NA),
	TMB = as.numeric(NA),
	TCR_Shannon = as.numeric(NA),
	TCR_Richness = as.numeric(NA),
	TCR_Evenness = as.numeric(NA),
	BCR_Shannon = as.numeric(NA),
	BCR_Richness = as.numeric(NA),
	BCR_Evenness = as.numeric(NA),	
	PFS.time = as.integer(NA),
	PFS.event = as.integer(NA),
	OS.time = as.integer(NA),
	OS.event = as.integer(NA),
	Responder = as.integer(cli_chen[,"Response"])
)
rownames(dat_chen) <- dat_chen$patientID
#> all(colnames(gex_chen) == rownames(dat_chen))
#[1] TRUE
save(gex_chen, file="./RData/gex_chen.RData")
save(dat_chen, file="./RData/dat_chen.RData")





###
#
# Cherry-pick key genes from the clinic or other sources (even if mutation status is used in clinic rather than gex)
# Note aliases! Like PDL1 <-> CD274
#
###

keyGenes <- c(
	"CD274", "PDL1", # PD-L1
	"PDCD1", "CD279"
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
	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7542663/
	# Antigen processing and presentation machinery
	# Score itself: "The APMS (sum of the log2 z-scores for each gene)"
	## Omit and include only the antigen presentation machinery score
	#"B2M",
	#"CALR",
	#"NLRC5",
	#"PSMB9",
	#"PSME1",
	#"PSME3",
	#"RFX5",
	#"HSP90AB1",
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
	"PDK-1"
	"AKT",
	# Mikko NGS syöpäpaneeli 1 lisäyksiä
	"PIK3CA",
	"KIT",
	"NRAS",
	"PDGFRA",
	# Other interesting? Cytokines, chemokines
	"CLCL10",
	"CSCL11",
	
	
)

###
#
# Lung IO specific "scores", and filtering out genes not present
#
###

# Interferon gamma response score
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5531419/
# Might be too melanoma specific, albeit covering a multitude of other cancers too
# Calculation
# "After performance of quantile normalization, a log10 transformation was applied, 
# and signature scores were calculated by averaging of the included genes for the IFN-gamma (6-gene) and expanded immune (18-gene) signatures."
# ... although ...
# "Logistic regression modeling was used to conduct the hypothesis testing associated with best overall response (BOR), and a Cox model was used for testing of PFS and OS."
IFNGscore1 <- function(gex, columns=TRUE){
	genes <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")
	# > all(c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG") %in% rownames(gex_synthetic_refseq105_genes_tpm))
	# [1] TRUE
	
}
# "The final 18-gene profile was derived through a cross-validated penalized regression modeling strategy in a large cohort of pembrolizumab-treated patients across 9 different tumor types."
IFNGscore2 <- function(gex, columns=TRUE){
	genes <- c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB")
	#> all(c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB") %in% rownames(gex_synthetic_refseq105_genes_tpm))
	#[1] TRUE
	
}

# Antigen processing and presentation machinery
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7542663/
# Score itself: "The APMS (sum of the log2 z-scores for each gene)"
APMscore <- function(gex, columns=TRUE){
	genes <- c("B2M", "CALR", "NLRC5", "PSMB9", "PSME1", "PSME3", "RFX5", "HSP90AB1")
	#> all(c("B2M", "CALR", "NLRC5", "PSMB9", "PSME1", "PSME3", "RFX5", "HSP90AB1") %in% rownames(gex_synthetic_refseq105_genes_tpm))
	#[1] TRUE
	
}

# Tumor inflammatory signature (160 genes)
# https://ascopubs.org/doi/abs/10.1200/JCO.2020.38.5_suppl.47
# "TIS, algorithmically defined as the mean mRNA expression of the 160 genes was developed with each tumor assigned into a weak, moderate or strong inflammation group"
# ... using Damotte et al. instead
# Inflammation was observed differently in varying PD-L1 expressed tumors
# "Strongly inflamed tumors presented with improved ORR to ICI in NSCLC"
# Possibly differing criteria for "response" from what DREAM uses: "clinical benefit was defined as complete or partial RECIST response while stable and progressive disease were defined as lack of clinical benefit."
TISscore <- function(gex, columns=TRUE){
	# Taken from Fig 1 panel d in the publ.
	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6829827/#MOESM2
	genes <- c("CD276", "HLA-DQA1", "CD274", "IDO1", "HLA-DRB1", "HLA-E", "CMKLR1", "PDCD1LG2", "PSMB10", "LAG3", "CXCL9", "STAT1", "CD8A", "CCL5", "NKG7", "TIGIT", "CD27", "CXCR6")
	#> all(c("CD276", "HLA-DQA1", "CD274", "IDO1", "HLA-DRB1", "HLA-E", "CMKLR1", "PDCD1LG2", "PSMB10", "LAG3", "CXCL9", "STAT1", "CD8A", "CCL5", "NKG7", "TIGIT", "CD27", "CXCR6") %in% rownames(gex_synthetic_refseq105_genes_tpm))
	#[1] TRUE
	
}

###
#
# Derive simple sample GSVA results to refine the genes to pathway-level information
#
###

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
library(GSVA)

# Read chosen set of GMTs; namely: hallmarks, oncogenic, and immunologic pathways
# Downloaded from e.g. https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C7
gmt_h <- GSEABase::getGmt(".\\MSigDB\\h.all.v7.2.symbols.gmt")
gmt_c6 <- GSEABase::getGmt(".\\MSigDB\\c6.all.v7.2.symbols.gmt")
gmt_c7 <- GSEABase::getGmt(".\\MSigDB\\c7.all.v7.2.symbols.gmt")

# Synthetic
gmt_synthetic <- rbind(
	GSVA::gsva(as.matrix(gex_synthetic_refseq105_genes_tpm), gmt_h),	# Hallmarks
	GSVA::gsva(as.matrix(gex_synthetic_refseq105_genes_tpm), gmt_c6)#,	# Oncogenic
	#GSVA::gsva(as.matrix(gex_synthetic_refseq105_genes_tpm), gmt_c7)	# Immunology
)
# TCGA
#gmt_luad <- rbind(
#	GSVA::gsva(gex_luad, gmt_h),	# Hallmarks
#	GSVA::gsva(gex_luad, gmt_c6),	# Oncogenic
#	GSVA::gsva(gex_luad, gmt_c7)	# Immunology
#)
#gmt_lusc <- rbind(
#	GSVA::gsva(gex_lusc, gmt_h),	# Hallmarks
#	GSVA::gsva(gex_lusc, gmt_c6),	# Oncogenic
#	GSVA::gsva(gex_lusc, gmt_c7)	# Immunology
#)
## Combined (gex transposed)
gmt_tcga <- rbind(
	GSVA::gsva(gex_tcga, gmt_h),	# Hallmarks
	GSVA::gsva(gex_tcga, gmt_c6)#,	# Oncogenic
	#GSVA::gsva(gex_tcga, gmt_c7)	# Immunology
)

# TIDEs
gmt_lauss <- rbind(
	GSVA::gsva(gex_lauss, gmt_h),	# Hallmarks
	GSVA::gsva(gex_lauss, gmt_c6)#,	# Oncogenic
	#GSVA::gsva(gex_lauss, gmt_c7)	# Immunology
)
gmt_kim <- rbind(
	GSVA::gsva(gex_kim, gmt_h),	# Hallmarks
	GSVA::gsva(gex_kim, gmt_c6)#,	# Oncogenic
	#GSVA::gsva(gex_kim, gmt_c7)	# Immunology
)
gmt_chen <- rbind(
	GSVA::gsva(gex_chen, gmt_h),	# Hallmarks
	GSVA::gsva(gex_chen, gmt_c6)#,	# Oncogenic
	#GSVA::gsva(gex_chen, gmt_c7)	# Immunology
)



###
#
# Immune cell deconvolution, using immunedeconv-package
#
###

#install.packages("remotes")
#remotes::install_github("icbi-lab/immunedeconv")
library(immunedeconv)


## Synthetic
## Updated: Separate idc xcell and mcp
# Renaming function that transforms the data table to an understandable numeric matrix
rename_idc <- function(x, method="xcell"){
	tmp <- gsub(" ", "_", gsub("\\+|\\-", "", c(x[,1])[[1]]))
	x <- as.matrix(x[,-1])
	rownames(x) <- paste(method, "_", tmp, sep="")
	x
}

idc_xce_synthetic <- rename_idc(immunedeconv::deconvolute(gex_synthetic_refseq105_genes_tpm, method="xcell"), method="xce")
idc_mcp_synthetic <- rename_idc(immunedeconv::deconvolute(gex_synthetic_refseq105_genes_tpm, method="mcp_counter"), method="mcp")

#idc_synthetic <- rbind(
#	immunedeconv::deconvolute(gex_synthetic_refseq105_genes_tpm, method="xcell"),
#	immunedeconv::deconvolute(gex_synthetic_refseq105_genes_tpm, method="mcp_counter")
#)

## TIDEs
idc_xce_lauss <- rename_idc(immunedeconv::deconvolute(gex_lauss, method="xcell"), method="xce")
idc_mcp_lauss <- rename_idc(immunedeconv::deconvolute(gex_lauss, method="mcp_counter"), method="mcp")
# Lauss
#idc_lauss <- rbind(
#	immunedeconv::deconvolute(gex_lauss, method="xcell"),
#	immunedeconv::deconvolute(gex_lauss, method="mcp_counter")
#)
#tmp <- idc_lauss[,1]
#idc_lauss <- as.matrix(idc_lauss[,-1])
#rownames(idc_lauss) <- c(tmp)[[1]]
# Kim
idc_xce_kim <- rename_idc(immunedeconv::deconvolute(gex_kim, method="xcell"), method="xce")
idc_mcp_kim <- rename_idc(immunedeconv::deconvolute(gex_kim, method="mcp_counter"), method="mcp")
#idc_kim <- rbind(
#	immunedeconv::deconvolute(gex_kim, method="xcell"),
#	immunedeconv::deconvolute(gex_kim, method="mcp_counter")
#)
#tmp <- idc_kim[,1]
#idc_kim <- as.matrix(idc_kim[,-1])
#rownames(idc_kim) <- c(tmp)[[1]]
# Chen
# - unable to run, key genes missing -

#rm(tmp)

###
#
# DeMixT
#
###

# https://github.com/wwylab/DeMixT
# devtools::install_github("wwylab/DeMixT")
# or rather...
# https://bioconductor.org/packages/release/bioc/html/DeMixT.html
library(DeMixT)


###
#
# ESTIMATE
#
###

# https://www.nature.com/articles/ncomms3612
# https://bioinformatics.mdanderson.org/public-software/estimate/
# http://r-forge.r-project.org/R/?group_id=2237
# install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)
#library(estimate)

# Nevermind, requires .gct input

# Preliminary testing whether any of the gmt / gene sets are found by lasso
if(FALSE){
	library(glmnet)
	library(survival)
	# LASSO testing

	### Lauss

	set.seed(1)
	# PFS
	lasso_sub1_lauss <- glmnet(x=t(gmt_lauss), y=survival::Surv(time=cli_lauss[,"PFS"], event=cli_lauss[,"PFS.Event"]), family="cox")
	lassocv_sub1_lauss <- cv.glmnet(x=t(gmt_lauss), y=survival::Surv(time=cli_lauss[,"PFS"], event=cli_lauss[,"PFS.Event"]), family="cox", nfolds=5)
	plot(lassocv_sub1_lauss)
	# No non-zeros...
	# OS
	lasso_sub2_lauss <- glmnet(x=t(gmt_lauss), y=survival::Surv(time=cli_lauss[,"OS"], event=cli_lauss[,"OS.Event"]), family="cox")
	lassocv_sub2_lauss <- cv.glmnet(x=t(gmt_lauss), y=survival::Surv(time=cli_lauss[,"OS"], event=cli_lauss[,"OS.Event"]), family="cox", nfolds=5)
	plot(lassocv_sub2_lauss)
	# No non-zeros...
	# Response
	lasso_sub3_lauss <- glmnet(x=t(gmt_lauss), y=cli_lauss[,"Response"], family="binomial")
	lassocv_sub3_lauss <- cv.glmnet(x=t(gmt_lauss), y=cli_lauss[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_lauss)

	### Kim

	set.seed(1)
	# Response
	lasso_sub3_kim <- glmnet(x=t(gmt_kim), y=cli_kim[,"Response"], family="binomial")
	lassocv_sub3_kim <- cv.glmnet(x=t(gmt_kim), y=cli_kim[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_kim)
	#> rownames(gmt_kim)[predict(lasso_sub3_kim, s=lassocv_sub3_kim$lambda.min, type="nonzero")[,1]]
	# [1] "IL15_UP.V1_DN"                                                            
	# [2] "GOLDRATH_EFF_VS_MEMORY_CD8_TCELL_DN"                                      
	# [3] "GSE17721_4H_VS_24H_POLYIC_BMDC_DN"                                        
	# [4] "GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_UP"                                    
	# [5] "GSE29615_CTRL_VS_DAY7_LAIV_FLU_VACCINE_PBMC_DN"                           
	# [6] "GSE36476_YOUNG_VS_OLD_DONOR_MEMORY_CD4_TCELL_40H_TSST_ACT_UP"             
	# [7] "GSE8515_IL1_VS_IL6_4H_STIM_MAC_DN"                                        
	# [8] "GSE22601_DOUBLE_NEGATIVE_VS_IMMATURE_CD4_SP_THYMOCYTE_UP"                 
	# [9] "GSE2585_CD80_HIGH_VS_LOW_MTEC_UP"                                         
	#[10] "GSE3920_IFNA_VS_IFNG_TREATED_FIBROBLAST_DN"                               
	#[11] "GSE6875_WT_VS_FOXP3_KO_TREG_DN"                                           
	#[12] "GSE7831_CPG_VS_INFLUENZA_STIM_PDC_4H_UP"                                  
	#[13] "GSE12484_HEALTHY_VS_PERIDONTITIS_NEUTROPHILS_DN"                          
	#[14] "GSE7459_UNTREATED_VS_IL6_TREATED_ACT_CD4_TCELL_DN"                        
	#[15] "GSE19374_UNINF_VS_LISTERIA_INFECTED_MACROPHAGE_UP"                        
	#[16] "GSE22432_CONVENTIONAL_CDC_VS_PLASMACYTOID_PDC_DN"                         
	#[17] "GSE34156_NOD2_LIGAND_VS_NOD2_AND_TLR1_TLR2_LIGAND_24H_TREATED_MONOCYTE_UP"
	#[18] "GSE11961_MARGINAL_ZONE_BCELL_VS_MEMORY_BCELL_DAY7_DN"                     
	#[19] "GSE11961_FOLLICULAR_BCELL_VS_MARGINAL_ZONE_BCELL_DN"                      
	#[20] "GSE42724_MEMORY_VS_B1_BCELL_UP"                                           
	#[21] "GSE46606_IRF4_KO_VS_WT_CD40L_IL2_IL5_1DAY_STIMULATED_BCELL_UP"
	#> rownames(gmt_kim)[predict(lasso_sub3_kim, s=lassocv_sub3_kim$lambda.1se, type="nonzero")[,1]]
	#[1] "IL21_UP.V1_DN"                                                "GSE10325_MYELOID_VS_LUPUS_MYELOID_UP"                        
	#[3] "GSE36476_YOUNG_VS_OLD_DONOR_MEMORY_CD4_TCELL_40H_TSST_ACT_UP" "GSE2585_CD80_HIGH_VS_LOW_MTEC_UP"                            
	#[5] "GSE41176_UNSTIM_VS_ANTI_IGM_STIM_TAK1_KO_BCELL_3H_DN"

	### Chen

	set.seed(1)
	# Response
	lasso_sub3_chen <- glmnet(x=t(gmt_chen), y=cli_chen[,"Response"], family="binomial")
	## Of note:
	#Warning message:
	#In lognet(x, is.sparse, ix, jx, y, weights, offset, alpha, nobs,  :
	#  one multinomial or binomial class has fewer than 8  observations; dangerous ground
	lassocv_sub3_chen <- cv.glmnet(x=t(gmt_chen), y=cli_chen[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_chen)
	#rownames(gmt_chen)[predict(lasso_sub3_chen, s=lassocv_sub3_chen$lambda.min, type="nonzero")[,1]]
	#> rownames(gmt_chen)[predict(lasso_sub3_chen, s=lassocv_sub3_chen$lambda.min, type="nonzero")[,1]]
	#[1] "GSE22601_IMMATURE_CD4_SINGLE_POSITIVE_VS_DOUBLE_POSITIVE_THYMOCYTE_UP"
	#[2] "GSE37301_HEMATOPOIETIC_STEM_CELL_VS_CD4_TCELL_DN"


}






# Preliminary testing of whether glmnet grabs anything
if(FALSE){

	### Lauss

	set.seed(1)
	# PFS
	lasso_sub1_lauss <- glmnet(x=t(idc_lauss), y=survival::Surv(time=cli_lauss[,"PFS"], event=cli_lauss[,"PFS.Event"]), family="cox")
	lassocv_sub1_lauss <- cv.glmnet(x=t(idc_lauss), y=survival::Surv(time=cli_lauss[,"PFS"], event=cli_lauss[,"PFS.Event"]), family="cox", nfolds=5)
	plot(lassocv_sub1_lauss)
	#> rownames(idc_lauss)[predict(lasso_sub1_lauss, s=lassocv_sub1_lauss$lambda.min, type="nonzero")[,1]]
	# [1] "Myeloid dendritic cell activated" "T cell CD4+ (non-regulatory)"     "T cell CD4+ central memory"      
	# [4] "T cell CD4+ effector memory"      "T cell CD8+ effector memory"      "Common lymphoid progenitor"      
	# [7] "Common myeloid progenitor"        "Cancer associated fibroblast"     "Hematopoietic stem cell"         
	#[10] "Macrophage M2"                    "Mast cell"                        "Neutrophil"                      
	#[13] "NK cell"                          "T cell NK"                        "B cell plasma"                   
	#[16] "T cell CD4+ Th1"                  "T cell CD4+ Th2"                  "T cell regulatory (Tregs)"       
	#[19] "stroma score"                     "cytotoxicity score"               "NK cell"                         
	#[22] "Myeloid dendritic cell"           "Neutrophil"
	# --> All cell types except couple...
	# OS
	lasso_sub2_lauss <- glmnet(x=t(idc_lauss), y=survival::Surv(time=cli_lauss[,"OS"], event=cli_lauss[,"OS.Event"]), family="cox")
	lassocv_sub2_lauss <- cv.glmnet(x=t(idc_lauss), y=survival::Surv(time=cli_lauss[,"OS"], event=cli_lauss[,"OS.Event"]), family="cox", nfolds=5)
	plot(lassocv_sub2_lauss)
	#> rownames(idc_lauss)[predict(lasso_sub2_lauss, s=lassocv_sub2_lauss$lambda.min, type="nonzero")[,1]]
	# [1] "Myeloid dendritic cell activated" "T cell CD4+ (non-regulatory)"     "T cell CD4+ effector memory"     
	# [4] "T cell CD8+ naive"                "T cell CD8+ effector memory"      "Class-switched memory B cell"    
	# [7] "Common lymphoid progenitor"       "Cancer associated fibroblast"     "Hematopoietic stem cell"         
	#[10] "Monocyte"                         "Neutrophil"                       "NK cell"                         
	#[13] "T cell NK"                        "T cell gamma delta"               "T cell CD4+ Th1"                 
	#[16] "T cell regulatory (Tregs)"        "stroma score"                     "T cell"                          
	#[19] "NK cell"                          "Myeloid dendritic cell"           "Neutrophil"                      
	#[22] "Endothelial cell"                 "Cancer associated fibroblast"    
	#> rownames(idc_lauss)[predict(lasso_sub2_lauss, s=lassocv_sub2_lauss$lambda.1se, type="nonzero")[,1]]
	# [1] "Myeloid dendritic cell activated" "T cell CD4+ effector memory"      "T cell CD8+ naive"               
	# [4] "T cell CD8+ effector memory"      "Common lymphoid progenitor"       "Granulocyte-monocyte progenitor" 
	# [7] "Hematopoietic stem cell"          "Macrophage M2"                    "T cell CD4+ Th1"                 
	#[10] "T cell regulatory (Tregs)"        "cytotoxicity score"               "NK cell"
	# --> Highly informative for both PFS/OS
	# Response
	lasso_sub3_lauss <- glmnet(x=t(idc_lauss), y=cli_lauss[,"Response"], family="binomial")
	lassocv_sub3_lauss <- cv.glmnet(x=t(idc_lauss), y=cli_lauss[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_lauss)
	# No non-zeros...

	### Kim

	set.seed(1)
	# Response
	lasso_sub3_kim <- glmnet(x=t(idc_kim), y=cli_kim[,"Response"], family="binomial")
	lassocv_sub3_kim <- cv.glmnet(x=t(idc_kim), y=cli_kim[,"Response"], family="binomial", nfolds=5)
	plot(lassocv_sub3_kim)
	#> rownames(idc_kim)[predict(lasso_sub3_kim, s=lassocv_sub3_kim$lambda.min, type="nonzero")[,1]]
	#[1] "T cell CD4+ central memory"   "Macrophage"                   "T cell CD4+ Th1"              "Endothelial cell"            
	#[5] "Cancer associated fibroblast"
	#> rownames(idc_kim)[predict(lasso_sub3_kim, s=lassocv_sub3_kim$lambda.1se, type="nonzero")[,1]]
	#[1] "Macrophage"                   "T cell CD4+ Th1"              "Endothelial cell"             "Cancer associated fibroblast"

}


# Only locally, ignored by .gitignore
#save.image("temp.RData")

