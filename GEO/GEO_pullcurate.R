###
#
# Pull & curate suitable datasets from GEO
#
###

# A lot of the may be utilizable best from: https://github.com/JasonACarter/IMPRES_Correspondence

# Studies suggested by the DREAM organizers available in GEO:

# GSE115821 - Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115821
# Auslander et al.
# Cancer(s): Melanoma/Neuroblastoma (?)
# Treatment(s): N=31 treated with anti-PD-1, N=10 treated with anti-CTLA-4
# N(s): 37, for some patients multiple samples at different time points or different biopsies at same time point
# Response types: Binary (?)
# Year: 2018
# Notes: 
# - Identified 15 pairwise transcriptomics relations between immune checkpoint genes; AUC=0.83; validated on 6 published existing cohorts at the time (2018) (n total 297)
# - Criticized at https://www.nature.com/articles/s41591-019-0671-4

# GSE121810 - Neoadjuvant anti-PD-1 immunotherapy promotes a survival benefit with intratumoral and systemic immune responses in recurrent glioblastoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121810
# Cloughesy et al.
# Cancer(s): Glioblastoma (recurrent, surgically resectable)
# Treatment(s): Neoadjuvant Pembrolizumab vs. neoadjuvant post-surgical PD-1 blockade alone
# N(s): 29 (?)
# Response types: Survival
# Year: 2019
# Notes: 
# - Recurrent disease
# - "This dataset contains the transcriptomes of recurrent glioblastoma with either neoadjuvant (1 dose) or adjuvant pembrolizumab treatment"
# - "This trial was registered with ClinicalTrials.gov under the identifier NCT02852655 (https://clinicaltrials.gov/ct2/show/NCT02852655)."

# GSE78220 - mRNA expressions in pre-treatment melanomas undergoing anti-PD-1 checkpoint inhibition therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220
# Hugo et al.
# Cancer(s): Melanoma
# Treatment(s): All with anti-PD-1, some responded others didn't
# N(s): 28 (?)
# Response types: Response vs. no response (binary)
# Year: 2016
# Notes:
# - "Mutations in cell adhesion genes and the DNA repair gene BRCA2 were enriched in responding tumors, and a high mutational load associated with improved survival."
# - "Innately resistant tumors displayed frequent transcriptomic up-expression of genes that enriched for mesenchymal transition, cell adhesion, ECM organization, wound-healing and angiogenesis."
# - "he transcriptomes of innate resistance also enriched for signatures indicating up-regulation of these processes. Notably, MAPK-targeted therapy (MAPKi) induced similar signatures in melanoma, suggesting that a form of MAPKi resistance mediates cross-resistance to anti-PD-1 therapy. Co-enrichment of IPRIM (Innate anti-PD-1 Resistance Induced by MAPKi) signatures defined a transcriptomic subset across advanced cancers, suggesting that attenuating processes underlying these signatures may augment anti-PD1 responses."

# GSE52562 - Gene expression profiling of tumor biopsies before and after pidilizumab therapy in patients with relapsed follicular lymphoma grade 1 or grade 2.
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52562
# Westin et al.
# Cancer(s): Follicular lymphoma (grade 1 or 2)
# Treatment(s): Pidilizumab
# N(s): N=8 pairs for pre vs. post-treatment, N=10 additional pre-treatment samples
# Response types: Progression free survival
# Year: 2014
# Notes: 
# - Comparison before and after pidilizumab

# GSE67501 - Expression data in human renal cell carcinoma samples from patients who did or did not respond to anti-PD-1 (Nivolumab) immunotherapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67501
# Ascierto et al. 2016
# Cancer(s): Renal cell carcinoma
# Treatment(s): Anti-PD-1 (Nivolumab)
# N(s): N=4 responders, N=7 non-responders
# Response types: Binary (Response vs. no response)
# Year: 2016
# Notes:
# - Focus on positive tumor expression of PD-L1
# - Microarray expression data

# GSE79691 - Transcriptional mechanisms of resistance to anti-PD-1 therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79691
# Ascierto et al. 2017
# Cancer(s): Melanoma, skin metastases
# Treatment(s): Nivolumab
# N(s): response (N=6), no response (N=4)
# Response types: Binary
# Year: 2017
# Notes: 
# - All samples from a single patient
# - Patient included in Ascierto et al. 2016?

## Other identified relevant studies

# GSE91061 - Molecular portraits of tumor mutational and micro-environmental sculpting by immune checkpoint blockade therapy
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061
# Riaz et al. 
# Cancer(s): Advanced melanoma
# Treatment(s): Ipilimumab, Nivolumab
# N(s): N=58 on-treatment and N=51 pre-treatment, 64 patients
# Response types:
# Year:
# Notes
# - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5685550/
# - Ipilimumab-naive
# - Abstract: "The mechanisms by which immune checkpoint blockade modulates tumor evolution during therapy are unclear. We assessed genomic changes in tumors from 68 patients with advanced melanoma, who progressed on ipilimumab or were ipilimumab-naive, before and after nivolumab initiation (CA209-038 study). Tumors were analyzed by whole-exome, transcriptome, and/or T cell receptor (TCR) sequencing. In responding patients, mutation and neoantigen load were reduced from baseline, and analysis of intratumoral heterogeneity during therapy demonstrated differential clonal evolution within tumors and putative selection against neoantigenic mutations on-therapy. Transcriptome analyses before and during nivolumab therapy revealed increases in distinct immune cell subsets, activation of specific transcriptional networks, and upregulation of immune checkpoint genes that were more pronounced in patients with response. Temporal changes in intratumoral TCR repertoire revealed expansion of T cell clones in the setting of neoantigen loss. Comprehensive genomic profiling data in this study provide insight into nivolumab's mechanism of action."

# GSE93157 - Programmed death 1 receptor blockade and immune-related gene expression profiling in non-small cell lung carcinoma, head and neck squamous cell carcinoma and melanoma
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93157
# Prat et al.
# Cancer(s): Non-small cell lung carcinoma, heda and neck squamous cell carcinoma, melanoma
# Treatment(s): Anti-PD1 (pembrolizumab or nivolumab)
# N(s): 
# Response types:
# Year:
# Notes
# - Has non-small cell lung carcinoma
# - https://cancerres.aacrjournals.org/content/77/13/3540.long
# - Focus heavily on immune cell decomposition

## Other

# Template
# Name
# URL
# Cite
# Cancer(s):
# Treatment(s):
# N(s):
# Response types:
# Year:
# Notes
# - 1
# - 2
# - 3

