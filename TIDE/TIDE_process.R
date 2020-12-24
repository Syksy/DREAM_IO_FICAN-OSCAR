##
#
# Processing datasets from TIDE database
#
##

setwd("C:\\Users\\Daniel\\DREAM_2020_IO\\TIDE")

# Lauss et al., Nat Commun 2017
# PFS, OS and 0/1 Response available
gex_lauss <- read.table(".\\Lauss2017_ACT_Melanoma_RNASeq\\ICB.Lauss2017_ACT_Melanoma.self_subtract", header=TRUE)
cli_lauss <- read.table(".\\Lauss2017_ACT_Melanoma_RNASeq\\ICB.Lauss2017_ACT_Melanoma.clinical", header=TRUE)

# Kim et al., Nat Medicine 2018
# (Pembrolizumab)
# Response 0/1 available

gex_kim <- read.table(".\\Kim2018_PD1_Gastric_RNASeq\\ICB.Kim2018_Pembrolizumab_Gastric.self_subtract", header=TRUE)
cli_kim <- read.table(".\\Kim2018_PD1_Gastric_RNASeq\\ICB.Kim2018_Pembrolizumab_Gastric.clinical", header=TRUE, nrows=57)
# cli_kim missing some responses in the end

# Chen et al., Cancer Discov 2016

gex_chen <- read.table(".\\Chen2016_PD1_Melanoma_Nanostring_Ipi.Prog\\ICB.Chen2016_PD1_Melanoma_Ipi.Prog.self_subtract", header=TRUE)
cli_chen <- read.table(".\\Chen2016_PD1_Melanoma_Nanostring_Ipi.Prog\\ICB.Chen2016_PD1_Melanoma_Ipi.Prog.clinical", header=TRUE)

