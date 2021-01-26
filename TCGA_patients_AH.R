# TCGA-datan kuratointia : LUSC ja LUAD chemopotilaat

lusc <- curatedTCGAData::curatedTCGAData(diseaseCode = "LUSC",assays=c("RNASeq2*","mRNA*"),dry.run = FALSE)
luad <- curatedTCGAData::curatedTCGAData(diseaseCode = "LUAD",assays=c("RNASeq2*","mRNA*"),dry.run = FALSE)

colnames(assay(luad.RNA@ExperimentList@listData$`LUAD_RNASeqGene-20160128`))

col.lusc <- as.matrix(lusc@colData)
col.luad <- as.matrix(luad@colData)

# Haetaan hoitotiedot
treatlusc <- col.lusc[,c(which(grepl("therapy",colnames(col.lusc))),which(grepl("treatment",colnames(col.lusc))))]
treatluad <- col.luad[,c(which(grepl("therapy",colnames(col.luad))),which(grepl("treatment",colnames(col.luad))))]

# Mitä eri hoitoja on tarjolla
treatypes.lusc <-unique(c(treatlusc[,"patient.drugs.drug.therapy_types.therapy_type"],treatlusc[,"patient.drugs.drug.2.therapy_types.therapy_type"],
         treatlusc[,"patient.drugs.drug.3.therapy_types.therapy_type"],treatlusc[,"patient.drugs.drug.4.therapy_types.therapy_type"],
         treatlusc[,"patient.drugs.drug.5.therapy_types.therapy_type"],treatlusc[,"patient.drugs.drug.6.therapy_types.therapy_type"],
         treatlusc[,"patient.drugs.drug.7.therapy_types.therapy_type"],treatlusc[,"patient.drugs.drug.8.therapy_types.therapy_type"],
         treatlusc[,"patient.drugs.drug.9.therapy_types.therapy_type"],treatlusc[,"patient.drugs.drug.10.therapy_types.therapy_type"],
         treatlusc[,"patient.drugs.drug.11.therapy_types.therapy_type"],treatlusc[,"patient.drugs.drug.12.therapy_types.therapy_type"]))
treatypes.luad <- unique(c(treatluad[,"patient.drugs.drug.therapy_types.therapy_type"],treatluad[,"patient.drugs.drug.2.therapy_types.therapy_type"],
         treatluad[,"patient.drugs.drug.3.therapy_types.therapy_type"],treatluad[,"patient.drugs.drug.4.therapy_types.therapy_type"],
         treatluad[,"patient.drugs.drug.5.therapy_types.therapy_type"],treatluad[,"patient.drugs.drug.6.therapy_types.therapy_type"],
         treatluad[,"patient.drugs.drug.7.therapy_types.therapy_type"],treatluad[,"patient.drugs.drug.8.therapy_types.therapy_type"],
         treatluad[,"patient.drugs.drug.9.therapy_types.therapy_type"],treatluad[,"patient.drugs.drug.10.therapy_types.therapy_type"]))

# Sädetystä saavat  
radiation.lusc <-sort(unique(c(which(col.lusc[,"radiation_therapy"]=="yes"),which(col.lusc[,"patient.follow_ups.follow_up.2.additional_radiation_therapy"]=="yes"),
                          which(col.lusc[,"patient.follow_ups.follow_up.2.radiation_therapy"]=="yes"),which(col.lusc[,"patient.follow_ups.follow_up.3.radiation_therapy"]=="yes"),
                          which(col.lusc[,"patient.follow_ups.follow_up.3.additional_radiation_therapy"]=="yes"),which(col.lusc[,"patient.follow_ups.follow_up.4.additional_radiation_therapy"]=="yes"),
                          which(col.lusc[,"patient.follow_ups.follow_up.3.radiation_therapy"]=="yes"),which(col.lusc[,"patient.follow_ups.follow_up.radiation_therapy"]=="yes"),
                          which(col.lusc[,"patient.follow_ups.follow_up.additional_radiation_therapy"]=="yes"),which(col.lusc[,"patient.new_tumor_events.new_tumor_event.additional_radiation_therapy"]=="yes"),
                          which(col.lusc[,"patient.radiation_therapy"]=="yes"),which(col.lusc[,"patient.radiations.radiation.2.radiation_treatment_ongoing"]=="yes"),
                          which(col.lusc[,"patient.radiations.radiation.radiation_treatment_ongoing"]=="yes"))))
radiation.luad <-sort(unique(c(which(col.luad[,"radiation_therapy"]=="yes"),which(col.luad[,"patient.follow_ups.follow_up.2.additional_radiation_therapy"]=="yes"),
                               which(col.luad[,"patient.follow_ups.follow_up.2.radiation_therapy"]=="yes"),which(col.luad[,"patient.follow_ups.follow_up.3.radiation_therapy"]=="yes"),
                               which(col.luad[,"patient.follow_ups.follow_up.3.additional_radiation_therapy"]=="yes"),which(col.luad[,"patient.follow_ups.follow_up.4.additional_radiation_therapy"]=="yes"),
                               which(col.luad[,"patient.follow_ups.follow_up.3.radiation_therapy"]=="yes"),which(col.luad[,"patient.follow_ups.follow_up.radiation_therapy"]=="yes"),
                               which(col.luad[,"patient.follow_ups.follow_up.additional_radiation_therapy"]=="yes"),which(col.luad[,"patient.new_tumor_events.new_tumor_event.additional_radiation_therapy"]=="yes"),
                               which(col.luad[,"patient.radiation_therapy"]=="yes"),which(col.luad[,"patient.radiations.radiation.2.radiation_treatment_ongoing"]=="yes"),
                               which(col.luad[,"patient.radiations.radiation.radiation_treatment_ongoing"]=="yes"))))

chemot.lusc <- sort(unique(c(which(treatlusc[,"patient.drugs.drug.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.11.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatlusc[,"patient.drugs.drug.12.therapy_types.therapy_type"]=="chemotherapy"))))

chemot.luad <- sort(unique(c(which(treatluad[,"patient.drugs.drug.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatluad[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatluad[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatluad[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatluad[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatluad[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatluad[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatluad[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatluad[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="chemotherapy"),
                        which(treatluad[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="chemotherapy"))))

target.lusc <- sort(unique(c(which(treatlusc[,"patient.drugs.drug.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.11.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.drugs.drug.12.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.follow_ups.follow_up.targeted_molecular_therapy"]=="yes"),
                             which(treatlusc[,"patient.follow_ups.follow_up.2.targeted_molecular_therapy"]=="yes"),
                             which(treatlusc[,"patient.follow_ups.follow_up.3.targeted_molecular_therapy"]=="yes"),
                             which(treatlusc[,"patient.follow_ups.follow_up.4.targeted_molecular_therapy"]=="yes"),
                             which(treatlusc[,"patient.targeted_molecular_therapy"]=="yes"))))

target.luad <- sort(unique(c(which(treatluad[,"patient.drugs.drug.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatluad[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatluad[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatluad[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatluad[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatluad[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatluad[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatluad[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatluad[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatluad[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="targeted molecular therapy"),
                             which(treatlusc[,"patient.follow_ups.follow_up.targeted_molecular_therapy"]=="yes"),
                             which(treatlusc[,"patient.follow_ups.follow_up.2.targeted_molecular_therapy"]=="yes"),
                             which(treatlusc[,"patient.follow_ups.follow_up.3.targeted_molecular_therapy"]=="yes"),
                             which(treatlusc[,"patient.follow_ups.follow_up.4.targeted_molecular_therapy"]=="yes"),
                             which(treatlusc[,"patient.targeted_molecular_therapy"]=="yes"))))

immuno.lusc <- sort(unique(c(which(treatlusc[,"patient.drugs.drug.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.11.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatlusc[,"patient.drugs.drug.12.therapy_types.therapy_type"]=="immunotherapy"))))
immuno.luad <- sort(unique(c(which(treatluad[,"patient.drugs.drug.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatluad[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatluad[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatluad[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatluad[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatluad[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatluad[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatluad[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatluad[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="immunotherapy"),
                             which(treatluad[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="immunotherapy"))))

muu.lusc <- sort(unique(c(which(treatlusc[,"patient.drugs.drug.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.11.therapy_types.therapy_type"]=="vaccine"),
                             which(treatlusc[,"patient.drugs.drug.12.therapy_types.therapy_type"]=="vaccine"),
                          which(treatlusc[,"patient.drugs.drug.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.11.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.12.therapy_types.therapy_type"]=="ancillary"),
                          which(treatlusc[,"patient.drugs.drug.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.11.therapy_types.therapy_type"]=="other, specify"),
                          which(treatlusc[,"patient.drugs.drug.12.therapy_types.therapy_type"]=="other, specify"))))
muu.luad <- sort(unique(c(which(treatluad[,"patient.drugs.drug.therapy_types.therapy_type"]=="vaccine"),
                             which(treatluad[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="vaccine"),
                             which(treatluad[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="vaccine"),
                             which(treatluad[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="vaccine"),
                             which(treatluad[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="vaccine"),
                             which(treatluad[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="vaccine"),
                             which(treatluad[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="vaccine"),
                             which(treatluad[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="vaccine"),
                             which(treatluad[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="vaccine"),
                             which(treatluad[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="vaccine"),
                          which(treatluad[,"patient.drugs.drug.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="ancillary"),
                          which(treatluad[,"patient.drugs.drug.therapy_types.therapy_type"]=="other, specify"),
                          which(treatluad[,"patient.drugs.drug.2.therapy_types.therapy_type"]=="other, specify"),
                          which(treatluad[,"patient.drugs.drug.3.therapy_types.therapy_type"]=="other, specify"),
                          which(treatluad[,"patient.drugs.drug.4.therapy_types.therapy_type"]=="other, specify"),
                          which(treatluad[,"patient.drugs.drug.5.therapy_types.therapy_type"]=="other, specify"),
                          which(treatluad[,"patient.drugs.drug.6.therapy_types.therapy_type"]=="other, specify"),
                          which(treatluad[,"patient.drugs.drug.7.therapy_types.therapy_type"]=="other, specify"),
                          which(treatluad[,"patient.drugs.drug.8.therapy_types.therapy_type"]=="other, specify"),
                          which(treatluad[,"patient.drugs.drug.9.therapy_types.therapy_type"]=="other, specify"),
                          which(treatluad[,"patient.drugs.drug.10.therapy_types.therapy_type"]=="other, specify"))))

## Noukitaan vain chemoa saavat
ei.targ.lusc <- chemot.lusc[which(!chemot.lusc%in%target.lusc)]
eimuuta.lusc <- ei.targ.lusc[which(!ei.targ.lusc%in%muu.lusc)]

ei.targ.luad <- chemot.luad[which(!chemot.luad%in%target.luad)]
ei.immuno.luad <- ei.targ.luad[which(!ei.targ.luad%in%immunot.luad)]
eimuuta.luad <- ei.immuno.luad[which(!ei.immuno.luad%in%muu.luad)]

### HUOM! Näemmä radiation unohdettu! Siitä johtunee ero.


## otetaan kaikki clindata näistä
lusc.clin <- col.lusc[eimuuta.lusc,]
luad.clin <- col.luad[eimuuta.luad,]

## Otetaan RNASeq2GeneNorma
### ONGELMA! Potilaskoodit erilaiset clin ja RNA Seq, jossa pidemmät koodit. Lyhyemmillä clin koodeilla useampi vaste, koska
### pidemmässä loppu on erilainen. Mikä valitaan?
lusc.idt <- which(grepl(paste(rownames(col.lusc)[eimuuta.lusc],collapse="-01A|"),colnames(assay(lusc@ExperimentList@listData$`LUSC_RNASeq2GeneNorm-20160128`))))
luad.idt <- which(grepl(paste(rownames(col.luad)[eimuuta.luad],collapse="-01A|"),colnames(assay(luad@ExperimentList@listData$`LUAD_RNASeq2GeneNorm-20160128`))))
lusc.RNAS2 <- assay(lusc@ExperimentList@listData$`LUSC_RNASeq2GeneNorm-20160128`)[,lusc.idt]
luad.RNAS2 <- assay(luad@ExperimentList@listData$`LUAD_RNASeq2GeneNorm-20160128`)[,luad.idt]

## Otetaan mRNAArray LUAD
luad.idt.mRNA <- which(grepl(paste(rownames(col.luad)[eimuuta.luad],collapse="|"),colnames(assay(luad@ExperimentList@listData$`LUAD_mRNAArray-20160128`))))
luad.mRNAA <- assay(luad@ExperimentList@listData$`LUAD_mRNAArray_TX_g4502a-20160128`)[,luad.idt.mRNA]

## Otetaan mRNAArray_TX
lusc.idt.TX <- which(grepl(paste(rownames(col.lusc)[eimuuta.lusc],collapse="|"),colnames(assay(lusc@ExperimentList@listData$`LUSC_mRNAArray_TX_g4502a-20160128`))))
lusc.mRNAA.TX <- assay(lusc@ExperimentList@listData$`LUSC_mRNAArray_TX_g4502a-20160128`)[,lusc.idt.TX]

## Otetaan mRNAArray_TX_ht_hg
lusc.idt.TXht <- which(grepl(paste(rownames(col.lusc)[eimuuta.lusc],collapse="|"),colnames(assay(lusc@ExperimentList@listData$`LUSC_mRNAArray_TX_ht_hg_u133a-20160128`))))
lusc.mRNAA.TXht <- assay(lusc@ExperimentList@listData$`LUSC_mRNAArray_TX_ht_hg_u133a-20160128`)[,lusc.idt.TXht]

## Otetaan mRNAArray_huex
lusc.idt.huex <- which(grepl(paste(rownames(col.lusc)[eimuuta.lusc],collapse="|"),colnames(assay(lusc@ExperimentList@listData$`LUSC_mRNAArray_huex-20160128`))))
lusc.mRNAA.huex <- assay(lusc@ExperimentList@listData$`LUSC_mRNAArray_huex-20160128`)[,lusc.idt.huex]


# OSCAR-testi
lusc.patients <- rownames(col.lusc)[eimuuta.lusc]
lusc.colid <- which(grepl(paste(rownames(col.lusc)[eimuuta.lusc],collapse="|"),rownames(col.lusc)))
lusc.y <- cbind(as.numeric(col.lusc[lusc.colid,"days_to_death"]),as.numeric(col.lusc[lusc.colid,"vital_status"]),as.numeric(col.lusc[lusc.colid,"days_to_last_followup"]))
lusc.y <- t(apply(lusc.y,MARGIN=1,FUN=function(z){
  if(is.na(z[1])){
    return(c(z[3],z[2]))
  }else{
    return(c(z[1],z[2]))
  }}))
tmp<-oscar(t(as.matrix(lusc.RNAS2[3:19,])),y=as.matrix(lusc.y),family="cox")
visu(tmp.lusc, y=c("goodness"),legend=F)
feat(tmp.lusc,k=6)


luad.patients <- rownames(col.luad)[eimuuta.luad]
luad.colid <- which(grepl(paste(rownames(col.luad)[eimuuta.luad],collapse="|"),rownames(col.luad)))
luad.y <- cbind(as.numeric(col.luad[luad.colid,"days_to_death"]),as.numeric(col.luad[luad.colid,"vital_status"]),as.numeric(col.luad[luad.colid,"days_to_last_followup"]))
luad.y <- t(apply(luad.y,MARGIN=1,FUN=function(z){
  if(is.na(z[1])){
    return(c(z[3],z[2]))
  }else{
    return(c(z[1],z[2]))
  }}))
tmp.luad<-oscar(t(as.matrix(luad.RNAS2[1:30,])),y=as.matrix(luad.y),family="cox")
visu(tmp.luad, y=c("goodness"),legend=F)
feat(tmp.luad,k=6)
## CHECK, ettei nollasarakkeita x-matriisissa!
cv.luad <- cv.oscar(tmp.luad)
             

#Yhteis
x.all <- cbind(c(rep(1,107),rep(0,19)),rbind(t(as.matrix(luad.RNAS2[1:30,])),t(as.matrix(lusc.RNAS2[1:30,]))))
colnames(x.all)[1]<-"luad"
tmp.all <- oscar(x=x.all,y=rbind(luad.y,lusc.y),family="cox")
cv.all <- cv.oscar(tmp.all, fold = 5)

