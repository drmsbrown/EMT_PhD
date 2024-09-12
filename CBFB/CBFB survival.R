
BiocManager::install("biocLite.R")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  #barcode = listSamples, 
                  legacy = TRUE)


# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
BRCAMatrix <- assay(BRCARnaseqSE,"subtype") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")

BRCARnaseqSE@colData@listData[["paper_BRCA_Subtype_PAM50"]]

BRCAMatrix_CBFb <- as.data.frame(BRCAMatrix[which(rownames(BRCAMatrix)=="CBFB|865"),])
BRCAMatrix_CBFb$barcode <- rownames(BRCAMatrix_CBFb)
BRCAMatrix_CBFb$ID <- BRCARnaseqSE$patient
BRCAMatrix_CBFb$subtype <- BRCARnaseqSE@colData@listData[["paper_BRCA_Subtype_PAM50"]]


library(TCGAbiolinks)
# Survival Analysis SA

clinical_patient_Cancer <- GDCquery_clinic("TCGA-BRCA","clinical")


#dataBRCAcomplete <- log2(BRCA_rnaseqv2)
#dataBRCAcomplete <- log2(dataNorm)

dataBRCA_CBFb <- as.data.frame(log2(BRCAMatrix_CBFb$`BRCAMatrix[which(rownames(BRCAMatrix) == "CBFB|865"), ]`))
dataBRCA_CBFb$raw <- BRCAMatrix_CBFb$`BRCAMatrix[which(rownames(BRCAMatrix) == "CBFB|865"), ]`
dataBRCA_CBFb$ID <- BRCAMatrix_CBFb$ID
dataBRCA_CBFb$barcode <- BRCAMatrix_CBFb$barcode
dataBRCA_CBFb$subtype <- BRCAMatrix_CBFb$subtype
colnames(dataBRCA_CBFb) <- c("log2", "raw", "ID", "barcode", "subtype")

library(tidyverse)
#find duplicates and only keep one
dataBRCA_CBFb_rm <- dataBRCA_CBFb %>% distinct(ID, .keep_all = TRUE)
  
  
dataBRCA_CBFb_rm <- dataBRCA_CBFb_rm[order(dataBRCA_CBFb_rm$ID),]
clinical_patient_Cancer <- clinical_patient_Cancer[order(clinical_patient_Cancer$submitter_id),]


inboth <- clinical_patient_Cancer$submitter_id %in% dataBRCA_CBFb_rm$ID
clinical_patient_Cancer <- clinical_patient_Cancer[inboth,]
clinical_patient_Cancer$subtype <- dataBRCA_CBFb_rm$subtype
namecheck <- cbind(as.data.frame(clinical_patient_Cancer$submitter_id), as.data.frame(dataBRCA_CBFb_rm$ID))


tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient_Cancer,BRCAMatrix,
                                  Genelist = c("ESR1"), Survresult = TRUE,ThreshTop=0.69,ThreshDown=0.33)


tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient_Cancer[which(clinical_patient_Cancer$subtype=="LumB"),],BRCAMatrix,
                                  Genelist = c("CBFB|865"), Survresult = TRUE,ThreshTop=0.69,ThreshDown=0.33)


inboth <- cfu$submitter_id %in% dataBRCA_CBFb_rm$ID
dataBRCA_CBFb_rm <- dataBRCA_CBFb_rm[inboth,]

clin_exp <- cbind(cfu, dataBRCA_CBFb_rm)
clin_exp$vital_status2 <- factor(clin_exp$vital_status, labels = c(1,0))
clin_exp$exp_quart <- cut(clin_exp$log2, c(0,10.877,11.837,13.7), 
                          labels = c("Low","Mid", "High"))

clin_exp$exp_quart <- relevel(factor(clin_exp$exp_quart), ref = "Mid")
clin_exp <- clin_exp[-which(is.na(clin_exp$vital_status))]


library(ggplot2)
library(survival)
library(survminer)

fit <- survfit(Surv(days_to_death, vital_status2) ~ exp_quart, data = clin_exp)
p1 <- ggsurvplot(fit, pval=TRUE, pval.method = FALSE, conf.int = FALSE, risk.table = "absolute", legend = "none") + 
  ggtitle("Overall survival by Hormone Status")


cfu








cfu <- clinical_patient_Cancer[clinical_patient_Cancer[, "bcr_patient_barcode"] %in% substr(colnames(BRCAMatrix), 1, 12), ]
if ("days_to_last_followup" %in% colnames(cfu)){
  colnames(cfu)[grep("days_to_last_followup", colnames(cfu))] <-"days_to_last_follow_up"
  cfu <- as.data.frame(subset(cfu,select = c(
    "bcr_patient_barcode",
    "days_to_death",
    "days_to_last_follow_up",
    "vital_status")))
}
# Set alive death to inf
if (length(grep("alive", cfu$vital_status, ignore.case = TRUE)) > 0){
  cfu[grep("alive", cfu$vital_status, ignore.case = TRUE), "days_to_death"] <-"-Inf"
}
  

# Set dead follow up to inf
if (length(grep("dead", cfu$vital_status, ignore.case = TRUE)) > 0){
  cfu[grep("dead", cfu$vital_status, ignore.case = TRUE), "days_to_last_follow_up"] <-"-Inf"
}
  

cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]

cfu$days_to_death <- as.numeric(as.character(cfu$days_to_death))
cfu$days_to_last_follow_up <-
  as.numeric(as.character(cfu$days_to_last_follow_up))
rownames(cfu) <- cfu[, "bcr_patient_barcode"] #mod1

cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]

cfu_complete <- cfu

