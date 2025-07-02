#####################################
#CHTN stage III analysis
#12/03/21
######################################
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(Rphenograph)
library(heatmaply)
library(gcookbook)
library(viridis)
rm(list=ls())
#calculate mean counts
dir1 <- "/Users/msb/Google Drive/Dartmouth/Raman Lab/TSA/CHTN/211202 pheno/"
dir2 <- "/Users/msb/Google Drive/Dartmouth/Raman Lab/TSA/CHTN/211202 pheno/figures/"
dir3 <- "/Users/msb/Google Drive/Dartmouth/Raman Lab/TSA/CHTN/"
final_dir <- "/Users/msb/Google Drive/Dartmouth/Raman Lab/TSA/CHTN/Final figures/"


colsEMT <- brewer.pal(7, "RdYlBu")
colsPar <- brewer.pal(10, "PiYG")
EMT_col <- c(colsEMT[7], colsEMT[6], colsEMT[5], colsEMT[3], colsEMT[2], colsEMT[1], colsPar[2])

cols_spectral <- brewer.pal(11, "Spectral")
cols_set1 <- brewer.pal(9, "Set1")
pretty_grad <- c("#003f5c", "#2f4b7c", "#665191", "#a05195", "#d45087", "#f95d6a","#ff7c43", "#ffa600")
blues <- brewer.pal(9, "Blues")

filepath <- list.files(path = dir1, pattern = "cell_seg_data.txt")
filepath <- filepath[-153]

file_names <- c()
for (i in 1:length(filepath)) {
  TMA <- strsplit(filepath, "CHTN")[[i]][2]
  TMA <- strsplit(TMA, "_")[[1]][1]
  tumor <- strsplit(filepath, ",")[[i]][2]
  image <- strsplit(filepath, "]")[[i]][1]
  image <- strsplit(image,",")[[1]][3]
  file_names <- as.data.frame(rbind(file_names, cbind(TMA, MSI = paste0(tumor, "_", image), ID = paste0(TMA, "_", tumor, "_", image))))
}

#read in patient data
pat_dat_all <- read.csv(paste0(dir3, "CHTN_StageIII_data_new.csv"), stringsAsFactors = F)
summary(pat_dat_all)
head(pat_dat_all)
ID <- c()
for (i in 1:dim(pat_dat_all)[1]){
  x <- paste0(pat_dat_all$CASESET[i], "_", pat_dat_all$MSI[i])
  ID <- as.data.frame(rbind(ID, x))
}
Therapy <- as.data.frame(rep("NA", dim(pat_dat_all)[1]))
colnames(Therapy) <- c("Therapy")
for (i in 1:dim(pat_dat_all)[1]){
  if (is.na(pat_dat_all$CHEMO_THERAPY_[i]) && is.na(pat_dat_all$RADIATION_THERAPY_[i]) && is.na(pat_dat_all$HORMONE_THERAPY_[i]) && is.na(pat_dat_all$OTHER_THERAPY_[i])){
    Therapy[i,] <- NA;
  }else if (pat_dat_all$OTHER_THERAPY_[i]=="0=Yes"){
    Therapy[i,] <- "Other";
  }else if (pat_dat_all$CHEMO_THERAPY_[i]=="1=Yes" && pat_dat_all$RADIATION_THERAPY_[i]=="0=No" && pat_dat_all$HORMONE_THERAPY_[i]== "0=No"){
    Therapy[i,] <- "Chemo";
  }else if (pat_dat_all$CHEMO_THERAPY_[i]=="0=No" && pat_dat_all$RADIATION_THERAPY_[i]=="1=Yes" && pat_dat_all$HORMONE_THERAPY_[i]== "0=No"){
    Therapy[i,] <- "Rad";
  }else if (pat_dat_all$CHEMO_THERAPY_[i]=="0=No" && pat_dat_all$RADIATION_THERAPY_[i]=="0=No" && pat_dat_all$HORMONE_THERAPY_[i]== "1=Yes"){
    Therapy[i,] <- "Hormone";
  }else (Therapy[i,] <- "Combined")
}
pat_dat_all$Therapy <- Therapy
pat_dat <- as.data.frame(cbind(ID = ID, MSI = pat_dat_all$MSI, TMA=pat_dat_all$CASESET, Row = pat_dat_all$ROW, Column = pat_dat_all$COLUMN, Sample = pat_dat_all$Sample, 
                               Case = pat_dat_all$CASE_IDENTIFIER, Age = pat_dat_all$AGE_AT_DIAGNOSIS, Therapy = Therapy, Status_at_death = pat_dat_all$CANCER_STATUS_AT_DEATH_,
                               Hist = pat_dat_all$MOST_PROMINENT_HISTOLOGICAL_TYP_, CauseofDeath = pat_dat_all$CANCER_STATUS_AT_DEATH_,
                               Score = pat_dat_all$TOTAL_SCORE_, Status = pat_dat_all$VITAL_STATUS_, Overall_Survival = pat_dat_all$FU_MONTHS_OS, 
                               Recurrence_site = pat_dat_all$TYPE_SITE_OF_FIRST_NON_BREAST_R_, RFS = pat_dat_all$FU_MONTHS_RFS, 
                               Met_RFS = pat_dat_all$FU_MONTHS_OUTSIDE_IPS_RFS, Local_RFS = pat_dat_all$FU_MONTHS_IPS_RFS,
                               ER = pat_dat_all$ER_SCORE, PR = pat_dat_all$PR_SCORE, HER2 = pat_dat_all$HER2_SCORE))
missing <- pat_dat[which(is.na(pat_dat$MSI==TRUE)),]
pat_dat <- pat_dat[-which(is.na(pat_dat$MSI==TRUE)),]

#reorder
pat_dat_18 <- pat_dat[which(pat_dat$TMA==18),]
pat_dat_18 <- pat_dat_18[order(pat_dat_18$MSI),]
pat_dat_19 <- pat_dat[which(pat_dat$TMA==19),]
pat_dat_19 <- pat_dat_19[order(pat_dat_19$MSI),]
pat_dat <- rbind(pat_dat_18, pat_dat_19)

#tumor_dat <- pat_dat[which(pat_dat$Sample=="Patient"),]
is.na(pat_dat) <- "NA"
#define subtype
#edit patient dat - [33,] has no ER PR HER2 data
pat_dat$ER[33] <- "2"
pat_dat$PR[33] <- "2"
pat_dat$HER2[33] <- "2"
subtype <- as.data.frame(rep("NA", dim(pat_dat)[1]))
colnames(subtype) <- c("Subtype")
for (i in 1:dim(pat_dat)[1]){
  if (is.na(pat_dat$ER[i]) && is.na(pat_dat$PR[i]) && is.na(pat_dat$HER2[i])){
    subtype[i,] <- NA;
  }else if (pat_dat$ER[i]=="0" && pat_dat$PR[i]=="0" && pat_dat$HER2[i]== "0"){
    subtype[i,] <- "TNBC";
  }else if (pat_dat$ER[i]=="1" && pat_dat$PR[i]=="1" && pat_dat$HER2[i]== "0"){
    subtype[i,] <- "HR+";
  }else if (pat_dat$ER[i]=="1" && pat_dat$PR[i]=="1" && pat_dat$HER2[i]== "1"){
    subtype[i,] <- "HER2+";
  }else if (pat_dat$ER[i]=="1" && pat_dat$PR[i]=="0" && pat_dat$HER2[i]== "0"){
    subtype[i,] <- "HR+";
  }else if (pat_dat$ER[i]=="0" && pat_dat$PR[i]=="1" && pat_dat$HER2[i]== "0"){
    subtype[i,] <- "HR+";
  }else if (pat_dat$ER[i]=="0" && pat_dat$PR[i]=="0" && pat_dat$HER2[i]== "1"){
    subtype[i,] <- "HER2+";
  }else if (pat_dat$ER[i]=="1" && pat_dat$PR[i]=="0" && pat_dat$HER2[i]== "1"){
      subtype[i,] <- "HER2+";
  }else if (pat_dat$ER[i]=="0" && pat_dat$PR[i]=="1" && pat_dat$HER2[i]== "1"){
    subtype[i,] <- "HER2+";
  }else if (pat_dat$ER[i]=="2" && pat_dat$PR[i]=="2" && pat_dat$HER2[i]== "2"){
    subtype[i,] <- "Unknown+";
  }else (subtype[i,] <- NA)
}
#define subtype HR+ or HR- regardless of HER2 status
subtype2 <- as.data.frame(rep("NA", dim(pat_dat)[1]))
colnames(subtype2) <- c("Subtype2")
for (i in 1:dim(pat_dat)[1]){
  if (is.na(pat_dat$ER[i]) && is.na(pat_dat$PR[i]) && is.na(pat_dat$HER2[i])){
    subtype2[i,] <- NA;
  }else if (pat_dat$ER[i]=="0" && pat_dat$PR[i]=="0"){
    subtype2[i,] <- "HR-";
  }else if (pat_dat$ER[i]=="1" && pat_dat$PR[i]=="1"){
    subtype2[i,] <- "HR+";
  }else if (pat_dat$ER[i]=="1" && pat_dat$PR[i]=="0"){
    subtype2[i,] <- "HR+";
  }else if (pat_dat$ER[i]=="0" && pat_dat$PR[i]=="1"){
    subtype2[i,] <- "HR+";
  }else if (pat_dat$ER[i]=="2" && pat_dat$PR[i]=="2" && pat_dat$HER2[i]== "2"){
    subtype2[i,] <- "Unknown+";
  }else (subtype2[i,] <- NA)
}
Met_status <- as.data.frame(rep("NA", dim(pat_dat)[1]))
colnames(Met_status) <- c("Met_status")
for (i in 1:dim(pat_dat)[1]){
  if (is.na(pat_dat$Recurrence_site[i])){
    Met_status[i,] <- NA;
  }else if (pat_dat$Recurrence_site[i]=="0=None"){
    Met_status[i,] <- "None";
  }else (Met_status[i,] <- "Metastasis")
}
pat_dat <- rbind(pat_dat_18, pat_dat_19)
pat_dat$Subtype <- subtype[,1]
pat_dat$Subtype2 <- subtype2[,1]
pat_dat$Met_status <- Met_status[,1]

####################################################
#HET & EMT SCORES
#Het <- read.csv(paste0(dir1, "Het_test_res.csv"), stringsAsFactors = F)
#Het <- Het[order(Het$name),]

Het2 <- read.csv(paste0(dir1, "Het_scores.csv"), stringsAsFactors = F)
Het2 <- Het2[order(Het2$name),]

#EMT <- read.csv(paste0(dir1, "scores_new.csv"), stringsAsFactors = F)
#EMT <- EMT[order(EMT$name),]

EMT2 <- read.csv(paste0(dir1, "EMT_scores.csv"), stringsAsFactors = F)
EMT2 <- EMT2[order(EMT2$name),]
EMT_new_weights <- read.csv(paste0(dir1, "EMT_scores_new.csv"), stringsAsFactors = F)
EMT_new_weights <- EMT_new_weights[order(EMT_new_weights$name),]

Scores <- as.data.frame(cbind(sample = Het2$name, Het_score = Het2$y_pred, EMT_score = EMT_new_weights$score))
Scores$EMT_score <- as.numeric(Scores$EMT_score)
Scores$Het_score <- factor(Scores$Het_score, levels = c("high", "mid", "low"))

namecheck <- as.data.frame(cbind(pat_dat$MSI, Scores$sample))

ggplot(pat_dat_tum, aes(x=V1, fill = Het_score)) +
  geom_bar() + 
  scale_fill_manual(name="Het Score", values = viridis(3)) +
  labs(x="Clone", y="Percent of Images", title = "Evol Tumor Heterogeneity Scores") +
  #scale_fill_viridis()+
  #facet_wrap(~V1) +
  theme_minimal() +
  theme(title =element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
#plot het and EMT scores
ppi=300
png(paste0(dir2, "EMTbyHet_dotplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(Scores, aes(x = Het_score, y= EMT_score, fill = Het_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 0.015, color = NA) +
  geom_boxplot(width=0.4, color = "gray24", fill = NA, size =0.4) +
  scale_fill_manual(name="Het Score", 
                    values = viridis(3)) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x="Heterogeneity Score", y="EMT Score", title = "CHTN Score Distribution") +
  #facet_wrap(~clone) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()


###########################
#EMT and het scores against prognosis

#annotate Patient data with Scores
pat_dat <- cbind(pat_dat, EMT_score = Scores$EMT_score, Het_score = Scores$Het_score)
#write.csv(pat_dat, paste0(dir3, "Patient_data_wScores_new.csv"))
pat_dat <- read.csv(paste0(dir3, "Patient_data_wScores_new.csv"), stringsAsFactors = F)

#remove bad cores
bad_ID <- c("18_1_B", "18_1_O", "18_2_D", "18_12_M", "18_1_C", "18_2_F", "18_2_R", "18_3_O", "18_5_Q", "18_6_C", "18_6_J","18_10_P",
            "18_7_B", "18_7_P", "18_9_O", "18_11_I", "18_11_J", "19_2_I", "19_3_C", "19_5_C", "19_5_R", "19_7_I", "19_9_C", "19_9_Q", "19_10_Q")
bad_cores <- c()
for (i in 1:length(bad_ID)){
  bad_cores <- rbind(bad_cores, which(pat_dat$V1==bad_ID[i]))
}
pat_dat <- pat_dat[-bad_cores,]
pat_dat_tum <- pat_dat[which(pat_dat$Sample=="Patient"),]

pat_dat_tum$Het_score <- factor(pat_dat_tum$Het_score, levels = c("high", "mid", "low"))

pat_dat_tum$EMT_score_disc <- cut(pat_dat_tum$EMT_score, c(0,0.29,0.69,1), 
                                  labels = c("Epi","Int", "Mes"))
pat_dat_tum$EMT_score_disc <- factor(pat_dat_tum$EMT_score_disc, 
                                     levels = c("Epi", "Int", "Mes"),
                                     labels = c("Epithelial", "Intermediate", "Mesenchymal"))


ppi=300
png(paste0(final_dir, "EMTbyHet_barplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(pat_dat_tum, aes(x = Het_score, fill = EMT_score_disc)) +
  geom_bar(position="fill")+
  scale_fill_manual(name="EMT Score", values = c("Epithelial" = EMT_col[1], "Intermediate" = "goldenrod2", "Mesenchymal" = EMT_col[6])) +
  #scale_fill_manual(name="EMT Score", values = c("Epi" = blues[4], "Int" = blues[6], "Mes" = blues[8])) +
  labs(x="Heterogeneity Score", y="Proportion of cells", title = "EMT and Heterogeneity Score Correlation") +
  #facet_wrap(~Subtype2) +
  theme_classic() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

ppi=300
png(paste0(final_dir, "EMTbyHet_dotplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(pat_dat_tum, aes(x = Het_score, y= EMT_score, fill = Het_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 0.015, color = NA) +
  geom_boxplot(width=0.4, color = "gray24", fill = NA, size =0.4) +
  scale_fill_manual(name="Het Score", 
                    values = c("high" = c("#FC4E07"), "mid" = c("#E7B800"), "low" = c("#0097a2"))) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x="Heterogeneity Score", y="EMT Score", title = "EMT Score Distribution By Heterogeneity Score") +
  facet_wrap(~Subtype2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

ppi=300
png(paste0(dir2, "Recurrence_EMTbyHet_dotplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(pat_dat_tum, aes(x = Het_score, y= EMT_score, fill = Het_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 0.015, color = NA) +
  geom_boxplot(width=0.4, color = "gray24", fill = NA, size =0.4) +
  scale_fill_manual(name="Het Score", 
                    values = c("high" = cols_set1[1], "mid" = cols_set1[3], "low" = cols_set1[2])) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x="Heterogeneity Score", y="EMT Score", title = "CHTN EMT Score Distribution By Recurrence") +
  facet_wrap(~Met_status) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

pat_dat_tum$EMT_score_disc <- cut(pat_dat_tum$EMT_score, c(0,0.29,0.69,1), 
                                  labels = c("Epi","Int", "Mes"))

ggplot(pat_dat_tum, aes(x = Het_score, y= EMT_score, fill = Het_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 0.015, color = NA) +
  #geom_boxplot(width=0.4, color = "gray24", fill = NA, size =0.4) +
  scale_fill_manual(name="Het Score", values = viridis(3)) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x="Heterogeneity Score", y="EMT Score", title = "CHTN EMT Score Distribution By Recurrence") +
  #facet_wrap(~Met_status) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))

ggplot(pat_dat_tum, aes(x = Het_score, y= EMT_score, fill = Het_score)) 
  
EMT_het <- matrix(c(2, 3, 8, 13, 46, 23, 14, 14, 1), ncol=3, byrow=TRUE)
rownames(EMT_het) <- c('High','Mid','Low ')
colnames(EMT_het) <- c('Epithelial','Intermediate','Mesenchymal')
#EMT_het <- as.table(EMT_het)

fisher.test(EMT_het)
library(vcd)

ppi=300
png(paste0(final_dir, "EMT_Het_correlation_mosaic.png"), height = 5*ppi, width = 5*ppi, res=ppi )
mosaicplot(EMT_het, shade = TRUE, xlab = "Heterogeneity  Score", ylab = "EMT Score", border = "black", main = "Score Correlation")
dev.off()

png(paste0(final_dir, "EMT_Het_correlation_balloon.png"), height = 4*ppi, width = 5.5*ppi, res=ppi )
ggballoonplot(EMT_het, size = "value", fill = "value", size.range = c(5, 25), show.label = TRUE, font.label = c(11, "white")) + scale_fill_viridis_c(option = "D") + 
  labs(x= "EMT Score", y = "Heterogeneity Score", main = "EMT and Heterogeneity Score Correlation") + theme_bw()
dev.off()
##########################################################################################################################
#Survival analysis
##########################################################################################################################

library(ggplot2)
library(survival)
library(survminer)
library(ggfortify)
library(cowplot)
library(patchwork)

scale_colour_discrete <- scale_colour_colorblind
?plot.survfit
?survfit

#colorblind pallete
c("#0097a2", "#E7B800", "#FC4E07")

pat_dat_tum <- pat_dat_tum[-18,]
pat_dat_tum$Status <- as.numeric(factor(pat_dat_tum$Status, labels = c(0,1))) #makes it 2 and 1 for some reason but it still works??
pat_dat_tum$Overall_Survival <- as.numeric(pat_dat_tum$Overall_Survival)
pat_dat_tum$RFS <- as.numeric(pat_dat_tum$RFS)
pat_dat_tum$Met_RFS <- as.numeric(pat_dat_tum$Met_RFS)
#set mid as reference
pat_dat_tum$Het_score = relevel(factor(pat_dat_tum$Het_score), ref = "mid")
pat_dat_tum$Subtype2 = relevel(factor(pat_dat_tum$Subtype2), ref = "HR+")
pat_dat_tum$Met_status <- factor(pat_dat_tum$Met_status, labels = c(1,0))
pat_dat_tum$Therapy <- factor(pat_dat_tum$Met_status, labels = c(1,0))

hist(pat_dat_tum$EMT_score)
Scores$count <- rep(1, length(Scores))
ggplot(Scores, aes(x = count, y= EMT_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 0.015, color = NA) +
  geom_boxplot(width=0.4, color = "gray24", fill = NA, size =0.4)

#split EMT score into quartiles 
summary(pat_dat_tum$EMT_score)
pat_dat_tum$EMT_score_disc <- cut(pat_dat_tum$EMT_score, c(0,0.29,0.69,1), 
                                  labels = c("Epi","Int", "Mes"))



#pat_dat_tum <- pat_dat_tum[which(pat_dat_tum$TMA=="19"),]

#pat_dat_tum_sub <- pat_dat_tum[-which(pat_dat_tum$Recurrence_site=="4=Never Disease Free"),]


#basic survival 
Subtype_fit <- survfit(Surv(Overall_Survival, Status) ~ Subtype2, data = pat_dat_tum)
p1 <- ggsurvplot(Subtype_fit, pval=TRUE, pval.method = FALSE, conf.int = FALSE, risk.table = "absolute", legend = "none", legend.labs = c("HR pos", "HR neg")) + 
  ggtitle("Overall survival by Hormone Status")
Therapy_fit <- survfit(Surv(Overall_Survival, Status) ~ Therapy, data = pat_dat_tum)
p2 <- ggsurvplot(Therapy_fit, pval=TRUE, pval.method = TRUE, conf.int = FALSE, risk.table = "absolute", legend = "none", legend.labs = c("Chemo", "Combined", "Hormone", "Rad")) + 
  ggtitle("Overall Survival by Therapy")
Therapy_fitHRpos <- survfit(Surv(Overall_Survival, Status) ~ Therapy, data = pat_dat_tum[which(pat_dat_tum$Subtype2=="HR+"),])
Therapy_fitHRneg <- survfit(Surv(Overall_Survival, Status) ~ Therapy, data = pat_dat_tum[which(pat_dat_tum$Subtype2=="HR-"),])
p3 <- ggsurvplot(Therapy_fitHRpos, pval=TRUE, pval.method = TRUE, conf.int = FALSE, risk.table = "absolute", legend = "none") + ggtitle("All survival byTherapy HR+")
p4 <- ggsurvplot(Therapy_fitHRneg, pval=TRUE, pval.method = TRUE, conf.int = FALSE, risk.table = "absolute", legend = "none") + ggtitle("All survival by Therapy HR-")

ppi=300
png(paste0(dir2, "Overall_survival_byTherapy.png"), height = 6*ppi, width = 5*ppi, res = ppi)
p2
dev.off()

###########################
#by Het score
###########################

#Overall Survival 


Het_fit <- survfit(Surv(Overall_Survival, Status) ~ Het_score, data = pat_dat_tum)
Het_fit_2 <- survfit(Surv(Overall_Survival, Status) ~ Het_score + EMT_score_disc, data = pat_dat_tum)
EMT_fit <- survfit(Surv(Overall_Survival, Status) ~ EMT_score_disc, data = pat_dat_tum)
Het_fit_HRpos <- survfit(Surv(Overall_Survival, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype2=="HR+"),])
Het_fit_HRneg <- survfit(Surv(Overall_Survival, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype2=="HR-"),])
#Het_fit_TNBC <- survfit(Surv(Overall_Survival, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype=="TNBC"),])
#Het_fit_HER2 <- survfit(Surv(Overall_Survival, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype=="HER2+"),])
#Het_fit_HR <- survfit(Surv(Overall_Survival, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype=="HR+"),])

c("#4575B4","goldenrod2", "#D73027")
p1 <- ggsurvplot(Subtype_fit, pval=TRUE, pval.method = TRUE, conf.int = FALSE, risk.table = "nrisk_cumevents", legend = "none") + ggtitle("All survival by Subtype")
p2 <- ggsurvplot(Het_fit, pval=FALSE, pval.method = FALSE, conf.int = FALSE, #legend = "right",
                 risk.table = "absolute", legend.labs = c("Mid", "High", "Low"), 
                 legend.title="Heterogeneity \n\ score", palette = c("grey27", "#FC4E07", "#0097a2" )) +
  ggtitle("Overall Survival")
p3 <- ggsurvplot(Het_fit_2, pval=FALSE, pval.method = FALSE, conf.int = FALSE, #legend = "right",
                 risk.table = "absolute", legend.labs = c("Mid + Epi", "Mid + Int", "Mid + Mes", 
                                                          "High + Epi", "High + Int", "High + Mes", 
                                                          "Low + Epi", "Low + Int", "Low + Mes"), 
                 legend.title="Combination score  \n\ (Heterogeneity + EMT)") + #, palette = c("grey27", "#FC4E07", "#0097a2" )) +
  ggtitle("Overall Survival")
p2_EMT <- ggsurvplot(Het_fit, pval=FALSE, pval.method = FALSE, conf.int = FALSE, #legend = "right",
                     risk.table = "absolute", legend.labs = c("Epithelial", "Intermediate", "Mesenchymal"), 
                     legend.title="EMT \n\ score", palette = c("#4575B4","goldenrod2", "#D73027")) +
  ggtitle("Overall Survival by EMT Score")
#p3 <- ggsurvplot(Het_fit_TNBC, pval=TRUE, pval.method = TRUE, conf.int = FALSE, risk.table = "nrisk_cumevents", legend = "none") +ggtitle("TNBC survival by Het score")
#p4 <- ggsurvplot(Het_fit_HER2, pval=TRUE, conf.int = FALSE, risk.table = "nrisk_cumevents", legend = "none") +ggtitle("HER2+ survival by Het score")
#p5 <- ggsurvplot(Het_fit_HR, pval=TRUE, conf.int = FALSE, risk.table = "nrisk_cumevents", legend = "none") +ggtitle("HR+ survival by Het score")
p3 <- ggsurvplot(Het_fit_HRpos, pval=FALSE, pval.method = FALSE, conf.int = FALSE, 
                 risk.table = "absolute", legend.labs = c("Mid", "High", "Low"), 
                 legend.title="Het \n\ score", palette = c("grey27", "#FC4E07", "#0097a2" )) +
  ggtitle("Overall Survival in HR positive disease")
p4 <- ggsurvplot(Het_fit_HRneg, pval=FALSE, pval.method = FALSE, conf.int = FALSE, 
                 risk.table = "absolute", legend.labs = c("Mid", "High", "Low"), 
                 legend.title="Het \n\ score", palette = c("grey27", "#FC4E07", "#0097a2" ), pval.coord = c(60, 0.15)) +
  ggtitle("Overall Survival in HR negative disease")

ppi=300
png(paste0(final_dir, "Overall_survival_HetbyEMT.png"), height = 8.5*ppi, width = 8*ppi, res = ppi)
p3
dev.off()

##########################
#Recurrance Free Survival

Subtype_fit <- survfit(Surv(RFS, Status) ~ Subtype2, data = pat_dat_tum)
Het_fit <- survfit(Surv(RFS, Status) ~ Het_score, data = pat_dat_tum)
Het_fit_HRpos <- survfit(Surv(RFS, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype2=="HR+"),])
Het_fit_HRneg <- survfit(Surv(RFS, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype2=="HR-"),])
#Het_fit_TNBC <- survfit(Surv(RFS, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype=="TNBC"),])
#Het_fit_HER2 <- survfit(Surv(RFS, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype=="HER2+"),])
#Het_fit_HR <- survfit(Surv(RFS, Status) ~ Het_score, data = pat_dat_tum[which(pat_dat_tum$Subtype=="HR+"),])

p1 <- ggsurvplot(Subtype_fit, pval=TRUE, pval.method = FALSE, conf.int = FALSE, risk.table = "absolute", legend = "none", legend.labs = c("HR pos", "HR neg")) + 
  ggtitle("RFS by Hormone Status")
p2 <- ggsurvplot(Het_fit, pval=TRUE, pval.method = FALSE, conf.int = FALSE, #legend = "right",
                 risk.table = "absolute", legend.labs = c("Mid", "High", "Low"), 
                 legend.title="Heterogeneity \n\ score", palette = c("grey27", "#FC4E07", "#0097a2" )) +
  ggtitle("Relapse Free Survival")
p3 <- ggsurvplot(Het_fit_HRpos, pval=FALSE, pval.method = FALSE, conf.int = FALSE, 
                 risk.table = "absolute", legend.labs = c("Mid", "High", "Low"), 
                 legend.title="Het \n\ score", palette = c("grey27", "#FC4E07", "#0097a2" )) +
  ggtitle("Relapse Free Survival in HR positive disease")
p4 <- ggsurvplot(Het_fit_HRneg, pval=TRUE, pval.method = FALSE, conf.int = FALSE, 
                 risk.table = "absolute", legend.labs = c("Mid", "High", "Low"), 
                 legend.title="Het \n\ score", palette = c("grey27", "#FC4E07", "#0097a2" ), pval.coord = c(60, 0.15)) +
  ggtitle("Relapse Free Survival in HR negative disease")


png(paste0(dir2, "RFS_by_HetScore.png"), height = 6*ppi, width = 5*ppi, res = ppi)
p2
dev.off()



#plot het scores by subtype
coxph(Surv(Overall_Survival, Status) ~ Subtype2, data = pat_dat_tum)
univ_Het_fit <- coxph(Surv(Overall_Survival, Status) ~ Age +Subtype2 +  Het_score, data = pat_dat_tum)
univ_EMT_fit <- coxph(Surv(Overall_Survival, Status) ~  Age + Subtype2 +EMT_score_disc, data = pat_dat_tum)
coxph(Surv(Overall_Survival, Status) ~ Het_score + EMT_score_disc, data = pat_dat_tum)
coxph(Surv(Overall_Survival, Status) ~ Subtype2 + Het_score, data = pat_dat_tum)
coxph(Surv(Overall_Survival, Status) ~ Subtype2 + Het_score + EMT_score_disc +Age, data = pat_dat_tum)
coxph(Surv(Overall_Survival, Status) ~ Subtype2 + Het_score*EMT_score_disc, data = pat_dat_tum)

summary(univ_Het_fit)
summary(univ_EMT_fit)

fit_basic <- coxph(Surv(Overall_Survival, Status) ~ Subtype2 + Age, data = pat_dat_tum)
summary(fit_basic)

fit <- coxph(Surv(Overall_Survival, Status) ~  Age + Subtype2+ Het_score + EMT_score_disc  , data = pat_dat_tum)
summary(fit)  #output provides HR CIs
confint(fit)  #coefficient CIs
exp(confint(fit))  #Also HR CIs
summary$conf.int

ppi=300
png(paste0(dir2, "Forest_fill_HR_model.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggforest(fit, data = NULL,
         main = "Hazard ratio",
         cpositions = c(0.02, 0.22, 0.4),
         fontsize = 0.7,
         refLabel = "reference",
         noDigits = 2)
dev.off()


fit_RFS <- coxph(Surv(RFS, Status) ~ Age + Subtype2 + Het_score + EMT_score_disc, data = pat_dat_tum)
summary(fit_RFS)
#hormone positive
fitHRpos <- coxph(Surv(RFS, Status) ~ Age + Het_score + EMT_score_disc, data = pat_dat_tum[which(pat_dat_tum$Subtype2=="HR+"),])
summary(fitHRpos)
#confint(fitHRpos)

#hormone negative
fitHRneg <- coxph(Surv(RFS, Status) ~ Age  + Het_score + EMT_score_disc, data = pat_dat_tum[which(pat_dat_tum$Subtype2=="HR-"),])
summary(fitHRneg)

fit_intx <- coxph(Surv(Overall_Survival, Status) ~ Age + Subtype2 + Het_score*EMT_score_disc, data = pat_dat_tum)
summary(fit_intx)
ggforest(fit_intx)
fitHRneg_intx <- coxph(Surv(Overall_Survival, Status) ~ Age + Het_score*EMT_score_disc, data = pat_dat_tum[which(pat_dat_tum$Subtype2=="HR-"),])
summary(fitHRneg_intx)

##########################
#Forest Plots
##########################
forest_df <- pat_dat_tum
colnames(forest_df) <- c(colnames(pat_dat_tum[1:23]), "Subtype_old", "Subtype", "Met_status", "EMT_score", "Het_Score", "EMT_Score")

fit <- coxph(Surv(Overall_Survival, Status) ~ Age + Subtype + Het_Score + EMT_Score, data = forest_df)
png(paste0(final_dir, "Overall_surv_forest.png"), height = 6*ppi, width = 6*ppi, res = ppi)
ggforest(fit, main = "Multivariate Cox Proportional Hazard Model",)
dev.off()

fit_HRneg <- coxph(Surv(Overall_Survival, Status) ~ Age + Het_Score + EMT_Score, data = forest_df[which(forest_df$Subtype=="HR-"),])
png(paste0(final_dir, "Overall_surv_forest_HR-.png"), height = 6*ppi, width = 6*ppi, res = ppi)
ggforest(fit_HRneg, main = "Multivariate Cox Proportional Hazard Model \n\ HR negative disease",)
dev.off()

fit_HRpos <- coxph(Surv(Overall_Survival, Status) ~ Age + Het_Score + EMT_Score, data = forest_df[which(forest_df$Subtype=="HR+"),])
png(paste0(final_dir, "Overall_surv_forest_HR+.png"), height = 6*ppi, width = 6*ppi, res = ppi)
ggforest(fit_HRpos, main = "Multivariate Cox Proportional Hazard Model \n\ HR positive disease",)
dev.off()


Het_HR_labels <- c("High Het", "Mid Het (ref)", "Low Het")
Het_HR <- data.frame(yAxis = length(Het_HR_labels):1,boxOdds = c(3.91, 1, 1.43),
                              boxCILow = c(2, 0, 0.69),
                              boxCIHigh = c(7.6, 0, 2.38))


