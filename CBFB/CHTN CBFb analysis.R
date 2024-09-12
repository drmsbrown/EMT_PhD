#####################################
#CHTN stage III CBFb stain analysis
#4/17/22
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
final_dir <- "/Users/msb/Google Drive/Dartmouth/Raman Lab/TSA/CHTN/Final figures/CBFb/"


pat_dat <- read.csv(paste0(dir3, "Patient_data_wScores&CBFb.csv"), stringsAsFactors = F)
pat_dat <- pat_dat[1:35]

#remove bad cores
bad_ID <- c("18_1_B", "18_1_O", "18_2_D", "18_12_M", "18_1_C", "18_2_F", "18_2_R", "18_3_O", "18_5_Q", "18_6_C", "18_6_J","18_10_P",
            "18_7_B", "18_7_P", "18_9_O", "18_11_I", "18_11_J", "19_2_I", "19_3_C", "19_5_C", "19_5_R", "19_7_I", "19_9_C", "19_9_Q", "19_10_Q")
bad_cores <- c()
for (i in 1:length(bad_ID)){
  bad_cores <- rbind(bad_cores, which(pat_dat$ID==bad_ID[i]))
}
pat_dat <- pat_dat[-bad_cores,]
pat_dat_tum <- pat_dat[which(pat_dat$Sample=="Patient"),]

pat_dat_tum$Het_score <- factor(pat_dat_tum$Het_score, levels = c("high", "mid", "low"))

pat_dat_tum$EMT_score_disc <- cut(pat_dat_tum$EMT_score, c(0,0.29,0.69,1), 
                                  labels = c("Epi","Int", "Mes"))
pat_dat_tum$EMT_score_disc <- factor(pat_dat_tum$EMT_score_disc, 
                                     levels = c("Epi", "Int", "Mes"),
                                     labels = c("Epithelial", "Intermediate", "Mesenchymal"))

pat_dat_tum <- pat_dat_tum[-which(pat_dat_tum$CBFb_stain=="No tumor"),]
pat_dat_tum <- pat_dat_tum[-17,]

CBFb <- as.data.frame(rep("NA", dim(pat_dat_tum)[1]))
colnames(CBFb) <- c("CBFb")
for (i in 1:dim(pat_dat_tum)[1]){
  if (is.na(pat_dat_tum$CBFb_stain_comp[i])){
    CBFb[i,] <- NA;
  }else if (pat_dat_tum$CBFb_stain_comp[i]=="negative"){
    CBFb[i,] <- 0;
  }else (CBFb[i,] <- 1)
}
pat_dat_tum$CBFb <- CBFb[,1]


table(pat_dat_tum$CBFb)

#Does CBFb yes or no correlate with EMT or Het score?
fisher.test(pat_dat_tum$EMT_score_disc, pat_dat_tum$CBFb)
fisher.test(pat_dat_tum$Het_score, pat_dat_tum$CBFb)
#NOPE

#what about nuclear or cytoplasmic staining?
fisher.test(pat_dat_tum$EMT_score_disc, pat_dat_tum$CBFb_stain_comp)
fisher.test(pat_dat_tum$Het_score, pat_dat_tum$CBFb_stain_comp)
#NOPE

#what about H_score?
ggplot(pat_dat_tum, aes(x = EMT_score, y = H_score)) +
  geom_jitter()

ppi=300
png(paste0(final_dir, "CBFbHscore_EMTscoredotplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(pat_dat_tum, aes(x = EMT_score_disc, y = H_score, fill = H_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 2, color = "grey", fill = "grey") +
  geom_boxplot(width=0.4, color = "gray24", fill = NA, size =0.4) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x="EMT score", y="H_score", title = "H_score of CBFb staining by EMT score") +
  #facet_wrap(~Subtype2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

aov <- aov(H_score ~ EMT_score_disc, pat_dat_tum)
summary(aov)
#NOT SIGNIFICANT
aov <- aov(H_score ~ EMT_score_disc, pat_dat_tum[which(pat_dat_tum$CBFb_stain_comp=="cytoplasm"),])
summary(aov)
aov <- aov(H_score ~ EMT_score_disc, pat_dat_tum[which(pat_dat_tum$CBFb_stain_comp=="nuclear"),])
summary(aov)
#extra not significant for cytoplasmic vs nuclear staining


#het score?
ppi=300
png(paste0(final_dir, "CBFbHscore_Hetscoredotplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(pat_dat_tum, aes(x = Het_score, y = H_score, fill = H_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 2, color = "grey", fill = "grey") +
  geom_boxplot(width=0.4, color = "gray24", fill = NA, size =0.4) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x="Het score", y="H_score", title = "H_score of CBFb staining by Het score") +
  #facet_wrap(~Subtype2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

aov <- aov(H_score ~ Het_score, pat_dat_tum)
summary(aov)
#Also not significant

#combined EMT and het??
aov <- aov(H_score ~ Het_score + EMT_score_disc, pat_dat_tum)
summary(aov)
#STILL NO


pat_dat_tum$CBFb_stain_perc <- as.numeric(pat_dat_tum$CBFb_stain_perc)
aov <- aov(CBFb_stain_perc ~ Het_score, pat_dat_tum)
summary(aov)

pat_dat_tum$CBFb_stain_perc_disc <- cut(pat_dat_tum$CBFb_stain_perc, c(0,24,49,74, 100), 
                                  labels = c("1st","2nd", "3rd", "4th"))
fisher.test(pat_dat_tum$CBFb_stain_perc_disc, pat_dat_tum$EMT_score_disc)

#Check for differences in subtype
HR_pos <- pat_dat_tum[which(pat_dat_tum$Subtype2=="HR+"),]
HR_neg <- pat_dat_tum[which(pat_dat_tum$Subtype2=="HR-"),]

hist(HR_pos$H_score)
hist(HR_neg$H_score)

ppi=300
png(paste0(final_dir, "CBFbHscore_hist_bysubtype.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(pat_dat_tum, aes(x = H_log)) +
  geom_histogram(binwidth = 0.125)+
  #facet_wrap(~Subtype2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

ppi=300
png(paste0(final_dir, "CBFbHscore_dotplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(pat_dat_tum, aes(x = Subtype2, y = H_score, fill = H_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 2, color = "grey", fill = "grey") +
  geom_boxplot(width=0.4, color = "gray24", fill = NA, size =0.4) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x="Hormone Status", y="H_score", title = "H_score of CBFb staining by Hormone Status") +
  #facet_wrap(~Subtype2) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()
#t.test(HR_pos$H_score,HR_neg$H_score)
wilcox.test(HR_pos$H_score,HR_neg$H_score)


#Survival analysis
library(ggplot2)
library(survival)
library(survminer)
library(ggfortify)
library(cowplot)
library(patchwork)


pat_dat_tum$Status <- as.numeric(factor(pat_dat_tum$Status, labels = c(0,1))) #makes it 2 and 1 for some reason but it still works??
pat_dat_tum$Overall_Survival <- as.numeric(pat_dat_tum$Overall_Survival)
pat_dat_tum$RFS <- as.numeric(pat_dat_tum$RFS)
pat_dat_tum$Met_RFS <- as.numeric(pat_dat_tum$Met_RFS)
#set mid as reference
pat_dat_tum$Het_score = relevel(factor(pat_dat_tum$Het_score), ref = "mid")
pat_dat_tum$Subtype2 = relevel(factor(pat_dat_tum$Subtype2), ref = "HR+")
pat_dat_tum$Met_status <- factor(pat_dat_tum$Met_status, labels = c(1,0))
pat_dat_tum$Therapy <- factor(pat_dat_tum$Met_status, labels = c(1,0))
pat_dat_tum$CBFb_stain_comp = relevel(factor(pat_dat_tum$CBFb_stain_comp), ref = "negative")



fit <- survfit(Surv(Overall_Survival, Status) ~ CBFb, data = pat_dat_tum)
p1 <- ggsurvplot(fit, pval=TRUE, pval.method = FALSE, conf.int = FALSE, risk.table = "absolute", 
                 legend.labs = c("CBFb neg", "CBFb pos")) + 
  ggtitle("Overall survival by CBFb expression")
fit <- survfit(Surv(Met_RFS, Status) ~ CBFb, data = pat_dat_tum)
p1 <- ggsurvplot(fit, pval=TRUE, pval.method = FALSE, conf.int = FALSE, risk.table = "absolute", 
                 legend.labs = c("No", "Yes")) + 
  ggtitle("Metastatic RFS by CBFb expression")
Comp_fit <- survfit(Surv(Overall_Survival, Status) ~ CBFb_stain_comp, data = pat_dat_tum)
p2 <- ggsurvplot(Comp_fit, pval=TRUE, pval.method = FALSE, conf.int = FALSE, risk.table = "absolute", 
                 legend.labs = c("Negative", "Cytoplasm", "Nuclear", "Combined")) + 
  ggtitle("Overall survival by CBFb compartment")

ppi=300
png(paste0(final_dir, "Overall_survival_byexpression.png"), height = 6*ppi, width = 7*ppi, res = ppi)
p1
dev.off()


fit <- coxph(Surv(Overall_Survival, Status) ~ Age + Subtype2 + CBFb, data = pat_dat_tum)
png(paste0(final_dir, "Overall_surv_forest.png"), height = 6*ppi, width = 6*ppi, res = ppi)
ggforest(fit, main = "Multivariate Cox Proportional Hazard Model",)
dev.off()

fit <- coxph(Surv(Overall_Survival, Status) ~ Age + Subtype2 + H_score, data = pat_dat_tum)

