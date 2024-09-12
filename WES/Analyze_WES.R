#################################
#Analyze mutations from WES
#Author: msbdith Brown
#4/20/20
#################################
rm(list=ls())
library("GenVisR")
library("vcfR")
library(vcd)
library("GenomicRanges")
library("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library("MutationalPatterns")
library(ggplot2)
library(RColorBrewer)


dir_in <- "/Users/msb/mountpt/Raman_Shoyu/MSB/SUM149_WES/F20FTSUSAT0349_HUMwpkE/VCF/"

dir1 <- "~/Google Drive/Dartmouth/Raman Lab/WES/files/"
dir2 <- "~/Google Drive/Dartmouth/Raman Lab/WES/figures/"

indel_gr <- readRDS(paste0(dir1, "indel_granges.RDS"))
snp_gr <- readRDS(paste0(dir1, "snp_granges.RDS"))


colsEMT <- brewer.pal(7, "RdYlBu")
colsPar <- brewer.pal(10, "PiYG")
EMT_col <- c(colsEMT[7], colsEMT[6], colsEMT[5], colsEMT[3], colsEMT[2], colsEMT[1], colsPar[2])

################################
#Test the random distribution of SNP or InDels 
#between clones to determine if they have similar or different genetic backgrounds
##############################
#convert  INDEL files to GRanges objects
#get files for indel csv
get_filepath <- function(x){
  paste0("/Users/msb/mountpt/Raman_Shoyu/MSB/SUM149_WES/F20FTSUSAT0349_HUMwpkE/",
         x, "/result_variation/indel/", x, ".indel.annot.csv")
}
indel_filepath <- c(get_filepath("149P"), get_filepath("sc1"), get_filepath("sc3"), 
                    get_filepath("sc4"), get_filepath("sc6"), 
                    get_filepath("sc8"), get_filepath("sc14"))

#convert csv to data frame for granges object
indel_df <- lapply(indel_filepath, function(x) as.data.frame(read.csv(x)))
#generate granges objects
indel_gr <- lapply(indel_df, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
names(indel_gr) <- c("P","sc1","sc3","sc4","sc6","sc8","sc14")
saveRDS(indel_gr, paste0(dir1, "indel_granges.RDS"))

#create consensus INDEL Granges object
indel_consensus <- union(indel_gr$P, indel_gr$sc1)
indel_consensus <- union(indel_consensus, indel_gr$sc3)
indel_consensus <- union(indel_consensus, indel_gr$sc4)
indel_consensus <- union(indel_consensus, indel_gr$sc6)
indel_consensus <- union(indel_consensus, indel_gr$sc8)
indel_consensus <- union(indel_consensus, indel_gr$sc14)
indel_consensus_df <- as.data.frame(indel_consensus)

#find overlaps and annotate indel_consensus file
overlaps <- findOverlaps(indel_gr$sc4, indel_consensus)
indicies <- subjectHits(overlaps)
indel_consensus_df[,"sc4"] <- 0
indel_consensus_df[,"sc4"][indicies] <- 1
indel_consensus_df[,"sc4"] <- factor(indel_consensus_df[,"sc4"], levels=c("1", "0"))
overlaps <- findOverlaps(indel_gr$sc6, indel_consensus)
indicies <- subjectHits(overlaps)
indel_consensus_df[,"sc6"] <- 0
indel_consensus_df[,"sc6"][indicies] <- 1
indel_consensus_df[,"sc6"] <- factor(indel_consensus_df[,"sc6"], levels=c("1", "0"))
overlaps <- findOverlaps(indel_gr$sc3, indel_consensus)
indicies <- subjectHits(overlaps)
indel_consensus_df[,"sc3"] <- 0
indel_consensus_df[,"sc3"][indicies] <- 1
indel_consensus_df[,"sc3"] <- factor(indel_consensus_df[,"sc3"], levels=c("1", "0"))
overlaps <- findOverlaps(indel_gr$sc8, indel_consensus)
indicies <- subjectHits(overlaps)
indel_consensus_df[,"sc8"] <- 0
indel_consensus_df[,"sc8"][indicies] <- 1
indel_consensus_df[,"sc8"] <- factor(indel_consensus_df[,"sc8"], levels=c("1", "0"))
overlaps <- findOverlaps(indel_gr$sc14, indel_consensus)
indicies <- subjectHits(overlaps)
indel_consensus_df[,"sc14"] <- 0
indel_consensus_df[,"sc14"][indicies] <- 1
indel_consensus_df[,"sc14"] <- factor(indel_consensus_df[,"sc14"], levels=c("1", "0"))
overlaps <- findOverlaps(indel_gr$sc1, indel_consensus)
indicies <- subjectHits(overlaps)
indel_consensus_df[,"sc1"] <- 0
indel_consensus_df[,"sc1"][indicies] <- 1
indel_consensus_df[,"sc1"] <- factor(indel_consensus_df[,"sc1"], levels=c("1", "0"))
overlaps <- findOverlaps(indel_gr$P, indel_consensus)
indicies <- subjectHits(overlaps)
indel_consensus_df[,"P"] <- 0
indel_consensus_df[,"P"][indicies] <- 1
indel_consensus_df[,"P"] <- factor(indel_consensus_df[,"P"], levels=c("1", "0"))

#write.csv(indel_consensus_df, file = paste0(dir1, "INDEL_consensus_overlaps.csv"))
indel_consensus_df <- read.csv(file = paste0(dir1, "INDEL_consensus_overlaps.csv"), stringsAsFactors = F)

#contingency test for SNP overlaps between clones
indel_contingency_test <- function(clone1, clone2, name1, name2, color){
  indel_table <- table(indel_consensus_df[,paste0(clone1)], indel_consensus_df[,paste0(clone2)])
  colnames(indel_table) <- c(paste("in", name1), paste("not in", name1))
  rownames(indel_table) <- c(paste("in", name2), paste("not in", name2))
  test <- fisher.test(indel_table)
  print(indel_table)
  print(test)
  mosaicplot(indel_table, xlab = name2, ylab = name1, main = paste(name1, "vs", name2, "Indels"), 
             color = c(EMT_col[1], color), type = c("deviance"), cex.axis = test$p.value)
  ppi=300
  png(paste0(dir2, paste(name1, "vs", name2, "Indels"), ".png"), width=5*ppi, height=5*ppi, res=ppi)
  mosaicplot(indel_table, xlab = name2, ylab = name1, main = paste(name1, "vs", name2, "Indels"), 
             color = c(EMT_col[1], color), type = c("deviance"))
  dev.off()
}

#SMALL p-value means that they share a lot of the same SNPs and are NOT randomly distributed = good thing
indel_contingency_test("sc4", "sc1", "E", "M2", EMT_col[6])
indel_contingency_test("sc4", "sc3", "E", "EM2", EMT_col[3])
indel_contingency_test("sc4", "sc6", "E", "EM1", EMT_col[2])
indel_contingency_test("sc4", "sc8", "E", "EM3", EMT_col[4])
indel_contingency_test("sc4", "sc14", "E", "M1", EMT_col[5])
indel_contingency_test("sc4", "P", "E", "Parental", EMT_col[7])

#plot odds ratios
odds_labels <- c("E vs EM1", "E vs EM2", "E vs EM3", "E vs M1", "E vs M2", "E vs P")
indel_odds_ratios <- data.frame(yAxis = length(odds_labels):1,boxOdds = c(23.64283, 22.142, 17.0957, 20.61282, 0.6612739, 22.39046),
                          boxCILow = c(22.35531, 20.9403, 16.19622, 19.49573, 0.6336216, 21.19359),
                          boxCIHigh = c(25.01407, 23.3987, 18.04925, 21.76195, 0.6900567, 23.67352))
ppi=300
png(paste0(dir2, "INDEL_oddsratio.png"), width=5*ppi, height=7*ppi, res=ppi)
p <- ggplot(indel_odds_ratios, aes(x = boxOdds, y = yAxis))
p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = length(odds_labels):1, labels = odds_labels) +
  ylab("") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="F.test p < 2.2e-16", size = 3.5, hjust = 0) +
  xlab("Odds ratio") + ggtitle("Overlap of Indels in Clonal Populations \n\ Compared to Epithelial")
dev.off()



#########################################################
#convert  SNP files to GRanges objects
#get files for snp csv
get_filepath <- function(x){
  paste0("/Users/msb/mountpt/Raman_Shoyu/MSB/SUM149_WES/F20FTSUSAT0349_HUMwpkE/",
         x, "/result_variation/snp/", x, ".snp.annot.csv")
}
snp_filepath <- c(get_filepath("149P"), get_filepath("sc1"), get_filepath("sc3"), 
                    get_filepath("sc4"), get_filepath("sc6"), 
                    get_filepath("sc8"), get_filepath("sc14"))
#convert csv to data frame for granges object
snp_df <- lapply(snp_filepath, function(x) as.data.frame(read.csv(x)))
#generate granges objects
snp_gr <- lapply(snp_df, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
names(snp_gr) <- c("P","sc1","sc3","sc4","sc6","sc8","sc14")
saveRDS(snp_gr, paste0(dir1, "snp_granges.RDS"))

#create consensus snp Granges object for analysis
snp_consensus <- union(snp_gr$P, snp_gr$sc1)
snp_consensus <- union(snp_consensus, snp_gr$sc3)
snp_consensus <- union(snp_consensus, snp_gr$sc4)
snp_consensus <- union(snp_consensus, snp_gr$sc6)
snp_consensus <- union(snp_consensus, snp_gr$sc8)
snp_consensus <- union(snp_consensus, snp_gr$sc14)
snp_consensus_df <- as.data.frame(snp_consensus)

overlaps <- findOverlaps(snp_gr$sc4, snp_consensus)
indicies <- subjectHits(overlaps)
snp_consensus_df[,"sc4"] <- 0
snp_consensus_df[,"sc4"][indicies] <- 1
snp_consensus_df[,"sc4"] <- factor(snp_consensus_df[,"sc4"], levels=c("1", "0"))
overlaps <- findOverlaps(snp_gr$sc6, snp_consensus)
indicies <- subjectHits(overlaps)
snp_consensus_df[,"sc6"] <- 0
snp_consensus_df[,"sc6"][indicies] <- 1
snp_consensus_df[,"sc6"] <- factor(snp_consensus_df[,"sc6"], levels=c("1", "0"))
overlaps <- findOverlaps(snp_gr$sc3, snp_consensus)
indicies <- subjectHits(overlaps)
snp_consensus_df[,"sc3"] <- 0
snp_consensus_df[,"sc3"][indicies] <- 1
snp_consensus_df[,"sc3"] <- factor(snp_consensus_df[,"sc3"], levels=c("1", "0"))
overlaps <- findOverlaps(snp_gr$sc8, snp_consensus)
indicies <- subjectHits(overlaps)
snp_consensus_df[,"sc8"] <- 0
snp_consensus_df[,"sc8"][indicies] <- 1
snp_consensus_df[,"sc8"] <- factor(snp_consensus_df[,"sc8"], levels=c("1", "0"))
overlaps <- findOverlaps(snp_gr$sc14, snp_consensus)
indicies <- subjectHits(overlaps)
snp_consensus_df[,"sc14"] <- 0
snp_consensus_df[,"sc14"][indicies] <- 1
snp_consensus_df[,"sc14"] <- factor(snp_consensus_df[,"sc14"], levels=c("1", "0"))
overlaps <- findOverlaps(snp_gr$sc1, snp_consensus)
indicies <- subjectHits(overlaps)
snp_consensus_df[,"sc1"] <- 0
snp_consensus_df[,"sc1"][indicies] <- 1
snp_consensus_df[,"sc1"] <- factor(snp_consensus_df[,"sc1"], levels=c("1", "0"))
overlaps <- findOverlaps(snp_gr$P, snp_consensus)
indicies <- subjectHits(overlaps)
snp_consensus_df[,"P"] <- 0
snp_consensus_df[,"P"][indicies] <- 1
snp_consensus_df[,"P"] <- factor(snp_consensus_df[,"P"], levels=c("1", "0"))

write.csv(snp_consensus_df, file = paste0(dir1, "SNP_consensus_overlaps.csv"))
snp_consensus_df <- read.csv(file = paste0(dir1, "SNP_consensus_overlaps.csv"), stringsAsFactors = F)
#contingency test for SNP overlaps between clones
snp_contingency_test <- function(clone1, clone2, name1, name2, color){
  SNP_table <- table(snp_consensus_df[,paste0(clone1)], snp_consensus_df[,paste0(clone2)])
  test <- fisher.test(SNP_table)
  colnames(SNP_table) <- c(paste("in", name1), paste("not in", name1))
  rownames(SNP_table) <- c(paste("in", name2), paste("not in", name2))
  print(SNP_table)
  print(test)
  mosaicplot(SNP_table)
  mosaicplot(SNP_table, xlab = name2, ylab = name1, main = paste(name1, "vs", name2, "SNPs"), 
             color = c(EMT_col[1], color), type = c("deviance"))
  ppi=300
  png(paste0(dir2, paste(name1, "vs", name2, "SNPs"), ".png"), width=5*ppi, height=5*ppi, res=ppi)
  mosaicplot(SNP_table, xlab = name2, ylab = name1, main = paste(name1, "vs", name2, "SNPs"), 
             color = c(EMT_col[1], color), type = c("deviance"))
  dev.off()
}
#SMALL p-value means that they share a lot of the same SNPs and are NOT randomly distributed = good thing
snp_contingency_test("sc4", "sc6", "E", "EM1", EMT_col[2])
snp_contingency_test("sc4", "sc3", "E", "EM2", EMT_col[3])
snp_contingency_test("sc4", "sc8", "E", "EM3", EMT_col[4])
snp_contingency_test("sc4", "sc14", "E", "M1", EMT_col[5])
snp_contingency_test("sc4", "sc1", "E", "M2", EMT_col[6])
snp_contingency_test("sc4", "P", "E", "Parental", EMT_col[7])

odds_labels <- c("E vs EM1", "E vs EM2", "E vs EM3", "E vs M1", "E vs M2", "E vs P")
SNP_odds_ratios <- data.frame(yAxis = length(odds_labels):1,boxOdds = c(136.847, 101.5619, 117.0523, 101.5692, 0.2265217, 139.6725),
                                boxCILow = c(132.5325, 98.06287,112.1309, 98.22278, 0.2216535, 135.1296),
                                boxCIHigh = c(141.5104, 104.71845, 121.2198, 105.40130, 0.2314949, 145.5372))
ppi=300
png(paste0(dir2, "SNP_oddsratio.png"), width=5*ppi, height=7*ppi, res=ppi)
p <- ggplot(SNP_odds_ratios, aes(x = boxOdds, y = yAxis))
p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = length(odds_labels):1, labels = odds_labels) +
  ylab("") +
  annotate(geom = "text", y =1.1, x = 3.5, label ="F.test p < 2.2e-16", size = 3.5, hjust = 0) +
  xlab("Odds ratio") + ggtitle("Overlap of SNPs in Clonal Populations \n\ Compared to Epithelial")
dev.off()
#######################################################
#Explore the data
indel_consensus_df <- read.csv(paste0(dir1, "INDEL_consensus_overlaps.csv"), stringsAsFactors = F)
snp_consensus_df <- read.csv(paste0(dir1, "SNP_consensus_overlaps.csv"), stringsAsFactors = F)

library(GenomicFeatures)
hg19.refseq.db <- makeTxDbFromUCSC(genome="hg19", table="refGene")
refseq.genes<- genes(hg19.refseq.db)
gns = granges(Homo.sapiens, column="SYMBOL")

#which indels are in sc1 but not others?
sc1_indels <- indel_consensus_df[which(indel_consensus_df$sc1=="1"),]
sc4_no <- sc1_indels[which(sc1_indels$sc4=="0"),]
sc4_sc6_no <- sc4_no[which(sc4_no$sc6=="0"),]
sc4_sc6_sc8_no <- sc4_sc6_no[which(sc4_sc6_no$sc8=="0"),]
sc4_sc6_sc8_sc14_no <- sc4_sc6_sc8_no[which(sc4_sc6_sc8_no$sc14=="0"),]
only_sc1_indels <- sc4_sc6_sc8_sc14_no[which(sc4_sc6_sc8_sc14_no$P=="0"),]
#convert to Granges
sc1_indels_gr <- makeGRangesFromDataFrame(only_sc1_indels, keep.extra.columns = F)

#which indels are in sc1 but not others?
sc1_snps <- snp_consensus_df[which(snp_consensus_df$sc1=="1"),]
sc4_no <- sc1_snps[which(sc1_snps$sc4=="0"),]
sc4_sc6_no <- sc4_no[which(sc4_no$sc6=="0"),]
sc4_sc6_sc8_no <- sc4_sc6_no[which(sc4_sc6_no$sc8=="0"),]
sc4_sc6_sc8_sc14_no <- sc4_sc6_sc8_no[which(sc4_sc6_sc8_no$sc14=="0"),]
only_sc1_snps <- sc4_sc6_sc8_sc14_no[which(sc4_sc6_sc8_sc14_no$P=="0"),]
#convert to Granges
sc1_snps_gr <- makeGRangesFromDataFrame(only_sc1_snps, keep.extra.columns = F)

#annotate Granges objects to generate gene list of sc1 unique indels and snps





