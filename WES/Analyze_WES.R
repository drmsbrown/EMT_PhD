#################################
#Analyze mutations from WES
#Author: Meredith Brown
#4/20/20
#################################
rm(list=ls())
library("GenVisR")
library("vcfR")
library("GenomicRanges")
library("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library("MutationalPatterns")


dir_in <- "/Users/msb/mountpt/Raman_Shoyu/MSB/SUM149_WES/F20FTSUSAT0349_HUMwpkE/VCF/"

dir1 <- "~/Google Drive/Dartmouth/Raman Lab/WES/files/"
dir2 <- "~/Google Drive/Dartmouth/Raman Lab/WES/figures/"

indel_gr <- readRDS(paste0(dir1, "indel_granges.RDS"))
snp_gr <- readRDS(paste0(dir1, "snp_granges.RDS"))

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

#contingency test for SNP overlaps between clones
indel_contingency_test <- function(clone1, clone2, name1, name2){
  #clone1 <- indel_gr$sc4
  #clone2 <- indel_gr$sc1
  #name1 <- "sc4"
  #name2 <- "sc1"
  overlaps1 <- findOverlaps(clone1, indel_consensus)
  indicies1 <- subjectHits(overlaps1)
  overlaps2 <- findOverlaps(clone2, indel_consensus)
  indicies2 <- subjectHits(overlaps2)
  indel_consensus_df[,paste0(name1)] <- 0
  indel_consensus_df[,paste0(name1)][indicies1] <- 1
  indel_consensus_df[,paste0(name1)] <- factor(indel_consensus_df[,paste0(name1)], levels=c("1", "0"))
  indel_consensus_df[,paste0(name2)] <- 0
  indel_consensus_df[,paste0(name2)][indicies2] <- 1
  indel_consensus_df[,paste0(name2)] <- factor(indel_consensus_df[,paste0(name2)], levels=c("1", "0"))
  indel_table <- table(indel_consensus_df[,paste0(name1)], indel_consensus_df[,paste0(name2)])
  test <- fisher.test(indel_table)
  print(indel_table)
  print(test)
  mosaicplot(indel_table, xlab = name1, ylab = name2)
}
#SMALL p-value means that they share a lot of the same SNPs and are NOT randomly distributed = good thing
indel_contingency_test(indel_gr$sc4, indel_gr$sc1, "sc4", "sc1")
indel_contingency_test(indel_gr$sc4, indel_gr$sc3, "sc4", "sc3")
indel_contingency_test(indel_gr$sc4, indel_gr$sc6, "sc4", "sc6")
indel_contingency_test(indel_gr$sc4, indel_gr$sc8, "sc4", "sc8")
indel_contingency_test(indel_gr$sc4, indel_gr$sc14, "sc4", "sc14")

indel_contingency_test(indel_gr$sc1, indel_gr$sc14, "sc1", "P")


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

#contingency test for SNP overlaps between clones
snp_contingency_test <- function(clone1, clone2, name1, name2){
  #clone1 <- snp_gr$sc4
  #clone2 <- snp_gr$sc1
  #name1 <- "sc4"
  #name2 <- "sc1"
  overlaps1 <- findOverlaps(clone1, snp_consensus)
  indicies1 <- subjectHits(overlaps1)
  overlaps2 <- findOverlaps(clone2, snp_consensus)
  indicies2 <- subjectHits(overlaps2)
  snp_consensus_df[,paste0(name1)] <- 0
  snp_consensus_df[,paste0(name1)][indicies1] <- 1
  snp_consensus_df[,paste0(name2)] <- 0
  snp_consensus_df[,paste0(name2)][indicies2] <- 1
  SNP_table <- table(snp_consensus_df[,paste0(name1)], snp_consensus_df[,paste0(name2)])
  test <- fisher.test(SNP_table)
  print(SNP_table)
  print(test)
  mosaicplot(SNP_table)
}
#SMALL p-value means that they share a lot of the same SNPs and are NOT randomly distributed = good thing
snp_contingency_test(snp_gr$sc4, snp_gr$sc1, "sc4", "sc1")
snp_contingency_test(snp_gr$sc4, snp_gr$sc3, "sc4", "sc3")
snp_contingency_test(snp_gr$sc4, snp_gr$sc6, "sc4", "sc6")
snp_contingency_test(snp_gr$sc4, snp_gr$sc8, "sc4", "sc8")
snp_contingency_test(snp_gr$sc4, snp_gr$sc14, "sc4", "sc14")
snp_contingency_test(snp_gr$sc4, snp_gr$P, "sc4", "P")


#######################################################
#Explore the data

#which indels are in sc1 but not others?
#Fill indel_consensus data frame
overlaps1 <- findOverlaps(indel_gr$sc4, indel_consensus)
indicies1 <- subjectHits(overlaps1)
indel_consensus_df[,"sc4"] <- 0
indel_consensus_df[,"sc4"][indicies1] <- 1
print(indel_consensus_df)
overlaps1 <- findOverlaps(indel_gr$sc3, indel_consensus)
indicies1 <- subjectHits(overlaps1)
indel_consensus_df[,"sc3"] <- 0
indel_consensus_df[,"sc3"][indicies1] <- 1
print(indel_consensus_df)
overlaps1 <- findOverlaps(indel_gr$sc6, indel_consensus)
indicies1 <- subjectHits(overlaps1)
indel_consensus_df[,"sc6"] <- 0
indel_consensus_df[,"sc6"][indicies1] <- 1
print(indel_consensus_df)
overlaps1 <- findOverlaps(indel_gr$sc8, indel_consensus)
indicies1 <- subjectHits(overlaps1)
indel_consensus_df[,"sc8"] <- 0
indel_consensus_df[,"sc8"][indicies1] <- 1
print(indel_consensus_df)
overlaps1 <- findOverlaps(indel_gr$sc14, indel_consensus)
indicies1 <- subjectHits(overlaps1)
indel_consensus_df[,"sc14"] <- 0
indel_consensus_df[,"sc14"][indicies1] <- 1
print(indel_consensus_df)
overlaps1 <- findOverlaps(indel_gr$P, indel_consensus)
indicies1 <- subjectHits(overlaps1)
indel_consensus_df[,"P"] <- 0
indel_consensus_df[,"P"][indicies1] <- 1
print(indel_consensus_df)
overlaps1 <- findOverlaps(indel_gr$sc1, indel_consensus)
indicies1 <- subjectHits(overlaps1)
indel_consensus_df[,"sc1"] <- 0
indel_consensus_df[,"sc1"][indicies1] <- 1
print(indel_consensus_df)

write.table(indel_consensus_df, paste0(dir1, "indel_consensus_df.csv"), sep = ",")

#try to do this as a loop...?
clone_list <- c(indel_gr$sc4, indel_gr$sc1, indel_gr$sc3, indel_gr$sc6, indel_gr$sc8, indel_gr$sc14, indel_gr$P)
name_list <- c("sc4","sc1","sc3","sc6","sc8","sc14","P")
fill_indel_consensus <- function(clone_list, name_list){
  for (c in clone_list){
    overlaps <- findOverlaps(clone_list[,c+1], indel_consensus)
    indicies <- subjectHits(overlaps)
    print(overlaps)
    print(indicies)
  }
  for (i in name_list){
    indel_consensus_df[,paste0(name_list[i+1])] <- 0
    indel_consensus_df[,paste0(name_list[i+1])][indicies] <- 1
    print(indel_consensus_df)
  }
  print(indel_consensus_df)
}
fill_indel_consensus(clone_list[1:7], name_list[1:7])
for (i in clone_list[1:7]){
  overlaps <- findOverlaps(clone_list[,i+1], indel_consensus)
  indicies <- subjectHits(overlaps)
  print(overlaps)
  print(indicies)
}
for (i in name_list){
  indel_consensus_df[,paste0(name_list[i+1])] <- 0
  indel_consensus_df[,paste0(name_list[i+1])][indicies] <- 1
  print(indel_consensus_df)
}
