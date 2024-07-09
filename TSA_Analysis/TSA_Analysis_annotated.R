#############################################################################
#Validate Phenotyping and Analyze TSA Output
#Author: Meredith Brown
#7/21/2021
#############################################################################

#Install and load in necessary libraries
#install.packages("ggplot2")
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(heatmaply)
rm(list = ls())

#set input and output directories
#input: Where are the files?
dir1 <- "/Users/yourname/example_dir/files_dir"
#output: for figures 
dir2 <- "/Users/yourname/example_dir/output_dir/"

#Get the files and filenames for TMAs
filepath <- list.files(path = dir1, pattern = "cell_seg_data.txt")
file_names <- c()
#Strips the TMA output file name to get the TMA coordinates (i.e. A1)
for (i in 1:length(filepath)) {
  TMA <- strsplit(filepath, "1050")[[i]][2]
  TMA <- strsplit(TMA, " ")[[1]][1]
  tumor <- strsplit(filepath, ",")[[i]][2]
  tumor <- paste0(TMA, tumor)
  image <- strsplit(filepath, "]")[[i]][1]
  image <- strsplit(image,",")[[1]][3]
  file_names <- as.data.frame(rbind(file_names, cbind(tumor, image, paste0(tumor, "_", image))))
}

colnames(file_names) <- c("tumor","image","ID")

#Annotate for tumor ID. Must be in order of file_names. Replace number (6) with number of replicates
print(file_names)
file_names$ID <- c(rep("sample1", 6), rep("sample2", 6), rep("sample3", 6), rep("sample4", 6), rep("sample5", 6), rep("sample6", 6))

#import, read, and compile the Mean marker data from Entire Cell
TMA_pheno <- c()
for (i in 1:length(filepath)){
  dat <- read.delim(paste0(dir1, filepath[i]), na.strings = "#N/A")
  dat <- dat[which(dat$Tissue.Category=="Tumor"),]
  dat_sub <- as.data.frame(cbind(Sample = dat$Sample.Name, Tissue = dat$Tissue.Category, Cell = dat$Cell.ID, Phenotype = dat$Phenotype,
                                 Cell_X_pos = dat$Cell.X.Position, Cell_Y_pos = dat$Cell.Y.Position))
  dat_sub <- cbind(dat_sub, OPAL520 = dat$Entire.Cell.K14..Opal.520..Total..Normalized.Counts..Total.Weighting.,
                   OPAL540 = dat$Entire.Cell.K8..Opal.540..Total..Normalized.Counts..Total.Weighting.,
                   OPAL570 = dat$Entire.Cell.ZEB1..Opal.570..Total..Normalized.Counts..Total.Weighting., 
                   OPAL620 = dat$Entire.Cell.Snail..Opal.620..Total..Normalized.Counts..Total.Weighting.,
                   OPAL650 = dat$Entire.Cell.Ecad..Opal.650..Total..Normalized.Counts..Total.Weighting.,
                   OPAL690 = dat$Entire.Cell.Vim..Opal.690..Total..Normalized.Counts..Total.Weighting.,
                   DAPI = dat$Entire.Cell.DAPI..DAPI..Total..Normalized.Counts..Total.Weighting.)
  dat_sub[is.na(dat_sub)] <-  0
  dat_sub$tumor <- rep(file_names$tumor[i], nrow(dat_sub))
  dat_sub$image <- rep(file_names$image[i], nrow(dat_sub))
  dat_sub$ID <- rep(file_names$ID[i], nrow(dat_sub))
  TMA_pheno <- rbind(TMA_pheno, dat_sub)
}

#write as a CSV file
write.csv(TMA_pheno, paste0(dir1, "TMA_evo_entirecell.csv"))

#add count column
TMA_pheno$count <- rep(1, nrow(TMA_pheno))
#remove cells that have not been phenotyped
no_pheno <- which(TMA_pheno$Phenotype=="")
TMA_pheno_rm <- TMA_pheno[-no_pheno,]

#rank phenotypes from epithelial to mesenchymal
TMA_pheno_rm$Phenotype <- factor(TMA_pheno_rm$Phenotype, 
                                 levels = c("Ecad only","Ecad+K8", "K14",  "Trip", "Snail", "Vimentin", "Vim+Zeb1"), 
                                 labels = c("E-cad only", "K8 & E-cad", "K14", "Triple+","Snail", "Vim only","Vim & ZEB1"))
#reorder TMA tumors
TMA_pheno_rm$tumor <- factor(TMA_pheno_rm$tumor, 
                             levels = c("A1", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A11",
                                        "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10"))
#plot phenotypes by sample
ppi=300
png(paste0(dir2, "Phenotypes_per_tumor.png"), height = 7*ppi, width = 6*ppi, res = ppi)
ggplot(TMA_pheno_rm, aes(fill=Phenotype, y =count, x=tumor)) + 
  geom_bar(position="fill", stat="identity") + ggtitle("Patient Tumor Phenotypes") + ylab("Percent of cells")
dev.off()

#Compile markers for normalization
OPAL520 <- as.data.frame(percentize(TMA_pheno_rm$OPAL520))
names(OPAL520) <- c("count")
OPAL520 <- cbind(OPAL520,  tumor= TMA_pheno_rm$tumor, image = TMA_pheno_rm$image, OPAL = rep('K14', nrow(TMA_pheno_rm)), phenotype = TMA_pheno_rm$Phenotype)
OPAL540 <- as.data.frame(percentize(TMA_pheno_rm$OPAL540))
names(OPAL540) <- c("count")
OPAL540 <- cbind(OPAL540, tumor = TMA_pheno_rm$tumor, image = TMA_pheno_rm$image, OPAL = rep('K8', nrow(TMA_pheno_rm)), phenotype = TMA_pheno_rm$Phenotype)
OPAL570 <- as.data.frame(percentize(TMA_pheno_rm$OPAL570))
names(OPAL570) <- c("count")
OPAL570 <- cbind(OPAL570, tumor = TMA_pheno_rm$tumor, image = TMA_pheno_rm$image, OPAL = rep('ZEB1', nrow(TMA_pheno_rm)), phenotype = TMA_pheno_rm$Phenotype)
OPAL620 <- as.data.frame(percentize(TMA_pheno_rm$OPAL620))
names(OPAL620) <- c("count")
OPAL620 <- cbind(OPAL620, tumor = TMA_pheno_rm$tumor, image = TMA_pheno_rm$image, OPAL = rep('Snail', nrow(TMA_pheno_rm)), phenotype = TMA_pheno_rm$Phenotype)
OPAL650 <- as.data.frame(percentize(TMA_pheno_rm$OPAL650))
names(OPAL650) <- c("count")
OPAL650 <- cbind(OPAL650, tumor = TMA_pheno_rm$tumor, image = TMA_pheno_rm$image, OPAL = rep('E-cadherin', nrow(TMA_pheno_rm)), phenotype = TMA_pheno_rm$Phenotype)
OPAL690 <- as.data.frame(percentize(TMA_pheno_rm$OPAL690))
names(OPAL690) <- c("count")
OPAL690 <- cbind(OPAL690, tumor = TMA_pheno_rm$tumor, image = TMA_pheno_rm$image, OPAL = rep('Vimentin', nrow(TMA_pheno_rm)), phenotype = TMA_pheno_rm$Phenotype)

TMA_pheno_rcfg <- rbind(OPAL520, OPAL540, OPAL570, OPAL650, OPAL690)
TMA_pheno_rcfg$OPAL <- factor(TMA_pheno_rcfg$OPAL, 
                              levels = c("E-cadherin", "K8", "K14", "Snail", "Vimentin", "ZEB1"))

#Check for correct phenotyping
ppi=300
png(paste0(dir2,"TMA_violin_OPAL_by_pheno.png"), height = 5*ppi, width = 6*ppi, res = ppi)
ggplot(TMA_pheno_rcfg, aes(x= OPAL, y=count, color = OPAL)) + 
  geom_violin(trim = F) + geom_boxplot(width=0.1, color = "grey") + 
  facet_wrap(~phenotype) + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
dev.off()


