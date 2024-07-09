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
dir1 <- "/Users/yourname/example_dir/files_dir/"
#output: for figures 
dir2 <- "/Users/yourname/example_dir/output_dir/"

#Get the files and filenames for TMAs
#This will recognize all cell segmentation txt files from the inForm output
filepath <- list.files(path = dir1, pattern = "cell_seg_data.txt")
file_names <- c()


###################################################################### IF THIS IS A TMA ########################################################################

#Strips the TMA output file name to get the TMA coordinates (i.e. A1)
#You will need the TMA number or identifier, and the coordinates for identification. TMA name will depend on your inForm file name but the coordinates 
#are always the 2nd and 3rd spaces in brackets i.e. Core[1,11,A] = image A11
#e.g. filename: "1050 #4_Core[1,11,A]_[6762,61721]_cell_seg_data.txt"
#               TMA = 1050 #4
#               TMA_row = A11
#               TMA_column = A
for (i in 1:length(filepath)) {
  TMA <- strsplit(filepath, "_Core")[[i]][1]
  TMA_row <- strsplit(filepath, ",")[[i]][2]
  TMA_column <- strsplit(filepath, "]")[[i]][1]
  TMA_column <- strsplit(TMA_column,",")[[1]][3]
  file_names <- as.data.frame(rbind(file_names, cbind(TMA, TMA_row, TMA_column, ID = paste0(TMA_row, "_", TMA_column),  group = NA, sample = NA)))
}
colnames(file_names) <- c("TMA","row","column", "ID", "group", "sample")

write.csv(file_names, paste0(dir1, "File_names.csv"))

######################### Go into your file directory and annotate "File_names.csv" with the correct sample names and save ##################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###################################################################### IF THIS IS NOT A TMA ########################################################################
for (i in 1:length(filepath)) {
  tumor <- strsplit(filepath, "_Core")[[i]][1]
  coordinates <- strsplit(filepath, "_")[[i]][2]
  file_names <- as.data.frame(rbind(file_names, cbind(TMA, tumor, coordinates, group = NA, sample = NA)))
}
colnames(file_names) <- c("tumor","coordinates", "group", "sample")

write.csv(file_names, paste0(dir1, "File_names.csv"))

######################### Go into your file directory and annotate "File_names.csv" with the correct sample names and save ##################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Read in your annotated file names with all the samples labeled
file_names <- read.csv(paste0(dir1, "File_names.csv"),stringsAsFactors = F)


######################### Check mean counts to confirm good overall staining. Can be done BEFORE phenotyping  ##########################################
TSA_mean <- c()
for (i in 1:length(filepath)){
  dat <- read.delim(paste0(dir1, filepath[i]), na.strings = "#N/A")
  dat <- dat[which(dat$Tissue.Category=="Tumor"),]
  dat_sub <- as.data.frame(cbind(Sample = dat$Sample.Name, Tissue = dat$Tissue.Category, Cell = dat$Cell.ID, Phenotype = dat$Phenotype,
                                 Cell_X_pos = dat$Cell.X.Position, Cell_Y_pos = dat$Cell.Y.Position))
  dat_sub <- cbind(dat_sub, OPAL520 = dat$Entire.Cell.Opal.520.Mean..Normalized.Counts..Total.Weighting.,
                   OPAL540 = dat$Entire.Cell.Opal.540.Mean..Normalized.Counts..Total.Weighting.,
                   OPAL570 = dat$Entire.Cell.Opal.570.Mean..Normalized.Counts..Total.Weighting., 
                   OPAL620 = dat$Entire.Cell.Opal.620.Mean..Normalized.Counts..Total.Weighting.,
                   OPAL650 = dat$Entire.Cell.Opal.650.Mean..Normalized.Counts..Total.Weighting.,
                   OPAL690 = dat$Entire.Cell.Opal.690.Mean..Normalized.Counts..Total.Weighting.,
                   DAPI =    dat$Entire.Cell.DAPI.Mean..Normalized.Counts..Total.Weighting.)
  dat_sub[is.na(dat_sub)] <-  0
  dat_sub$sample <- rep(file_names$sample[i], nrow(dat_sub))
  dat_sub$group <- rep(file_names$group[i], nrow(dat_sub))
  TSA_mean <- rbind(TSA_mean, dat_sub)
}

Plot_mean_counts <- function(TSA, dir, png, title){
  OPAL520 <- as.data.frame(TSA_mean$OPAL520)
  names(OPAL520) <- c("count")
  OPAL520 <- cbind(OPAL520,  batch= TSA_mean$ID, OPAL = rep('K14', nrow(TSA_mean)))
  OPAL540 <- as.data.frame(TSA_mean$OPAL540)
  names(OPAL540) <- c("count")
  OPAL540 <- cbind(OPAL540, batch= TSA_mean$ID, OPAL = rep('K8', nrow(TSA_mean)))
  OPAL570 <- as.data.frame(TSA_mean$OPAL570)
  names(OPAL570) <- c("count")
  OPAL570 <- cbind(OPAL570, batch= TSA_mean$ID, OPAL = rep('ZEB1', nrow(TSA_mean)))
  OPAL620 <- as.data.frame(TSA_mean$OPAL620)
  names(OPAL620) <- c("count")
  OPAL620 <- cbind(OPAL620, batch= TSA_mean$ID, OPAL = rep('Snail', nrow(TSA_mean)))
  OPAL650 <- as.data.frame(TSA_mean$OPAL650)
  names(OPAL650) <- c("count")
  OPAL650 <- cbind(OPAL650, batch= TSA_mean$ID, OPAL = rep('E-cadherin', nrow(TSA_mean)))
  OPAL690 <- as.data.frame(TSA_mean$OPAL690)
  names(OPAL690) <- c("count")
  OPAL690 <- cbind(OPAL690, batch= TSA_mean$ID, OPAL = rep('Vimentin', nrow(TSA_mean)))
  
  TSA_mean_rcfg_raw <- rbind(OPAL520, OPAL540, OPAL570, OPAL620, OPAL650, OPAL690)
  
  ppi=300
  p <- ggplot(TSA_mean_rcfg_raw, aes(x=OPAL, y = count)) + 
    geom_boxplot() + ggtitle(title)
  png(paste0(dir2, png), height = 7*ppi, width = 8*ppi, res = ppi )
  p
  dev.off()
  
  print(p)
}

#output a boxplot to check the raw mean counts of each OPAL. Can be found in your figures directory
Plot_mean_counts(TSA_mean, dir2, "CHTN_mean_counts.png", "CHTN mean counts")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################### Import and extract OPAL raw counts for analysis ##################################################

#First check your Column Names for Entire Cell and ensure they match for the input below
test <- read.delim(paste0(dir1, filepath[1]), na.strings = "#N/A")
#these may say something like Entire.Cell.Opal.520..Total or Entire.Cell.K8.Opal.520..Total. #### Periods matter! ###

#import, read, and compile the total marker data from Entire Cell. 
#Check that these inputs match your column names. ALSO check for Tissue Category name - may be "Tumor" or "tumor" 
OPAL_total <- c()
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
  dat_sub$group <- rep(file_names$group[i], nrow(dat_sub))
  dat_sub$sample <- rep(file_names$sample[i], nrow(dat_sub))
  OPAL_total <- rbind(OPAL_total, dat_sub)
}

#write as a CSV file 
write.csv(OPAL_total, paste0(dir1, "Raw_entirecell.csv"))

#add count column
OPAL_total$count <- rep(1, nrow(OPAL_total))
#remove cells that have not been phenotyped
no_pheno <- which(OPAL_total$Phenotype=="")
OPAL_total_rm <- OPAL_total[-no_pheno,]

################## Check your phenotype names ##############
unique(OPAL_total_rm$Phenotype)

#rank phenotypes from epithelial to mesenchymal. "levels" must match your CURRENT phenotype names and be in E-M order. 
# You can re-name phenotypes in "labels" section. Just be sure to keep the order
OPAL_total_rm$Phenotype <- factor(OPAL_total_rm$Phenotype, 
                                  levels = c("Ecad only","K8+ecad", "K14",  "Trip+", "K8+vim","Snail+", "Vim only", "Vim+ZEB1"), 
                                  labels = c("E-cad only", "K8 & E-cad", "KRT", "Triple+", "KRT & Vim","Snail", "Vim only","Vim & ZEB1"))
#reorder tumors or samples 
OPAL_total_rm$group <- factor(OPAL_total_rm$group, 
                             levels = c("group 1", "group 2", "group 3"))

#plot phenotypes by sample or group
ppi=300
png(paste0(dir2, "Phenotypes_per_tumor.png"), height = 7*ppi, width = 6*ppi, res = ppi)
ggplot(OPAL_total_rm, aes(fill=Phenotype, y =count, x=group)) + 
  geom_bar(position="fill", stat="identity") + ggtitle("Patient Tumor Phenotypes") + ylab("Percent of cells")
dev.off()

######################################### Tidy data and normalize markers to a percentile 0-1 for violin plots #######################################

#Compile markers for normalization
OPAL520 <- as.data.frame(percentize(OPAL_total_rm$OPAL520))
names(OPAL520) <- c("count")
OPAL520 <- cbind(OPAL520,  group= OPAL_total_rm$group, sample = OPAL_total_rm$sample, OPAL = rep('K14', nrow(OPAL_total_rm)), phenotype = OPAL_total_rm$Phenotype)
OPAL540 <- as.data.frame(percentize(OPAL_total_rm$OPAL540))
names(OPAL540) <- c("count")
OPAL540 <- cbind(OPAL540, group = OPAL_total_rm$group, sample = OPAL_total_rm$sample, OPAL = rep('K8', nrow(OPAL_total_rm)), phenotype = OPAL_total_rm$Phenotype)
OPAL570 <- as.data.frame(percentize(OPAL_total_rm$OPAL570))
names(OPAL570) <- c("count")
OPAL570 <- cbind(OPAL570, group = OPAL_total_rm$group, sample = OPAL_total_rm$sample, OPAL = rep('ZEB1', nrow(OPAL_total_rm)), phenotype = OPAL_total_rm$Phenotype)
OPAL620 <- as.data.frame(percentize(OPAL_total_rm$OPAL620))
names(OPAL620) <- c("count")
OPAL620 <- cbind(OPAL620, group = OPAL_total_rm$group, sample = OPAL_total_rm$sample, OPAL = rep('Snail', nrow(OPAL_total_rm)), phenotype = OPAL_total_rm$Phenotype)
OPAL650 <- as.data.frame(percentize(OPAL_total_rm$OPAL650))
names(OPAL650) <- c("count")
OPAL650 <- cbind(OPAL650, group = OPAL_total_rm$group, sample = OPAL_total_rm$sample, OPAL = rep('E-cadherin', nrow(OPAL_total_rm)), phenotype = OPAL_total_rm$Phenotype)
OPAL690 <- as.data.frame(percentize(OPAL_total_rm$OPAL690))
names(OPAL690) <- c("count")
OPAL690 <- cbind(OPAL690, group = OPAL_total_rm$group, sample = OPAL_total_rm$sample, OPAL = rep('Vimentin', nrow(OPAL_total_rm)), phenotype = OPAL_total_rm$Phenotype)

OPAL_total_rcfg <- rbind(OPAL520, OPAL540, OPAL570, OPAL650, OPAL690)
OPAL_total_rcfg$OPAL <- factor(OPAL_total_rcfg$OPAL, 
                              levels = c("E-cadherin", "K8", "K14", "Snail", "Vimentin", "ZEB1"))

################################ Check that your phenotypes have been assigned correctly i.e. "Ecad only" has high ecad and nothing else ##############

ppi=300
png(paste0(dir2,"TMA_violin_OPAL_by_pheno.png"), height = 5*ppi, width = 6*ppi, res = ppi)
ggplot(OPAL_total_rcfg, aes(x= OPAL, y=count, color = OPAL)) + 
  geom_violin(trim = F) + geom_boxplot(width=0.1, color = "grey") + 
  facet_wrap(~phenotype) + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
dev.off()

#see EMT&Het code for analysis of EMT and Heterogeneity scores
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


