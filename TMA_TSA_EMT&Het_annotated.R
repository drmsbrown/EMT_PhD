#############################################################################
#Analyze and Plot EMT and Heterogeneity Scores
#Author: Meredith Brown
#7/21/2021
#############################################################################

#Install and load in necessary libraries
#install.packages("ggplot2")
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(heatmaply)
library(gcookbook)
library(viridis)
rm(list = ls())

#set input and output directories
#input: Where are the files?
dir1 <- "/Users/yourname/example_dir/files_dir/"
#output: for figures 
dir2 <- "/Users/yourname/example_dir/output_dir/"

#color palletes
colsEMT <- brewer.pal(7, "RdYlBu")
colsPar <- brewer.pal(10, "PiYG")
EMT_col <- c(colsEMT[7], colsEMT[6], colsEMT[5], colsEMT[3], colsEMT[2], colsEMT[1], colsPar[2])


#Read in your annotated file names with all the samples labeled
file_names <- read.csv(paste0(dir1, "File_names.csv"),stringsAsFactors = F)


##################################################### Import HET & EMT scores from python readout #####################################################
#Use python code to generate heterogeneity and EMT scores for each file. Save these scores into the main file directory

Het <- read.csv(paste0(dir1, "Het_score.csv"), stringsAsFactors = F)
Het <- Het[order(Het$name),]

EMT <- read.csv(paste0(dir1, "EMT_score.csv"), stringsAsFactors = F)
EMT <- EMT[order(EMT$name),]


Scores <- as.data.frame(cbind(sample = Het$name, Het_score = Het$y_pred, EMT_score = EMT$score))
Scores$EMT_score <- as.numeric(Scores$EMT_score)
Scores$Het_score <- factor(Scores$Het_score, levels = c("high", "mid", "low"))

#check that the filenames and Het/EMT scores match
namecheck <- as.data.frame(cbind(file_names$sample, Scores$sample))

Scores <- cbind(Scores, group = file_names$group, sample = file_names$sample)

#Create a tertile EMT score
Scores$EMT_score_disc <- cut(pat_dat_tum$EMT_score, c(0,0.29,0.69,1), 
                                  labels = c("Epithelial","Intermediate", "Mesenchymal"))


#################### Plot Heterogeneity scores by group #########################
#heterogeneity score by group. ## out facet_wrap for no subset
ppi=300
png(paste0(dir2, "Het_barplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(Scores, aes(x=sample, fill = Het_score)) +
  geom_bar() + 
  scale_fill_manual(name="Het Score", values = viridis(3)) +
  labs(x="Clone", y="Percent of Images", title = "Tumor Heterogeneity Scores") +
  facet_wrap(~group) +
  theme_minimal() +
  theme(title =element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

#EMT score by group
ppi=300
png(paste0(dir2, "EMT_dotplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(Scores, aes(x=sample, fill = EMT_score)) +
  geom_bar() + 
  labs(x="Clone", y="Percent of Images", title = "Tumor Heterogeneity Scores") +
  facet_wrap(~group) +
  theme_minimal() +
  theme(title =element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

#EMT boxplot tertiles
ppi=300
png(paste0(dir2, "EMT_barplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(Scores, aes(x=sample, fill = EMT_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 0.015, color = NA) + scale_fill_manual(name="EMT Score", values = c("Epithelial" = EMT_col[1], "Intermediate" = "goldenrod2", "Mesenchymal" = EMT_col[6])) +
  labs(x="Clone", y="Percent of Images", title = "Tumor Heterogeneity Scores") +
  facet_wrap(~group) +
  theme_minimal() +
  theme(title =element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

#dotplot of EMT score stratified by Het score
ppi=300
png(paste0(dir2, "EMTbyHet_dotplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
ggplot(Scores, aes(x = Het_score, y= EMT_score, fill = Het_score)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all", width = 0.2, binwidth = 0.015, color = NA) +
  geom_boxplot(width=0.4, color = "gray24", fill = NA, size =0.4) +
  scale_fill_manual(name="Het Score", 
                    values = c("high" = c("#FC4E07"), "mid" = c("#E7B800"), "low" = c("#0097a2"))) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x="Heterogeneity Score", y="EMT Score", title = "EMT Score Distribution By Heterogeneity Score") +
  facet_wrap(~group) +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.5)), 
        #title =element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1)),
        axis.title.y.left = element_text(size=rel(1)))
dev.off()

ppi=300
png(paste0(dir2, "EMTbyHet_barplot.png"), height = 6*ppi, width = 5*ppi, res = ppi)
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


