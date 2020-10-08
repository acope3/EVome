
rm(list = ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")
library("tidyverse")
library("plyr")
library("qdap")
library("stringr")
library("DataCombine")
library("ggplot2")
library("gplots")

getwd()
setwd("/Users/Inna/Documents/Barr lab/EV proteomics/Combined peptide identification")

raw_peptides<-read.csv("raw_peptide_table.csv")

#Subset species-specific data into separate data frames

CAEEL_raw_peptides <- subset (raw_peptides, 
                              (Species == "CAEEL" | 
                                Species == "GAG_POL_GFP" | 
                                protein == "MS5322-tagRFP") &
                                model == 1)

ECOLI_raw_peptides <- subset (raw_peptides, 
                              ((Species == "E.coli" &
                                protein != "MS5322-tagRFP") &
                                 model == 1))

#Extract peptide sequences from all the runs for each species separately 
#in order to remap them against correcponding proteomes later
CAEEL_identified_peptides <- c(CAEEL_raw_peptides$sequence)
ECOLI_identified_peptides <- c(ECOLI_raw_peptides$sequence)

N_CAEEL_total_peptides <- length (CAEEL_identified_peptides)
N_ECOLI_total_peptides <- length (ECOLI_identified_peptides)

#Calculate peptide frequencies for each species
CAEEL_peptide_frequencies <- data.frame(table(CAEEL_identified_peptides))
ECOLI_peptide_frequencies <- data.frame(table(ECOLI_identified_peptides))

N_CAEEL_peptide_variety <- length (CAEEL_peptide_frequencies$CAEEL_identified_peptides)
N_ECOLI_peptide_variety <- length (ECOLI_peptide_frequencies$ECOLI_identified_peptides)

#Extract vector of observed unique peptides
CAEEL_tofind <- CAEEL_peptide_frequencies$CAEEL_identified_peptides
ECOLI_tofind <- ECOLI_peptide_frequencies$ECOLI_identified_peptides

#Import species proteomes sourced from BioMart
CAEEL_proteome_info<-read.csv("CAEEL_Proteome_info.csv")
ECOLI_proteome_info<-read.csv("ECOLI_Proteome_info.csv")


#here is a simple function to append one (or more) item to a list/ 
#It will be needed in the for loops later
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

#Generate empty lists to use in the loops below
CAEEL_L<-list()
ECOLI_L<-list()

#Locate peptide sequences in the C.elegans proteome
for (i in 1:length(CAEEL_tofind)) { 
  
  temp_vector <- grep(paste(CAEEL_tofind[i]), CAEEL_proteome_info$peptide)
  
  CAEEL_L<-lappend(CAEEL_L, temp_vector)
  names(CAEEL_L)[i] <- paste0(CAEEL_tofind[i])
  
}
head(CAEEL_L)


#Locate peptide sequences in the E. coli proteome
for (i in 1:length(ECOLI_tofind)) { 
  
  temp_vector <- grep(paste(ECOLI_tofind[i]), ECOLI_proteome_info$Sequence)
  
  ECOLI_L<-lappend(ECOLI_L, temp_vector)
  names(ECOLI_L)[i] <- paste0(ECOLI_tofind[i])
  
}
head(ECOLI_L)

#generate empty lists to use in the loops below
CAEEL_L2<-list()
ECOLI_L2<-list()

#subsetting C. elegans gene names for each peptide from the C.elegans proteome file 
for (i in 1:length(CAEEL_L)) { 
  
  x<-CAEEL_L[[i]]
  temp_vector <- CAEEL_proteome_info$external_gene_name[x]
  
  # temp_vector<-do.call(c, temp_vector) #concatenate several vectors into one vector
  
  temp_vector<-unique(temp_vector)
  
  CAEEL_L2<-lappend(CAEEL_L2, temp_vector)
  
  
}
head(CAEEL_L2,20)

CAEEL_all_possible_proteins<-unique(unlist(CAEEL_L2))
write.csv(CAEEL_all_possible_proteins, "CAEEL_all_possible_proteins.csv")
N_CAEEL_all_possible_proteins <- length (CAEEL_all_possible_proteins)

#subsetting E. coli gene names from the E. coli proteome file
for (i in 1:length(ECOLI_L)) { 
  
  x<-ECOLI_L[[i]]
  temp_vector <- ECOLI_proteome_info$Gene_name[x]
  
  # temp_vector<-do.call(c, temp_vector) #concatenate several vectors into one vector
  
  temp_vector<-unique(temp_vector)
  
  ECOLI_L2<-lappend(ECOLI_L2, temp_vector)
  
}
head(ECOLI_L2,20)

ECOLI_all_possible_proteins<-unique(unlist(ECOLI_L2))
write.csv(ECOLI_all_possible_proteins, "ECOLI_all_possible_proteins.csv")
N_ECOLI_all_possible_proteins <- length (ECOLI_all_possible_proteins)


#generate empty vector
CAEEL_name_strings<-c()
ECOLI_name_strings<-c()

#populate vector with strings of proteins to which peptides were mapped
for (i in (1:length(CAEEL_L2))) {
  
  g<-CAEEL_L2[[i]]
  m<-paste(as.character(g), collapse=", ")
  CAEEL_name_strings<-append(CAEEL_name_strings, m)
}

head(CAEEL_name_strings)
head(CAEEL_peptide_frequencies)
length(CAEEL_name_strings)
length(CAEEL_peptide_frequencies$Freq)

for (i in (1:length(ECOLI_L2))) {
  
  g<-ECOLI_L2[[i]]
  m<-paste(as.character(g), collapse=", ")
  ECOLI_name_strings<-append(ECOLI_name_strings, m)
}

head(ECOLI_name_strings)
head(ECOLI_peptide_frequencies)
length(ECOLI_name_strings)
length(ECOLI_peptide_frequencies$Freq)

#Compile a table of all identified peptides with gene names
CAEEL_peptides_remapped <- data.frame(Peptide_sequence=CAEEL_peptide_frequencies$CAEEL_identified_peptides, 
                                      Peptide_count=CAEEL_peptide_frequencies$Freq,
                                      Gene_names=CAEEL_name_strings)

ECOLI_peptides_remapped <- data.frame(Peptide_sequence=ECOLI_peptide_frequencies$ECOLI_identified_peptides, 
                                Peptide_count=ECOLI_peptide_frequencies$Freq,
                                Gene_names=ECOLI_name_strings)

#Remove empty peptides that were not assigned to a protein in Biomart 
# (all of them (27 in total) were checked manually, are not evident E.coli proteins)
ECOLI_notmapped <- ECOLI_peptides_remapped[ECOLI_peptides_remapped$Gene_names=="", ]
ECOLI_peptides_remapped<-ECOLI_peptides_remapped[!(ECOLI_peptides_remapped$Gene_names==""), ]

#Remove these peptides from the raw database
ECOLI_raw_peptides2 <- ECOLI_raw_peptides[ ! ECOLI_raw_peptides$sequence %in% ECOLI_notmapped$Peptide_sequence, ]


write.csv(CAEEL_peptides_remapped, file="CAEEL_peptides_remapped.csv")
write.csv(ECOLI_peptides_remapped, file="ECOLI_peptides_remapped.csv")
write.csv(ECOLI_raw_peptides2, file="ECOLI_raw_peptides2.csv")


#Subsetting only peptides that are mapped to one protein
CAEEL_uniquely_mapped <- grepl.sub(data = CAEEL_peptides_remapped, pattern = ",", 
                                   Var = "Gene_names", keep.found=FALSE)
head(CAEEL_uniquely_mapped,30)
length(CAEEL_uniquely_mapped$Peptide_sequence)
N_CAEEL_uniquely_mapping_peptides <- length(CAEEL_uniquely_mapped$Peptide_sequence)
N_CAEEL_unique_proteins <- length(unique(CAEEL_uniquely_mapped$Gene_names))
write.csv(CAEEL_uniquely_mapped, "CAEEL_uniquely_mapped.csv")

#subsetting non-unique peptides that map to a few proteins of isoforms
CAEEL_nonuniquely_mapped <- grepl.sub(data = CAEEL_peptides_remapped, pattern = ",", 
                                   Var = "Gene_names", keep.found=TRUE)
write.csv(CAEEL_nonuniquely_mapped, "CAEEL_nonuniquely_mapped.csv")



ECOLI_uniquely_mapped <- grepl.sub(data = ECOLI_peptides_remapped, pattern = ",", 
                                   Var = "Gene_names", keep.found=FALSE)
head(ECOLI_uniquely_mapped,30)
length(ECOLI_uniquely_mapped$Peptide_sequence)
N_ECOLI_uniquely_mapping_peptides <- length(ECOLI_uniquely_mapped$Peptide_sequence)
N_ECOLI_unique_proteins <- length(unique(ECOLI_uniquely_mapped$Gene_names))

write.csv(ECOLI_uniquely_mapped, "ECOLI_uniquely_mapped.csv")
#subsetting non-unique peptides that map to a few proteins of isoforms
ECOLI_nonuniquely_mapped <- grepl.sub(data = ECOLI_peptides_remapped, pattern = ",", 
                                      Var = "Gene_names", keep.found=TRUE)
write.csv(ECOLI_nonuniquely_mapped, "ECOLI_nonuniquely_mapped.csv")

# Calculates sum of identified peptide for each protein
CAEEL_genes_with_unique_peptide_counts <- aggregate(CAEEL_uniquely_mapped$Peptide_count, 
                                      by = list(Category=CAEEL_uniquely_mapped$Gene_names), 
                                      FUN=sum)
head(CAEEL_genes_with_unique_peptide_counts)

write.csv(CAEEL_genes_with_unique_peptide_counts, "CAEEL_genes_with_unique_peptide_counts.csv")

ECOLI_genes_with_unique_peptide_counts <- aggregate(ECOLI_uniquely_mapped$Peptide_count, 
                                                    by = list(Category=ECOLI_uniquely_mapped$Gene_names), 
                                                    FUN=sum)
head(ECOLI_genes_with_unique_peptide_counts)

write.csv(ECOLI_genes_with_unique_peptide_counts, "ECOLI_genes_with_unique_peptide_counts.csv")

#Figuring out genes that were identified only by nonuniquely mapped peptides

CAAEL_nonunique <- unlist(strsplit(CAEEL_nonuniquely_mapped$Gene_names, split=", "))
head(CAAEL_nonunique,5)
length(unique(CAAEL_nonunique))
CAAEL_nonunique<-unique(CAAEL_nonunique)
write.csv(CAAEL_nonunique, "CAAEL_nonunique.csv")


#---------------------------------------
#Processing of individual samples

# Create an empty list of vectors
CAEEL_peptides_by_sample <- list()
ECOLI_peptides_by_sample <- list()

# Populate the list of vectors with peptide sequences identified in each sample

CAEEL_peptides_by_sample[[1]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "light")
CAEEL_peptides_by_sample[[2]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "heavy")
CAEEL_peptides_by_sample[[3]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "1")
CAEEL_peptides_by_sample[[4]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "2")
CAEEL_peptides_by_sample[[5]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "3")
CAEEL_peptides_by_sample[[6]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "4")
CAEEL_peptides_by_sample[[7]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "5")
CAEEL_peptides_by_sample[[8]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "6")
CAEEL_peptides_by_sample[[9]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "7")
CAEEL_peptides_by_sample[[10]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "8")
CAEEL_peptides_by_sample[[11]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "9")
CAEEL_peptides_by_sample[[12]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "10")
CAEEL_peptides_by_sample[[13]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "11")
CAEEL_peptides_by_sample[[14]]<-subset(CAEEL_raw_peptides$sequence, CAEEL_raw_peptides$sample == "upper" |
                                         CAEEL_raw_peptides$sample == "lower")

names(CAEEL_peptides_by_sample) <- c("Light_01", "Heavy_01", "Light_02", "Heavy_02",
                                     "Light_03", "Heavy_03", "Light_04", "Heavy_04",
                                     "Light_05", "Heavy_05", "Light_06", "Light_07",
                                     "Heavy_06", "Mixed")

ECOLI_peptides_by_sample[[1]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "light")
ECOLI_peptides_by_sample[[2]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "heavy")
ECOLI_peptides_by_sample[[3]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "1")
ECOLI_peptides_by_sample[[4]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "2")
ECOLI_peptides_by_sample[[5]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "3")
ECOLI_peptides_by_sample[[6]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "4")
ECOLI_peptides_by_sample[[7]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "5")
ECOLI_peptides_by_sample[[8]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "6")
ECOLI_peptides_by_sample[[9]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "7")
ECOLI_peptides_by_sample[[10]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "8")
ECOLI_peptides_by_sample[[11]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "9")
ECOLI_peptides_by_sample[[12]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "10")
ECOLI_peptides_by_sample[[13]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "11")
ECOLI_peptides_by_sample[[14]]<-subset(ECOLI_raw_peptides$sequence, ECOLI_raw_peptides$sample == "upper" |
                                         ECOLI_raw_peptides$sample == "lower")

names(ECOLI_peptides_by_sample) <- c("Light_01", "Heavy_01", "Light_02", "Heavy_02",
                                     "Light_03", "Heavy_03", "Light_04", "Heavy_04",
                                     "Light_05", "Heavy_05", "Light_06", "Light_07",
                                     "Heavy_06", "Mixed")

CAEEL_df_frequencies <- list()
ECOLI_df_frequencies <- list()
#For each vector in the list calculate frequencies of each peptide, 
#generate 11 data frames, put them into one list
for (i in (1:length(CAEEL_peptides_by_sample))) {
  
  num <- str_pad(i, 2, pad="0")
  Samplename <- paste("CAEEL_peptides_Sample_", num, sep = '')
  assign (Samplename, data.frame(table(CAEEL_peptides_by_sample[i])))
  CAEEL_df_frequencies[[i]] <- get(Samplename)
  
}

for (i in (1:length(ECOLI_peptides_by_sample))) {
  
  num <- str_pad(i, 2, pad="0")
  Samplename <- paste("ECOLI_peptides_Sample_", num, sep = '')
  assign (Samplename, data.frame(table(ECOLI_peptides_by_sample[i])))
  ECOLI_df_frequencies[[i]] <- get(Samplename)

}

#For each peptide in the list of dataframes find gene name (in the Peptides_remapped) that it maps to

CAEEL_df_all_proteins <- list()
ECOLI_df_all_proteins <- list()

for (i in (1:length(CAEEL_df_frequencies))) {
  
  num <- str_pad(i, 2, pad="0")
  Samplename <- paste("CAEEL_df_all_proteins_Sample_", num, sep = '')
  assign (Samplename, merge(CAEEL_df_frequencies[i], CAEEL_peptides_remapped[ , c("Peptide_sequence", "Gene_names")], 
                            by.x="Var1", by.y="Peptide_sequence",
                            all.x=TRUE, all.y=FALSE))
  CAEEL_df_all_proteins[[i]] <- get(Samplename)
write.csv(CAEEL_df_all_proteins[[i]], paste("CAEEL_Sample_", num, "all_peptides.csv"))
}

for (i in (1:length(ECOLI_df_frequencies))) {
  
  num <- str_pad(i, 2, pad="0")
  Samplename <- paste("ECOLI_df_all_proteins_Sample_", num, sep = '')
  assign (Samplename, merge(ECOLI_df_frequencies[i], ECOLI_peptides_remapped[ , c("Peptide_sequence", "Gene_names")], 
                            by.x="Var1", by.y="Peptide_sequence",
                            all.x=TRUE, all.y=FALSE))
  ECOLI_df_all_proteins[[i]] <- get(Samplename)
  write.csv(ECOLI_df_all_proteins[[i]], paste("ECOLI_Sample_", num, "all_peptides.csv"))
}

#Subsetting peptides that map to only one protein
CAEEL_df_unique_proteins <- list()
ECOLI_df_unique_proteins <-list()
for (i in (1:length(CAEEL_df_all_proteins))) {
  
  num <- str_pad(i, 2, pad="0")
  Samplename <- paste("CAEEL_df_unique_proteins_", num, sep = '')
  assign (Samplename, grepl.sub(data = CAEEL_df_all_proteins[[i]], pattern = ",", 
                                Var = "Gene_names", keep.found=FALSE))
  CAEEL_df_unique_proteins[[i]] <- get(Samplename)
  write.csv(CAEEL_df_unique_proteins[[i]], paste("CAEEL_Sample_", num, "unique_peptides.csv"))
}

for (i in (1:length(ECOLI_df_all_proteins))) {
  
  num <- str_pad(i, 2, pad="0")
  Samplename <- paste("ECOLI_df_unique_proteins_", num, sep = '')
  assign (Samplename, grepl.sub(data = ECOLI_df_all_proteins[[i]], pattern = ",", 
                                Var = "Gene_names", keep.found=FALSE))
  ECOLI_df_unique_proteins[[i]] <- get(Samplename)
  write.csv(ECOLI_df_unique_proteins[[i]], paste("ECOLI_Sample_", num, "unique_peptides.csv"))
}

#Generate summary tables for each sample with sums of the peptide by gene name

CAEEL_df_summaries_by_gene_uniquepeptides <- list()
ECOLI_df_summaries_by_gene_uniquepeptides <- list()

for (i in (1:length(CAEEL_df_unique_proteins))) {
  
  num <- str_pad(i, 2, pad="0")
  Samplename <- paste("CAEEL_df_bygene_unique_peptide_counts_", num, sep = '')
  assign (Samplename, aggregate(CAEEL_df_unique_proteins[[i]]$Freq, 
                                by = list(CAEEL_df_unique_proteins[[i]]$Gene_names), 
                                FUN=sum))
  CAEEL_df_summaries_by_gene_uniquepeptides[[i]] <- get(Samplename)
  write.csv(CAEEL_df_summaries_by_gene_uniquepeptides[[i]], paste("CAEEL_Sample_", num, "unique_protein_summary.csv"))
}

for (i in (1:length(ECOLI_df_unique_proteins))) {
  
  num <- str_pad(i, 2, pad="0")
  Samplename <- paste("ECOLI_df_bygene_unique_peptide_counts_", num, sep = '')
  assign (Samplename, aggregate(ECOLI_df_unique_proteins[[i]]$Freq, 
                                by = list(ECOLI_df_unique_proteins[[i]]$Gene_names), 
                                FUN=sum))
  ECOLI_df_summaries_by_gene_uniquepeptides[[i]] <- get(Samplename)
  write.csv(ECOLI_df_summaries_by_gene_uniquepeptides[[i]], paste("ECOLI_Sample_", num, "unique_protein_summary.csv"))
}

#Combine dataframes with peptide counts for each sample into one data frame
CAEEL_gene_vector<-CAEEL_genes_with_unique_peptide_counts$Category
CAEEL_summary <- data.frame(matrix(ncol=1,nrow=length(CAEEL_genes_with_unique_peptide_counts$Category), 
                                   dimnames=list(NULL, c("Gene_name"))))
CAEEL_summary$Gene_name <- CAEEL_gene_vector
colnames(CAEEL_summary)<-c("Gene_name")

for (i in (1:length(CAEEL_df_summaries_by_gene_uniquepeptides))) {
  
  CAEEL_summary<-merge(CAEEL_summary, CAEEL_df_summaries_by_gene_uniquepeptides[i],  
                            by.x="Gene_name", by.y="Group.1", 
                            all.x=TRUE, all.y=TRUE)
   num <- str_pad(i, 2, pad="0")
  colnames(CAEEL_summary)[i+1] <- paste("Sample_", num, sep = '')
 
}

ECOLI_gene_vector<-ECOLI_genes_with_unique_peptide_counts$Category
ECOLI_summary <- data.frame(matrix(ncol=1,nrow=length(ECOLI_genes_with_unique_peptide_counts$Category), dimnames=list(NULL, c("Gene_name"))))
ECOLI_summary$Gene_name <- ECOLI_gene_vector
colnames(ECOLI_summary)<-c("Gene_name")

for (i in (1:length(ECOLI_df_summaries_by_gene_uniquepeptides))) {
  
  ECOLI_summary<-merge(ECOLI_summary, ECOLI_df_summaries_by_gene_uniquepeptides[i],  
                       by.x="Gene_name", by.y="Group.1", 
                       all.x=TRUE, all.y=TRUE)
  num <- str_pad(i, 2, pad="0")
  colnames(ECOLI_summary)[i+1] <- paste("ECOLI_Sample_", num, sep = '')
  
}

#Replace NA peptide counts with 0 values (i.e. no such peptides were identified in the given sample)

CAEEL_summary[, 2:15][is.na(CAEEL_summary[, 2:15])] <- 0
ECOLI_summary[, 2:15][is.na(ECOLI_summary[, 2:15])] <- 0
names(CAEEL_summary) <- c("Gene_name", "Light_01", "Heavy_01", "Light_02", "Heavy_02",
                                     "Light_03", "Heavy_03", "Light_04", "Heavy_04",
                                     "Light_05", "Heavy_05", "Light_06", "Light_07",
                                     "Heavy_06", "Mixed")
names(ECOLI_summary) <- c("Gene_name", "Light_01", "Heavy_01", "Light_02", "Heavy_02",
                          "Light_03", "Heavy_03", "Light_04", "Heavy_04",
                          "Light_05", "Heavy_05", "Light_06", "Light_07",
                          "Heavy_06", "Mixed")

write.csv(CAEEL_summary, "CAEEL_Proteomics_summary.csv")
write.csv(ECOLI_summary, "ECOLI_Proteomics_summary.csv")

CAEEL_summary_with_uniprot <- merge (CAEEL_summary, DB_uniprot_cds_genes[,c("external_gene_name","uniprot" )],
                                     by.x = "Gene_name", by.y = "external_gene_name",
                                     all.x=TRUE, all.y=FALSE)
head(CAEEL_summary_with_uniprot,40)
write.csv(CAEEL_summary_with_uniprot, "CAEEL_summary_with_uniprot.csv")

help(ConstructGenesTree)
GetNamesTaxa()
library("UniprotR")
GetProteomeFasta(UP000001940,
                 directorypath = "/Users/Inna/Documents/Barr lab/EV proteomics/Combined peptide identification")

Uniprot_proteome <- read.delim("uniprot-proteome_UP000001940.tab", stringsAsFactors = FALSE, quote = "", sep = "\t")
#---------------------------------------------------------------

#Principal Component Analysis
library("ggfortify")
tCAEEL_summary <- as.data.frame(t(CAEEL_summary[2:15]))

tCAEEL_summary$sample_type <- c("Light", "Heavy", "Light","Heavy",
                        "Light", "Heavy", "Light", 
                        "Heavy", "Light", "Heavy", 
                        "Light", "Light", "Heavy", "Mixed")

CAEEL_pca_res <- prcomp(tCAEEL_summary[1:2888], scale. = TRUE)
p1<-autoplot(CAEEL_pca_res, data = tCAEEL_summary, label=TRUE, colour = "sample_type", label.size = 3)



tECOLI_summary <- as.data.frame(t(ECOLI_summary[2:15]))


tECOLI_summary$sample_type <- c("Light", "Heavy", "Light","Heavy",
                                "Light", "Heavy", "Light", 
                                "Heavy", "Light", "Heavy", 
                                "Light", "Light", "Heavy", "Mixed")

ECOLI_pca_res <- prcomp(tECOLI_summary[1:1495], scale. = TRUE)
p2<-autoplot(ECOLI_pca_res, data = tECOLI_summary, label=TRUE, colour = "sample_type", label.size = 3)

ggsave("CAEEL_PCA_plot.pdf", p1, width=4, height=4, 
       device = "pdf", units = "in")
ggsave("ECOLI_PCA_plot.pdf", p2, width=4, height=4, 
       device = "pdf", units = "in")

#Correlation plots
library("PerformanceAnalytics")

chart.Correlation(log10(CAEEL_summary[2:15]+1),
                  method="pearson",
                  histogram=TRUE,
                  pch=16)

#--------------------------------------------
#Saturation test attempt for C.elegans

library("gplots")

CAEEL_venn_vector_list <- list(CAEEL_df_bygene_unique_peptide_counts_05$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_04$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_03$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_06$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_10$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_08$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_01$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_12$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_02$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_07$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_09$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_13$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_14$Group.1,
                      CAEEL_df_bygene_unique_peptide_counts_11$Group.1)

CAEEL_tmp <- CAEEL_venn_vector_list[[1]]
CAEEL_venn_values <- c()
CAEEL_saturation_plot_data <- c() 
for (i in (1:length(CAEEL_venn_vector_list)-1)) {
  
  CAEEL_venn_values <- venn(list(unique(CAEEL_tmp), CAEEL_venn_vector_list[[i+1]]), 
              show.plot=FALSE)
  CAEEL_saturation_plot_data <- c(CAEEL_saturation_plot_data, CAEEL_venn_values[2] 
                            )
  CAEEL_tmp <- c(CAEEL_tmp, CAEEL_venn_vector_list[[i+1]])
  
}
CAEEL_saturation_plot_data

CAEEL_saturation_cumulative <- c(length(CAEEL_venn_vector_list[[1]]))

CAEEL_saturation_cumulative


CAEEL_df_saturation <- data.frame(c(1:14), CAEEL_saturation_cumulative)
colnames(CAEEL_df_saturation)<-c("Consequtive_sample", "Cumulative_newly_discovered_proteins")
CAEEL_df_saturation

#Plotting saturation graphs for C. elegans
library ("ggplot2")
library("ggrepel")

p1 <-  ggplot(CAEEL_df_saturation, aes(x=Consequtive_sample, y=Cumulative_newly_discovered_proteins))+ 
  geom_point(color="black", shape=16, size=2)+
  scale_x_continuous(limits=c(1,14), breaks=seq(1, 14, by=1))+
  ylim(0,2500)+
  ylab("Cumulative newly discovered proteins")+
  xlab("Consequtive samples")+
  theme_minimal()+
  theme(axis.ticks.y.left = element_line(color="black"),
        axis.ticks.x = element_line(color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title.y = element_text(size = 10, color="black", 
                                    margin = margin(t=0, r=10, b=0, l=0)), 
        axis.title.x = element_text(size = 10, color="black", 
                                    margin = margin(t=10, r=0, b=0, l=0)),
        axis.text.y = element_text(size = 10, color="black"), 
        axis.text.x=element_text(size=10, 
                                 color="black"))

ggsave("CAEEL_saturation_cumulative_descending.pdf", p1, width=4, height=4, 
       device = "pdf", units = "in")

#Saturation test attempt for E. coli

library("gplots")

ECOLI_venn_vector_list <- list(ECOLI_df_bygene_unique_peptide_counts_02$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_04$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_06$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_08$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_10$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_13$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_14$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_01$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_03$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_05$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_07$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_09$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_11$Group.1,
                               ECOLI_df_bygene_unique_peptide_counts_12$Group.1)

ECOLI_tmp <- ECOLI_venn_vector_list[[1]]
ECOLI_venn_values <- c()
ECOLI_saturation_plot_data <- c() 
for (i in (1:length(ECOLI_venn_vector_list)-1)) {
  
  ECOLI_venn_values <- venn(list(unique(ECOLI_tmp), ECOLI_venn_vector_list[[i+1]]), 
                            show.plot=FALSE)
  ECOLI_saturation_plot_data <- c(ECOLI_saturation_plot_data, ECOLI_venn_values[2] 
  )
  ECOLI_tmp <- c(ECOLI_tmp, ECOLI_venn_vector_list[[i+1]])
  
}
ECOLI_saturation_plot_data

ECOLI_saturation_cumulative <- c(length(ECOLI_venn_vector_list[[1]]))

for (i in (1:(length(ECOLI_saturation_plot_data)-1))) {
  
  n <- sum (ECOLI_saturation_plot_data[i+1], ECOLI_saturation_cumulative[i])
  ECOLI_saturation_cumulative <- c(ECOLI_saturation_cumulative, n)
  
}
ECOLI_saturation_cumulative


ECOLI_df_saturation <- data.frame(c(1:14), ECOLI_saturation_cumulative)
colnames(ECOLI_df_saturation)<-c("Consequtive_sample", "Cumulative_newly_discovered_proteins")
ECOLI_df_saturation

#Plotting saturation graphs for E. coli
library ("ggplot2")
library("ggrepel")

p1 <-  ggplot(ECOLI_df_saturation, aes(x=Consequtive_sample, y=Cumulative_newly_discovered_proteins))+ 
  geom_point(color="black", shape=16, size=2)+
  scale_x_continuous(limits=c(1,14), breaks=seq(1, 14, by=1))+
  ylim(0,1500)+
  ylab("Cumulative newly discovered proteins")+
  xlab("Consequtive samples")+
  theme_minimal()+
  theme(axis.ticks.y.left = element_line(color="black"),
        axis.ticks.x = element_line(color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        axis.title.y = element_text(size = 10, color="black", 
                                    margin = margin(t=0, r=10, b=0, l=0)), 
        axis.title.x = element_text(size = 10, color="black", 
                                    margin = margin(t=10, r=0, b=0, l=0)),
        axis.text.y = element_text(size = 10, color="black"), 
        axis.text.x=element_text(size=10, 
                                 color="black"))

ggsave("ECOLI_saturation_cumulative_heavy_first.pdf", p1, width=4, height=4, 
       device = "pdf", units = "in")


  
library("VennDiagram")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
library("limma")

vennDiagram(vennCounts(cbind(CAEEL_df_unique_proteins_01$Gene_names, 
                            CAEEL_df_unique_proteins_02$Gene_names)))

for (i in (1:(length(CAEEL_saturation_plot_data)-1))) {
  
  n <- sum (CAEEL_saturation_plot_data[i+1], CAEEL_saturation_cumulative[i])
  CAEEL_saturation_cumulative <- c(CAEEL_saturation_cumulative, n)
  
}
p<-calculate.overlap(list(unique(CAEEL_df_unique_proteins_01$Gene_names), 
                          unique(CAEEL_df_unique_proteins_02$Gene_names)))

grid.arrange(p1, p2, nrow = 1)
