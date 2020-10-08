rm(list = ls())

library("tidyverse")
library("plyr")
library("qdap")
library("stringr")

getwd()
setwd("/Users/Inna/Documents/Barr lab/EV proteomics/Combined Enrichment analysis")

#import identified gene table
Merged_uniquely_mapped<-read.csv("CAEEL_Proteomics_summary.csv", stringsAsFactors = FALSE)
head(Merged_uniquely_mapped)


#extract identified genes into a vector
total_genes_identified <- unique(Merged_uniquely_mapped$Gene_name)
length(total_genes_identified)

#Source add'l info about identified genes info Biomart
library("biomaRt")
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useMart("ensembl",dataset="celegans_gene_ensembl")
output<-c("external_gene_name", "wormbase_gene", "description", "chromosome_name")
input<-c("external_gene_name")
tada_info<-getBM(
  attributes=output,
  filters=input,
  values=total_genes_identified,
  mart=ensembl)
head(tada_info)

#Remove info from square brackets in the description
tada.info <- as_tibble(tada_info)
tada_info <- mutate(tada_info,
                         description = rm_between(tada_info$description, "[", "]"))
tada.info <- as.data.frame(tada_info)
head(tada_info)

#Compare vectors "total_genes_identified" and tada_info$external_gene_name
length(setdiff(total_genes_identified, tada_info$external_gene_name))
length(setdiff(tada_info$external_gene_name, total_genes_identified))
length(intersect(total_genes_identified, tada_info$external_gene_name))
p<-which(duplicated(tada_info$external_gene_name))
p

CAEEL_proteome <- read.csv("CAEEL_Proteome_Biomart_corrected.csv", stringsAsFactors = FALSE)

#Compile table of identified genes
Summary_Merged_uniquely_mapped <- merge(Merged_uniquely_mapped, CAEEL_proteome[,c("external_gene_name", "wormbase_gene")], 
                                        by.x="Gene_name", by.y="external_gene_name",
                                        all.x=TRUE, all.y=FALSE)

Summary_Merged_uniquely_mapped <- unique(Summary_Merged_uniquely_mapped)
head(Summary_Merged_uniquely_mapped)

#Reorder and rename some columns
Summary_Merged_uniquely_mapped$X <- NULL

#Import Cao et al tables
TableS3_Cao<-read.csv("Cao_TableS3.csv", stringsAsFactors = FALSE)
TableS4_Cao<-read.csv("Cao_TableS4.csv",stringsAsFactors = FALSE)
TableS5_Cao<-read.csv("Cao_TableS5.csv", stringsAsFactors = FALSE)

#Calculate enrichment scores as tissue-specific TPM / average TPM across all the tissues 
# (6 - number of tissue types) to calculate average in the denominator
head(TableS3_Cao)
TableS3_Cao$NeuroEnrich <- TableS3_Cao$Neurons*6/rowSums(TableS3_Cao[,3:9])
TableS3_Cao$GliaEnrich <- TableS3_Cao$Glia*6/rowSums(TableS3_Cao[,3:9])
TableS3_Cao$HypodermisEnrich <- TableS3_Cao$Hypodermis*6/rowSums(TableS3_Cao[,3:9])
TableS3_Cao$PharynxEnrich <- TableS3_Cao$Pharynx*6/rowSums(TableS3_Cao[,3:9])
TableS3_Cao$Body_wall_muscleEnrich <- TableS3_Cao$Body_wall_muscle*6/rowSums(TableS3_Cao[,3:9])
TableS3_Cao$IntestineEnrich <- TableS3_Cao$Intestine*6/rowSums(TableS3_Cao[,3:9])
TableS3_Cao$GonadEnrich <- TableS3_Cao$Gonad*6/rowSums(TableS3_Cao[,3:9])

#Replace NA with 0
TableS3_Cao[is.na(TableS3_Cao)] <- 0
head(TableS3_Cao)

#Calculate enrichment scores as Ciliated_sensory_neurons TPM / average TPM across all the non-neuronal tissues
head(TableS4_Cao)
colnames(TableS4_Cao)
TableS4_Cao$CiliatedEnrich <- TableS4_Cao$Ciliated_sensory_neuron*27/rowSums(TableS4_Cao[,3:29])

#Replace NA values with 0
which(is.na(TableS4_Cao$CiliatedEnrich))
TableS4_Cao[is.na(TableS4_Cao)] <- 0
which(is.infinite(TableS4_Cao$CiliatedEnrich))

#Calculate enrichment scores as Cholinergic (cluster 15) TPM / average TPM across all the neurons
head(TableS5_Cao)
ncol(TableS5_Cao)
TableS5_Cao$IL2_Enrich <- TableS5_Cao$Cholinergic_.15.*40/rowSums(TableS5_Cao[,3:42])

#Replace NA values with 0
which(is.na(TableS5_Cao$IL2_Enrich))
TableS5_Cao[is.na(TableS5_Cao)] <- 0
which(is.infinite(TableS5_Cao$IL2_Enrich))

#Extract tissue enrichment scores into data frame Summary_Merged_uniquely_mapped
head(Summary_Merged_uniquely_mapped)
WBGeneIDs_tofind <- Summary_Merged_uniquely_mapped$wormbase_gene

Scores_df <- merge(Summary_Merged_uniquely_mapped, TableS3_Cao, 
                   by.x = "wormbase_gene", by.y="gene_id",
                   all.x = TRUE, all.y=FALSE)
nrow(Scores_df)

#Drop column "symbol"
Scores_df<-within(Scores_df, rm("symbol"))
head(Scores_df)

#Extract ciliated neurons enrichment score
Scores_df <- merge(Scores_df, TableS4_Cao, 
                   by.x = "wormbase_gene", by.y="gene_id",
                   all.x = TRUE, all.y=FALSE)
#Drop column "symbol" for the second time
head(Scores_df)
Scores_df<-within(Scores_df, rm("symbol"))
head(Scores_df)

#Extract Il2_Enrich score from TableS5_Cao
Scores_df <- merge(Scores_df, TableS5_Cao, 
                   by.x = "wormbase_gene", by.y="gene_id",
                   all.x = TRUE, all.y=FALSE)

#Drop column "symbol" for the third time
head(Scores_df)
Scores_df<-within(Scores_df, rm("symbol"))
head(Scores_df)

#Reorder columns
#colnames(Scores_df)
#reorder_vector <- c(1,2,3,4,5,6,7,8,16,17,18,19,20,21,22,
#                    50,91,9,10,11,12,13,14,15,23,24,25,26,
 #                   27,28,29,30,31,32,33,34,35,36,37,38,39,
  #                  40,41,42,43,44,45,46,47,48,49,51,52,53,
   #                 54,55,56,57,58,59,60,61,62,63,64,65,66,
    #                67,68,69,70,71,72,73,74,75,76,77,78,79,
     #               80,81,82,83,84,85,86,87,88,89,90)
#Scores_df <- Scores_df[reorder_vector]

write.csv(Scores_df, "Scores giant table.csv")
Scores_df <-read.csv("Scores giant table.csv")


#GENERATE PLOTS with ENRICHMENT SCORES

library ("ggplot2")
library("ggrepel")


#Neuronal and Ciliated_over_other_tissues
p1<-ggplot(Scores_df, aes(x=NeuroEnrich, y=CiliatedEnrich), alpha=0.5)+ 
  geom_point(color="grey", shape=17, size=0.8)+
  #Present in neurons and ciliated neurons
  geom_point(data = Scores_df[(Scores_df$NeuroEnrich >5) & (Scores_df$CiliatedEnrich>5), ],
             aes(x=NeuroEnrich, y=CiliatedEnrich), color="#601A4A", size=1.6)+
  geom_text_repel(data = Scores_df[(Scores_df$NeuroEnrich >5) & (Scores_df$CiliatedEnrich>5), ],
                  aes(x=NeuroEnrich, y=CiliatedEnrich, label=Gene_name), size=2)+
  #Present mostly in ciliated neurons
  geom_point(data = Scores_df[(Scores_df$NeuroEnrich <5) & (Scores_df$CiliatedEnrich>5), ],
             aes(x=NeuroEnrich, y=CiliatedEnrich), color="#EE442F", size=1.6)+
  geom_text_repel(data = Scores_df[(Scores_df$NeuroEnrich <5) & (Scores_df$CiliatedEnrich>5), ],
                  aes(x=NeuroEnrich, y=CiliatedEnrich, label=Gene_name), size=2)+
  
  #scale_y_continuous(limits=c(0,28), breaks=c(0,4,8,12,16,20,24,28))+
  xlim(0,7)+
  ylab("Enrichment in ciliated neurons over other tissues")+
  xlab("Enrichment in neurons over other tissues")+
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

ggsave("NeuroEnrich_over_CiliatedEnrich.pdf", p1, width=4, height=4, 
       device = "pdf", units = "in")

# error of 4 genes being excluded is related to the fact that they were not idenitified in Cao et al 
# these genes are Y41C4A.32, K04A8.21, F53E2.2, Y49F6B.20


colnames(Scores_df)
#IL_Enriched over CiliatedEnrich
p2<-ggplot(Scores_df, aes(x=IL2_Enrich, y=CiliatedEnrich), alpha=0.5)+ 
  geom_point(color="grey", shape=17, size=0.8)+
  #Present in IL2 neurons  ciliated neurons
  geom_point(data = Scores_df[(Scores_df$IL2_Enrich >5) & (Scores_df$CiliatedEnrich >5), ],
             aes(x=IL2_Enrich, y=CiliatedEnrich), color="#601A4A", size=1.6)+
  geom_text_repel(data = Scores_df[(Scores_df$IL2_Enrich >5) & (Scores_df$CiliatedEnrich >5), ],
                  aes(x=IL2_Enrich, y=CiliatedEnrich, label=Gene_name), size=2)+
  #Present mostly in ciliated neurons
  geom_point(data = Scores_df[(Scores_df$IL2_Enrich >5) & (Scores_df$CiliatedEnrich <5), ],
             aes(x=IL2_Enrich, y=CiliatedEnrich), color="#EE442F", size=1.6)+
  geom_text_repel(data = Scores_df[(Scores_df$IL2_Enrich >10) & (Scores_df$CiliatedEnrich <5), ],
                  aes(x=IL2_Enrich, y=CiliatedEnrich, label=Gene_name), size=2)+
  #Present mostly in ciliated neurons but excluded from IL2s (EVNs)
  geom_point(data = Scores_df[(Scores_df$IL2_Enrich <5) & (Scores_df$CiliatedEnrich >5), ],
             aes(x=IL2_Enrich, y=CiliatedEnrich), color="green", size=1.6)+
  geom_text_repel(data = Scores_df[(Scores_df$IL2_Enrich <5) & (Scores_df$CiliatedEnrich >5), ],
                  aes(x=IL2_Enrich, y=CiliatedEnrich, label=Gene_name), size=2)+
  
  #scale_y_continuous(limits=c(0,28), breaks=c(0,4,8,12,16,20,24,28))+
  theme_minimal()+
  xlim(0,43)+
  ylab("Enrichment in ciliated neurons over other tissues")+
  xlab("Enrichment in IL2 over other neuronal types")+
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

ggsave("IL2_Enrich_over_CiliatedEnrich_v2.pdf", p2, width=10, height=10, 
       device = "pdf", units = "in")

