rm(list = ls())

library("dplyr")
getwd()
setwd("/Users/Inna/Documents/Barr lab/EV proteomics/Combined Enrichment analysis/Domain Enrichment")

#Import doamin info for the whole proteome (sourced from Ensembl Biomart)
ElegansGenomeDomainInfo <- read.csv("ElegansGenomeInterProInfo.csv")

#Kill empty InterProID entries
ElegansGenomeDomainInfo <- ElegansGenomeDomainInfo[-which(ElegansGenomeDomainInfo$Interpro.ID == ""), ]

#Extarct domain IDs
GlobalDomainList <- unique(ElegansGenomeDomainInfo$Interpro.ID)

#Extract unique genes to calculate the size of the proteome
GlobalGeneList <- unique(ElegansGenomeDomainInfo$Gene.stable.ID)
length(GlobalGeneList)

#Calculate frequency of occurence for each domain
a <- rle(sort(ElegansGenomeDomainInfo$Interpro.ID))
BackgroundFrequencies <- data.frame(InterProDomain=a$values, BackgroundOccurence=a$lengths)

#Import doamin info for the EVome (sourced from Ensembl Biomart)
EVomeDomainInfo <- read.csv("EVomeDomainInfo.csv")
#Kill empty InterProID entries
EVomeDomainInfo <- EVomeDomainInfo[-which(EVomeDomainInfo$Interpro.ID == ""), ]

#Extract unique EV proteins to calculate the size of the EVome with assigned domains
EVomeDomainProteinList <- unique(EVomeDomainInfo$Gene.stable.ID)
length(EVomeDomainProteinList)

#Calculate frequency of occurence for each domain
b <- rle(sort(EVomeDomainInfo$Interpro.ID))
EVomeFrequencies <- data.frame(InterProDomain=b$values, EVomeOccurence=b$lengths)

#Generate data frame of EVome domain IDs with EVome and Background frequencies
Summary <- merge(EVomeFrequencies, BackgroundFrequencies, 
                                        by = "InterProDomain", 
                                        all.x=TRUE, all.y=FALSE)
Summary$EVomeFrequency <- Summary$EVomeOccurence/2888
Summary$ProteomeFrequency <- Summary$BackgroundOccurence/20191
Summary$Enrichment <- Summary$EVomeFrequency/Summary$ProteomeFrequency

#Calculate p-value using hypergeometric distribution model
#EVome size = 2,888 proteins
#Proteome size = 14,904 proteins
#http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
#phyper(success-in-sample, success-in-bkgd, failure-in-bkgd, sample-size)

Summary$p_value <- 1-phyper(Summary$EVomeOccurence-1, 
                                     Summary$BackgroundOccurence,
                            20191-Summary$BackgroundOccurence,
                            2888)

#Holm-Bonferroni method to calculate p-adjusted
Summary$padj <- p.adjust(Summary$p_value, method="holm") 
Summary$neglog10_padj <- -log10(Summary$padj) 


library("dplyr")

SpC_byDomain <- EVomeDomainInfo %>%
  group_by(Interpro.ID) %>%
  summarise(SumCounts = sum(SpectralCounts))

Summary <- merge(Summary, SpC_byDomain, 
                 by.x="InterProDomain", by.y="Interpro.ID",
                 all.x=TRUE, all.y=FALSE)

Summary <- merge(Summary, EVomeDomainInfo[,c("Interpro.ID","Interpro.Description")], 
                 by.x="InterProDomain", by.y="Interpro.ID",
                 all.x=TRUE, all.y=FALSE)
Summary <- unique(Summary)

write.csv(Summary, "SummaryDomainEnrichment.csv")





