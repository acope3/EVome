library("tidyverse")
library("plyr")
library("qdap")
library("stringr")
library("DataCombine")
library("ggplot2")
library("gplots")

install.packages("Rcpp")
library("Rcpp")
library("devtools")
install_github("js229/Vennerable")
library("Vennerable")
library("gridExtra")
library("grid"))
library("ggplotify")
library("gridGraphics")
library("VennDiagram")
library("ggrepel")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("graph")

getwd()
setwd("/Users/Inna/Documents/Barr lab/EV proteomics/Combined Enrichment analysis")

MaleEnriched <- read.csv("Kim_et_al_comparison.csv")
MaleEnriched$Tissue <- as.factor(MaleEnriched$Tissue)

p1<-ggplot(MaleEnriched, aes(x=log2(foldChange), y=log2(Peptide.counts), color=Tissue), alpha=0.5)+ 
  geom_point(shape=17, size=1.5)+
  #scale_y_continuous(limits=c(0,28), breaks=c(0,4,8,12,16,20,24,28))+
  theme_minimal()+
  scale_color_manual(values = c("Glia" = "#04F60B", "Neurons" = "#D50167", 
                                "Reproductive system" = "grey", "NoSpecificTissue" = "00A0D7"))+
  
  #Neurons highlighting
  geom_point(data = MaleEnriched[(MaleEnriched$Tissue == "Neurons"), ],
             aes(x=log2(foldChange), y=log2(Peptide.counts)), color="#D50167", size=3)+
  geom_text_repel(data = MaleEnriched[(MaleEnriched$Tissue == "Neurons"), ],
                  aes(x=log2(foldChange), y=log2(Peptide.counts), label=Gene_name), size=5)+
  #xlim(0,43)+
  ylab("log2 (Peptide counts)")+
  xlab("log2 (Enrichment in males over hermaphrodites)")+
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


p2<-ggplot(MaleEnriched, aes(x=log2(foldChange), y=log2(Peptide.counts), color=Tissue), alpha=0.5)+ 
  geom_point(shape=17, size=1.5)+
  #scale_y_continuous(limits=c(0,28), breaks=c(0,4,8,12,16,20,24,28))+
  theme_minimal()+
  scale_color_manual(values = c("Glia" = "#04F60B", "Neurons" = "#D50167", 
                                "Reproductive system" = "grey", "NoSpecificTissue" = "00A0D7"))+
  
  #Glia highlighting
  geom_point(data = MaleEnriched[(MaleEnriched$Tissue == "Glia"), ],
             aes(x=log2(foldChange), y=log2(Peptide.counts)), color="#04F60B", size=3)+
  geom_text_repel(data = MaleEnriched[(MaleEnriched$Tissue == "Glia"), ],
                  aes(x=log2(foldChange), y=log2(Peptide.counts), label=Gene_name), size=5)+
  #xlim(0,43)+
  ylab("log2 (Peptide counts)")+
  xlab("log2 (Enrichment in males over hermaphrodites)")+
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

p3<-ggplot(MaleEnriched, aes(x=log2(foldChange), y=log2(Peptide.counts), color=Tissue), alpha=0.5)+ 
  geom_point(shape=17, size=1.5)+
  #scale_y_continuous(limits=c(0,28), breaks=c(0,4,8,12,16,20,24,28))+
  theme_minimal()+
  scale_color_manual(values = c("Glia" = "#04F60B", "Neurons" = "#D50167", 
                                "Reproductive system" = "grey", "NoSpecificTissue" = "00A0D7"))+
  
  #Glia highlighting
  geom_point(data = MaleEnriched[(MaleEnriched$Tissue == "NoSpecificTissue"), ],
             aes(x=log2(foldChange), y=log2(Peptide.counts)), color="#00A0D7", size=3)+
  geom_text_repel(data = MaleEnriched[(MaleEnriched$Tissue == "NoSpecificTissue"), ],
                  aes(x=log2(foldChange), y=log2(Peptide.counts), label=Gene_name), size=5)+
  #xlim(0,43)+
  ylab("log2 (Peptide counts)")+
  xlab("log2 (Enrichment in males over hermaphrodites)")+
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

ggsave("Male-specific Neuronal.pdf", p1, width=10, height=10, 
       device = "pdf", units = "in")

ggsave("Male-specific Glial.pdf", p2, width=10, height=10, 
       device = "pdf", units = "in")

ggsave("Male-specific NoSpecificTissue.pdf", p3, width=10, height=10, 
       device = "pdf", units = "in")
