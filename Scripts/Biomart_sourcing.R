
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library("biomaRt")
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useMart("ensembl",dataset="celegans_gene_ensembl")
filters=listFilters(ensembl)
attributes = listAttributes(ensembl)
attributes
list

output<-c("wormbase_cds","uniprotswissprot", "uniprotsptrembl", "wormbase_gene","external_gene_name", "description")
input<-c("transcript_biotype")
Tada<-getBM(
  attributes=output,
  filters=input,
  values="protein_coding",
  mart=ensembl)

head(Tada, 40)



seq = getSequence(id = Tada$wormbase_cds, 
                  type = "wormbase_cds", 
                  seqType = "peptide", 
                  mart = ensembl)

head(seq)

tada <- merge(Tada, seq , by="wormbase_cds")
head(tada,1)

output<-c("uniparc", "wormbase_gene","wormbase_cds")
input<-c("transcript_biotype")
Tada<-getBM(
  attributes=output,
  filters=input,
  values="protein_coding",
  mart=ensembl)

write.csv(tada, file="Proteome_info.csv")
write.csv(Tada, file="Uniprot_IDs.csv")

head(CAEEL_proteome_info)
head(Tada)

Proteome_uniprot_IDs <- merge(Tada, CAEEL_proteome_info , by="wormbase_cds",
              all.x=FALSE, all.y=TRUE)
              
write.csv(Proteome_uniprot_IDs, file="Proteome_uniprot_IDs.csv")

install.packages("UniprotR")
library("UniprotR")

#Matching Haiyan's uniprot IDs to the Wormbase transcript CDSs
cds_vector <- read.csv("Uniprot_identified.csv")


ConvertID(cds_vector$Uniprot, ID_from = "ACC" , ID_to = "WORMBASE_TRS_ID", 
          directorypath = "/Users/Inna/Documents/Barr lab/EV proteomics/Combined peptide identification")

DB_IDs<-read.csv("Database identifiers.csv")
DB_IDs_transcript_vector<-unlist(strsplit(DB_IDs$To.WORMBASE_TRS_ID, split=", "))


#Locate cds_ID in Uniprot identifiers
DB_L<-list()
for (i in 1:length(DB_IDs_transcript_vector)) { 
  
  temp_vector <- grep(paste(DB_IDs_transcript_vector[i]), DB_IDs$To.WORMBASE_TRS_ID)
  
  DB_L<-lappend(DB_L, temp_vector)
  names(DB_L)[i] <- paste0(DB_IDs_transcript_vector[i])
  
}
head(DB_L)

#subsetting Uniprot names for each wormbase cds 
DB_L2<-list()
for (i in 1:length(DB_L)) { 
  
  x<-DB_L[[i]]
  temp_vector <- DB_IDs$From_UniProtKB_AC_ID[x]
  
  # temp_vector<-do.call(c, temp_vector) #concatenate several vectors into one vector
  
  temp_vector<-unique(temp_vector)
  
  DB_L2<-lappend(DB_L2, temp_vector)
  
  
}

DB_IDs_final <- do.call(rbind,mapply(cbind, DB_IDs_transcript_vector, DB_L2))
DB_IDs_final<-as.data.frame(DB_IDs_final)
names(DB_IDs_final)<-c("wormbase_cds", "uniprot")
head(DB_IDs_final)

write.csv(DB_IDs_final, "wormbase cds-to-uniprot.csv")

DB_uniprot_cds_genes <- merge(DB_IDs_final, CAEEL_proteome_info, by="wormbase_cds",
                              all.x=TRUE, all.y=FALSE)
DB_uniprot_cds_genes$X <- NULL
head(DB_uniprot_cds_genes)
