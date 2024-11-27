# remove the list
rm(list=ls())

# call the library
library(dplyr)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db) 


#setwd("C://Users/mahmu052/Documents/UMN/PROJ_MACHINE/RNA_SEQ_RUN/data_01222024/simrun_sen_unix/")
# read the cell age gene
dat_sen <- read.csv("./senMayo.csv")

# Take the gene column
dat_sen_gene <- as.data.frame(dat_sen$Gene.human.)

# Change column name and unique gene
colnames(dat_sen_gene)[1]<- "gene"
dat_uniq_sen_gene <- as.data.frame(dat_sen_gene$gene %>%unique())
colnames(dat_uniq_sen_gene)[1]<- "gene"

#############################
#####Annotationdbi###########
#############################
# Example with a list of gene symbols
gene_symbols <- dat_uniq_sen_gene$gene
ensembl_ids <- select(org.Hs.eg.db, keys = gene_symbols, columns = "ENSEMBL", keytype = "SYMBOL")
ensemble_ids_uniq <- unique(ensembl_ids$ENSEMBL)

# read the final data
dat_final <- read.csv("./final_dat.csv")
dat_final<-dat_final[,-1]


# Common genes between cellage and final data
fin_gene <- colnames(dat_final[,-c(1,2)])
sen_gene <- ensemble_ids_uniq
common_genes <- intersect(sen_gene,fin_gene)


# Take the first part
dat_part1<- dat_final[,c(1,2)]

# Take the second part
dat_part2 <- dat_final[,common_genes]

# Create the final dataset
dat_sen_final <- cbind(dat_part1, dat_part2)

# save the final data set for cell age
write.csv(dat_sen_final,"./senmayo_final_dat.csv")

