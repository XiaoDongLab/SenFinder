###################################
###########Shamsed Mahmud##########
########22 January, 2024 ##########
###################################

#####################################
######## Direct Create the data set##
#####################################


###############################
### Call The library###########
###############################

#clear the R environment
rm(list=ls())

#library(tidyverse)
#library(biomaRt)
#call library
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)

###############################
## Read the GSE56293 dataset###
###############################

###############
#read GSE56293#
###############
dat_gse56293=read.table("./GSE56293/SRR1210368.genes.results",sep="\t",header=TRUE)
dat=as.data.frame(dat_gse56293$gene_id)
colnames(dat)="gene_id"

for(x in 368:373){
  
  new_d=read.table(paste0("./GSE56293/SRR1210",x,".genes.results"), sep="\t", header=TRUE)
  new_d=new_d[,c(1,6)]
  colnames(new_d)[2]<-paste0("SRR1210",x)
  dat=merge(dat,new_d,by="gene_id")
}
# remove the file
rm(dat_gse56293)
rm(new_d)

# Assign to dat2
dat2 <- dat

#remove the dat file
rm(dat)

########################
###read GSE 63577 1#####
########################
dat_gse63577=read.table("./GSE63577/GSE63577_1/SRR1660534.genes.results",sep="\t",header=TRUE)
dat=as.data.frame(dat_gse63577$gene_id)
colnames(dat)="gene_id"

for(x in 534:560){
  
  new_d=read.table(paste0("./GSE63577/GSE63577_1/SRR1660",x,".genes.results"), sep="\t", header=TRUE)
  new_d=new_d[,c(1,6)]
  colnames(new_d)[2]<-paste0("SRR1660",x)
  dat=merge(dat,new_d,by="gene_id")
}
# remove the file
rm(dat_gse63577)
rm(new_d)

# Assign to dat2
dat4 <- dat

#remove the dat file
rm(dat)


########################
###read GSE 63577 2#####
########################
dat_gse63577=read.table("./GSE63577/GSE63577_2/SRR2751110.genes.results",sep="\t",header=TRUE)
dat=as.data.frame(dat_gse63577$gene_id)
colnames(dat)="gene_id"

for(x in 110:127){
  
  new_d=read.table(paste0("./GSE63577/GSE63577_2/SRR2751",x,".genes.results"), sep="\t", header=TRUE)
  new_d=new_d[,c(1,6)]
  colnames(new_d)[2]<-paste0("SRR2751",x)
  dat=merge(dat,new_d,by="gene_id")
}
# remove the file
rm(dat_gse63577)
rm(new_d)

# Assign to dat2
dat5 <- dat

#remove the dat file
rm(dat)

########################
###read GSE 63577 3#####
########################
dat_gse63577=read.table("./GSE63577/GSE63577_3/SRR1560895.genes.results",sep="\t",header=TRUE)
dat=as.data.frame(dat_gse63577$gene_id)
colnames(dat)="gene_id"

for(x in 895:897){
  
  new_d=read.table(paste0("./GSE63577/GSE63577_3/SRR1560",x,".genes.results"), sep="\t", header=TRUE)
  new_d=new_d[,c(1,6)]
  colnames(new_d)[2]<-paste0("SRR1560",x)
  dat=merge(dat,new_d,by="gene_id")
}
# remove the file
rm(dat_gse63577)
rm(new_d)

# Assign to dat2
dat8 <- dat

#remove the dat file
rm(dat)

##################
#read GSE64553####
##################

dat_gse64553=read.table("./GSE64553/SRR1736309.genes.results",sep="\t",header=TRUE)
dat=as.data.frame(dat_gse64553$gene_id)
colnames(dat)="gene_id"

for(x in 309:368){
  
  new_d=read.table(paste0("./GSE64553/SRR1736",x,".genes.results"), sep="\t", header=TRUE)
  new_d=new_d[,c(1,6)]
  colnames(new_d)[2]<-paste0("SRR1736",x)
  dat=merge(dat,new_d,by="gene_id")
}
# remove the file
rm(dat_gse64553)
rm(new_d)

# Assign to dat2
dat6 <- dat

#remove the dat file
rm(dat)


##################
#read GSE EMTB ###
##################

dat_emtb=read.table("./GSEALLFILES_UTH/ERR1805188.genes.results",sep="\t",header=TRUE)
dat=as.data.frame(dat_emtb$gene_id)
colnames(dat)="gene_id"

for(x in 188:265){
  
  
    new_d=read.table(paste0("./GSEALLFILES_UTH/ERR1805",x,".genes.results"), sep="\t", header=TRUE)
    new_d=new_d[,c(1,6)]
    colnames(new_d)[2]<-paste0("ERR1805",x)
    dat=merge(dat,new_d,by="gene_id")
  
}
# remove the file
rm(dat_emtb)
rm(new_d)

# Assign to dat2
dat7 <- dat

#remove the dat file
rm(dat)

###############################################
############create the dataset#################
###############################################

#Get the first gene coloumn to get the common genes
col2<-as.data.frame(dat2[,1])
colnames(col2)[1]<-"Gene"

col4<-as.data.frame(dat4[,1])
colnames(col4)[1]<-"Gene"

col5<-as.data.frame(dat5[,1])
colnames(col5)[1]<-"Gene"

col6<-as.data.frame(dat6[,1])
colnames(col6)[1]<-"Gene"

col7<-as.data.frame(dat7[,1])
colnames(col7)[1]<-"Gene"

col8<-as.data.frame(dat8[,1])
colnames(col8)[1]<-"Gene"

##Check number of genes overlapping check less number here
#common genes
col<-as.data.frame(Reduce(intersect,list(col2$Gene,col4$Gene,col5$Gene,col6$Gene,col7$Gene,col8$Gene)))
colnames(col)<-"gene_id"

############make left join############
dat2_n<-merge(col,dat2,all.x= T)
dat4_n<-merge(col,dat4,all.x= T)
dat5_n<-merge(col,dat5,all.x= T)
dat6_n<-merge(col,dat6,all.x= T)
dat7_n<-merge(col,dat7,all.x= T)
dat8_n<-merge(col,dat8,all.x= T)

##########cbind all ag files###########
dat<-join_all(list(dat2_n,dat4_n,dat5_n,dat6_n,dat7_n,dat8_n),by="gene_id",type="inner")

#transform data
dat_t <-setNames(data.frame(t(dat[,-1])),dat[,1])
temp_list<-colnames(dat)
dat_t$sample_id=temp_list[-1]


# save the dataset
write.csv(dat_t,"./dat_set_wodv.csv")

# read created dataset
rm(list=ls())

# read the dataset
dat <- read.csv("./dat_set_wodv.csv")
colnames(dat)[1]="sample_id"

# read phenotype dataset
datn <- read.csv("./all_sample_new_orig11012023.csv")
datn<-datn[,c(2,3)]

#merge dataset and delete old senes colulmn
findat <- cbind(datn,dat)
findat <- findat[,-2]

# write the final data
write.csv(findat,"./check_findata.csv")


# Create final data
# remove data new_senescence with level 3
prework_dat <- findat[findat$new_senescence!=3,]

# remove data new_senescence with level 2
prework_lev_2 <- prework_dat[prework_dat$new_senescence==2,]

# save the data with level 2
write.csv(prework_lev_2,"./check_unknown_dat.csv")

# remove level 2 from the dataset
final_dat <- prework_dat[prework_dat$new_senescence!=2,]

# save the final data
write.csv(final_dat,"./final_dat.csv")


