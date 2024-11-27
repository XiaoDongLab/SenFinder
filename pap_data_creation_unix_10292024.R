# remove the list
rm(list=ls())

# call the library
library(dplyr)
library(tidyr)


#############################
#### Pap data creation ######
#############################
# read pap gene
dat_pap <- read.csv("./Machine_project_genelist.csv")
dat_pap_gene <- dat_pap$Ensemble.gene.id


# read the final data
dat_final <- read.csv("./final_dat.csv")
dat_final<-dat_final[,-1]


# Common genes between cellage and final data
fin_gene <- colnames(dat_final[,-c(1,2)])
pap_gene <- dat_pap_gene
common_genes <- intersect(fin_gene,pap_gene)


# Take the first part
dat_part1<- dat_final[,c(1,2)]

# Take the second part
dat_part2 <- dat_final[,common_genes]

# Create the final dataset
dat_pap_final <- cbind(dat_part1, dat_part2)

# save the final data set for cell age
write.csv(dat_pap_final,"./pap_final_dat_10212024.csv")


######################################################################
######## Comparing heatmap of 114 samples and 4 methods prediction####
######################################################################

library(tidyverse)

set.seed(123)

dat_pap_final <- read.csv("./pap_final_dat_10212024.csv")
# Exclude the senescence status and sample number for PCA
pca_data <- dat_pap_final[, 3:22]

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Extract PCA scores (first two principal components)
pca_scores <- as.data.frame(pca_result$x[, 1:2])

# Add the senescence status and sample number back to the PCA scores
pca_scores$Senescence_Status <- dat_pap_final$new_senescence
pca_scores$Sample_Number <- as.factor(dat_pap_final$sample_no)

pca_scores <- pca_scores %>%
  mutate(Senescence_Status = recode(Senescence_Status, `0` = "Non-sen", `1` = "Sen"))

pdf("./pca_plot_sample_10292024.pdf")
# Plot the PCA result and add axis labels
ggplot(pca_scores, aes(x = PC1, y = PC2, color = as.factor(Senescence_Status), shape = Sample_Number)) +
  geom_point(size = 4) +
  labs(title = "PCA of Samples with Senescence Status",
       x = "Principal Component 1 (PC1)",  # X-axis label
       y = "Principal Component 2 (PC2)",  # Y-axis label
       color = "Senescence Status",
       shape = "Sample Number") +
  theme_minimal()+
  scale_color_manual(values = c("Sen" = "red", "Non-sen" = "green")) +
  theme(  # Add axis lines back
    axis.line = element_line(color = "black"),  # X and Y axis lines
    axis.ticks = element_line(color = "black"), # Add axis ticks
    panel.grid = element_blank()                # Remove grid lines for clarity
  )
dev.off()


######################################################
######### TSNE PLOT ##################################
######################################################
library(Rtsne)

# take tsne data
tsne_data <- dat_pap_final[, 3:22]

# Perform t-SNE
set.seed(123)  # Set seed for reproducibility
tsne_result <- Rtsne(tsne_data, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)

# Extract t-SNE coordinates (2D)
tsne_scores <- as.data.frame(tsne_result$Y)
colnames(tsne_scores) <- c("Dim1", "Dim2")

# Add the senescence status and sample number back to the t-SNE scores
tsne_scores$Senescence_Status <- dat_pap_final$new_senescence
tsne_scores$Sample_Number <- as.factor(dat_pap_final$sample_no)


pdf("./tsne_plot_sample.pdf")
# Plot the t-SNE result and color by senescence status and shape by sample number
ggplot(tsne_scores, aes(x = Dim1, y = Dim2, color = as.factor(Senescence_Status), shape = Sample_Number)) +
  geom_point(size = 4) +
  labs(title = "t-SNE of Samples with Senescence Status",
       x = "t-SNE Dimension 1",  # X-axis label
       y = "t-SNE Dimension 2",  # Y-axis label
       color = "Senescence Status",
       shape = "Sample Number") +
  theme_minimal() +
  theme(  # Add axis lines back
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank()
  )
dev.off()

######################################################
######### UMAP PLOT ##################################
######################################################


# Load necessary packages
library(umap)

# Exclude the senescence status and sample number from UMAP input data
umap_data <- dat_pap_final[, 3:22]

# Perform UMAP
umap_result <- umap(umap_data)

# Extract UMAP coordinates (2D)
umap_scores <- as.data.frame(umap_result$layout)
colnames(umap_scores) <- c("Dim1", "Dim2")

# Add the senescence status and sample number back to the UMAP scores
umap_scores$Senescence_Status <- dat_pap_final$new_senescence
umap_scores$Sample_Number <- as.factor(dat_pap_final$sample_no)

pdf("./umap_plot_sample.pdf")
# Plot the UMAP result and color by senescence status and shape by sample number
ggplot(umap_scores, aes(x = Dim1, y = Dim2, color = as.factor(Senescence_Status), shape = Sample_Number)) +
  geom_point(size = 4) +
  labs(title = "UMAP of Samples with Senescence Status",
       x = "UMAP Dimension 1",  # X-axis label
       y = "UMAP Dimension 2",  # Y-axis label
       color = "Senescence Status",
       shape = "Sample Number") +
  theme_minimal() +
  theme(  # Add axis lines back
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank()
  )

dev.off()



#############################################
####### Heatmap for 114 samples #############
#############################################

# load package
library(pheatmap)
library(tidyr)
library(dplyr)
library("dendextend")


# read dataset
dat_pap_final <- read.csv("./pap_final_dat_10212024_mod.csv")

# the dataset is
heatmap_dat <- dat_pap_final[,4:23]
heatmap_dat_f <- t(heatmap_dat)
colnames(heatmap_dat_f)<-dat_pap_final$sample_name

# log10 expression of data
heatmap_dat_log10 <- log10(heatmap_dat_f+1)



# Sample status
samp_stat <- dat_pap_final$new_senescence
sen_stat <- ifelse(samp_stat==1,"Sen","Non-sen")
my_sample_col <-  data.frame(sample=sen_stat)
row.names(my_sample_col) <- dat_pap_final$sample_name

# Define colors for each group
group_colors <- c("Non-sen" = "green", "Sen" = "red")
my_sample_col$Sample <- factor(my_sample_col$sample, levels = unique(my_sample_col$sample))
my_sample_col <- data.frame(my_sample_col[,-1])
row.names(my_sample_col) <- colnames(heatmap_dat_log10)
colnames(my_sample_col)[1] <-"Sample" 

custom_colors <- colorRampPalette(c("deepskyblue4", "white", "orangered3"))(50)
breaks <- seq(-2, 2, length.out = 51)

pdf(file="./allsample_heatmap_log10scale.pdf",width=7,height=8)
pheatmap(t(heatmap_dat_log10),
         annotation_row = my_sample_col,
         cutree_cols = 2,
         breaks = breaks,
         color=custom_colors,
         scale = "column",
         annotation_colors = list(Sample=group_colors),
         fontsize_row = 5, 
         cluster_rows=FALSE)

dev.off()

pdf(file="./allsample_heatmap_log10scale_row.pdf",width=7,height=8)
pheatmap(t(heatmap_dat_log10),
         annotation_row = my_sample_col,
         cutree_cols = 2,
         breaks = breaks,
         color=custom_colors,
         scale = "column",
         annotation_colors = list(Sample=group_colors),
         fontsize_row = 4, 
         cluster_rows=TRUE)

dev.off()

