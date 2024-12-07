#clear the R environment
rm(list=ls())

# Call library
library(dplyr)
library(tidyr)
library(ggplot2)
library(lattice)
library(caret)
library(e1071)


# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the required arguments are provided
if (length(args) != 2) {
  stop("Two arguments are required: an RDS file and a CSV file.")
}

# Assign arguments to variables
rds_file <- args[1]
csv_file <- args[2]

# Load the RDS file
svm_model <- readRDS(rds_file)

# Load the CSV file
new_data <- read.csv(csv_file)

# modify data
features <- new_data[, -which(names(new_data) == "Sample_id")] 
sample_labels <- new_data$Sample_id

# predict the data

# Generate predictions
predicted_labels <- predict(svm_model, newdata = features)

# Convert X1 to Yes and X0 to No
converted_labels <- ifelse(predicted_labels == "X1", "Yes", "No")

# Create a dataframe
df <- data.frame(Column1 = sample_labels, Column2 = converted_labels)

# Modify column names
colnames(df) <- c("Sample_id", "Senescence_status")

# Save the updated dataset
write.csv(df, "./predicted_results.csv", row.names = FALSE)




