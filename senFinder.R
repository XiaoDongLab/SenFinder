# Call the library

#clear the R environment
rm(list=ls())

# Call library
library(dplyr)
library(tidyr)
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
features <- new_data[, -which(names(new_data) == "senescence")] 
actual_labels <- new_data$senescence

# predict the data

# Generate predictions
predicted_labels <- predict(svm_model, newdata = features)


# Confusion Matrix
conf_matrix <- table(Predicted = predicted_labels, Actual = actual_labels)

# Calculate Accuracy
accuracy <- (conf_matrix[1, 1] + conf_matrix[2, 2]) / sum(conf_matrix)


colnames(conf_matrix) <- c("Predicted_Neg", "Predicted_Pos")
rownames(conf_matrix) <- c("Actual_Neg", "Actual_Pos")



# Output file
output_file <- "classification_metrics.txt"

# Redirect output to file
sink(output_file)

# Print results
cat("Accuracy:", accuracy, "\n")
cat("Confusion Matrix:\n")
print(conf_matrix)

# Close the connection
sink()

# Add predictions to the original data
new_data$Predicted_SenescenceStatus <- predicted_labels

# Save the updated dataset
write.csv(new_data, "predicted_results.csv", row.names = FALSE)




