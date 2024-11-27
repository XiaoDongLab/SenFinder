
#################################################
############ Xiao Comments 10262024##############
#################################################

####### Call library#############

#clear the R environment
rm(list=ls())

# Call library
library(tidyverse)
library(caret)
library(e1071)
library(caret)
library(mRMRe)
library(randomForest)
library(rpart)
library(rrcovHD)
library(iml)

# read the paper 20 gene set data
dat <- read.csv("./pap_final_dat_10212024.csv")

# remove the sample column
dat <- dat[,-2]

# Train with svm with full data
# Train an SVM model using the full dataset

set.seed(123)

levels(dat$new_senescence) <- make.names(levels(dat$new_senescence))
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE)

svm_model <- train(as.factor(new_senescence) ~ ., data = dat, method = "svmLinear", trControl = train_control)



# subset sen and non sen data
dat$new_senescence <- as.factor(dat$new_senescence)

# Split your data into senescent and non-senescent samples
sen_dat <- dat[dat$new_senescence == 1, ]
# now remove the sen status column
sen_dat <- sen_dat[,-1]

non_sen_dat <- dat[dat$new_senescence == 0, ]
# now remove the non-sen status column
non_sen_dat <- non_sen_dat[,-1]

##############################################
### Now create senthetic data and predict ####
##############################################

# Initialize storage for sen counts
sen_counts <- data.frame(weight = numeric(), sen_count = integer())

# Loop over the range of proportions
for (w in seq(0.05, 1, by = 0.05)) {
  # Count sen predictions for this weight
  sen_count_for_weight <- 0
  
  for (i in 1:1000) {
    # Randomly select one sample from each group
    non_sen_sample <- non_sen_dat[sample(1:nrow(non_sen_dat), 1), ]
    sen_sample <- sen_dat[sample(1:nrow(sen_dat), 1), ]
    
    # Calculate the synthetic sample
    synthetic_sample <- (1-w) * non_sen_sample + w * sen_sample
    
    
    # Predict using the SVM model
    prediction <- predict(svm_model, as.data.frame(synthetic_sample))
    
    # Check if prediction is "sen" and increment counter if true
    if (prediction == "X1") {
      sen_count_for_weight <- sen_count_for_weight + 1
    }
  }
  
  # Store the sen count for the current weight
  sen_counts <- rbind(sen_counts, data.frame(weight = w, sen_count = sen_count_for_weight))
}

# Check the sen counts across weights
print(sen_counts)

write.csv(sen_counts,"svm_sen_counts.csv")

## Draw the plot
# Calculate the percentage of "sen" counts for each weight
sen_counts$sen_percentage <- (sen_counts$sen_count / 1000) 


pdf("svm_weight_proportion.pdf")
# Create the plot
ggplot(sen_counts, aes(x = weight, y = sen_percentage)) +
  geom_line(color = "blue",size=1.5) +                     # Line for percentage of sen counts
  geom_point(color = "blue", size=2.5) +                    # Points at each increment
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "red",size=1.5) +  # 50% reference line at 0.5
  labs(x = "Weight (from 0.0 to 1 by 0.05)", 
       y = "Proportion of Sen Counts", 
       title = "Sen Count Proportion vs. Weight") +
  theme_minimal() +
  theme(
    axis.line = element_line(size = 1.5, color = "black"),  # Thicker axis lines
    panel.grid = element_blank(),                           # Remove grid lines
    axis.ticks = element_line(size = 1.5)                   # Thicker ticks on axes
  ) +
  scale_x_continuous(breaks = seq(0.0, 1, by = 0.05), limits = c(0.0, 1)) +  # X-axis ticks from 0.0 to 1
  scale_y_continuous(breaks = seq(0.0, 1, by = 0.05), limits = c(0.0, 1))    # Y-axis ticks from 0.0 to 1
dev.off()

#######################################
###############RF######################


# subset sen and non sen data
dat$new_senescence <- as.factor(dat$new_senescence)

# Split your data into senescent and non-senescent samples
sen_dat <- dat[dat$new_senescence == 1, ]
# now remove the sen status column
sen_dat <- sen_dat[,-1]

non_sen_dat <- dat[dat$new_senescence == 0, ]
# now remove the non-sen status column
non_sen_dat <- non_sen_dat[,-1]

# Train model
set.seed(123)
dat$new_senescence <- as.factor(dat$new_senescence)
levels(dat$new_senescence) <- make.names(levels(dat$new_senescence))
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE)
rf_model <- train(as.factor(new_senescence) ~ ., data = dat, method = "rf", trControl = train_control)

##############################################
### Now create senthetic data and predict ####
##############################################

# Initialize storage for sen counts
sen_counts <- data.frame(weight = numeric(), sen_count = integer())

# Loop over the range of proportions
for (w in seq(0.05, 1, by = 0.05)) {
  # Count sen predictions for this weight
  sen_count_for_weight <- 0
  
  for (i in 1:1000) {
    # Randomly select one sample from each group
    non_sen_sample <- non_sen_dat[sample(1:nrow(non_sen_dat), 1), ]
    sen_sample <- sen_dat[sample(1:nrow(sen_dat), 1), ]
    
    # Calculate the synthetic sample
    synthetic_sample <- (1-w) * non_sen_sample + w * sen_sample
    
    
    # Predict using the RF model
    prediction <- predict(rf_model, as.data.frame(synthetic_sample))
    
    # Check if prediction is "sen" and increment counter if true
    if (prediction == "X1") {
      sen_count_for_weight <- sen_count_for_weight + 1
    }
  }
  
  # Store the sen count for the current weight
  sen_counts <- rbind(sen_counts, data.frame(weight = w, sen_count = sen_count_for_weight))
}

# Check the sen counts across weights
print(sen_counts)

write.csv(sen_counts,"rf_sen_counts.csv")

## Draw the plot
# Calculate the percentage of "sen" counts for each weight
sen_counts$sen_percentage <- (sen_counts$sen_count / 1000) 


pdf("rf_weight_proportion.pdf")
# Create the plot
ggplot(sen_counts, aes(x = weight, y = sen_percentage)) +
  geom_line(color = "blue",size=1.5) +                     # Line for percentage of sen counts
  geom_point(color = "blue", size=2.5) +                    # Points at each increment
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "red",size=1.5) +  # 50% reference line at 0.5
  labs(x = "Weight (from 0.0 to 1 by 0.05)", 
       y = "Proportion of Sen Counts", 
       title = "Sen Count Proportion vs. Weight") +
  theme_minimal() +
  theme(
    axis.line = element_line(size = 1.5, color = "black"),  # Thicker axis lines
    panel.grid = element_blank(),                           # Remove grid lines
    axis.ticks = element_line(size = 1.5)                   # Thicker ticks on axes
  ) +
  scale_x_continuous(breaks = seq(0.0, 1, by = 0.05), limits = c(0.0, 1)) +  # X-axis ticks from 0.0 to 1
  scale_y_continuous(breaks = seq(0.0, 1, by = 0.05), limits = c(0.0, 1))    # Y-axis ticks from 0.0 to 1
dev.off()

###############################
###### Desicion Tree ##########
###############################


# subset sen and non sen data
dat$new_senescence <- as.factor(dat$new_senescence)

# Split your data into senescent and non-senescent samples
sen_dat <- dat[dat$new_senescence == 1, ]
# now remove the sen status column
sen_dat <- sen_dat[,-1]

non_sen_dat <- dat[dat$new_senescence == 0, ]
# now remove the non-sen status column
non_sen_dat <- non_sen_dat[,-1]

# Train model
set.seed(123)
#dat$new_senescence <- as.factor(dat$new_senescence)
levels(dat$new_senescence) <- make.names(levels(dat$new_senescence))
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE)
dt_model <- train(as.factor(new_senescence) ~ ., data = dat, method = "rpart", trControl = train_control)

##############################################
### Now create senthetic data and predict ####
##############################################

# Initialize storage for sen counts
sen_counts <- data.frame(weight = numeric(), sen_count = integer())

# Loop over the range of proportions
for (w in seq(0.05, 1, by = 0.05)) {
  # Count sen predictions for this weight
  sen_count_for_weight <- 0
  
  for (i in 1:1000) {
    # Randomly select one sample from each group
    non_sen_sample <- non_sen_dat[sample(1:nrow(non_sen_dat), 1), ]
    sen_sample <- sen_dat[sample(1:nrow(sen_dat), 1), ]
    
    # Calculate the synthetic sample
    synthetic_sample <- (1-w) * non_sen_sample + w * sen_sample
    
    
    # Predict using the RF model
    prediction <- predict(dt_model, as.data.frame(synthetic_sample))
    
    # Check if prediction is "sen" and increment counter if true
    if (prediction == "X1") {
      sen_count_for_weight <- sen_count_for_weight + 1
    }
  }
  
  # Store the sen count for the current weight
  sen_counts <- rbind(sen_counts, data.frame(weight = w, sen_count = sen_count_for_weight))
}

# Check the sen counts across weights
print(sen_counts)

write.csv(sen_counts,"dt_sen_counts.csv")

## Draw the plot
# Calculate the percentage of "sen" counts for each weight
sen_counts$sen_percentage <- (sen_counts$sen_count / 1000) 


pdf("dt_weight_proportion.pdf")
# Create the plot
ggplot(sen_counts, aes(x = weight, y = sen_percentage)) +
  geom_line(color = "blue",size=1.5) +                     # Line for percentage of sen counts
  geom_point(color = "blue", size=2.5) +                    # Points at each increment
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "red",size=1.5) +  # 50% reference line at 0.5
  labs(x = "Weight (from 0.0 to 1 by 0.05)", 
       y = "Proportion of Sen Counts", 
       title = "Sen Count Proportion vs. Weight") +
  theme_minimal() +
  theme(
    axis.line = element_line(size = 1.5, color = "black"),  # Thicker axis lines
    panel.grid = element_blank(),                           # Remove grid lines
    axis.ticks = element_line(size = 1.5)                   # Thicker ticks on axes
  ) +
  scale_x_continuous(breaks = seq(0.0, 1, by = 0.05), limits = c(0.0, 1)) +  # X-axis ticks from 0.0 to 1
  scale_y_continuous(breaks = seq(0.0, 1, by = 0.05), limits = c(0.0, 1))    # Y-axis ticks from 0.0 to 1
dev.off()

#############################################
########   RSimca ###########################
#############################################

# subset sen and non sen data
dat$new_senescence <- as.factor(dat$new_senescence)

# Split your data into senescent and non-senescent samples
sen_dat <- dat[dat$new_senescence == 1, ]
# now remove the sen status column
sen_dat <- sen_dat[,-1]

non_sen_dat <- dat[dat$new_senescence == 0, ]
# now remove the non-sen status column
non_sen_dat <- non_sen_dat[,-1]



# Set the seed for reproducibility
set.seed(123)

#dat$new_senescence <- as.factor(dat$new_senescence)
levels(dat$new_senescence) <- make.names(levels(dat$new_senescence))

# Create 10-fold cross-validation
folds <- createFolds(dat$new_senescence, k = 10, list = TRUE)

# Initialize storage for predictions
fold_models <- list()

# Loop through each fold for cross-validation and store models
for (i in seq_along(folds)) {
  # Split data into training and test for the current fold
  train_index <- folds[[i]]
  train_data <- dat[train_index, ]
  
  # Train the RSIMCA model on the training set
  rsimca_model <- RSimca(
    as.matrix(train_data[, -which(names(train_data) == "new_senescence")]), 
    grouping = train_data$new_senescence
  )
  
  # Store the model for later use
  fold_models[[i]] <- rsimca_model
}

# Initialize storage for sen counts
sen_counts <- data.frame(weight = numeric(), sen_count = integer())

# Loop over the range of proportions from 0.05 to 1 with increments of 0.05
for (w in seq(0.05, 1, by = 0.05)) {
  sen_count_for_weight <- 0
  
  for (i in 1:1000) {
    # Randomly select one sample from each group
    non_sen_sample <- non_sen_dat[sample(1:nrow(non_sen_dat), 1), ]
    sen_sample <- sen_dat[sample(1:nrow(sen_dat), 1), ]
    
    # Calculate the synthetic sample
    synthetic_sample <- (1-w) * non_sen_sample + w * sen_sample
   
    
    # Use each fold model to predict and take majority voting for the class
    predictions <- sapply(fold_models, function(model) {
      pred <- predict(model, as.matrix(synthetic_sample))
      pred@classification
    })
    
    # Majority vote on the class predictions from the 10 models
    majority_class <- ifelse(mean(predictions == "X1") > 0.5, "sen", "non_sen")
    
    # Count as "sen" if the majority class is "sen"
    if (majority_class == "sen") {
      sen_count_for_weight <- sen_count_for_weight + 1
    }
  }
  
  # Store the sen count for the current weight
  sen_counts <- rbind(sen_counts, data.frame(weight = w, sen_count = sen_count_for_weight))
}

# Check the sen counts across weights
print(sen_counts)

write.csv(sen_counts,"RSimca_sen_counts.csv")

## Draw the plot
# Calculate the percentage of "sen" counts for each weight
sen_counts$sen_percentage <- (sen_counts$sen_count / 1000) 


pdf("rsimca_weight_proportion.pdf")
# Create the plot
ggplot(sen_counts, aes(x = weight, y = sen_percentage)) +
  geom_line(color = "blue",size=1.5) +                     # Line for percentage of sen counts
  geom_point(color = "blue", size=2.5) +                    # Points at each increment
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "red",size=1.5) +  # 50% reference line at 0.5
  labs(x = "Weight (from 0.0 to 1 by 0.05)", 
       y = "Proportion of Sen Counts", 
       title = "Sen Count Proportion vs. Weight") +
  theme_minimal() +
  theme(
    axis.line = element_line(size = 1.5, color = "black"),  # Thicker axis lines
    panel.grid = element_blank(),                           # Remove grid lines
    axis.ticks = element_line(size = 1.5)                   # Thicker ticks on axes
  ) +
  scale_x_continuous(breaks = seq(0.0, 1, by = 0.05), limits = c(0.0, 1)) +  # X-axis ticks from 0.0 to 1
  scale_y_continuous(breaks = seq(0.0, 1, by = 0.05), limits = c(0.0, 1))    # Y-axis ticks from 0.0 to 1
dev.off()


###############################################
###### Feature importance #####################
###############################################

###############################
##### Use iml packages#########
###############################

###############################################################
#### Train SVM with 10-fold cross validation###################
######### with 114 samples and calculate feature impportance###
###############################################################


#######################
##### SVM #############
#######################
set.seed(123)

# subset sen and non sen data
dat$new_senescence <- as.factor(dat$new_senescence)
levels(dat$new_senescence) <- make.names(levels(dat$new_senescence))
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE)
svm_model <- train(as.factor(new_senescence) ~ ., data = dat, method = "svmLinear", trControl = train_control)

# Create a predictor object for iml
predictor <- Predictor$new(
  model = svm_model,
  data = dat[, -which(names(dat) == "new_senescence")],  # Exclude target variable
  y = dat$new_senescence
)


# Calculate feature importance using permutation method
feature_importance <- FeatureImp$new(predictor, loss = "ce")  # "ce" is for classification error

# Plot feature importance
feature_importance$plot()

importance_data <- feature_importance$results

write.csv(importance_data,"./svm_model_feature_importance.csv")

#################################
##########Random Forest #########
#################################
# Load necessary libraries
library(randomForest)
library(iml)

# Step 1: Train the Random Forest Model
# Train model
set.seed(123)
# Ensure target variable is a factor
dat$new_senescence <- as.factor(dat$new_senescence)
levels(dat$new_senescence) <- make.names(levels(dat$new_senescence))

# Define 10-fold cross-validation control
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE)

# Train the Random Forest model using caret
rf_model <- train(
  new_senescence ~ ., 
  data = dat, 
  method = "rf", 
  trControl = train_control
)

# Step 2: Create a Predictor object for iml
predictor_rf <- Predictor$new(
  model = rf_model,
  data = dat[, -which(names(dat) == "new_senescence")],  # Exclude target variable
  y = dat$new_senescence  # Target variable
)

# Step 3: Calculate feature importance using permutation method
feature_importance_rf <- FeatureImp$new(predictor_rf, loss = "ce")  # "ce" for classification error

# Step 4: Plot and view feature importance
feature_importance_rf$plot()
importance_results <- feature_importance_rf$results  # Access raw importance results
print(importance_results)

write.csv(importance_results,"./rf_model_feature_importance.csv")


#######################
#### Desicion Tree ####
#######################

# Train model
set.seed(123)
dat$new_senescence <- as.factor(dat$new_senescence)
levels(dat$new_senescence) <- make.names(levels(dat$new_senescence))
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE)
dt_model <- train(as.factor(new_senescence) ~ ., data = dat, method = "rpart", trControl = train_control)



# Step 2: Create a Predictor object for iml
predictor_dt <- Predictor$new(
  model = dt_model,
  data = dat[, -which(names(dat) == "new_senescence")],  # Exclude target variable
  y = dat$new_senescence  # Target variable
)

# Step 3: Calculate feature importance using permutation method
feature_importance_dt <- FeatureImp$new(predictor_dt, loss = "ce")  # "ce" for classification error

# Step 4: Plot and view feature importance
feature_importance_dt$plot()
importance_results <- feature_importance_dt$results  # Access raw importance results
print(importance_results)

write.csv(importance_results,"./dt_model_feature_importance.csv")

##############################
########RSimca ###############
##############################

# Set the seed for reproducibility
set.seed(123)

# subset sen and non sen data
dat$new_senescence <- as.factor(dat$new_senescence)

#dat$new_senescence <- as.factor(dat$new_senescence)
levels(dat$new_senescence) <- make.names(levels(dat$new_senescence))

# Create 10-fold cross-validation
folds <- createFolds(dat$new_senescence, k = 10, list = TRUE)


custom_predict <- function(model, newdata) {
  # Predict with RSimca and get classifications
  predictions <- predict(model, as.matrix(newdata))
  
  # Convert classifications to factor format, as required by iml
  return(factor(predictions@classification, levels = c("X1", "X0")))
}

# Initialize storage for feature importances
feature_importances <- data.frame()

# Loop through each fold for cross-validation, store models, and calculate feature importance
for (i in seq_along(folds)) {
  # Split data into training set for the current fold
  train_index <- folds[[i]]
  train_data <- dat[train_index, ]
  
  # Train the RSIMCA model on the training set
  rsimca_model <- RSimca(
    as.matrix(train_data[, -which(names(train_data) == "new_senescence")]), 
    grouping = train_data$new_senescence
  )
  
  # Use iml package for feature importance with lower and upper bounds
  predictor <- Predictor$new(
    model = rsimca_model, 
    data = train_data[, -which(names(train_data) == "new_senescence")], 
    y = train_data$new_senescence,
    predict.function = custom_predict  # Use custom prediction function
  )
  
  # Calculate feature importance with permutations
  importance <- FeatureImp$new(predictor, loss = "ce")  # "ce" for classification error, adjust as needed
  
  # Extract results
  importance_results <- data.frame(
    fold = i,
    feature = importance$results$feature,
    importance = importance$results$importance,
    lower_bound = importance$results$importance.05,
    upper_bound = importance$results$importance.95,
    permutation_error= importance$results$permutation.error
  )
  
  # Store importance results for this fold
  feature_importances <- rbind(feature_importances, importance_results)
}

# Print the feature importance summary with bounds
print(feature_importances)
write.csv(feature_importances,"./rsimca_model_feature_importance_seg.csv")

# Calculate average importance and bounds across folds
library(dplyr)
detach(package:plyr)

# Check that 'feature' is indeed a column and is grouped correctly
feature_summary <- feature_importances %>%
  group_by(feature_importances$feature) %>%
  summarize(
    mean_importance = mean(importance, na.rm = TRUE),
    mean_lower_bound = mean(lower_bound, na.rm = TRUE),
    mean_upper_bound = mean(upper_bound, na.rm = TRUE),
    mean_permutation_error = mean(permutation_error, na.rm = TRUE)
  
  )%>%
  arrange(desc(mean_importance))


# Print the feature importance summary with bounds
print(feature_summary)


write.csv(feature_summary,"./rsimca_model_feature_importance.csv")









