# Load required libraries
# Ensure you have these installed: install.packages(c("psych", "mnormt", "OptimalCutpoints", "dplyr"))
# install.packages("dplyr")
library(psych)
library(OptimalCutpoints)
library(dplyr)


# --- Data Loading and Preparation ---
# For this script to run, you must load your DATASET1, DATASET2, and DATASET3 data frames first.
# The data sets are loaded using the Rstudio environment


# --- Data set formats--- 
# DATASET1 = ID: int, civil_status: int, cancer_family: int,
# prostdtfirst: Date, metastasis: int, date_met: Date, birthdate: Date,
# age_metastasis: num, age_PCD: num

# DATASET2 = ID: num, monocyte: num, urea: num, lymphocyte: num, thrombocyte: num, albumin: num,
# creatinine: num, globulin: num, psa: num, coag_factor: num

# DATASET3 = ID: num, code: chr (e.g. "DE119", "CM819"), date_diag: Date
# ------------------------

FULL_DATA <- MERGED_DATA


# --- TASK 1: Measures of test performance based on PSA ---

# Create predicted metastasis variables based on PSA cutoffs.
FULL_DATA$predicted_met_20 <- ifelse(FULL_DATA$psa >= 20, 1, 0)
FULL_DATA$predicted_met_100 <- ifelse(FULL_DATA$psa >= 100, 1, 0)

# Create factors for the confusion matrices.
actual_status <- factor(FULL_DATA$metastasis, levels = c(0, 1), labels = c("Actual: No Met", "Actual: Met"))
predicted_status_20 <- factor(FULL_DATA$predicted_met_20, levels = c(0, 1), labels = c("Predicted: No Met", "Predicted: Met"))
predicted_status_100 <- factor(FULL_DATA$predicted_met_100, levels = c(0, 1), labels = c("Predicted: No Met", "Predicted: Met"))

# Generate the confusion matrix for the 20 ng/ml cutoff
cat("\n\n--- Confusion Matrix (PSA Cutoff >= 20 ng/ml) ---\n")
conf_matrix_20 <- table(predicted_status_20, actual_status)
print(conf_matrix_20)

# Generate the confusion matrix for the 100 ng/ml cutoff
cat("\n\n--- Confusion Matrix (PSA Cutoff >= 100 ng/ml) ---\n")
conf_matrix_100 <- table(predicted_status_100, actual_status)
print(conf_matrix_100)

# Function to calculate performance metrics
calculate_metrics <- function(conf_matrix) {
  TP <- conf_matrix[2, 2]; TN <- conf_matrix[1, 1]
  FP <- conf_matrix[2, 1]; FN <- conf_matrix[1, 2]
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  accuracy    <- ifelse((TP + TN + FP + FN) > 0, (TP + TN) / (TP + TN + FP + FN), 0)
  return(c(Sensitivity = sensitivity, Specificity = specificity, Accuracy = accuracy))
}

metrics_20 <- calculate_metrics(conf_matrix_20)
metrics_100 <- calculate_metrics(conf_matrix_100)

performance_summary <- data.frame(
  "Prostate cancer metastasis" = c("Sensitivity", "Specificity", "Accuracy"),
  "PSA > 20 ng/ml" = sprintf("%.3f (%.1f%%)", metrics_20, metrics_20 * 100),
  "PSA > 100 ng/ml" = sprintf("%.3f (%.1f%%)", metrics_100, metrics_100 * 100),
  check.names = FALSE
)

cat("\n\n--- Performance Metrics Summary ---\n")
print(performance_summary, row.names = FALSE)


# --- TASK 3: Using the ROC curve to find the optimal cutoff ---

cat("\n\n--- Optimal Cutoff Analysis (Task 3) ---\n")

# Method 1: Youden Index
# This method maximizes (sensitivity + specificity - 1)
cutpoint_youden <- optimal.cutpoints(X = "psa", status = "metastasis", tag.healthy = 0,
                                     methods = "Youden", data = FULL_DATA,
                                     ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
cat("\n--- Youden Index Method Summary ---\n")
# Show's best cut off is at 23% sensititivty = True positives, and Specificity = True Negatives
# AOC seems to be 0.67, which means it is a poor indicator (0.5 being random and 1 being perfect)
summary(cutpoint_youden)
plot(cutpoint_youden)


# Method 2: Maximize Product of Sensitivity and Specificity
# This can be a good alternative when costs of misclassification are unknown.
cutpoint_max_spse <- optimal.cutpoints(X = "psa", status = "metastasis", tag.healthy = 0,
                                       methods = "MaxSpSe", data = FULL_DATA,
                                       ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
cat("\n--- Max Product of Sp & Se Method Summary ---\n")
# Show's best cut off at 17 ng/ml
summary(cutpoint_max_spse)
plot(cutpoint_max_spse)


# --- TASK 4: Principal Component Analysis ---

cat("\n\n--- Principal Component Analysis (Task 4) ---\n")

# Define the biomarker columns to be used for PCA
biomarker_cols <- c("monocyte", "urea", "lymphocyte", "thrombocyte", "albumin", "creatinine", "globulin", "psa", "coag_factor", "erythrocyte")

# The prcomp function cannot handle missing values (NA).
# We must first create a new data frame that excludes any rows with NA in the selected biomarker columns.
# This ensures that the data used for PCA is complete and that subsequent steps (like binding columns) are aligned.
FULL_DATA_PCA <- FULL_DATA[complete.cases(FULL_DATA[, biomarker_cols]), ]

# Perform PCA on the cleaned biomarker data, ensuring to center and scale the variables.
# We select the columns directly from our new, clean data frame.
pca_model <- prcomp(FULL_DATA_PCA[, biomarker_cols], center = TRUE, scale. = TRUE)

# Print the summary of the PCA model
cat("\n--- PCA Model Summary ---\n")
print(summary(pca_model))

# Create a scree plot to visualize the variance explained by each component
cat("\n--- Scree Plot ---\n")
screeplot(pca_model, type = "l", main = "Scree Plot of PCA Components")
abline(h = 1, col = "red", lty = 2)
legend("topright", legend = "Eigenvalue = 1", col = "red", lty = 2)


# Calculate and display Eigenvalues to help decide number of PCs to keep
eigenvalues <- pca_model$sdev^2
cat("\n--- Eigenvalues ---\n")
print(eigenvalues)
cat("\nBased on the Eigenvalue > 1 criterion, we would keep the components with a value greater than 1.\n")


# Inspect the component loadings (rotation matrix)
cat("\n--- Component Loadings (Rotation) ---\n")
print(pca_model$rotation)

# Let's say we decide to keep the first 3 Principal Components based on the scree plot and eigenvalues
num_pcs_to_keep <- 3
components <- pca_model$x[, 1:num_pcs_to_keep]
print(components)
# Generate a new dataset by combining the PCs with the cleaned original dataset.
# This works because both data frames now have the same number of rows.
FINAL_DATASET <- bind_cols(FULL_DATA_PCA, as.data.frame(components))

# View the first few rows of the new combined dataset
print(head(FINAL_DATASET))

# Export the final dataset to a CSV file
#write.csv(FINAL_DATASET, "prostate_cancer_with_pcs.csv", row.names = FALSE)
#cat("\nFinal dataset with Principal Components has been saved to 'prostate_cancer_with_pcs.csv'\n")

