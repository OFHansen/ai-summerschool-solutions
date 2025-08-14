# This document is designed to be run after Task 1 and 2 documents,
# DO NOT run tasks 1 and 2 after this file, without disabling the
# MASS package, it creates some issues when running Task 1 and 2.
# Thanks for listening to my Ted talk.
#
#
#
# Module 1.3: Logistic Regression Analysis for Prostate Cancer Metastasis Prediction
# This script requires that Task1.R and Task2.R have been run first to create the necessary datasets

# Load required libraries
library(dplyr)
library(car)        # for VIF calculation
library(MASS)       # for stepAIC
library(caret)      # for confusion matrix and cross-validation
library(pROC)       # for ROC curves
library(ggplot2)    # for plotting
library(gridExtra)  # for arranging multiple plots
library(tidyr)      # for data reshaping

# Ensure we have the required datasets from previous tasks
# MERGED_DATA should be available from Task1.R
# FINAL_DATASET should be available from Task2.R (with PCA components)

cat("=== MODULE 1.3: LOGISTIC REGRESSION ANALYSIS ===\n\n")

# First, let's create a binary urinary retention variable from the merged data
# Based on the disease mapping from Task1.R, we should have "Urinary.Retention" column
if(!"Urinary.Retention" %in% names(MERGED_DATA)) {
  # If not available, create a dummy variable for demonstration
  set.seed(123)
  MERGED_DATA$Urinary.Retention <- sample(c(0,1), nrow(MERGED_DATA), replace = TRUE, prob = c(0.8, 0.2))
  cat("Note: Created dummy Urinary.Retention variable for demonstration\n\n")
}

# =============================================================================
# =============================================================================
# TASK 1: Simple Logistic Regression Model (30 min)
# =============================================================================
# =============================================================================

cat("--- TASK 1: Simple Logistic Regression Model ---\n\n")

# Simple model: metastasis ~ urinary retention
simple_model <- glm(metastasis ~ Urinary.Retention, data = MERGED_DATA, family = "binomial")

cat("Simple Model Summary:\n")
print(summary(simple_model))

# Interpret the effect
coef_urin_ret <- coef(simple_model)["Urinary.Retention"]
odds_ratio <- exp(coef_urin_ret)
cat("\n--- Interpretation ---\n")
cat("Coefficient for urinary retention:", round(coef_urin_ret, 4), "\n")
cat("Odds ratio:", round(odds_ratio, 4), "\n")
cat("Interpretation: Patients with urinary retention have", round(odds_ratio, 2), 
    "times the odds of having metastatic prostate cancer compared to those without.\n\n")

# Calculate probabilities for patients with and without urinary retention
prob_with_retention <- predict(simple_model, newdata = data.frame(Urinary.Retention = 1), type = "response")
prob_without_retention <- predict(simple_model, newdata = data.frame(Urinary.Retention = 0), type = "response")

cat("Predicted Probabilities:\n")
cat("- Patient WITH urinary retention:", round(prob_with_retention, 4), "\n")
cat("- Patient WITHOUT urinary retention:", round(prob_without_retention, 4), "\n\n")

# Plot 1: Simple logistic regression curve
cat("Creating plot for simple logistic regression...\n")
plot_data <- data.frame(
  Urinary.Retention = c(0, 1),
  probability = c(prob_without_retention, prob_with_retention)
)

p1 <- ggplot(MERGED_DATA, aes(x = factor(Urinary.Retention), y = metastasis)) +
  geom_jitter(alpha = 0.6, width = 0.2, height = 0.05) +
  geom_point(data = plot_data, aes(x = factor(Urinary.Retention), y = probability), 
             color = "red", size = 4, shape = 15) +
  geom_line(data = plot_data, aes(x = as.numeric(factor(Urinary.Retention)), y = probability), 
            color = "red", size = 1.2) +
  labs(title = "Simple Logistic Regression: Metastasis ~ Urinary Retention",
       x = "Urinary Retention (0=No, 1=Yes)", 
       y = "Probability of Metastasis",
       subtitle = "Red squares show predicted probabilities") +
  theme_minimal() +
  ylim(0, 1)

print(p1)

# =============================================================================
# =============================================================================
# TASK 2: Multiple Logistic Regression Models (30 min)
# =============================================================================
# =============================================================================

cat("--- TASK 2: Multiple Logistic Regression Models ---\n\n")

# Model 1: All variables except PCs and non-relevant variables
# Exclude: ID, prostdtfirst, date_met, birthdate, age_metastasis (since it's after outcome)
exclude_vars <- c("ID", "prostdtfirst", "date_met", "birthdate", "age_metastasis", "creatinine")

# Get all variable names except excluded ones
all_vars <- names(MERGED_DATA)[!names(MERGED_DATA) %in% c("metastasis", exclude_vars)]
# Remove any PC variables if they exist
all_vars <- all_vars[!grepl("^PC", all_vars)]

# Create formula
formula1 <- as.formula(paste("metastasis ~", paste(all_vars, collapse = " + ")))

cat("Model 1: All available variables\n")
cat("Formula:", deparse(formula1), "\n\n")

# Fit the model
full_model <- glm(formula1, data = MERGED_DATA, family = "binomial")

cat("Full Model Summary:\n")
print(summary(full_model))

# Check for collinearity using VIF
cat("\n--- Variance Inflation Factor (VIF) ---\n")
vif_values <- vif(full_model)
print(vif_values)

# Flag high VIF values (> 5 or 10)
high_vif <- vif_values[vif_values > 5]
if(length(high_vif) > 0) {
  cat("\nWarning: High VIF values detected (>5):\n")
  print(high_vif)
} else {
  cat("\nAll VIF values are acceptable (<5).\n")
}

# Model 2: Using Principal Components (if available from Task2.R)
if(exists("FINAL_DATASET") && any(grepl("^PC", names(FINAL_DATASET)))) {
  cat("\n\n--- Model 2: Using Principal Components ---\n")
  
  # Get PC variable names
  pc_vars <- names(FINAL_DATASET)[grepl("^PC", names(FINAL_DATASET))]
  
  # Get non-biochemical variables
  non_biochem_vars <- c("civil_status", "cancer_family", "age_PCD")
  if("Urinary.Retention" %in% names(FINAL_DATASET)) {
    non_biochem_vars <- c(non_biochem_vars, "Urinary.Retention")
  }
  
  # Include disease variables
  disease_vars <- names(FINAL_DATASET)[grepl("^(Type|Essential|Age|Hyper|Osteo)", names(FINAL_DATASET))]
  
  # Create formula with PCs
  pc_formula_vars <- c(non_biochem_vars, disease_vars, pc_vars)
  pc_formula_vars <- pc_formula_vars[pc_formula_vars %in% names(FINAL_DATASET)]
  
  pc_formula <- as.formula(paste("metastasis ~", paste(pc_formula_vars, collapse = " + ")))
  
  cat("PC Model Formula:", deparse(pc_formula), "\n\n")
  
  pc_model <- glm(pc_formula, data = FINAL_DATASET, family = "binomial")
  
  cat("PC Model Summary:\n")
  print(summary(pc_model))
  
  # Check VIF for PC model
  cat("\n--- VIF for PC Model ---\n")
  vif_pc <- vif(pc_model)
  print(vif_pc)
  
  cat("\n--- Advantages and Disadvantages of Using Principal Components ---\n")
  cat("Advantages:\n")
  cat("- Reduces multicollinearity among biochemical variables\n")
  cat("- Reduces dimensionality while retaining most variance\n")
  cat("- Can handle highly correlated predictors\n\n")
  cat("Disadvantages:\n")
  cat("- Loss of interpretability (PCs are linear combinations)\n")
  cat("- May lose some important variable-specific information\n")
  cat("- Difficult to explain clinical relevance to practitioners\n\n")
} else {
  cat("\n\nNote: Principal components not available. Skipping PC model.\n\n")
  pc_model <- NULL
}

# =============================================================================
# =============================================================================
# TASK 3: Model Selection and Improvement (15 min)
# =============================================================================
# =============================================================================

cat("--- TASK 3: Model Selection and Improvement ---\n\n")

# Use stepwise selection on the full model
cat("Performing stepwise variable selection...\n")
step_model <- stepAIC(full_model, direction = "both", trace = FALSE)

cat("Final Selected Model:\n")
print(summary(step_model))

# Get selected variables
selected_vars <- names(coef(step_model))[-1]  # Exclude intercept
cat("\nSelected variables:", paste(selected_vars, collapse = ", "), "\n\n")

# Calculate predicted probabilities
probabilities <- predict(step_model, newdata = MERGED_DATA, type = "response")

# Create predicted classes using 0.5 cutoff
predicted_classes <- ifelse(probabilities > 0.5, 1, 0)
predicted_classes <- as.factor(predicted_classes)
observed_classes <- as.factor(MERGED_DATA$metastasis)

# Create confusion matrix
cat("--- Confusion Matrix (Cutoff = 0.5) ---\n")
conf_matrix <- confusionMatrix(predicted_classes, observed_classes, positive = "1")
print(conf_matrix)

# Calculate performance metrics manually
TP <- conf_matrix$table[2,2]
TN <- conf_matrix$table[1,1]
FP <- conf_matrix$table[2,1]
FN <- conf_matrix$table[1,2]

accuracy <- (TP + TN) / (TP + TN + FP + FN)
sensitivity <- TP / (TP + FN)
specificity <- TN / (TN + FP)

cat("\n--- Performance Metrics ---\n")
cat("Accuracy:", round(accuracy, 4), "\n")
cat("Sensitivity:", round(sensitivity, 4), "\n")
cat("Specificity:", round(specificity, 4), "\n\n")

# Plot 2: ROC Curve for the selected model
cat("Creating ROC curve for the final selected model...\n")
roc_obj <- roc(MERGED_DATA$metastasis, probabilities, quiet = TRUE)
auc_value <- auc(roc_obj)

p2 <- ggroc(roc_obj) +
  labs(title = paste("ROC Curve for Final Logistic Regression Model"),
       subtitle = paste("AUC =", round(auc_value, 3)),
       x = "Specificity (1 - False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  theme_minimal() +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50")

print(p2)

# Plot 3: Predicted probabilities distribution
cat("Creating predicted probabilities distribution plot...\n")
prob_data <- data.frame(
  probability = probabilities,
  metastasis = factor(MERGED_DATA$metastasis, labels = c("No Metastasis", "Metastasis"))
)

p3 <- ggplot(prob_data, aes(x = probability, fill = metastasis)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Distribution of Predicted Probabilities",
       subtitle = "Red line shows 0.5 cutoff threshold",
       x = "Predicted Probability of Metastasis",
       y = "Count",
       fill = "Actual Status") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "coral"))

print(p3)

# =============================================================================
# =============================================================================
# TASK 4: Cross-Validation (45 min)
# =============================================================================
# =============================================================================

cat("--- TASK 4: 10-Fold Cross-Validation ---\n\n")
train.control <- trainControl(method = "cv", number = 10, savePredictions="final") 
set.seed(123)
model_for_test <- train(metastasis ~ psa + Type.2.Diabetes + Osteoporosis + Essential.Hypertension + Urinary.Retention, 
                        data=MERGED_DATA, method = "glm", family = "binomial", trControl = train.control)
#Inspecting the results for each folder
pred_fold <- model_for_test$pred

pred_fold$class <- ifelse(pred_fold$pred >= 0.5, 1, 0)
pred_fold$equal <- ifelse(pred_fold$class == pred_fold$obs, 1, 0)

#Calculating the accuracy per folder
library(dplyr)
eachfold <- pred_fold %>% group_by(Resample) %>% summarise_at(vars(equal),list(Accuracy = mean))
print(eachfold)
# Confusion matrix — make sure both are factors with the same levels
cm_trained_model <- confusionMatrix(factor(pred_fold$class, levels = c(0,1)),
                                    factor(pred_fold$obs, levels = c(0,1)))
print(cm_trained_model)
#accuracy = 0.71
#sens = 0.92
#spec = 0.34
#pred_fold <- model_for_test$pred 
#print(pred_fold)

cat("\n--- Cross-Validation Performance Metrics ---\n")
cat("CV Accuracy:", round(cv_conf_matrix$overall["Accuracy"], 4), "\n")
cat("CV Sensitivity:", round(cv_conf_matrix$byClass["Sensitivity"], 4), "\n")
cat("CV Specificity:", round(cv_conf_matrix$byClass["Specificity"], 4), "\n\n")

# Plot 4: Cross-validation accuracy across folds
cat("Creating cross-validation accuracy plot...\n")
p4 <- ggplot(eachfold, aes(x = Resample, y = Accuracy)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = avg_accuracy, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Cross-Validation Accuracy Across Folds",
       subtitle = paste("Average accuracy:", round(avg_accuracy, 3), "± ", round(sd_accuracy, 3)),
       x = "Fold", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

print(p4)

# =============================================================================
# =============================================================================
# TASK 5: Closing Remarks and Comparisons
# =============================================================================
# =============================================================================

cat("--- TASK 5: Closing Remarks ---\n\n")

# Compare with PSA-based classification
# Using PSA cutoffs from literature
MERGED_DATA$psa_pred_20 <- ifelse(MERGED_DATA$psa >= 20, 1, 0)
MERGED_DATA$psa_pred_100 <- ifelse(MERGED_DATA$psa >= 100, 1, 0)

# PSA >= 20 performance
psa20_accuracy <- mean(MERGED_DATA$psa_pred_20 == MERGED_DATA$metastasis, na.rm = TRUE)
psa20_sensitivity <- sum(MERGED_DATA$psa_pred_20 == 1 & MERGED_DATA$metastasis == 1, na.rm = TRUE) / 
  sum(MERGED_DATA$metastasis == 1, na.rm = TRUE)
psa20_specificity <- sum(MERGED_DATA$psa_pred_20 == 0 & MERGED_DATA$metastasis == 0, na.rm = TRUE) / 
  sum(MERGED_DATA$metastasis == 0, na.rm = TRUE)

# PSA >= 100 performance
psa100_accuracy <- mean(MERGED_DATA$psa_pred_100 == MERGED_DATA$metastasis, na.rm = TRUE)
psa100_sensitivity <- sum(MERGED_DATA$psa_pred_100 == 1 & MERGED_DATA$metastasis == 1, na.rm = TRUE) / 
  sum(MERGED_DATA$metastasis == 1, na.rm = TRUE)
psa100_specificity <- sum(MERGED_DATA$psa_pred_100 == 0 & MERGED_DATA$metastasis == 0, na.rm = TRUE) / 
  sum(MERGED_DATA$metastasis == 0, na.rm = TRUE)

cat("--- Performance Comparison ---\n")
comparison_table <- data.frame(
  Method = c("Logistic Regression (Training)", "Logistic Regression (CV)", "PSA >= 20", "PSA >= 100"),
  Accuracy = c(round(accuracy, 4), round(avg_accuracy, 4), round(psa20_accuracy, 4), round(psa100_accuracy, 4)),
  Sensitivity = c(round(sensitivity, 4), round(cv_conf_matrix$byClass["Sensitivity"], 4), 
                  round(psa20_sensitivity, 4), round(psa100_sensitivity, 4)),
  Specificity = c(round(specificity, 4), round(cv_conf_matrix$byClass["Specificity"], 4), 
                  round(psa20_specificity, 4), round(psa100_specificity, 4))
)
print(comparison_table)

# Plot 5: Performance comparison visualization
cat("Creating performance comparison plot...\n")
comparison_long <- comparison_table %>%
  tidyr::pivot_longer(cols = c(Accuracy, Sensitivity, Specificity), 
                      names_to = "Metric", values_to = "Value")

p5 <- ggplot(comparison_long, aes(x = Method, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  labs(title = "Performance Comparison: Different Prediction Methods",
       x = "Method", y = "Performance Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  ylim(0, 1)

print(p5)

# Plot 6: Coefficient plot for final model
cat("Creating coefficient plot for final model...\n")
coef_data <- data.frame(
  Variable = names(coef(step_model))[-1],  # Exclude intercept
  Coefficient = coef(step_model)[-1],
  SE = summary(step_model)$coefficients[-1, "Std. Error"]
)

coef_data$Lower <- coef_data$Coefficient - 1.96 * coef_data$SE
coef_data$Upper <- coef_data$Coefficient + 1.96 * coef_data$SE
coef_data$Significant <- abs(coef_data$Coefficient) / coef_data$SE > 1.96

p6 <- ggplot(coef_data, aes(x = reorder(Variable, Coefficient), y = Coefficient)) +
  geom_point(aes(color = Significant), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Significant), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Coefficients of Final Logistic Regression Model",
       subtitle = "Error bars show 95% confidence intervals",
       x = "Variables", y = "Log Odds (Coefficient)") +
  theme_minimal() +
  coord_flip() +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red"),
                     labels = c("Non-significant", "Significant (p<0.05)"))

print(p6)

cat("Creating coefficient plot for final model...\n")
coef_data <- data.frame(
  Variable = names(coef(pc_model))[-1],  # Exclude intercept
  Coefficient = coef(pc_model)[-1],
  SE = summary(pc_model)$coefficients[-1, "Std. Error"]
)

coef_data$Lower <- coef_data$Coefficient - 1.96 * coef_data$SE
coef_data$Upper <- coef_data$Coefficient + 1.96 * coef_data$SE
coef_data$Significant <- abs(coef_data$Coefficient) / coef_data$SE > 1.96

p7 <- ggplot(coef_data, aes(x = reorder(Variable, Coefficient), y = Coefficient)) +
  geom_point(aes(color = Significant), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Significant), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Coefficients of Final Logistic Regression Model",
       subtitle = "Error bars show 95% confidence intervals",
       x = "Variables", y = "Log Odds (Coefficient)") +
  theme_minimal() +
  coord_flip() +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red"),
                     labels = c("Non-significant", "Significant (p<0.05)"))

print(p7)

cat("\n--- Analysis of Final Model Variables ---\n")
final_coefficients <- summary(step_model)$coefficients
significant_vars <- rownames(final_coefficients)[final_coefficients[,"Pr(>|z|)"] < 0.05]
significant_vars <- significant_vars[significant_vars != "(Intercept)"]

if(length(significant_vars) > 0) {
  cat("Significant predictors in the final model:\n")
  for(var in significant_vars) {
    coef_val <- final_coefficients[var, "Estimate"]
    p_val <- final_coefficients[var, "Pr(>|z|)"]
    odds_ratio <- exp(coef_val)
    cat("-", var, ": OR =", round(odds_ratio, 3), ", p =", round(p_val, 4), "\n")
  }
} else {
  cat("No statistically significant predictors found.\n")
}

cat("\n--- Study Design Recommendations ---\n")
cat("For a future study design, consider:\n")
cat("1. Prospective data collection to ensure temporal relationships\n")
cat("2. Standardized follow-up periods for metastasis detection\n")
cat("3. Additional biomarkers beyond traditional biochemical tests\n")
cat("4. Imaging-based features as potential predictors\n")
cat("5. Larger sample size, especially for metastatic cases\n")
cat("6. External validation cohort from different centers\n")
cat("7. Cost-effectiveness analysis of the prediction model\n")
cat("8. Integration with clinical decision support systems\n\n")

cat("=== ANALYSIS COMPLETE ===\n")

# Save the final model for future use
# save(step_model, file = "final_metastasis_prediction_model.RData")
# cat("Final model saved to 'final_metastasis_prediction_model.RData'\n")