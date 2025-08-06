# Load required libraries
# Ensure you have these installed:
# install.packages(c("dplyr", "tidyr", "psych", "GGally"))
library(dplyr)
library(tidyr)
library(psych)
library(GGally)
library(readr)
DATASET1 <- read_csv("Exercise1.1/DATASET1.csv")
DATASET2 <- read_csv("Exercise1.1/DATASET2.csv")
DATASET3 <- read_csv("Exercise1.1/DATASET3.csv")

# --- Data Loading ---
# For this script to run, you must load your DATASET1, DATASET2, and DATASET3
# data frames first, as they would be in your RStudio environment.
# I will alias DATASET1 as MYDATA as requested.
MYDATA <- DATASET1

# --- TASK 1 (DATASET 1) - Basic Statistics & Summary Table ---

cat("--- Task 1: DATASET 1 Analysis ---\n\n")

# Inspect the data structure and get a summary
cat("--- Structure of MYDATA (DATASET 1) ---\n")
str(MYDATA)
cat("\n--- Summary of MYDATA (DATASET 1) ---\n")
print(summary(MYDATA))

# Check for missing data
cat("\n--- Missing Values in MYDATA (DATASET 1) ---\n")
print(sapply(MYDATA, function(x) sum(is.na(x))))

# Ensure date columns are in Date format
MYDATA$prostdtfirst <- as.Date(MYDATA$prostdtfirst)
MYDATA$date_met <- as.Date(MYDATA$date_met)
MYDATA$birthdate <- as.Date(MYDATA$birthdate)

# Calculate Age at PC Diagnosis
MYDATA$age_PCD <- as.numeric(difftime(MYDATA$prostdtfirst, MYDATA$birthdate, units = "days")) / 365.25

# Calculate Age at Metastasis (for those with metastasis = 1)
MYDATA$age_metastasis <- ifelse(MYDATA$metastasis == 1,
                                as.numeric(difftime(MYDATA$date_met, MYDATA$birthdate, units = "days")) / 365.25,
                                NA)

# Create the summary table
age_pcd_stats <- describe(MYDATA$age_PCD)[c("mean", "sd")]
age_met_stats <- describe(MYDATA$age_metastasis)[c("mean", "sd")]

summary_table_1 <- data.frame(
  Characteristic = c("Age at PC diagnosis",
                     "Age at metastasis",
                     "Metastasis - Yes",
                     "Metastasis - No",
                     "Family cancer - Yes",
                     "Family cancer - No",
                     "Married/cohabiting",
                     "Single/divorced/widow"),
  Value = c(
    paste0(round(age_pcd_stats$mean, 1), " ± ", round(age_pcd_stats$sd, 1)),
    paste0(round(age_met_stats$mean, 1), " ± ", round(age_met_stats$sd, 1)),
    paste0(sum(MYDATA$metastasis == 1, na.rm=T), " (", round(mean(MYDATA$metastasis == 1, na.rm=T) * 100, 1), "%)"),
    paste0(sum(MYDATA$metastasis == 0, na.rm=T), " (", round(mean(MYDATA$metastasis == 0, na.rm=T) * 100, 1), "%)"),
    paste0(sum(MYDATA$cancer_family == 1, na.rm=T), " (", round(mean(MYDATA$cancer_family == 1, na.rm=T) * 100, 1), "%)"),
    paste0(sum(MYDATA$cancer_family == 0, na.rm=T), " (", round(mean(MYDATA$cancer_family == 0, na.rm=T) * 100, 1), "%)"),
    paste0(sum(MYDATA$civil_status == 1, na.rm=T), " (", round(mean(MYDATA$civil_status == 1, na.rm=T) * 100, 1), "%)"),
    paste0(sum(MYDATA$civil_status == 2, na.rm=T), " (", round(mean(MYDATA$civil_status == 2, na.rm=T) * 100, 1), "%)")
  )
)

cat("\n--- Summary Table for DATASET 1 based on example table ---\n")
print(summary_table_1)


# --- TASK 2 (DATASET 2) - Inspecting Biochemical Data ---

cat("\n\n--- Task 2: DATASET 2 Analysis ---\n\n")

# Check for missing data in DATASET 2
cat("--- Missing Values in DATASET 2 ---\n")
print(sapply(DATASET2, function(x) sum(is.na(x))))

# Create a correlation matrix plot
# Note: This plot can take a moment to generate if the dataset is large.
# We will exclude the non-numeric 'ID' column.
cat("\n--- Correlation Plot for DATASET 2 ---\n")
corr_plot <- ggcorr(DATASET2[, -1], hjust = 0.75, size = 4, layout.exp = 1)
print(corr_plot)


# --- TASK 3 (DATASET 3) - Preparing the Data ---
cat("\n\n--- Task 3: DATASET 3 Preparation ---\n\n")

# Map provided disease codes with values found from www.medinfo.dk
DATASET3_MAPPED <- DATASET3 %>%
  mutate(disease_name = case_when(
    code == "DE119" ~ "Type 2 Diabetes",
    code == "DM819" ~ "Osteoporosis",
    code == "DI109" ~ "Essential Hypertension",
    code == "DR339" ~ "Urinary Retention",
    code == "DE780" ~ "Hypercholesterolemia",
    code == "DH911" ~ "Age-related Hearing Loss",
    TRUE ~ "Other" # Default case
  ))

# Reshape the data from long to wide format using the new disease names
DATASET3_WIDE <- DATASET3_MAPPED %>%
  mutate(has_disease = 1) %>%
  select(ID, disease_name, has_disease) %>% # Keep only necessary columns
  pivot_wider(id_cols = ID,
              names_from = disease_name,
              values_from = has_disease,
              values_fill = 0) # Fill missing combinations with 0

# Clean column names to be valid R variable names (e.g., replace spaces with dots)
names(DATASET3_WIDE) <- make.names(names(DATASET3_WIDE))

cat("\n--- Reshaped DATASET 3 with Mapped Names (first 5 rows) ---\n")
print(head(DATASET3_WIDE, 5))

# --- TASK 4 (DATASET 1-3) - Merging the Data ---

cat("\n\n--- Task 4: Merging All Datasets ---\n\n")

# Merge the three datasets together using dplyr's join functions
# Start with MYDATA (DATASET 1) and join the others to it.

MERGED_DATA <- MYDATA %>%
  left_join(DATASET2, by = "ID") %>%
  left_join(DATASET3_WIDE, by = "ID")

# Replace NA values in the newly added disease columns with 0
# This assumes that if a patient is not in DATASET3_WIDE, they have none of the diseases listed.
disease_cols <- names(DATASET3_WIDE)[-1] # Get disease column names, excluding ID
MERGED_DATA <- MERGED_DATA %>%
  mutate(across(all_of(disease_cols), ~replace_na(., 0)))

# Impute missing values for biochemical data using the mean
# This is a simple imputation method as suggested by the exercise guide.
cat("\n--- Imputing Missing Values in Biochemical Data using Mean ---\n")
biochem_cols <- names(DATASET2)[-1] # Get biochem column names, excluding ID

MERGED_DATA <- MERGED_DATA %>%
  mutate(across(all_of(biochem_cols), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# Verify that there are no more missing values
cat("\n--- Missing Values in Final Merged Dataset ---\n")
print(sapply(MERGED_DATA, function(x) sum(is.na(x))))

cat("\n--- Final Merged Dataset (first 5 rows) ---\n")
print(head(MERGED_DATA, 5))

print(filter(MERGED_DATA, !is.na(age_metastasis)))