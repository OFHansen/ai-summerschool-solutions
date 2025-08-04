# Load required libraries (only psych needed, dplyr is optional)
# install.packages("psych") # run this first if you don't have psych
library(psych)

# Convert date columns to Date format
mydata$prostdtfirst <- as.Date(mydata$prostdtfirst)
mydata$date_met <- as.Date(mydata$date_met)
mydata$birthdate <- as.Date(mydata$birthdate)

# Calculate age at prostate cancer diagnosis (using base R)
mydata$age_PCD <- as.numeric(difftime(mydata$prostdtfirst, mydata$birthdate, units = "days")) / 365.25

# Calculate age at metastasis (only for patients with metastasis = 1)
mydata$age_metastasis <- ifelse(mydata$metastasis == 1,
                                as.numeric(difftime(mydata$date_met, mydata$birthdate, units = "days")) / 365.25,
                                NA)

# 1. Age at PC diagnosis - mean and SD
describe(mydata$age_PCD)[c("mean", "sd")]

# 2. Age at metastasis (for those with metastasis = 1)
patients_with_met <- subset(mydata, metastasis == 1)
if(nrow(patients_with_met) > 0) {
  describe(patients_with_met$age_metastasis)[c("mean", "sd")]
} else {
  print("No patients with metastasis")
}

# 3. Metastasis percentages
table(mydata$metastasis)
prop.table(table(mydata$metastasis)) * 100

# 4. Family cancer history
table(mydata$cancer_family)
prop.table(table(mydata$cancer_family)) * 100

# 5. Civil status
table(mydata$civil_status)
prop.table(table(mydata$civil_status)) * 100

# Summary table creation
summary_table <- data.frame(
  Characteristic = c("Age at PC diagnosis (mean ± SD)",
                     "Age at metastasis (mean ± SD)",
                     "Metastasis - Yes",
                     "Metastasis - No", 
                     "Family cancer - Yes",
                     "Family cancer - No",
                     "Married/cohabiting",
                     "Single/divorced/widow"),
  Value = c(
    paste0(round(describe(mydata$age_PCD)$mean, 1), " ± ", 
           round(describe(mydata$age_PCD)$sd, 1)),
    ifelse(nrow(patients_with_met) > 0,
           paste0(round(describe(patients_with_met$age_metastasis)$mean, 1), " ± ", 
                  round(describe(patients_with_met$age_metastasis)$sd, 1)),
           "No cases"),
    paste0(sum(mydata$metastasis == 1), " (", 
           round(mean(mydata$metastasis == 1) * 100, 1), "%)"),
    paste0(sum(mydata$metastasis == 0), " (", 
           round(mean(mydata$metastasis == 0) * 100, 1), "%)"),
    paste0(sum(mydata$cancer_family == 1), " (", 
           round(mean(mydata$cancer_family == 1) * 100, 1), "%)"),
    paste0(sum(mydata$cancer_family == 0), " (", 
           round(mean(mydata$cancer_family == 0) * 100, 1), "%)"),
    paste0(sum(mydata$civil_status == 1), " (", 
           round(mean(mydata$civil_status == 1) * 100, 1), "%)"),
    paste0(sum(mydata$civil_status == 2), " (", 
           round(mean(mydata$civil_status == 2) * 100, 1), "%)")
  )
)

# Display the summary table
print(summary_table)