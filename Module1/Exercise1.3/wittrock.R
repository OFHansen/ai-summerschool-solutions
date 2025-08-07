
print("TASK 1")
#TASK 1
# Make a simple model where you predict metastasis based only on the variable
# that indicates whether a patient has urinary retention or not.

# Tip: You can use the glm function in R with family = "binomial" for logistic regression.

# Check and discuss the results from the model using the summary function:
# summary(your.model.name)

#the p value of the model shows there is a very significant relation between urinary.retention and metastasis.
#it also has log odds of 0.8, which is very strong

model_urinary_retention <- glm(metastasis ~ Urinary.Retention, data=MERGED_DATA, family = "binomial")
print(summary_model_urinary_retention <- summary(model_urinary_retention))

# What is the estimated effect of urinary retention on the probability of having metastatic
# prostate cancer? How would you interpret this effect?

# Based on this very simple model, calculate the probability of metastasis:
# - For a man with urinary retention (e.g., urin_ret = 1)
# - For a man without urinary retention (e.g., urin_ret = 0)

print(task_1_model_urinary_retention_true <- predict(model_urinary_retention, data.frame(Urinary.Retention = 1), type="response")) # = 0.5255973 
print(task_1_model_urinary_retention_false <- predict(model_urinary_retention, data.frame(Urinary.Retention = 0), type="response")) # = 0.3083141
print("______________________________________________________________________________")
print("TASK 2")

library("car")
print("-----------------VIF old model no PCA-------------------------")

task_2_model <- glm(metastasis ~ age_PCD + age_metastasis + monocyte + urea + lymphocyte + thrombocyte +
                      erythrocyte + albumin + creatinine + globulin + psa + coag_factor + Type.2.Diabetes +
                      Osteoporosis + Essential.Hypertension + Urinary.Retention + Hypercholesterolemia +
                      Age.related.Hearing.Loss,
                    data=MERGED_DATA, family = "binomial")
print(vif(task_2_model))
#including (almost) everything leads to a model where a lot of the variables have a high VIF value, indicating multicollinearity, and also 
#R gives a warning about it not cenverging, which is not good.
print("-----------------VIF old model corrected no PCA-------------------")

task_2_model_corrected <- glm(metastasis ~ monocyte + urea + lymphocyte + thrombocyte +
                                erythrocyte + psa + Type.2.Diabetes +
                                Osteoporosis + Essential.Hypertension + Urinary.Retention + Hypercholesterolemia +
                                Age.related.Hearing.Loss,
                              data=MERGED_DATA, family = "binomial")
print(vif(task_2_model_corrected))
#after running vif() on the old model, i removed anything with a value of over 5. this created a model where 
#multicollinearity is eliminated, and the model should be valid
print(summary(task_2_model_corrected))

print("-----------------------VIF PCA model--------------------------")

new_task_2_model <- glm(metastasis ~ PC1+PC2+PC3,
                        data=FINAL_DATASET, family = "binomial")
print(vif(new_task_2_model))
#in this model the VIF values are in the range 1.1-1.5, so it is very certain multicollinearity has been eliminated,
#so it is a valid model
print("TASK 3")
print("--------------------------------------------------------------")
#backwards stepwise variable selection, until everything has p<0.10
selected_model <- glm(metastasis ~  psa + Type.2.Diabetes + Osteoporosis + Essential.Hypertension + Urinary.Retention,
                      data=MERGED_DATA, family = "binomial")
print(summary(selected_model))

probabilities <- predict(selected_model, newdata = MERGED_DATA, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "1", "0")
predicted.classes <- as.factor(predicted.classes)
observed.classes <- MERGED_DATA$metastasis
pred.results <- data.frame(observed.classes, probabilities, predicted.classes) 
#print(pred.results)
library("caret")
print(confusionMatrix(data = predicted.classes, as.factor(observed.classes)))
#output:
#           Reference
# Prediction   0   1
#           0 680 277
#           1  58 144
#accuracy = 0.711
#sensitivity = 0.92
#specificity = 0.34
#initial:monocyte + urea + lymphocyte + thrombocyte + erythrocyte + psa + Type.2.Diabetes +
#Osteoporosis + Essential.Hypertension + Urinary.Retention + Hypercholesterolemia + Age.related.Hearing.Loss,

#order of removing variables: monocyte,age.related.hearing.loss, urea, hypercholesterolemia, thrombycyte, erythrocyte

#TASK 4
print("-----------------------TASK 4------------------------")

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

#TASK 5
print("_______________________TASK 5_________________________")
print("Trained model confusion matrix")
print(cm_trained_model)
print("PSA model confusion matrix")
probabilities <- predict(new_task_2_model, newdata = FINAL_DATASET, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "1", "0")
predicted.classes <- as.factor(predicted.classes)
observed.classes <- FINAL_DATASET$metastasis
pred.results <- data.frame(observed.classes, probabilities, predicted.classes) 
cm_psa_model <- confusionMatrix(data = predicted.classes, as.factor(observed.classes))
print(cm_psa_model)
#model comparison:
#model    accuracy   sensitivity   specificity
#trained    0.71         0.92         0.34
#PSA        0.67          0.98        0.08 
#looking strictly at the accuracy, the trained model is slightly better.
#the sensitivity of the PSA model is very high, but on the other hand it has very low specificity, 
#so it will have many false positives. This is good for catching all possible cases of metastasis, 
#but will "waste" a lot of time on those that don't actually have it, and may overcrowd hospitals
#and make wait times longer for those who do have it, which is bad. Because it uses PCs, it is
#impossible to tell which factors are indicating the metastasis, which is something both patients
#and medical professionals are interested in.
#The trained model has a lower sensitivity, but also significantly higher specificity. This means
#that slightly fewer cases will be found, but much less false positives occur, which reduces the
#strain on hospitals. It is also possible to determine exactly which variables are in play
#when classifying, which is information patients and medical professionals are interested in.
#I do not know which one of these two is the best. In a real life scenario I would ask for the
#opinion of the medical professionals who will use the model in practice.

#variable selection:
#I am not a medical professional, so can not speak for the exact reason that specific variables have
#on classification. My general assumption about these, is that all the removed ones are merely indicators
#or closely related to those that were kept. Similarly to producing the PCs, those that with the highest degree of 
#multicollinearity are removed.

#new study discussion
#I would try to collect a lot of different data, that is not related to each other, 
#to eliminate as much multicollinearity as possible.
#instead of binary classification, i would seperate it into more categories like "no risk", "low risk", "medium risk",
#and "high risk". This way, low specificity is not as much as a problem, as they would in most cases be classified as
#"low risk" instead of a simple "yes", so the hospital is aware they don't have to spend as many resources on this person,
#as they do on one in the medium or high risk group. 

