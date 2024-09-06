##Prerequisites
library(tidyverse)
library(lavaan)

##Load tibble 
sample <- read_csv("sample.csv")
subsample <- read_csv("subsample.csv")

##Create character vectors
cardiometabolic <- c("std_bmi", "std_waist_circumference","std_sbp", "std_dbp","std_pulse_pressure","log_crp", 
                     "std_hba1c","std_hdl","std_ldl", "std_cholesterol","log_triglycerides","PC1", "PC2", "PC3")

label_cardiometabolic <- c("Body mass index","Waist circumference","Systolic blood pressure","Diastolic blood pressure","Pulse pressure",
                           "C-reactive protein","Glycated hemoglobin","High-density lipoprotein cholesterol","Low-density lipoprotein cholesterol",
                           "Total cholesterol","Triglycerides","Principal component 1", "Principal component 2", "Principal component 3")

cognitive_tests <- c("std_numeric_memory", "std_fluid_intelligence", "log_trail_making_test_b", "std_matrix_test", 
                     "std_symbol_digit_substitution", "std_tower_rearranging", "std_paired_associate_learning", "log_pairs_matching")

liver_fat <- c("log_liver_fat", "nafld", "mafld", "masld", "metald")

label_liver_fat <- c("Liver fat", "NAFLD", "MAFLD", "MASLD", "METALD")

prefix <- c("A effect", "Direct effect", "B effect", "Indirect effect", "Total effect")

variables_wmh <- paste(prefix,rep(label_cardiometabolic, each = 5), sep=" ")
variables_cognition <- paste(prefix,rep(label_liver_fat, each = 5), sep=" ")

mediation_names <- c("Variable","Standardized regression coefficient","Standard error", "Z-score", "P-value", "Lower confidence interval", "Upper confidence interval")

##Create mediation models
model_wmh_mediation <- lapply(cardiometabolic, function(x) c(paste("log_liver_fat~a*", x, 
  "+ std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + bristol + current_smoker + former_smoker + std_alcohol_consumption \n log_wmh~c*", x, 
  "+ std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + bristol + current_smoker + former_smoker + std_alcohol_consumption \n log_wmh~b*log_liver_fat \n ab := a*b \n total := c + (a*b)")))


model_wmh_mediation_by_sex <- 
  lapply(cardiometabolic, function(x) c(paste("log_liver_fat~a*", x, 
  "+ std_age_imaging + std_age_squared_imaging + reading + newcastle + bristol + current_smoker + former_smoker + std_alcohol_consumption \n log_wmh~c*", x,
  "+ std_age_imaging + std_age_squared_imaging + reading + newcastle + bristol + current_smoker + former_smoker + std_alcohol_consumption \n log_wmh~b*log_liver_fat \n ab := a*b \n total := c + (a*b)")))     
  

##No participants were assessed at Bristol therefore removed
model_cognition_mediation <- 
  lapply(liver_fat, function(x) c(paste("log_wmh~a*", x, 
  "+ std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n cognition~c*", x,
  "+ std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n cognition~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))  
  
model_cognition_mediation_by_sex <- lapply(liver_fat, function(x) c(paste("log_wmh~a*", x, 
  "+ std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n cognition~c*", x,
  "+ std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n cognition~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))  

##WMH mediation
wmh_mediation <- lapply(model_wmh_mediation, function(x) sem(x, se = "bootstrap", bootstrap = 10000, sample))
results_wmh_mediation <- lapply(1:14, function(x) filter(standardizedsolution(wmh_mediation[[x]]), label %in% c("a", "b", "c", "ab", "total")))
fit_measures_wmh_mediation <- lapply(1:14, function(x) fitmeasures(wmh_mediation[[x]], c("chisq", "pvalue", "cfi", "tli", "rmsea", "srmr")))
results_wmh_mediation <- tibble(bind_rows(results_wmh_mediation))
results_wmh_mediation %>% select(est.std, se, z, pvalue, ci.lower, ci.upper)
results_wmh_mediation <- tibble(bind_rows(results_wmh_mediation))
results_wmh_mediation <- bind_cols(variables_wmh, results_wmh_mediation)
results_wmh_mediation <- results_wmh_mediation %>% select(...1, est.std, se, z, pvalue, ci.lower, ci.upper)
names(results_wmh_mediation) <- mediation_names

men_wmh_mediation <- lapply(model_wmh_mediation_by_sex, function(x) sem(x, se = "bootstrap", bootstrap = 10, subset(sample, sex == 1)))
men_results_wmh_mediation <- lapply(1:14, function(x) filter(standardizedsolution(men_wmh_mediation[[x]]), label %in% c("a", "b", "c", "ab", "total")))
fit_measures_men_wmh_mediation <- lapply(1:14, function(x) fitmeasures(men_wmh_mediation[[x]], c("chisq", "pvalue", "cfi", "tli", "rmsea", "srmr")))
men_results_wmh_mediation <- tibble(bind_rows(men_results_wmh_mediation))
men_results_wmh_mediation %>% select(est.std, se, z, pvalue, ci.lower, ci.upper)
men_results_wmh_mediation <- tibble(bind_rows(men_results_wmh_mediation))
men_results_wmh_mediation <- bind_cols(paste("Men",variables_wmh), men_results_wmh_mediation)
men_results_wmh_mediation <- men_results_wmh_mediation %>% select(...1, est.std, se, z, pvalue, ci.lower, ci.upper)
names(men_results_wmh_mediation) <- mediation_names

women_wmh_mediation <- lapply(model_wmh_mediation_by_sex, function(x) sem(x, se = "bootstrap", bootstrap = 10000, subset(sample, sex == 0)))
women_results_wmh_mediation <- lapply(1:14, function(x) filter(standardizedsolution(women_wmh_mediation[[x]]), label %in% c("a", "b", "c", "ab", "total")))
fit_measures_women_wmh_mediation <- lapply(1:14, function(x) fitmeasures(women_wmh_mediation[[x]], c("chisq", "pvalue", "cfi", "tli", "rmsea", "srmr")))
women_results_wmh_mediation <- tibble(bind_rows(women_results_wmh_mediation))
women_results_wmh_mediation %>% select(est.std, se, z, pvalue, ci.lower, ci.upper)
women_results_wmh_mediation <- tibble(bind_rows(women_results_wmh_mediation))
women_results_wmh_mediation <- bind_cols(paste("Women",variables_wmh), women_results_wmh_mediation)
women_results_wmh_mediation <- women_results_wmh_mediation %>% select(...1, est.std, se, z, pvalue, ci.lower, ci.upper)
names(women_results_wmh_mediation) <- mediation_names

results_wmh_mediation <- bind_rows(results_wmh_mediation, men_results_wmh_mediation, women_results_wmh_mediation)
write_csv(results_wmh_mediation, "wmh_mediation.csv")

##Cognition mediation
cognition_mediation <- lapply(model_cognition_mediation, function(x) sem(x, se = "bootstrap", bootstrap = 10000, subsample))
results_cognition_mediation <- lapply(1:5, function(x) filter(standardizedsolution(cognition_mediation[[x]]), label %in% c("a", "b", "c", "ab", "total")))
fit_measures_cognition_mediation <- lapply(1:5, function(x) fitmeasures(cognition_mediation[[x]], c("chisq", "pvalue", "cfi", "tli", "rmsea", "srmr")))
results_cognition_mediation <- tibble(bind_rows(results_cognition_mediation))
results_cognition_mediation <- bind_cols(variables_cognition, results_cognition_mediation)
results_cognition_mediation <- results_cognition_mediation %>% select(...1, est.std, se, z, pvalue, ci.lower, ci.upper)
names(results_cognition_mediation) <- mediation_names

men_cognition_mediation <- lapply(model_cognition_mediation_by_sex, function(x) sem(x, se = "bootstrap", bootstrap = 10000, subset(subsample, sex == 1)))
men_results_cognition_mediation <- lapply(1:5, function(x) filter(standardizedsolution(men_cognition_mediation[[x]]), label %in% c("a", "b", "c", "ab", "total")))
fit_measures_men_cognition_mediation <- lapply(1:5, function(x) fitmeasures(men_cognition_mediation[[x]], c("chisq", "pvalue", "cfi", "tli", "rmsea", "srmr")))
men_results_cognition_mediation <- tibble(bind_rows(men_results_cognition_mediation))
men_results_cognition_mediation %>% select(est.std, se, z, pvalue, ci.lower, ci.upper)
men_results_cognition_mediation <- tibble(bind_rows(men_results_cognition_mediation))
men_results_cognition_mediation <- bind_cols(paste("Men",variables_cognition), men_results_cognition_mediation)
men_results_cognition_mediation <- men_results_cognition_mediation %>% select(...1, est.std, se, z, pvalue, ci.lower, ci.upper)
names(men_results_cognition_mediation) <- mediation_names

women_cognition_mediation <- lapply(model_cognition_mediation_by_sex, function(x) sem(x, se = "bootstrap", bootstrap = 10000, subset(subsample, sex == 0)))
women_results_cognition_mediation <- lapply(1:5, function(x) filter(standardizedsolution(women_cognition_mediation[[x]]), label %in% c("a", "b", "c", "ab", "total")))
fit_measures_women_cognition_mediation <- lapply(1:5, function(x) fitmeasures(women_cognition_mediation[[x]], c("chisq", "pvalue", "cfi", "tli", "rmsea", "srmr")))
women_results_cognition_mediation <- tibble(bind_rows(women_results_cognition_mediation))
women_results_cognition_mediation %>% select(est.std, se, z, pvalue, ci.lower, ci.upper)
women_results_cognition_mediation <- tibble(bind_rows(women_results_cognition_mediation))
women_results_cognition_mediation <- bind_cols(paste("Women",variables_cognition), women_results_cognition_mediation)
women_results_cognition_mediation <- women_results_cognition_mediation %>% select(...1, est.std, se, z, pvalue, ci.lower, ci.upper)
names(women_results_cognition_mediation) <- mediation_names

results_cognition_mediation <- bind_rows(results_cognition_mediation, men_results_cognition_mediation, women_results_cognition_mediation)



##Mediation cognitive tests
variables_cognitive_tests <- c(paste(paste(prefix, rep(label_cognitive_tests, each = 5), sep = " "), rep(label_liver_fat, each=40), sep=" "),
                               paste(paste("Men", prefix, rep(label_cognitive_tests, each = 5), sep = " "), rep(label_liver_fat, each=40), sep=" "),
                               paste(paste("Women",prefix, rep(label_cognitive_tests, each = 5), sep = " "), rep(label_liver_fat, each=40), sep=" "))

model_cognitive_tests_mediation <- 
        c(c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*log_liver_fat + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*log_liver_fat + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))), 
          
          c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*nafld + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*nafld + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))), 
          
          c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*mafld + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*mafld + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))), 
          
          c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*masld + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*masld + + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))), 
          
          c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*metald + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*metald + + std_age_imaging + std_age_squared_imaging + sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)"))))
        )

model_cognitive_tests_mediation_by_sex <- 
        c(c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*log_liver_fat + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*log_liver_fat + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))), 
          
          c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*nafld + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*nafld + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))), 
          
          c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*mafld + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*mafld + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))), 
          
          c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*masld + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*masld + + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)")))), 
          
          c(lapply(cognitive_tests, function(x) c(paste("log_wmh~a*metald + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + std_icv \n",
                                                        x, "~c*metald + + std_age_imaging + std_age_squared_imaging + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education +  intermediate_education \n",
                                                        x, "~b*log_wmh + std_icv \n ab := a*b \n total := c + (a*b)"))))
        )

cognitive_test_mediation <- c(lapply(model_cognitive_tests_mediation, function(x) sem(x, se = "bootstrap", bootstrap = 10000, sample)),
                              lapply(model_cognitive_tests_mediation_by_sex, function(x) sem(x, se = "bootstrap", bootstrap = 10000, subset(sample, sex == 1))),
                              lapply(model_cognitive_tests_mediation_by_sex, function(x) sem(x, se = "bootstrap", bootstrap = 10000, subset(sample, sex == 0))))
results_cognitive_test_mediation <- lapply(1:length(cognitive_test_mediation), function(x) filter(standardizedsolution(cognitive_test_mediation[[x]]), label %in% c("a", "b", "c", "ab", "total")))
fit_measures_cognitive_test_mediation <- lapply(1:length(cognitive_test_mediation), function(x) fitmeasures(cognitive_test_mediation[[x]], c("chisq", "pvalue", "cfi", "tli", "rmsea", "srmr")))
results_cognitive_test_mediation <- tibble(bind_rows(results_cognitive_test_mediation))
results_cognitive_test_mediation <- bind_cols(variables_cognitive_tests, results_cognitive_test_mediation)
results_cognitive_test_mediation <- results_cognitive_test_mediation %>% select(...1, est.std, se, z, pvalue, ci.lower, ci.upper)
names(results_cognitive_test_mediation) <- mediation_names

write_csv(results_cognitive_test_mediation, "cognitive_tests_mediation.csv")

