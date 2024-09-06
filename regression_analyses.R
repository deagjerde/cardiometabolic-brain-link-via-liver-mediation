##Prerequisites
library(tidyverse)

##Load tibble 
sample <- read_csv("sample.csv")
subsample <- read_csv("subsample.csv")

##Create character vectors
cardiometabolic <- c("std_bmi", "std_waist_circumference","std_sbp", "std_dbp","std_pulse_pressure","log_crp", 
                     "std_hba1c","std_hdl","std_ldl", "std_cholesterol","log_triglycerides","PC1", "PC2", "PC3")

label_cardiometabolic <- c("Body mass index","Waist circumference","Systolic blood pressure","Diastolic blood pressure","Pulse pressure",
                           "C-reactive protein","Glycated hemoglobin","High-density lipoprotein cholesterol","Low-density lipoprotein cholesterol",
                           "Total cholesterol","Triglycerides","Principal component 1", "Principal component 2", "Principal component 3")

liver_fat <- c("log_liver_fat", "nafld", "mafld", "masld", "metald")

label_liver_fat <- c("Liver fat", "NAFLD", "MAFLD", "MASLD", "METALD")

cognitive_tests <- c("std_numeric_memory", "std_fluid_intelligence", "log_trail_making_test_b", "std_matrix_test", 
                     "std_symbol_digit_substitution", "std_tower_rearranging", "std_paired_associate_learning", "log_pairs_matching")

label_cognitive_tests <- c("Numeric memory", "Fluid intelligence", "Trail making test B", "Matrix test", "Symbol digit substitution", "Tower rearranging", "Paired associate learning", "Pairs matching")

regression_names <- c("Variable","Regression coefficient","Standard error", "T-statistic", "P-value", "Degrees of freedom")

##Create regression models
model_liver_regression <-  lapply(cardiometabolic, function(x) c(paste("log_liver_fat ~", x, 
  '+ poly(std_age_imaging, 2)*sex + reading + newcastle + bristol + current_smoker + former_smoker + std_alcohol_consumption'))) 

model_liver_regression_sex <-  lapply(cardiometabolic, function(x) c(paste("log_liver_fat ~", x, 
  '*sex + poly(std_age_imaging, 2)*sex + reading + newcastle + bristol + current_smoker + former_smoker + std_alcohol_consumption'))) 

model_wmh_regression <- lapply(c(cardiometabolic, liver_fat), function(x) c(paste("log_wmh ~", x,
  '+ poly(std_age_imaging, 2)*sex + reading + newcastle + bristol + current_smoker + former_smoker + std_alcohol_consumption')))

model_wmh_regression_sex <- lapply(c(cardiometabolic, liver_fat), function(x) c(paste("log_wmh ~", x,
  '*sex + poly(std_age_imaging, 2)*sex + reading + newcastle + bristol + current_smoker + former_smoker + std_alcohol_consumption')))   
                                                                       
model_cognition_regression <- lapply(c(cardiometabolic, liver_fat, "log_wmh"), function(x) c(paste("cognition ~", x, 
  '+ poly(std_age_imaging, 2)*sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education + intermediate_education')))

model_cognition_regression[[20]] <- paste0(model_cognition_regression[[20]], " + std_icv", collapse = " ")

model_cognition_regression_sex <- lapply(c(cardiometabolic, liver_fat, "log_wmh"), function(x) c(paste("cognition ~", x, 
  '*sex + poly(std_age_imaging, 2)*sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education + intermediate_education')))

##Edit model for WMH as it requires adjustment for ICV
model_cognition_regression_sex[[20]] <- paste0(model_cognition_regression_sex[[20]], " + std_icv", collapse = " ")

model_cognitive_tests_regression <- c(lapply(cognitive_tests, function(x) c(paste(x, "~", 
  paste0(c(cardiometabolic, liver_fat), '+ poly(std_age_imaging, 2)*sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education + intermediate_education')))),
  lapply(cognitive_tests, function(x) c(paste(x, '~ log_wmh + poly(std_age_imaging, 2)*sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education + intermediate_education + std_icv'))))

model_cognitive_tests_regression_sex <- c(lapply(cognitive_tests, function(x) c(paste(x, "~", 
  paste0(c(cardiometabolic, liver_fat), '*sex + poly(std_age_imaging, 2)*sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education + intermediate_education')))),
  lapply(cognitive_tests, function(x) c(paste(x, '~ log_wmh*sex + poly(std_age_imaging, 2)*sex + reading + newcastle + current_smoker + former_smoker + std_alcohol_consumption + higher_education + intermediate_education + std_icv'))))

##Liver fat regression
liver_regression <- lapply(model_liver_regression, function(x) lm(x, sample))
results_liver_regression <- lapply(1:length(liver_regression), function(x) c(summary(liver_regression[[x]])$coef[2,], summary(liver_regression[[x]])$df[2]))
results_liver_regression <- bind_rows(results_liver_regression)
results_liver_regression <- bind_cols(label_cardiometabolic, results_liver_regression)
names(results_liver_regression) <- regression_names
results_liver_regression <- results_liver_regression %>% mutate("Partial correlation coefficient" = `T-statistic`/sqrt(`T-statistic`*`T-statistic`+ `Degrees of freedom`))
results_liver_regression <- results_liver_regression %>% mutate("Fischer Z" = 0.5*log((1 + `Partial correlation coefficient`)/(1 - `Partial correlation coefficient`)))
results_liver_regression <- results_liver_regression %>% mutate("Lower confidence interval" = `Fischer Z` - 1.96*(1/sqrt(nrow(sample)-3)),
                                                                "Upper confidence interval" = `Fischer Z` + 1.96*(1/sqrt(nrow(sample)-3)))
results_liver_regression <- results_liver_regression %>% mutate("Lower confidence interval" = (exp(2*`Lower confidence interval`) - 1) / (exp(2*`Lower confidence interval`) + 1),
                                                                "Upper confidence interval" = (exp(2*`Upper confidence interval`) - 1) / (exp(2*`Upper confidence interval`) + 1)) %>%
                                                                  select(-c(`Fischer Z`, `Degrees of freedom`))
write_csv(results_liver_regression, "output/results_liver_regression.csv")

##liver fat regression sex interaction
liver_regression_sex <- lapply(model_liver_regression_sex, function(x) lm(x, sample))
results_liver_regression_sex <- lapply(1:length(liver_regression_sex), function(x) c(summary(liver_regression_sex[[x]])$coef[12,], summary(liver_regression_sex[[x]])$df[2]))
results_liver_regression_sex <- bind_rows(results_liver_regression_sex)
results_liver_regression_sex <- bind_cols(label_cardiometabolic, results_liver_regression_sex)
names(results_liver_regression_sex) <- regression_names
results_liver_regression_sex <- results_liver_regression_sex %>% mutate("Partial correlation coefficient" = `T-statistic`/sqrt(`T-statistic`*`T-statistic`+ `Degrees of freedom`))
results_liver_regression_sex <- results_liver_regression_sex %>% mutate("Fischer Z" = 0.5*log((1 + `Partial correlation coefficient`)/(1 - `Partial correlation coefficient`)))
results_liver_regression_sex <- results_liver_regression_sex %>% mutate("Lower confidence interval" = `Fischer Z` - 1.96*(1/sqrt(nrow(sample)-3)),
                                                                "Upper confidence interval" = `Fischer Z` + 1.96*(1/sqrt(nrow(sample)-3)))
results_liver_regression_sex <- results_liver_regression_sex %>% mutate("Lower confidence interval" = (exp(2*`Lower confidence interval`) - 1) / (exp(2*`Lower confidence interval`) + 1),
                                                                "Upper confidence interval" = (exp(2*`Upper confidence interval`) - 1) / (exp(2*`Upper confidence interval`) + 1)) %>%
  select(-c(`Fischer Z`, `Degrees of freedom`))
write_csv(results_liver_regression_sex, "output/results_liver_regression_sex.csv")

##WMH regression
wmh_regression <- lapply(model_wmh_regression, function(x) lm(x, sample))
results_wmh_regression <- lapply(1:length(wmh_regression), function(x) c(summary(wmh_regression[[x]])$coef[2,], summary(wmh_regression[[x]])$df[2]))
results_wmh_regression <- bind_rows(results_wmh_regression)
results_wmh_regression <- bind_cols(c(label_cardiometabolic, label_liver_fat), results_wmh_regression)
names(results_wmh_regression) <- regression_names
results_wmh_regression <- results_wmh_regression %>% mutate("Partial correlation coefficient" = `T-statistic`/sqrt(`T-statistic`*`T-statistic`+ `Degrees of freedom`))
results_wmh_regression <- results_wmh_regression %>% mutate("Fischer Z" = 0.5*log((1 + `Partial correlation coefficient`)/(1 - `Partial correlation coefficient`)))
results_wmh_regression <- results_wmh_regression %>% mutate("Lower confidence interval" = `Fischer Z` - 1.96*(1/sqrt(nrow(sample)-3)),
                                    "Upper confidence interval" = `Fischer Z` + 1.96*(1/sqrt(nrow(sample)-3)))
results_wmh_regression <- results_wmh_regression %>% mutate("Lower confidence interval" = (exp(2*`Lower confidence interval`) - 1) / (exp(2*`Lower confidence interval`) + 1),
                                          "Upper confidence interval" = (exp(2*`Upper confidence interval`) - 1) / (exp(2*`Upper confidence interval`) + 1)) %>%
  select(-c(`Fischer Z`, `Degrees of freedom`))
write_csv(results_wmh_regression, "output/results_wmh_regression.csv")

##WMH regression sex interaction
wmh_regression_sex <- lapply(model_wmh_regression_sex, function(x) lm(x, sample))
results_wmh_regression_sex <- lapply(1:length(wmh_regression_sex), function(x) c(summary(wmh_regression_sex[[x]])$coef[12,], summary(wmh_regression_sex[[x]])$df[2]))
results_wmh_regression_sex <- bind_rows(results_wmh_regression_sex)
results_wmh_regression_sex <- bind_cols(c(label_cardiometabolic, label_liver_fat), results_wmh_regression_sex)
names(results_wmh_regression_sex) <- regression_names
results_wmh_regression_sex <- results_wmh_regression_sex %>% mutate("Partial correlation coefficient" = `T-statistic`/sqrt(`T-statistic`*`T-statistic`+ `Degrees of freedom`))
results_wmh_regression_sex <- results_wmh_regression_sex %>% mutate("Fischer Z" = 0.5*log((1 + `Partial correlation coefficient`)/(1 - `Partial correlation coefficient`)))
results_wmh_regression_sex <- results_wmh_regression_sex %>% mutate("Lower confidence interval" = `Fischer Z` - 1.96*(1/sqrt(nrow(sample)-3)),
                                                            "Upper confidence interval" = `Fischer Z` + 1.96*(1/sqrt(nrow(sample)-3)))
results_wmh_regression_sex <- results_wmh_regression_sex %>% mutate("Lower confidence interval" = (exp(2*`Lower confidence interval`) - 1) / (exp(2*`Lower confidence interval`) + 1),
                                                            "Upper confidence interval" = (exp(2*`Upper confidence interval`) - 1) / (exp(2*`Upper confidence interval`) + 1)) %>%
  select(-c(`Fischer Z`, `Degrees of freedom`))
write_csv(results_wmh_regression_sex, "output/results_wmh_regression_sex.csv")

##Cognition regression
cognition_regression <- lapply(model_cognition_regression, function(x) lm(x, subsample))
results_cognition_regression <- lapply(1:length(cognition_regression), function(x) c(summary(cognition_regression[[x]])$coef[2,], summary(cognition_regression[[x]])$df[2]))
results_cognition_regression <- bind_rows(results_cognition_regression)
results_cognition_regression <- bind_cols(c(label_cardiometabolic, label_liver_fat, "WMH"), results_cognition_regression)
names(results_cognition_regression) <- regression_names
results_cognition_regression <- results_cognition_regression %>% mutate("Partial correlation coefficient" = `T-statistic`/sqrt(`T-statistic`*`T-statistic`+ `Degrees of freedom`))
results_cognition_regression <- results_cognition_regression %>% mutate("Fischer Z" = 0.5*log((1 + `Partial correlation coefficient`)/(1 - `Partial correlation coefficient`)))
results_cognition_regression <- results_cognition_regression %>% mutate("Lower confidence interval" = `Fischer Z` - 1.96*(1/sqrt(nrow(sample)-3)),
                                                            "Upper confidence interval" = `Fischer Z` + 1.96*(1/sqrt(nrow(sample)-3)))
results_cognition_regression <- results_cognition_regression %>% mutate("Lower confidence interval" = (exp(2*`Lower confidence interval`) - 1) / (exp(2*`Lower confidence interval`) + 1),
                                                            "Upper confidence interval" = (exp(2*`Upper confidence interval`) - 1) / (exp(2*`Upper confidence interval`) + 1)) %>%
  select(-c(`Fischer Z`, `Degrees of freedom`))

write_csv(results_cognition_regression, "output/results_cognition_regression.csv")

##Cognition regression sex interaction
cognition_regression_sex <- lapply(model_cognition_regression_sex, function(x) lm(x, subsample))
results_cognition_regression_sex <- lapply(1:length(cognition_regression_sex), function(x) c(summary(cognition_regression_sex[[x]])$coef[13,], summary(cognition_regression_sex[[x]])$df[2]))
results_cognition_regression_sex[[20]] <- c(summary(cognition_regression_sex[[20]])$coef[14,], summary(cognition_regression_sex[[20]])$df[2])
results_cognition_regression_sex <- bind_rows(results_cognition_regression_sex)
results_cognition_regression_sex <- bind_cols(c(label_cardiometabolic, label_liver_fat, "WMH"), results_cognition_regression_sex)
names(results_cognition_regression_sex) <- regression_names
results_cognition_regression_sex <- results_cognition_regression_sex %>% mutate("Partial correlation coefficient" = `T-statistic`/sqrt(`T-statistic`*`T-statistic`+ `Degrees of freedom`))
results_cognition_regression_sex <- results_cognition_regression_sex %>% mutate("Fischer Z" = 0.5*log((1 + `Partial correlation coefficient`)/(1 - `Partial correlation coefficient`)))
results_cognition_regression_sex <- results_cognition_regression_sex %>% mutate("Lower confidence interval" = `Fischer Z` - 1.96*(1/sqrt(nrow(sample)-3)),
                                                                        "Upper confidence interval" = `Fischer Z` + 1.96*(1/sqrt(nrow(sample)-3)))
results_cognition_regression_sex <- results_cognition_regression_sex %>% mutate("Lower confidence interval" = (exp(2*`Lower confidence interval`) - 1) / (exp(2*`Lower confidence interval`) + 1),
                                                                        "Upper confidence interval" = (exp(2*`Upper confidence interval`) - 1) / (exp(2*`Upper confidence interval`) + 1)) %>%
  select(-c(`Fischer Z`, `Degrees of freedom`))

write_csv(results_cognition_regression_sex, "output/results_cognition_regression_sex.csv")


##Individual cognitive tests regression
cognitive_tests_regression <- lapply(unlist(model_cognitive_tests_regression), function(x) lm(x, sample))
results_cognitive_tests_regression <- lapply(1:length(cognitive_tests_regression), function(x) c(summary(cognitive_tests_regression[[x]])$coef[2,], summary(cognitive_tests_regression[[x]])$df[2]))
results_cognitive_tests_regression <- bind_rows(results_cognitive_tests_regression)
results_cognitive_tests_regression <- bind_cols(c(paste(rep(label_cognitive_tests,each=19), c(label_cardiometabolic, label_liver_fat)), paste(label_cognitive_tests, "WMH")),
                                                  results_cognitive_tests_regression)
names(results_cognitive_tests_regression) <- c(regression_names)

sample_sizes <- sample %>% 
  select(all_of(cognitive_tests)) %>%
  summarise(across(everything(), ~ sum(!is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Sample size") %>% 
  select('Sample size')

results_cognitive_tests_regression$'Sample sizes' <- c(rep(unlist(sample_sizes), each = 19), unlist(sample_sizes))

results_cognitive_tests_regression <- results_cognitive_tests_regression %>% mutate("Partial correlation coefficient" = `T-statistic`/sqrt(`T-statistic`*`T-statistic`+ `Degrees of freedom`))
results_cognitive_tests_regression <- results_cognitive_tests_regression %>% mutate("Fischer Z" = 0.5*log((1 + `Partial correlation coefficient`)/(1 - `Partial correlation coefficient`)))
results_cognitive_tests_regression <- results_cognitive_tests_regression %>% mutate("Lower confidence interval" = `Fischer Z` - 1.96*(1/sqrt(`Sample sizes`-3)),
                                                                                    "Upper confidence interval" = `Fischer Z` + 1.96*(1/sqrt(`Sample sizes`-3)))
results_cognitive_tests_regression <- results_cognitive_tests_regression %>% 
  mutate("Lower confidence interval" = (exp(2*`Lower confidence interval`) - 1) / (exp(2*`Lower confidence interval`) + 1),
         "Upper confidence interval" = (exp(2*`Upper confidence interval`) - 1) / (exp(2*`Upper confidence interval`) + 1)) %>%
  select(-c(`Fischer Z`, `Degrees of freedom`, `Sample sizes`))

write_csv(results_cognitive_tests_regression, "output/results_cognitive_tests_regression.csv")

##Individual cognitive tests regression sex interaction
cognitive_tests_regression_sex <- lapply(unlist(model_cognitive_tests_regression_sex), function(x) lm(x, sample))
results_cognitive_tests_regression_sex <- lapply(1:length(cognitive_tests_regression_sex), function(x) c(summary(cognitive_tests_regression_sex[[x]])$coef[13,], summary(cognitive_tests_regression_sex[[x]])$df[2]))
results_cognitive_tests_regression_sex <- bind_rows(results_cognitive_tests_regression_sex)
results_cognitive_tests_regression_sex <- bind_cols(c(paste(rep(label_cognitive_tests,each=19), c(label_cardiometabolic, label_liver_fat)), paste(label_cognitive_tests, "WMH")),
                                                    results_cognitive_tests_regression_sex)
names(results_cognitive_tests_regression_sex) <- c(regression_names)
for(x in 153:160){
  results_cognitive_tests_regression_sex[x,2] <- summary(cognitive_tests_regression_sex[[x]])$coef[14,1]
  results_cognitive_tests_regression_sex[x,3] <- summary(cognitive_tests_regression_sex[[x]])$coef[14,2]
  results_cognitive_tests_regression_sex[x,4] <- summary(cognitive_tests_regression_sex[[x]])$coef[14,3]
  results_cognitive_tests_regression_sex[x,5] <- summary(cognitive_tests_regression_sex[[x]])$coef[14,4]
  results_cognitive_tests_regression_sex[x,6] <- summary(cognitive_tests_regression_sex[[x]])$df[2]
}

sample_sizes <- sample %>% 
  select(all_of(cognitive_tests)) %>%
  summarise(across(everything(), ~ sum(!is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Sample size") %>% 
  select('Sample size')

results_cognitive_tests_regression_sex$'Sample sizes' <- c(rep(unlist(sample_sizes), each = 19), unlist(sample_sizes))

results_cognitive_tests_regression_sex <- results_cognitive_tests_regression_sex %>% mutate("Partial correlation coefficient" = `T-statistic`/sqrt(`T-statistic`*`T-statistic`+ `Degrees of freedom`))
results_cognitive_tests_regression_sex <- results_cognitive_tests_regression_sex %>% mutate("Fischer Z" = 0.5*log((1 + `Partial correlation coefficient`)/(1 - `Partial correlation coefficient`)))
results_cognitive_tests_regression_sex <- results_cognitive_tests_regression_sex %>% mutate("Lower confidence interval" = `Fischer Z` - 1.96*(1/sqrt(`Sample sizes`-3)),
                                                                                    "Upper confidence interval" = `Fischer Z` + 1.96*(1/sqrt(`Sample sizes`-3)))
results_cognitive_tests_regression_sex <- results_cognitive_tests_regression_sex %>% 
  mutate("Lower confidence interval" = (exp(2*`Lower confidence interval`) - 1) / (exp(2*`Lower confidence interval`) + 1),
         "Upper confidence interval" = (exp(2*`Upper confidence interval`) - 1) / (exp(2*`Upper confidence interval`) + 1)) %>%
  select(-c(`Fischer Z`, `Degrees of freedom`, `Sample sizes`))

write_csv(results_cognitive_tests_regression_sex, "output/results_cognitive_tests_regression_sex.csv")

