##Prerequisites
library(tidyverse)

##Load tibble 
sample <- read_csv("sample.csv")
subsample <- read_csv("subsample.csv")

##T-test between men and women
names_t_test_sex <- c("Variables", "Men mean", "Men standard deviation", "Women mean", "Women standard deviation", "T-statistic", "P-value")
variables_t_test <- c("age_baseline", "age_imaging", "bmi", "waist_circumference", "sbp", "dbp", "pulse_pressure",
                      "crp", "hba1c", "hdl", "ldl", "cholesterol", "triglycerides", "liver_fat", "wmh", "alcohol_consumption")
labels_t_test <- c("Age baseline", "Age imaging", "Body mass index", "Waist circumference", "Systolic blood pressure",
                  "Diastolic blood pressure", "Pulse pressure", "C-reactive protein", "Glycated hemoglobin", "High-density lipoprotein cholesterol",
                  "Low-density lipoprotein cholesterol", "Total cholesterol", "Triglycerides", "Liver fat", "White matter hyperintensities",
                  "Alcohol consumption")
t_test_sex <- tibble(as.data.frame((matrix(data = NA, nrow = length(variables_t_test), ncol = length(names_t_test_sex)))))
names(t_test_sex) <- names_t_test_sex
t_test_sex <- t_test_sex %>% mutate(
  `Variables` = labels_t_test,
  `Men mean` = unlist(lapply(variables_t_test, function(x) sample %>% filter(sex == 1) %>% summarise(mean(sample[[x]],na.rm=T)))),
  `Women mean` = unlist(lapply(variables_t_test, function(x) sample %>% filter(sex == 0) %>% summarise(mean(sample[[x]],na.rm=T)))),
  `Men standard deviation` = unlist(lapply(variables_t_test, function(x) sample %>% filter(sex == 1) %>% summarise(sd(sample[[x]],na.rm=T)))),
  `Women standard deviation` = unlist(lapply(variables_t_test, function(x) sample %>% filter(sex == 0) %>% summarise(sd(sample[[x]],na.rm=T)))),
  `T-statistic` = unlist(lapply(variables_t_test, function(x) t.test(subset(sample, sex == 1)[[x]], subset(sample, sex == 0)[[x]], alternative = "two.sided", var.equal = FALSE)$statistic)),
  `P-value`= unlist(lapply(variables_t_test, function(x) t.test(subset(sample, sex == 1)[[x]], subset(sample, sex == 0)[[x]], alternative = "two.sided", var.equal = FALSE)$p.value)))

write_csv(t_test_sex, "output/t_test_sex.csv")

##T-test between total sample and subsample
names_t_test_samples <- c("Variables", "Total sample mean", "Total sample standard deviation", "Subsample mean", "Subsample standard deviation", "T-statistic", "P-value")
t_test_samples <- tibble(as.data.frame((matrix(data = NA, nrow = length(variables_t_test), ncol = length(names_t_test_samples)))))
names(t_test_samples) <- names_t_test_samples
t_test_samples <- t_test_samples %>% mutate(
  "Variables" = labels_t_test,
  "Total sample mean" = unlist(lapply(variables_t_test, function(x) mean(sample[[x]], na.rm=T))),
  "Subsample mean" = unlist(lapply(variables_t_test, function(x) mean(subsample[[x]], na.rm=T))),
  "Total sample standard deviation" = unlist(lapply(variables_t_test, function(x) sd(sample[[x]], na.rm=T))),
  "Subsample standard deviation" = unlist(lapply(variables_t_test, function(x) sd(subsample[[x]], na.rm=T))),
  "T-statistic" = unlist(lapply(variables_t_test, function(x) t.test(sample[[x]], subsample[[x]], alternative = "two.sided", var.equal = FALSE)$statistic)),
  "P-value" = unlist(lapply(variables_t_test, function(x) t.test(sample[[x]], subsample[[x]], alternative = "two.sided", var.equal = FALSE)$p.value)))

write_csv(t_test_samples, "output/t_test_samples.csv")

##Chi-squared between men and women
names_chisq_sex <- c("Variables", "Men number", "Men proportion", "Women number", "Women proportion", "Chi-squared", "P-value")
variables_chisq <- c("asian", "black", "mixed", "other", "white", "higher_education", "intermediate_education", "cheadle", "newcastle", "reading", "bristol",
                     "hypertension", "diabetes", "dyslipidemia", "steatosis", "nafld", "mafld", "masld", "metald", "current_smoker", "former_smoker")
labels_chisq <- c("Asian", "Black", "Mixed", "Other ethnicity", "White", "Higher education", "Intermediate education", "Cheadle", "Newcastle", "Reading", "Bristol",
                  "Hypertension", "Diabetes", "Dyslipidemia", "Steatosis", "NAFLD", "MAFLD", "MASLD", "METALD", "Current smoker", "Former smoker")
chisq_sex <- tibble(as.data.frame((matrix(data = NA, nrow = length(variables_chisq), ncol = length(names_chisq_sex)))))
names(chisq_sex) <- names_chisq_sex
chisq_sex <- chisq_sex %>% mutate(
  "Variables" = labels_chisq,
  "Men number" = unlist(lapply(variables_chisq, function(x) summarise(sample, sum(subset(sample, sex ==1)[[x]] == 1)))),
  "Men proportion" = unlist(lapply(variables_chisq, function(x) summarise(sample, sum(subset(sample, sex ==1)[[x]] == 1))/nrow(subset(sample, sex == 1)))),
  "Women number" = unlist(lapply(variables_chisq, function(x) summarise(sample, sum(subset(sample, sex ==0)[[x]] == 1)))),
  "Women proportion" = unlist(lapply(variables_chisq, function(x) summarise(sample, sum(subset(sample, sex ==0)[[x]] == 1))/nrow(subset(sample, sex == 0)))),
  "Chi-squared" = unlist(lapply(variables_chisq, function(x) chisq.test(table(sample$sex, sample[[x]]))$statistic)),
  "P-value" = unlist(lapply(variables_chisq, function(x) chisq.test(table(sample$sex, sample[[x]]))$p.value))
)

write_csv(chisq_sex, "output/chisquared_sex.csv")

##Chi-squared test between the total sample and subsample
names_chisq_samples <- c("Variables", "Total sample number", "Total sample proportion", "Subsample number", "Subsample proportion", "Chi-squared", "P-value")
chisq_samples <- tibble(as.data.frame((matrix(data = NA, nrow = length(variables_chisq)+1, ncol = length(names_chisq_samples)))))
names(chisq_samples) <- names_chisq_samples
chisq_samples <- chisq_samples %>% mutate(
  Variables = c("Sex",labels_chisq),
  "Total sample number" = unlist(lapply(c("sex",variables_chisq), function(x) sum(sample[[x]], na.rm=T))),
  "Total sample proportion" = unlist(lapply(c("sex",variables_chisq), function(x) sum(sample[[x]], na.rm=T)/nrow(sample))),
  "Subsample number" = unlist(lapply(c("sex",variables_chisq), function(x) sum(subsample[[x]], na.rm=T))),
  "Subsample proportion" = unlist(lapply(c("sex",variables_chisq), function(x) sum(subsample[[x]], na.rm=T)/nrow(subsample))),
  "Chi-squared" = unlist(lapply(c("sex",variables_chisq), function(x) chisq.test(matrix(data = c(sum(sample[[x]]==0), sum(sample[[x]]==1), sum(subsample[[x]]==0), sum(subsample[[x]]==1)), nrow = 2, ncol = 2))$statistic)),
  "P-value" = unlist(lapply(c("sex",variables_chisq), function(x) chisq.test(matrix(data = c(sum(sample[[x]]==0), sum(sample[[x]]==1), sum(subsample[[x]]==0), sum(subsample[[x]]==1)), nrow = 2, ncol = 2))$p.value))
)

write_csv(chisq_samples, "output/chisquared_samples.csv")