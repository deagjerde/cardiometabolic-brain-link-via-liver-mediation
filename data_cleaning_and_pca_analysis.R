library(tidyverse)
library(missMDA)

##Import UK Biobank data fields
df <- read_csv("imported_data_fields.csv")

##Change the names of the columns
df <- df %>% rename("id" = "eid", "trail_making_test_b" = "6350-2.0", "matrix_test" = "6373-2.0", "tower_rearranging" = "21004-2.0", "symbol_digit_substitution" = "23324-2.0",
              "wmh" = "25781-2.0", "icv" = "26521-2.0", "center" = "54-2.0", "pairs_matching" = "399-2.2", "numeric_memory" = "4282-2.0", "fluid_intelligence" = "20016-2.0", 
              "reaction_time" = "20023-2.0", "paired_associate_learning" = "20197-2.0", "liver_fat" = "24352-2.0", "pilot_liver_fat" = "22436-2.0", "sex" = "31-0.0",
              "waist_circumference" = "48-0.0", "date_baseline" = "53-0.0", "date_imaging" = "53-2.0", "manual_sbp_1" = "93-0.0", "manual_sbp_2" = "93-0.1",
              "manual_dbp_1" = "94-0.0", "manual_dbp_2" = "94-0.1", "alcohol_intake_frequency" = "1558-2.0", "weekly_red_wine" = "1568-2.0", "weekly_white_wine" = "1578-2.0", 
              "weekly_beer" = "1588-2.0", "weekly_spirits" = "1598-2.0","weekly_fortified_wine" = "1608-2.0", "automated_dbp_1" = "4079-0.0", 
              "automated_dbp_2" = "4079-0.1", "automated_sbp_1" = "4080-0.0", "automated_sbp_2" = "4080-0.1", "monthly_red_wine" =  "4407-2.0", 
              "monthly_white_wine" = "4418-2.0", "monthly_beer" = "4429-2.0", "monthly_spirits" = "4440-2.0", "monthly_fortified_wine" = "4451-2.0", 
              "monthly_other_alcohol" = "4462-2.0", "weekly_other_alcohol" = "5364-2.0", "education_1" = "6138-0.0", "education_2" = "6138-0.1", 
              "education_3" = "6138-0.2", "education_4" = "6138-0.3", "education_5" = "6138-0.4", "education_6" = "6138-0.5",
              "pilot_education_1" = "10722-0.0", "pilot_education_2" = "10722-0.1", "pilot_education_3" = "10722-0.2", "pilot_education_4" = "10722-0.3",
              "pilot_education_5" = "10722-0.4", "smoking_status" = "20116-2.0", "alcohol_drinker_status" = "20117-2.0", 
              "ethnicity" = "21000-0.0", "bmi" = "21001-0.0", "age_baseline" = "21003-0.0", "age_imaging" = "21003-2.0", "cholesterol" = "30690-0.0", 
              "crp" = "30710-0.0", "hba1c" = "30750-0.0", "hdl" = "30760-0.0", "ldl" = "30780-0.0", "triglycerides" = "30870-0.0")

names(df) <- gsub("41280", "date_diagnosis", names(df))
names(df) <- gsub("41270", "diagnosis", names(df))
names(df) <- gsub("20003", "medication", names(df))
names(df) <- gsub("-", "_", names(df))
names(df) <- gsub("\\.", "_", names(df))

##Ensure that dates are classified as dates and filter out participants who did not attend imaging assessment
df <- df %>% 
  mutate(across(starts_with("date"), as.Date)) %>% 
  filter(!is.na(date_imaging))

##Check if participants have a diagnosis of exclusion before imaging assessment
dates <- df %>% select(id, date_imaging, starts_with("date_diagnosis"))
diagnosis  <- df %>% select(id, starts_with("diagnosis"))
names(dates) <- c("id", "date_imaging", 1:(length(dates)-2))
names(diagnosis) <- c("id",1:(length(diagnosis)-1))

exclusion <- read_csv("diagnoses_of_exclusion.csv")
exclude <- list()
for(i in 1:(length(diagnosis)-1)){
  exclude[[i]] <- if_else(grepl(paste(paste0("\\<",exclusion$coding,"\\>"), 
  sep = " ", collapse = "|"), diagnosis[[i]]) & dates[[i]] <= dates$date_imaging, 1, 0)
  diagnosis[[paste0("exclude_",i)]] <- exclude[[i]]
}

diagnosis <- diagnosis %>% unite("exclude",starts_with("exclude"),sep=" ")
df$exclude <- if_else(grepl("1",diagnosis$exclude),1,0)

##Merge data fields that have pilot data fields
df <- df %>% mutate(
  liver_fat = coalesce(liver_fat, pilot_liver_fat),
  education_1 = coalesce(education_1, pilot_education_1),
  education_2 = coalesce(education_2, pilot_education_2),
  education_3 = coalesce(education_3, pilot_education_3),
  education_4 = coalesce(education_4, pilot_education_4),
  education_5 = coalesce(education_5, pilot_education_5))

##Create mean systolic and diastolic blood pressure from automated and manual values
df <- df %>% mutate(
  mean_automated_sbp = 
    case_when(!is.na(automated_sbp_1) & !is.na(automated_sbp_2) ~ (automated_sbp_1 + automated_sbp_2) / 2,
    !is.na(automated_sbp_1) ~ automated_sbp_1,
    TRUE ~ automated_sbp_2),
  mean_automated_dbp = 
    case_when(!is.na(automated_dbp_1) & !is.na(automated_dbp_2) ~ (automated_dbp_1 + automated_dbp_2) / 2,
              !is.na(automated_dbp_1) ~ automated_dbp_1,
              TRUE ~ automated_dbp_2),
  mean_manual_sbp = 
    case_when(!is.na(manual_sbp_1) & !is.na(manual_sbp_2) ~ (manual_sbp_1 + manual_sbp_2) / 2,
              !is.na(manual_sbp_1) ~ manual_sbp_1,
              TRUE ~ manual_sbp_2),
  mean_manual_dbp = 
    case_when(!is.na(manual_dbp_1) & !is.na(manual_dbp_2) ~ (manual_dbp_1 + manual_dbp_2) / 2,
              !is.na(manual_dbp_1) ~ manual_dbp_1,
              TRUE ~ manual_dbp_2))
           
df <- df %>%
  mutate(
  sbp = coalesce(mean_automated_sbp, mean_manual_sbp),
  dbp = coalesce(mean_automated_dbp, mean_manual_dbp),
  pulse_pressure = (sbp - dbp))

##Create variables for cardiometabolic medication use
lipid_lowering_drugs <- read_csv("lipid_lowering_drugs_used_in_the_uk_biobank.csv")
antihypertensives <- read_csv("antihypertensive_drugs_used_in_the_uk_biobank.csv")
antidiabetics <- read_csv("antidiabetics_used_in_the_uk_biobank.csv")

df <- df %>% unite("medication",starts_with("medication_0"), sep=" ")
df$lipid_lowering <- if_else(grepl(paste(lipid_lowering_drugs$coding, collapse = "|"), df$medication), 1, 0)
df$antihypertensive <- if_else(grepl(paste(antihypertensives$coding, collapse = "|"), df$medication), 1, 0)
df$antidiabetic <- if_else(grepl(paste(antidiabetics$coding, collapse = "|"), df$medication), 1, 0)

##Replace special values with NA
df <- df %>% mutate(across(.cols = c(grep("weekly|monthly|alcohol|smoking",colnames(df),value = T)), ~na_if(., -3))) %>%
  mutate(across(.cols = c(grep("weekly|monthly|alcohol|smoking",colnames(df),value = T)), ~na_if(., -1))) %>%
  mutate(trail_making_test_b = na_if(trail_making_test_b, 0)) %>%
  mutate(numeric_memory = na_if(numeric_memory, -1))


##Create variable for weekly alcohol unit consumption
##Converts drinks into UK alcohol units (8g alcohol)
df <- df %>% mutate(across(grep("wine",colnames(df)),~ . * 1.5)) %>%
  mutate(across(grep("beer",colnames(df)),~ . * 3)) %>%
  mutate(across(grep("monthly",colnames(df)),~ . / 4.29))

##Based on alcohol intake frequency participants responded for weekly or monthly consumption, merge these fields
df <- df %>% mutate(
  beer = coalesce(weekly_beer, monthly_beer),
  fortified_wine = coalesce(weekly_fortified_wine, monthly_fortified_wine),
  other_alcohol = coalesce(weekly_other_alcohol, monthly_other_alcohol),
  red_wine = coalesce(weekly_red_wine, monthly_red_wine),
  white_wine = coalesce(weekly_white_wine, monthly_white_wine),
  spirits = coalesce(weekly_spirits, monthly_spirits))

df$alcohol_consumption <- rowSums(df[, c("beer", "fortified_wine", "other_alcohol", "red_wine", "white_wine", "spirits")], na.rm = TRUE)

##Convert UK units into grams of alcohol
df$alcohol_consumption <- df$alcohol_consumption*8

##Create variable for participants who lack information on alcohol consumption or have discordant information
df <- df %>% mutate(
  missing_alcohol = if_else(is.na(df$alcohol_intake_frequency) | (df$alcohol_intake_frequency == 1 & df$alcohol_consumption == 0) | 
                              (df$alcohol_intake_frequency == 2 & df$alcohol_consumption == 0) | (df$alcohol_intake_frequency == 3 & df$alcohol_consumption == 0), 1, 0))

##Normalize and standardize continuous variables
df <- df %>%
  mutate(
  std_bmi = ((bmi - mean(bmi, na.rm=T)) / sd(bmi, na.rm=T)),
  std_waist_circumference = ((waist_circumference - mean(waist_circumference, na.rm=T)) / sd(waist_circumference, na.rm=T)),
  std_sbp = ((sbp - mean(sbp, na.rm=T)) / sd(sbp, na.rm=T)),
  std_dbp = ((dbp - mean(dbp, na.rm=T)) / sd(dbp, na.rm=T)),
  std_pulse_pressure = ((pulse_pressure - mean(pulse_pressure, na.rm=T)) / sd(pulse_pressure, na.rm=T)),
  std_hba1c = ((hba1c - mean(hba1c, na.rm=T)) / sd(hba1c, na.rm=T)),
  std_hdl = ((hdl - mean(hdl, na.rm=T)) / sd(hdl, na.rm=T)),
  std_ldl = ((ldl - mean(ldl, na.rm=T)) / sd(ldl, na.rm=T)),
  std_cholesterol = ((cholesterol - mean(cholesterol, na.rm=T)) / sd(cholesterol, na.rm=T)),
  std_icv = ((icv - mean(icv, na.rm=T)) / sd(icv, na.rm=T)),
  
  std_age_imaging = ((age_imaging - mean(age_imaging, na.rm=T)) / sd(age_imaging, na.rm=T)),
  std_alcohol_consumption = ((alcohol_consumption - mean(alcohol_consumption, na.rm=T)) / sd(alcohol_consumption, na.rm=T)),
         
  std_numeric_memory = ((numeric_memory - mean(numeric_memory, na.rm=T)) / sd(numeric_memory, na.rm=T)),
  std_fluid_intelligence = ((fluid_intelligence - mean(fluid_intelligence, na.rm=T)) / sd(fluid_intelligence, na.rm=T)),
  std_matrix_test = ((matrix_test - mean(matrix_test, na.rm=T)) / sd(matrix_test, na.rm=T)),
  std_symbol_digit_substitution = ((symbol_digit_substitution - mean(symbol_digit_substitution, na.rm=T)) / sd(symbol_digit_substitution, na.rm=T)),
  std_tower_rearranging = ((tower_rearranging - mean(tower_rearranging, na.rm=T)) / sd(tower_rearranging, na.rm=T)),
  std_paired_associate_learning = ((paired_associate_learning - mean(paired_associate_learning, na.rm=T)) / sd(paired_associate_learning, na.rm=T)))

df <- df %>%
  mutate(
  log_triglycerides = log(triglycerides + 1),
  log_crp = log(crp + 1),
  log_wmh = log(wmh + 1),
  log_liver_fat = log(liver_fat + 1),
  log_trail_making_test_b = log(trail_making_test_b + 1),
  log_reaction_time = log(reaction_time + 1),
  log_pairs_matching = log(pairs_matching + 1),
  std_age_squared_imaging = std_age_imaging * std_age_imaging)

##Create dummy variables for education, assessment center, ethnic background, smoking, and cardiometabolic/steatotic diagnoses
df <- df %>% unite("education", starts_with("education"), sep = " ", remove = F)

df <- df %>% mutate(
  higher_education = if_else(grepl("1", education), 1, 0),
  intermediate_education = if_else(grepl("2|3", education) & !grepl("1", education), 1, 0)
)

df <- df %>% mutate(
  cheadle = if_else(center == "11025", 1, 0),
  reading = if_else(center == "11026", 1, 0),
  newcastle = if_else(center == "11027", 1, 0),
  bristol = if_else(center == "11028", 1, 0)
)

df <- df %>% mutate(
  asian = if_else(ethnicity %in% c("3001", "3002", "3003", "3004", "3", "5"), 1, 0),
  black = if_else(ethnicity %in% c("4", "4001", "4002", "4003"), 1,0),
  mixed = if_else(ethnicity %in% c("2", "2001", "2002", "2003", "2004"), 1, 0),
  other = if_else(ethnicity %in% c("6"), 1, 0),
  white = if_else(ethnicity %in% c("1", "1001", "1002", "1003"), 1, 0))

df <- df %>% mutate(
  current_smoker = if_else(smoking_status == 2, 1, 0),
  former_smoker = if_else(smoking_status == 1, 1, 0))

df <- df %>% mutate(
  hypertension = if_else(antihypertensive == 1 | sbp >= 140 | dbp >= 90, 1, 0),
  diabetes = if_else(antidiabetic == 1 | hba1c >= 48, 1, 0),
  dyslipidemia = if_else(lipid_lowering == 1 | hdl < 1.03 | ldl > 4.13 | cholesterol >= 6.20 | triglycerides > 2.25, 1, 0)
)

df <- df %>% mutate(
  mafld_criteria = 
    if_else((waist_circumference  >= 102 & sex == 1) | (waist_circumference >= 88 & sex == 0) |
    (waist_circumference  >= 90 & sex == 1 & asian == 1) | (waist_circumference >= 80 & sex == 0 & asian == 1), 1, 0) +
    if_else(triglycerides >= 1.7, 1,0) +
    if_else((hdl < 1 & sex == 1) | (hdl < 1.3 & sex == 0), 1, 0) +
    if_else(lipid_lowering == 1, 1, 0) +
    if_else(hba1c >= 39 & diabetes == 0, 1, 0) +
    if_else(crp >= 2, 1,0) +
    if_else(sbp >= 130 | dbp >= 85 | antihypertensive == 1, 1, 0))

df <- df %>% mutate(
  masld_criteria = 
  if_else((waist_circumference >= 80 & sex == 0) | (waist_circumference >= 94 & sex == 1) | (waist_circumference >= 90 & sex == 1 & asian == 1), 1,0) +
  if_else(bmi >= 25 | (bmi >= 23 & asian == 1), 1, 0) +
  if_else(hba1c >= 39 & hba1c < 48, 1, 0) +
  if_else(sbp >= 130 | dbp >= 85, 1, 0) +
  if_else(triglycerides >= 1.7, 1,0) +
  if_else((hdl < 1 & sex == 1) | (hdl < 1.3 & sex == 0), 1, 0))


df <- df %>% mutate(
  steatosis = if_else(liver_fat >= 5, 1, 0),
  nafld = if_else((steatosis == 1 & sex == 1 & alcohol_consumption <= 210) | (steatosis == 1 & sex == 0 & alcohol_consumption <= 140), 1, 0),
  mafld = if_else((steatosis == 1 & bmi >= 25 & asian == 0) | (steatosis == 1 & bmi >= 23 & asian == 1) |
                    (steatosis == 1 & diabetes == 1)  | (steatosis == 1 & mafld_criteria >= 2), 1, 0),
  masld = if_else((steatosis == 1 & alcohol_consumption <= 210 & sex == 1 & masld_criteria >= 1) | (steatosis == 1 & alcohol_consumption <= 140 & sex == 0 & masld_criteria >= 1), 1, 0),
  metald = if_else((steatosis == 1 & alcohol_consumption > 210 & alcohol_consumption <= 420 & sex == 1 & masld_criteria >= 1) | 
                     (steatosis == 1 & alcohol_consumption > 140 & alcohol_consumption <= 350 & sex == 0 & masld_criteria >= 1), 1, 0))


##Create data set for pariticpants who are eligible (ie non-missing data on WMH, ICV, and liver fat)
eligible <- df %>% filter(!is.na(wmh), !is.na(icv), !is.na(liver_fat))
write_csv(eligible, "output/eligible.csv")

##Information on exclusion based on missing data and diagnoses
flow_chart <- tibble(
  eligible_inclusion = nrow(eligible),
  exclusion_diagnosis = sum(eligible$exclude == 1),
  age_imaging = sum(is.na(eligible$age_imaging)),
  education = sum(is.na(eligible$education)),
  sex = sum(is.na(eligible$sex)),
  center = sum(is.na(eligible$center)),
  bmi = sum(is.na(eligible$bmi)),
  waist_circumference = sum(is.na(eligible$waist_circumference)),
  sbp = sum(is.na(eligible$sbp)),
  dbp = sum(is.na(eligible$dbp)),
  crp = sum(is.na(eligible$crp)),
  hba1c = sum(is.na(eligible$hba1c)),
  hdl = sum(is.na(eligible$hdl)),
  ldl = sum(is.na(eligible$ldl)),
  cholesterol = sum(is.na(eligible$cholesterol)),
  triglycerides = sum(is.na(eligible$triglycerides)),
  alcohol_consumption = sum(eligible$missing_alcohol == 1),
  smoking_status = sum(is.na(eligible$smoking_status)))

write_csv(flow_chart,"output/flow_chart.csv")


##Create data set with participants who do not have diagnoses of exclusion and who are not missing data
sample <- eligible %>%
  filter(exclude == 0, !is.na(age_imaging), !is.na(education), !is.na(sex), missing_alcohol == 0, !is.na(smoking_status),
         !is.na(center), !is.na(sbp), !is.na(dbp), !is.na(bmi), !is.na(waist_circumference), !is.na(hba1c),
         !is.na(crp), !is.na(hdl), !is.na(ldl), !is.na(cholesterol), !is.na(triglycerides))

##Create columns without outliers for PCA analysis
cardiometabolic <- c("std_bmi", "std_waist_circumference", "std_sbp", "std_dbp", "std_pulse_pressure", "log_crp",
                     "std_hba1c", "std_hdl", "std_ldl", "std_cholesterol", "log_triglycerides")

cognitive_tests <-  c("std_numeric_memory", "std_fluid_intelligence", "log_trail_making_test_b", "std_matrix_test", "std_symbol_digit_substitution",
                      "std_tower_rearranging", "std_paired_associate_learning", "log_pairs_matching")

sample <- sample %>%
  mutate(across(
    all_of(c(cognitive_tests, cardiometabolic)),
    list(thresholded = ~ {
      mean_value <- mean(.x, na.rm = TRUE)
      sd_value <- sd(.x, na.rm = TRUE)
      upper_threshold <- mean_value + 3 * sd_value
      lower_threshold <- mean_value - 3 * sd_value
      .x_thresholded <- .x
      .x_thresholded[.x_thresholded > upper_threshold | .x_thresholded < lower_threshold] <- NA
      .x_thresholded
    }),
    .names = "{col}_thresholded"
  ))

missing_data <- sample %>% 
  select(all_of(c(cardiometabolic,paste0(cardiometabolic,"_thresholded"), cognitive_tests, paste0(cognitive_tests, "_thresholded")))) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "NAs")
  
write_csv(missing_data, "output/missing_data.csv")


#Check that no participants have more than 50% missing cardiometabolic data
sample %>%
  filter(if_any(all_of((paste0(cardiometabolic,"_thresholded"))), ~ !is.na(.))) %>%
  nrow()

sample %>%
  filter(if_any(all_of(paste0(cardiometabolic,"_thresholded")), ~ is.na(.x))) %>%
  nrow()

sample %>%
  mutate(missing = rowMeans(is.na(across(all_of(paste0(cardiometabolic,"_thresholded")))))) %>%
  filter(missing <= 0.5) %>%
  select(-missing) %>%
  nrow()

fraction_missing_cardiometabolic_data <- sample %>%
  summarise(across(all_of(paste0(cardiometabolic,"_thresholded")),
                   ~ sum(is.na(.)) / n(),
                   .names = "{col}"
  )) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Fraction NA")

write_csv(missing_data, "output/fraction_missing_cardiometabolic_data.csv")

#Estimate the number of components to include for the imputation analyses
cardiometabolic_component_number <- estim_ncpPCA(select(sample, all_of(paste0(cardiometabolic,"_thresholded"))), ncp.max = 5)

#Imputation. Run algorithm 10 times to avoid convergence to a local minimum
cardiometabolic_imputed <- imputePCA(select(sample, all_of(paste0(cardiometabolic,"_thresholded"))), ncp = cardiometabolic_component_number$ncp, 
                               scale = 10, method = "Regularized", 
                               nb.init = 10, seed = 1,
                               threshold = 1e-06, maxiter = 10000)

#Perform PCA on the imputed data
cardiometabolic_pca <- prcomp(cardiometabolic_imputed$completeObs) 
cardiometabolic_rotation <- as.data.frame(as.matrix(cardiometabolic_pca$rotation, ncol=12))
summary(cardiometabolic_pca)
cardiometabolic_rotation[,1] <- cardiometabolic_rotation[,1] * -1
write_csv(cardiometabolic_rotation, "output/cardiometabolic_pca_rotation.csv")
sample <- bind_cols(sample, cardiometabolic_pca$x[,c("PC1", "PC2", "PC3")])
sample$PC1 <- sample$PC1 *-1

#Create subsample for participants with less than 50% missing cognitive data
sample %>%
  filter(if_any(all_of((paste0(cognitive_tests,"_thresholded"))), ~ !is.na(.))) %>%
  nrow()

subsample <- sample %>%
  filter(if_any(all_of((paste0(cognitive_tests,"_thresholded"))), ~ !is.na(.)))

subsample %>%
  filter(if_any(all_of(paste0(cognitive_tests,"_thresholded")), ~ is.na(.x))) %>%
  nrow()

subsample <- subsample %>%
  mutate(missing = rowMeans(is.na(across(all_of(paste0(cognitive_tests,"_thresholded")))))) %>%
  filter(missing <= 0.5) %>%
  select(-missing)

fraction_missing_cognitive_data <- subsample %>%
  summarise(across(all_of(paste0(cognitive_tests,"_thresholded")),
    ~ sum(is.na(.)) / n(),
    .names = "{col}"
  )) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Fraction NA")
write_csv(fraction_missing_cognitive_data, "output/fraction_missing_cognitive_data.csv")

#Estimate the number of components to include for the imputation analyses
cognition_component_number <- estim_ncpPCA(select(subsample, all_of(paste0(cognitive_tests,"_thresholded"))), ncp.max = 5)

#Imputation. Run algorithm 10 times to avoid convergence to a local minimum
cognition_imputed <- imputePCA(select(subsample, all_of(paste0(cognitive_tests,"_thresholded"))), ncp = cognition_component_number$ncp, 
                               scale = 10, method = "Regularized", 
                               nb.init = 10, seed = 1,
                               threshold = 1e-06, maxiter = 10000)

#Perform PCA on the imputed data
cognition_pca <- prcomp(cognition_imputed$completeObs) 
cognition_rotation <- as.data.frame(as.matrix(cognition_pca$rotation, ncol=9))
summary(cognition_pca)
write_csv(cognition_rotation, "output/cognition_pca_rotation.csv")

subsample <- bind_cols(subsample, cognition_pca$x[,"PC1"])
subsample <- rename(subsample, "cognition" = "...960")

##Compute correlations of cardiometabolic variables
correlations_cardiometabolic <- tibble(as.data.frame(matrix(data = unlist(lapply(cardiometabolic, function(y) lapply(cardiometabolic, function(x) cor.test(sample[[y]], sample[[x]])$estimate))),  
               nrow = length(cardiometabolic), ncol = length(cardiometabolic), byrow = T)))

significance_cardiometabolic <- tibble(as.data.frame(matrix(data = unlist(lapply(cardiometabolic, function(y) lapply(cardiometabolic, function(x) cor.test(sample[[y]], sample[[x]])$p.value))),  
                                                         nrow = length(cardiometabolic), ncol = length(cardiometabolic), byrow = T)))

label_cardiometabolic <- c("Body mass index","Waist circumference","Systolic blood pressure","Diastolic blood pressure","Pulse pressure",
                           "C-reactive protein","Glycated hemoglobin","High-density lipoprotein cholesterol","Low-density lipoprotein cholesterol",
                           "Total cholesterol","Triglycerides")

names(correlations_cardiometabolic) <- label_cardiometabolic
correlations_cardiometabolic <- bind_cols(label_cardiometabolic,correlations_cardiometabolic) %>% rename("Variables" = "...1")

write_csv(correlations_cardiometabolic, "output/correlations_cardiometabolic.csv")

names(significance_cardiometabolic) <- label_cardiometabolic
significance_cardiometabolic <- bind_cols(label_cardiometabolic,significance_cardiometabolic) %>% rename("Variables" = "...1")

write_csv(significance_cardiometabolic, "output/significance_cardiometabolic.csv")

##Compute correlations of cognitive variables
correlations_cognitive_tests <- tibble(as.data.frame(matrix(data = unlist(lapply(cognitive_tests, function(y) lapply(cognitive_tests, function(x) cor.test(sample[[y]], sample[[x]])$estimate))),  
                                                            nrow = length(cognitive_tests), ncol = length(cognitive_tests), byrow = T)))

significance_cognitive_tests <- tibble(as.data.frame(matrix(data = unlist(lapply(cognitive_tests, function(y) lapply(cognitive_tests, function(x) cor.test(sample[[y]], sample[[x]])$p.value))),  
                                                            nrow = length(cognitive_tests), ncol = length(cognitive_tests), byrow = T)))

label_cognitive_tests <- c("Numeric memory", "Fluid intelligence", "Trail making test B", "Matrix test", "Symbol digit substitution", "Tower rearranging", "Paired associate learning", "Pairs matching")

names(correlations_cognitive_tests) <- label_cognitive_tests
correlations_cognitive_tests <- bind_cols(label_cognitive_tests,correlations_cognitive_tests) %>% rename("Variables" = "...1")

write_csv(correlations_cognitive_tests, "output/correlations_cognitive_tests.csv")

names(significance_cognitive_tests) <- label_cognitive_tests
significance_cognitive_tests <- bind_cols(label_cognitive_tests,significance_cognitive_tests) %>% rename("Variables" = "...1")

write_csv(significance_cognitive_tests, "output/significance_cognitive_tests.csv")

##Write finished clean files
write_csv(sample,"output/sample.csv")
write_csv(subsample,"output/subsample.csv")


