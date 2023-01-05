################################ General information ###################################
# Adverse outcomes due to RASi initation after AKI hospitalization
# Roemer Janse - 4 December 2020
# Code for outcome derivation of AKI hospitalization
### -------------------------------------------------------------------------------- ###

timestamp()

rm(list = ls())

pacman::p_load("dplyr", "tidyverse")

setwd("~/Datasets/SCREAM2")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sample.Rdata")
memory.limit(size = 60000)

# To allow a lab follow-up of at least a year, we exclude everyone with a discharge date after 31st of dec 2017 for all outcomes with
# lab components (i.e., hyperkalemia, ckd progression, and eskd). Thus, these outcomes have a cohort up until 31st of dec 2017 while
# the coded outcomes have a cohort up until the 31st of dec 2018.
sample_lab <- sample %>% dplyr::select(lopnr, cohort, index_dt, dis_dt, adm_creat, adm_egfr, censor_dt, censor_dt_lab) %>% 
    filter(dis_dt <= as.Date("2017-12-31"))

sample <- sample %>% dplyr::select(lopnr, cohort, index_dt, adm_creat, adm_egfr, censor_dt, censor_dt_lab)

outcome <- function(dataframe, codes){
    ev_diagn <- dataframe %>% filter(grepl(codes, diagnosis)) %>% dplyr::select(-diagnosis) #keep only the ICD codes we are interested in
    event <- sample %>% left_join(ev_diagn, "lopnr") %>% 
        mutate(event = ifelse(!(is.na(as.Date(bdat))) & as.Date(bdat) > index_dt & as.Date(bdat) <= censor_dt, 1, 0)) %>% 
        arrange(cohort, lopnr, desc(event), as.Date(bdat)) %>% group_by(cohort, lopnr) %>%  slice(1L) %>% ungroup() %>% 
        mutate(event_dt = as.Date(ifelse(event == 1, as.Date(bdat), censor_dt), origin = "1970-01-01"),
               time2event = as.numeric(event_dt - index_dt)) %>% 
        dplyr::select(lopnr, cohort, event, time2event) 
    return(event)
}

outcome_lab <- function(dataframe, codes){
    ev_diagn <- dataframe %>% filter(grepl(codes, diagnosis)) %>% dplyr::select(-diagnosis) #keep only the ICD codes we are interested in
    event <- sample_lab %>% left_join(ev_diagn, "lopnr") %>% 
        mutate(event = ifelse(!(is.na(as.Date(bdat))) & as.Date(bdat) > index_dt & as.Date(bdat) <= censor_dt, 1, 0)) %>% 
        arrange(cohort, lopnr, desc(event), as.Date(bdat)) %>% group_by(cohort, lopnr) %>%  slice(1L) %>% ungroup() %>% 
        mutate(event_dt = as.Date(ifelse(event == 1, as.Date(bdat), censor_dt), origin = "1970-01-01"),
               time2event = as.numeric(event_dt - index_dt)) %>% 
        dplyr::select(lopnr, cohort, event, time2event) 
    return(event)
}

##### 1. Mortality #####
load("death.Rdata")

death <- death %>% arrange(lopnr, dodsdat) %>% group_by(lopnr) %>% slice(1L) %>% ungroup()

oc_mort <- sample %>% left_join(death, "lopnr") %>% 
    mutate(event_mort = ifelse(!(is.na(as.Date(dodsdat))) & as.Date(dodsdat) > index_dt & as.Date(dodsdat) <= censor_dt, 1, 0)) %>% 
    arrange(cohort, lopnr, desc(event_mort), as.Date(dodsdat)) %>% group_by(cohort, lopnr) %>%  slice(1L) %>% ungroup() %>% 
    mutate(event_dt = as.Date(ifelse(event_mort == 1, as.Date(dodsdat), censor_dt), origin = "1970-01-01"),
           time2death = as.numeric(event_dt - index_dt)) %>% 
    dplyr::select(lopnr, cohort, event_mort, time2death) 

##### 2. Recurrent AKI #####
load("slv_diagnoses_long.Rdata")

oc_aki <- outcome(slv_diagnoses_long, "^N17") %>% rename(event_aki = event, time2aki = time2event)

##### 3. Heart failure #####
oc_hf <- outcome(slv_diagnoses_long, "^I110|^I130|^I132|^I50") %>% rename(event_hf = event, time2hf = time2event)

##### 4. End-Stage Kidney Disease #####
load("snr_rrt.Rdata")

load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/diagnoses.Rdata")

relevant_diags <- sample_lab %>% left_join(diagnoses, "lopnr") %>% dplyr::select(-(cohort:censor_dt))
rm(diagnoses)

## ESKD criteria (CKD G5 [ICD-10 codes], start KRT or eGFR < 15 twice consecutively)
# - ICD-10 codes N18.5 
oc_eskd_ICD <- outcome_lab(relevant_diags, "^N185") %>% rename(event_eskd_ICD = event, time2eskd_ICD = time2event)

# - Initiation of KRT
oc_eskd_KRT <- sample_lab %>% left_join(snr_rrt, "lopnr") %>% 
    mutate(krt_dt = as.Date(rrt_date, origin = "1970-01-01"), 
           event_eskd_KRT = ifelse(!(is.na(krt_dt)) & krt_dt > index_dt & krt_dt <= censor_dt & grepl("^HD|^PD|^TX", event_type), 1, 0)) %>% 
    arrange(cohort, lopnr, desc(event_eskd_KRT), krt_dt) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    mutate(event_dt = as.Date(ifelse(event_eskd_KRT == 1, krt_dt, censor_dt), origin = "1970-01-01"),
           time2eskd_KRT = as.numeric(event_dt - index_dt)) %>% 
    dplyr::select(lopnr, cohort, event_eskd_KRT, time2eskd_KRT) 

# - Outpatient eGFR < 15 ml/min/1.73m2
# Using linear regression (see review Kidney International Edouard Fu) & AJKD
# Need date of birth, sex and CKD-EPI formula for eGFR calculation
load("demo.Rdata")
load("lab_values.Rdata")

ckd_epi <- function(creatinine, age, female){
    k <- ifelse(female == 1, 62, 80)
    alpha <- ifelse(female == 1, -0.329, -0.411)
    return(ifelse(female == 1,
                  141 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.209)) * (0.993 ^ age) * 1.018,
                  141 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.209)) * (0.993 ^ age)))
}

labs <- filter(lab_values, (test == "crea" | test == "vencrea" | test == "artcrea") & ip == 0) %>% rename(creat = result)

placeholder <- sample_lab %>% dplyr::select(lopnr, cohort, index_dt, censor_dt, adm_egfr) %>% 
    left_join(labs, "lopnr") %>% 
    left_join(demo %>% dplyr::select(lopnr, dob, female), "lopnr") %>% # For date of birth
    filter(as.Date(datum) >= index_dt & as.Date(datum) <= censor_dt) %>% 
    arrange(cohort, lopnr, datum) %>% 
    mutate(age = round(as.numeric(as.Date(datum) - as.Date(dob)) / 365.25),
           fu_egfr = ckd_epi(creat, age, female)) %>% # Calculate age and eGFR
    group_by(cohort, lopnr, datum) %>% mutate(fu_egfr = mean(fu_egfr)) %>% slice(1L) %>% ungroup() # Calculate average eGFR on one day

placeholder2 <- placeholder %>% filter(fu_egfr < 15)

egfr15_outcome <- placeholder %>% filter(lopnr %in% placeholder2$lopnr) %>% # n = 
    group_by(cohort, lopnr) %>% mutate(month = (as.numeric(as.Date(datum) - index_dt) / 30), # Rescale to 0 months
                                       month = month - min(month)) %>% ungroup() %>% 
    group_by(cohort, lopnr) %>% mutate(max_month = max(month)) %>% ungroup() %>% 
    mutate(intercept = NA, beta = NA, solution = NA)

rm(placeholder,placeholder2)

# Make a loop, select patient one by one, run linear regression, solve linear regression, save output  
egfr15_outcome_c1 <- filter(egfr15_outcome, cohort == 1)
egfr15_outcome_c2 <- filter(egfr15_outcome, cohort == 2)

dat_c1 <- unique(egfr15_outcome_c1$lopnr)
dat_c2 <- unique(egfr15_outcome_c2$lopnr)

for(i in dat_c1){
    MyText <- as.numeric(i)
    subset <- egfr15_outcome_c1[egfr15_outcome_c1$lopnr == MyText, ]
    # linear regression
    lm <- lm(fu_egfr ~ month, data = subset)
    y = 15
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    egfr15_outcome_c1 <- egfr15_outcome_c1 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                      beta = ifelse(lopnr == MyText, b, beta),
                                                      solution = ifelse(lopnr == MyText, x, solution)) 
}

for(i in dat_c2){
    MyText <- as.numeric(i)
    subset <- egfr15_outcome_c2[egfr15_outcome_c2$lopnr == MyText, ]
    # linear regression
    lm <- lm(fu_egfr ~ month, data = subset)
    y = 15
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    egfr15_outcome_c2 <- egfr15_outcome_c2 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                      beta = ifelse(lopnr == MyText, b, beta),
                                                      solution = ifelse(lopnr == MyText, x, solution)) 
}

egfr15_outcome <- rbind(egfr15_outcome_c1, egfr15_outcome_c2)

# Plot check infinite solutions (i.e., slope == 0)
# test <- egfr15_outcome %>% filter(lopnr == "ID0005026")
# 
# model <- lm(egfr_fu ~ month, data = test)
# plot(test$month, test$egfr_fu)
# abline(model)

# Check whether eGFR 15 is crossed before the last eGFR measurement 
# If slope is negative AND eGFR is crossed between baseline and last measurement, then we have an event
# If slope is negative AND eGFR is crossed BEFORE baseline BUT baseline eGFR>15, then we also have an event.
# ... --> There are some individuals who have all measurements under 15, and therefore if you fit a linear line they cross the Y-axis 
# before baseline
# Of course, these individuals also have an event. We set their time2event to 1 day
oc_eskd_GFR <- egfr15_outcome %>% mutate(event_eskd_GFR = ifelse((beta < 0 & solution > 0 & solution < max_month)|
                                                                 (beta < 0 & solution < 0 & adm_egfr > 15), 1, 0)) %>% 
    arrange(cohort, lopnr) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% filter(event_eskd_GFR == 1) %>% 
    mutate(time2eskd_GFR = ifelse(solution > 0, (solution * 30), 1)) %>% 
    dplyr::select(lopnr, cohort, event_eskd_GFR, time2eskd_GFR)

rm(lm, subset, a, b, dat_c1, dat_c2, egfr15_outcome_c1, egfr15_outcome_c2, i, MyText, x, y, egfr15_outcome) 

placeholder <- sample_lab %>% dplyr::select(cohort, lopnr, index_dt, censor_dt_lab)

oc_eskd_GFR <- placeholder %>% left_join(oc_eskd_GFR, c("lopnr", "cohort")) %>% 
    replace_na(list(event_eskd_GFR = 0)) %>% 
    mutate(time2eskd_GFR = ifelse(event_eskd_GFR == 0, as.numeric(censor_dt_lab - index_dt), time2eskd_GFR))

# Create ESKD combined endpoint
oc_eskd <- sample_lab %>% dplyr::select(lopnr, cohort) %>% 
    left_join(oc_eskd_ICD, c("lopnr", "cohort")) %>% 
    left_join(oc_eskd_KRT, c("lopnr", "cohort")) %>% 
    left_join(oc_eskd_GFR, c("lopnr", "cohort")) %>% 
    mutate(event_eskd = pmax(event_eskd_ICD, event_eskd_KRT, event_eskd_GFR, na.rm = TRUE),
           time2eskd = pmin(time2eskd_ICD, time2eskd_KRT, time2eskd_GFR)) %>% 
    dplyr::select(lopnr, cohort, event_eskd, time2eskd)

##### 5. CKD Progression #####
# eGFR decline >30%
placeholder <- sample_lab %>% dplyr::select(lopnr, cohort, index_dt, censor_dt, adm_egfr) %>% 
    left_join(labs, "lopnr") %>% 
    left_join(demo %>% dplyr::select(lopnr, dob, female), "lopnr") %>% # For date of birth
    filter(as.Date(datum) >= index_dt & as.Date(datum) <= censor_dt) %>% 
    arrange(cohort, lopnr, datum) %>% 
    mutate(age = round(as.numeric(as.Date(datum) - as.Date(dob)) / 365.25),
           fu_egfr = ckd_epi(creat, age, female)) %>% # Calculate age and eGFR
    group_by(cohort, lopnr, datum) %>% mutate(fu_egfr = mean(fu_egfr)) %>% slice(1L) %>% ungroup() %>% # Calculate average eGFR on one day
    filter(!(is.na(adm_egfr))) %>% 
    mutate(egfr30decline = ifelse((fu_egfr / adm_egfr) <= 0.7, 1, 0),
           threshold = adm_egfr * 0.7)

placeholder2 <- placeholder %>% filter(egfr30decline == 1) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup()

egfr30decline_outcome <- placeholder %>% filter(lopnr %in% placeholder2$lopnr) %>% # n = 
    group_by(cohort, lopnr) %>% mutate(month = as.numeric(as.Date(datum) - index_dt)/30, # rescale to 0 months
                                       month = month - min(month)) %>% ungroup() %>% 
    group_by(cohort, lopnr) %>% mutate(max_month = max(month)) %>% ungroup() %>% 
    mutate(intercept = NA, beta = NA, solution = NA)

rm(placeholder,placeholder2)

# Make a loop, select patient one by one, run linear regression, solve linear regression, save output  
egfr30decline_outcome_c1 <- filter(egfr30decline_outcome, cohort == 1)
egfr30decline_outcome_c2 <- filter(egfr30decline_outcome, cohort == 2)

dat_c1 <- unique(egfr30decline_outcome_c1$lopnr)
dat_c2 <- unique(egfr30decline_outcome_c2$lopnr)

for(i in dat_c1){
    MyText <- as.numeric(i)
    subset <- egfr30decline_outcome_c1[egfr30decline_outcome_c1$lopnr == MyText, ]
    # linear regression
    lm <- lm(fu_egfr ~ month, data = subset)
    y = subset$threshold
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    egfr30decline_outcome_c1 <- egfr30decline_outcome_c1 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                                    beta = ifelse(lopnr == MyText, b, beta),
                                                                    solution = ifelse(lopnr == MyText, x, solution)) 
}

for(i in dat_c2){
    MyText <- as.numeric(i)
    subset <- egfr30decline_outcome_c2[egfr30decline_outcome_c2$lopnr == MyText, ]
    # linear regression
    lm <- lm(fu_egfr ~ month, data = subset)
    y = subset$threshold
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    egfr30decline_outcome_c2 <- egfr30decline_outcome_c2 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                                    beta = ifelse(lopnr == MyText, b, beta),
                                                                    solution = ifelse(lopnr == MyText, x, solution)) 
}

egfr30decline_outcome <- rbind(egfr30decline_outcome_c1, egfr30decline_outcome_c2)

# Check whether eGFR 30% decline is crossed before the last eGFR measurement 
# If slope is negative AND 30% eGFR decline is crossed between baseline and last measurement, then we have an event
# If slope is negative AND eGFR is crossed BEFORE baseline, then we also have an event.
# In this case, the time2event should be the first measurement that is below 30% eGFR decline
oc_ckd_decl <- egfr30decline_outcome %>% mutate(event_ckd_decl = ifelse((beta < 0 & solution > 0 & solution < max_month)|
                                                                        (beta < 0 & solution < 0), 1, 0)) %>% 
    filter(event_ckd_decl == 1) %>% group_by(cohort, lopnr) %>% 
    mutate(first_under30 = as.Date(min(as.Date(datum)), origin = "1970-01-01")) %>% slice(1L) %>% ungroup() %>%
    # If threshold is crossed before baseline, time2event is equal to first measurement that eGFR>30% is crossed
    mutate(time2ckd_decl = ifelse(solution > 0,solution * 30, as.numeric(first_under30 - index_dt))) %>% 
    dplyr::select(lopnr, cohort, event_ckd_decl,time2ckd_decl)

rm(lm, subset, a, b, egfr30decline_outcome_c1, egfr30decline_outcome_c2, dat_c1, dat_c2, i, MyText, x, y) 

placeholder <- sample_lab %>% dplyr::select(cohort, lopnr, index_dt, censor_dt_lab)

oc_ckd_decl <- placeholder %>% left_join(oc_ckd_decl, c("lopnr", "cohort")) %>% 
    replace_na(list(event_ckd_decl = 0)) %>% 
    mutate(time2ckd_decl = ifelse(event_ckd_decl == 0, as.numeric(censor_dt_lab - index_dt), time2ckd_decl))

# Doubling of baseline serum creatinine
placeholder <- sample_lab %>% dplyr::select(lopnr, cohort, index_dt, censor_dt, adm_creat) %>% 
    left_join(labs,"lopnr") %>% 
    filter(as.Date(datum) >= index_dt & as.Date(datum) <= censor_dt) %>% 
    arrange(cohort, lopnr, datum) %>% 
    group_by(cohort, lopnr, datum) %>% 
    mutate(creat = mean(creat)) %>% 
    slice(1L) %>% ungroup() # calculate average creatinine on one day

placeholder2 <- placeholder %>% mutate(creat_double_outcome = ifelse((creat / adm_creat) >= 2, 1, 0)) %>% 
    arrange(cohort, lopnr, desc(creat_double_outcome), datum) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    filter(!is.na(creat_double_outcome)) %>% # remove those without baseline creatinine measurement 
    filter(creat_double_outcome == 1)

creat_double_outcome <- placeholder %>% filter(lopnr %in% placeholder2$lopnr) %>% # n =
    group_by(cohort, lopnr) %>% mutate(month = as.numeric(as.Date(datum) - index_dt)/30, # rescale to 0 months
                                       month = month-min(month)) %>% ungroup() %>% 
    group_by(cohort, lopnr) %>% mutate(max_month = max(month)) %>% ungroup() %>% 
    mutate(intercept = NA, beta = NA, solution = NA,
           threshold = 2 * adm_creat)

rm(placeholder, placeholder2)

# Make a loop, select patient one by one, run linear regression, solve linear regression, save output  
creat_double_outcome_c1 <- filter(creat_double_outcome, cohort == 1)
creat_double_outcome_c2 <- filter(creat_double_outcome, cohort == 2)

dat_c1 <- unique(creat_double_outcome_c1$lopnr)
dat_c2 <- unique(creat_double_outcome_c2$lopnr)

for(i in dat_c1){
    MyText <- as.numeric(i)
    subset <- creat_double_outcome_c1[creat_double_outcome_c1$lopnr == MyText, ]
    # linear regression
    lm <- lm(creat ~ month, data = subset)
    y = subset$threshold
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    creat_double_outcome_c1 <- creat_double_outcome_c1 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                                  beta = ifelse(lopnr == MyText, b, beta),
                                                                  solution = ifelse(lopnr == MyText, x, solution)) 
}

for(i in dat_c2){
    MyText <- as.numeric(i)
    subset <- creat_double_outcome_c2[creat_double_outcome_c2$lopnr == MyText, ]
    # linear regression
    lm <- lm(creat ~ month, data = subset)
    y = subset$threshold
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    creat_double_outcome_c2 <- creat_double_outcome_c2 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                                  beta = ifelse(lopnr == MyText, b, beta),
                                                                  solution = ifelse(lopnr == MyText, x, solution)) 
}

creat_double_outcome <- rbind(creat_double_outcome_c1, creat_double_outcome_c2)

# Check whether doubling creatinine is crossed before the last creatinine measurement 
# If slope is positive AND doubling creatinine is crossed between baseline and last measurement, then we have an event
oc_ckd_dubl <- creat_double_outcome %>% mutate(event_ckd_dubl = ifelse((beta > 0 & solution > 0 & solution < max_month), 1, 0)) %>% 
    group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    filter(event_ckd_dubl == 1) %>% 
    mutate(time2ckd_dubl = ifelse(solution > 0, solution * 30, 1)) %>% 
    dplyr::select(cohort, lopnr, event_ckd_dubl, time2ckd_dubl) 

rm(lm, subset, a, b, creat_double_outcome_c1, creat_double_outcome_c2, dat_c1, dat_c2, i, MyText, x, y, creat_double_outcome) 

placeholder <- sample_lab %>% dplyr::select(cohort, lopnr, index_dt, censor_dt_lab)

oc_ckd_dubl <- placeholder %>% left_join(oc_ckd_dubl, c("lopnr", "cohort")) %>% 
    replace_na(list(event_ckd_dubl = 0)) %>% 
    mutate(time2ckd_dubl = ifelse(event_ckd_dubl == 0, as.numeric(censor_dt_lab - index_dt), time2ckd_dubl))

oc_ckd <- sample_lab %>% dplyr::select(lopnr, cohort) %>% 
    left_join(oc_ckd_decl, c("lopnr", "cohort")) %>% 
    left_join(oc_ckd_dubl, c("lopnr", "cohort")) %>% 
    mutate(event_ckd = pmax(event_ckd_decl, event_ckd_dubl),
           time2ckd = pmin(time2ckd_decl, time2ckd_dubl)) %>% 
    dplyr::select(lopnr, cohort, event_ckd, time2ckd)

##### 6. Hyperkalemia #####
# Derive all potassium measurements
pot <- lab_values %>% filter((test == "ab_potas" | test == "potas"| test == "other_potas") & ip == 0) %>%
    dplyr::select(lopnr, datum, result) %>%
    mutate(pot_dt = as.Date(datum, origin = "1970-01-01")) %>% 
    rename(pot = result)

oc_hk <- sample_lab %>% left_join(pot, "lopnr") %>% 
    mutate(event_hk_moderate = ifelse(pot >= 5.5 & pot_dt > index_dt & pot_dt <= censor_dt_lab, 1, 0)) %>% 
    replace_na(list(event_hk_moderate = 0)) %>% 
    mutate(event_dt_moderate = as.Date(ifelse(event_hk_moderate == 1, pot_dt, censor_dt_lab), origin = "1970-01-01"),
           time2hk_moderate = as.numeric(event_dt_moderate - index_dt)) %>% 
    arrange(cohort, lopnr, desc(event_hk_moderate), event_dt_moderate) %>% 
    group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    dplyr::select(lopnr, cohort, event_hk_moderate, time2hk_moderate)

##### 7. Major adverse cardiovascular event #####
# Defined as either cardiovascular death, non-fatal stroke or non-fatal myocardial infarction
dods <- death %>% filter(diag_no == "ULORSAK") %>% 
    mutate(bdat = as.Date(dodsdat, origin = "1970-01-01")) %>%
    rename(diagnosis = diagnos)

oc_mace_death <- outcome(dods, "^I") %>% rename(event_mace_death = event, time2mace_death = time2event)
oc_mace_mi <- outcome(slv_diagnoses_long, "^I21|^I22|^I23") %>% rename(event_mace_mi = event, time2mace_mi = time2event)
oc_mace_stroke <- outcome(slv_diagnoses_long, "^I63") %>% rename(event_mace_stroke = event, time2mace_stroke = time2event)

oc_mace <- sample %>% dplyr::select(lopnr, cohort) %>% 
    left_join(oc_mace_death, c("lopnr", "cohort")) %>% 
    left_join(oc_mace_mi, c("lopnr", "cohort")) %>% 
    left_join(oc_mace_stroke, c("lopnr", "cohort")) %>% 
    mutate(event_mace = pmax(event_mace_death, event_mace_mi, event_mace_stroke),
           time2mace = pmin(time2mace_death, time2mace_mi, time2mace_stroke),
           event = ifelse(time2mace == time2mace_mi, "MI",
                          ifelse(time2mace == time2mace_stroke, "Stroke",
                                 ifelse(time2mace == time2mace_death, "Death", NA)))) %>% 
    dplyr::select(lopnr, cohort, event_mace, time2mace, event)

##### 8. Composite of death and MACE #####
oc_comp <- sample %>% dplyr::select(lopnr, cohort) %>%
    left_join(oc_mort, c("lopnr", "cohort")) %>%
    left_join(oc_mace, c("lopnr", "cohort")) %>%
    mutate(event_comp = pmax(event_mort, event_mace),
           time2comp = pmin(time2death, time2mace),
           event_type = ifelse(time2comp == time2death & time2comp == time2mace, "Death",
                               ifelse(time2comp == time2mace & time2comp != time2death, event,
                                      ifelse(time2comp == time2death & time2comp != time2mace, "Death", NA))),
           event_type = ifelse(event_comp == 0, "No event", event_type)) %>%
    dplyr::select(lopnr, cohort, event_comp, time2comp, event_type)

##### 9. Merging outcomes to data #####
# Also calculate time until outcome
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sample.Rdata")

population <- sample %>% left_join(oc_mort, c("lopnr", "cohort")) %>% 
    left_join(oc_aki, c("lopnr", "cohort")) %>%
    left_join(oc_hf, c("lopnr", "cohort")) %>% 
    left_join(oc_eskd, c("lopnr", "cohort")) %>%
    left_join(oc_ckd, c("lopnr", "cohort")) %>% 
    left_join(oc_hk, c("lopnr", "cohort")) %>% 
    left_join(oc_mace, c("lopnr", "cohort")) %>% 
    left_join(oc_comp, c("lopnr", "cohort")) %>%
    mutate(lab_outcome = ifelse(!(is.na(event_ckd)), 1, 0),
           lab_outcome = ifelse(is.na(lab_outcome), 0, lab_outcome))

##### 10. Individual components of composite outcome ####
population <- population %>% mutate(event_comp_death = ifelse(event_type == "Death", 1, 0),
                                    time2comp_death = ifelse(event_comp_death == 1, time2comp,
                                                             as.numeric(censor_dt - index_dt)),
                                    event_comp_mi = ifelse(event_type == "MI", 1, 0),
                                    time2comp_mi = ifelse(event_comp_mi == 1, time2comp,
                                                          as.numeric(censor_dt - index_dt)),
                                    event_comp_stroke = ifelse(event_type == "Stroke", 1, 0),
                                    time2comp_stroke = ifelse(event_comp_stroke == 1, time2comp,
                                                              as.numeric(censor_dt - index_dt)))

save(population, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/population.Rdata")

rm(list = ls())

#######################################################
############ Sensitivity analysis outcomes ############ 
#######################################################
setwd("D:/Datasets/SCREAM2")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sa_landmark6months.Rdata")
memory.limit(size = 60000)

# To allow a lab follow-up of at least a year, we exclude everyone with a discharge date after 31st of dec 2017 for all outcomes with
# lab components (i.e., hyperkalemia, ckd progression, and eskd). Thus, these outcomes have a cohort up until 31st of dec 2017 while
# the coded outcomes have a cohort up until the 31st of dec 2018.
sample_lab <- sa_landmark6months %>% dplyr::select(lopnr, cohort, index_dt, dis_dt, adm_creat, adm_egfr, censor_dt, censor_dt_lab) %>% 
    filter(dis_dt <= as.Date("2017-12-31"))


sample <- sa_landmark6months %>% dplyr::select(lopnr, cohort, index_dt, adm_creat, adm_egfr, censor_dt, censor_dt_lab)

outcome <- function(dataframe, codes){
    ev_diagn <- dataframe %>% filter(grepl(codes, diagnosis)) %>% dplyr::select(-diagnosis) #keep only the ICD codes we are interested in
    event <- sample %>% left_join(ev_diagn, "lopnr") %>% 
        mutate(event = ifelse(!(is.na(as.Date(bdat))) & as.Date(bdat) > index_dt & as.Date(bdat) <= censor_dt, 1, 0)) %>% 
        arrange(cohort, lopnr, desc(event), as.Date(bdat)) %>% group_by(cohort, lopnr) %>%  slice(1L) %>% ungroup() %>% 
        mutate(event_dt = as.Date(ifelse(event == 1, as.Date(bdat), censor_dt), origin = "1970-01-01"),
               time2event = as.numeric(event_dt - index_dt)) %>% 
        dplyr::select(lopnr, cohort, event, time2event) 
    return(event)
}

outcome_lab <- function(dataframe, codes){
    ev_diagn <- dataframe %>% filter(grepl(codes, diagnosis)) %>% dplyr::select(-diagnosis) #keep only the ICD codes we are interested in
    event <- sample_lab %>% left_join(ev_diagn, "lopnr") %>% 
        mutate(event = ifelse(!(is.na(as.Date(bdat))) & as.Date(bdat) > index_dt & as.Date(bdat) <= censor_dt, 1, 0)) %>% 
        arrange(cohort, lopnr, desc(event), as.Date(bdat)) %>% group_by(cohort, lopnr) %>%  slice(1L) %>% ungroup() %>% 
        mutate(event_dt = as.Date(ifelse(event == 1, as.Date(bdat), censor_dt), origin = "1970-01-01"),
               time2event = as.numeric(event_dt - index_dt)) %>% 
        dplyr::select(lopnr, cohort, event, time2event) 
    return(event)
}

##### 1. Mortality #####
load("death.Rdata")

death <- death %>% arrange(lopnr, dodsdat) %>% group_by(lopnr) %>% slice(1L) %>% ungroup()

oc_mort <- sample %>% left_join(death, "lopnr") %>% 
    mutate(event_mort = ifelse(!(is.na(as.Date(dodsdat))) & as.Date(dodsdat) > index_dt & as.Date(dodsdat) <= censor_dt, 1, 0)) %>% 
    arrange(cohort, lopnr, desc(event_mort), as.Date(dodsdat)) %>% group_by(cohort, lopnr) %>%  slice(1L) %>% ungroup() %>% 
    mutate(event_dt = as.Date(ifelse(event_mort == 1, as.Date(dodsdat), censor_dt), origin = "1970-01-01"),
           time2death = as.numeric(event_dt - index_dt)) %>% 
    dplyr::select(lopnr, cohort, event_mort, time2death) 

##### 2. Recurrent AKI #####
load("slv_diagnoses_long.Rdata")

oc_aki <- outcome(slv_diagnoses_long, "^N17") %>% rename(event_aki = event, time2aki = time2event)

##### 3. Heart failure #####
oc_hf <- outcome(slv_diagnoses_long, "^I110|^I130|^I132|^I50") %>% rename(event_hf = event, time2hf = time2event)

##### 4. End-Stage Kidney Disease #####
load("snr_rrt.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/diagnoses.Rdata")

relevant_diags <- sample_lab %>% left_join(diagnoses, "lopnr") %>% dplyr::select(-(cohort:censor_dt))
rm(diagnoses)

## ESKD criteria (CKD G5 [ICD-10 codes], start KRT or eGFR < 15 twice consecutively)
# - ICD-10 codes N18.5 
oc_eskd_ICD <- outcome_lab(relevant_diags, "^N185") %>% rename(event_eskd_ICD = event, time2eskd_ICD = time2event)

# - Initiation of KRT
oc_eskd_KRT <- sample_lab %>% left_join(snr_rrt, "lopnr") %>% 
    mutate(krt_dt = as.Date(rrt_date, origin = "1970-01-01"), 
           event_eskd_KRT = ifelse(!(is.na(krt_dt)) & krt_dt > index_dt & krt_dt <= censor_dt & grepl("^HD|^PD|^TX", event_type), 1, 0)) %>% 
    arrange(cohort, lopnr, desc(event_eskd_KRT), krt_dt) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    mutate(event_dt = as.Date(ifelse(event_eskd_KRT == 1, krt_dt, censor_dt), origin = "1970-01-01"),
           time2eskd_KRT = as.numeric(event_dt - index_dt)) %>% 
    dplyr::select(lopnr, cohort, event_eskd_KRT, time2eskd_KRT) 

# - Outpatient eGFR < 15 ml/min/1.73m2
# Using linear regression (see review Kidney International Edouard Fu) & AJKD
# Need date of birth, sex and CKD-EPI formula for eGFR calculation
load("demo.Rdata")
load("lab_values.Rdata")

ckd_epi <- function(creatinine, age, female){
    k <- ifelse(female == 1, 62, 80)
    alpha <- ifelse(female == 1, -0.329, -0.411)
    return(ifelse(female == 1,
                  141 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.209)) * (0.993 ^ age) * 1.018,
                  141 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.209)) * (0.993 ^ age)))
}

labs <- filter(lab_values, (test == "crea" | test == "vencrea" | test == "artcrea") & ip == 0) %>% rename(creat = result)

placeholder <- sample_lab %>% dplyr::select(lopnr, cohort, index_dt, censor_dt, adm_egfr) %>% 
    left_join(labs, "lopnr") %>% 
    left_join(demo %>% dplyr::select(lopnr, dob, female), "lopnr") %>% # For date of birth
    filter(as.Date(datum) >= index_dt & as.Date(datum) <= censor_dt) %>% 
    arrange(cohort, lopnr, datum) %>% 
    mutate(age = round(as.numeric(as.Date(datum) - as.Date(dob)) / 365.25),
           fu_egfr = ckd_epi(creat, age, female)) %>% # Calculate age and eGFR
    group_by(cohort, lopnr, datum) %>% mutate(fu_egfr = mean(fu_egfr)) %>% slice(1L) %>% ungroup() # Calculate average eGFR on one day

placeholder2 <- placeholder %>% filter(fu_egfr < 15)

egfr15_outcome <- placeholder %>% filter(lopnr %in% placeholder2$lopnr) %>% # n = 
    group_by(cohort, lopnr) %>% mutate(month = (as.numeric(as.Date(datum) - index_dt) / 30), # Rescale to 0 months
                                       month = month - min(month)) %>% ungroup() %>% 
    group_by(cohort, lopnr) %>% mutate(max_month = max(month)) %>% ungroup() %>% 
    mutate(intercept = NA, beta = NA, solution = NA)

rm(placeholder,placeholder2)

# Make a loop, select patient one by one, run linear regression, solve linear regression, save output  
egfr15_outcome_c1 <- filter(egfr15_outcome, cohort == 1)
egfr15_outcome_c2 <- filter(egfr15_outcome, cohort == 2)

dat_c1 <- unique(egfr15_outcome_c1$lopnr)
dat_c2 <- unique(egfr15_outcome_c2$lopnr)

for(i in dat_c1){
    MyText <- as.numeric(i)
    subset <- egfr15_outcome_c1[egfr15_outcome_c1$lopnr == MyText, ]
    # linear regression
    lm <- lm(fu_egfr ~ month, data = subset)
    y = 15
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    egfr15_outcome_c1 <- egfr15_outcome_c1 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                      beta = ifelse(lopnr == MyText, b, beta),
                                                      solution = ifelse(lopnr == MyText, x, solution)) 
}

for(i in dat_c2){
    MyText <- as.numeric(i)
    subset <- egfr15_outcome_c2[egfr15_outcome_c2$lopnr == MyText, ]
    # linear regression
    lm <- lm(fu_egfr ~ month, data = subset)
    y = 15
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    egfr15_outcome_c2 <- egfr15_outcome_c2 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                      beta = ifelse(lopnr == MyText, b, beta),
                                                      solution = ifelse(lopnr == MyText, x, solution)) 
}

egfr15_outcome <- rbind(egfr15_outcome_c1, egfr15_outcome_c2)

# Plot check infinite solutions (i.e., slope == 0)
# test <- egfr15_outcome %>% filter(lopnr == "ID0005026")
# 
# model <- lm(egfr_fu ~ month, data = test)
# plot(test$month, test$egfr_fu)
# abline(model)

# Check whether eGFR 15 is crossed before the last eGFR measurement 
# If slope is negative AND eGFR is crossed between baseline and last measurement, then we have an event
# If slope is negative AND eGFR is crossed BEFORE baseline BUT baseline eGFR>15, then we also have an event.
# ... --> There are some individuals who have all measurements under 15, and therefore if you fit a linear line they cross the Y-axis 
# before baseline
# Of course, these individuals also have an event. We set their time2event to 1 day
oc_eskd_GFR <- egfr15_outcome %>% mutate(event_eskd_GFR = ifelse((beta < 0 & solution > 0 & solution < max_month)|
                                                                     (beta < 0 & solution < 0 & adm_egfr > 15), 1, 0)) %>% 
    arrange(cohort, lopnr) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% filter(event_eskd_GFR == 1) %>% 
    mutate(time2eskd_GFR = ifelse(solution > 0, (solution * 30), 1)) %>% 
    dplyr::select(lopnr, cohort, event_eskd_GFR, time2eskd_GFR)

rm(lm, subset, a, b, dat_c1, dat_c2, egfr15_outcome_c1, egfr15_outcome_c2, i, MyText, x, y, egfr15_outcome) 

placeholder <- sample_lab %>% dplyr::select(cohort, lopnr, index_dt, censor_dt_lab)

oc_eskd_GFR <- placeholder %>% left_join(oc_eskd_GFR, c("lopnr", "cohort")) %>% 
    replace_na(list(event_eskd_GFR = 0)) %>% 
    mutate(time2eskd_GFR = ifelse(event_eskd_GFR == 0, as.numeric(censor_dt_lab - index_dt), time2eskd_GFR))

# Create ESKD combined endpoint
oc_eskd <- sample_lab %>% dplyr::select(lopnr, cohort) %>% 
    left_join(oc_eskd_ICD, c("lopnr", "cohort")) %>% 
    left_join(oc_eskd_KRT, c("lopnr", "cohort")) %>% 
    left_join(oc_eskd_GFR, c("lopnr", "cohort")) %>% 
    mutate(event_eskd = pmax(event_eskd_ICD, event_eskd_KRT, event_eskd_GFR, na.rm = TRUE),
           time2eskd = pmin(time2eskd_ICD, time2eskd_KRT, time2eskd_GFR)) %>% 
    dplyr::select(lopnr, cohort, event_eskd, time2eskd)

##### 5. CKD Progression #####
# eGFR decline >30%
placeholder <- sample_lab %>% dplyr::select(lopnr, cohort, index_dt, censor_dt, adm_egfr) %>% 
    left_join(labs, "lopnr") %>% 
    left_join(demo %>% dplyr::select(lopnr, dob, female), "lopnr") %>% # For date of birth
    filter(as.Date(datum) >= index_dt & as.Date(datum) <= censor_dt) %>% 
    arrange(cohort, lopnr, datum) %>% 
    mutate(age = round(as.numeric(as.Date(datum) - as.Date(dob)) / 365.25),
           fu_egfr = ckd_epi(creat, age, female)) %>% # Calculate age and eGFR
    group_by(cohort, lopnr, datum) %>% mutate(fu_egfr = mean(fu_egfr)) %>% slice(1L) %>% ungroup() %>% # Calculate average eGFR on one day
    filter(!(is.na(adm_egfr))) %>% 
    mutate(egfr30decline = ifelse((fu_egfr / adm_egfr) <= 0.7, 1, 0),
           threshold = adm_egfr * 0.7)

placeholder2 <- placeholder %>% filter(egfr30decline == 1) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup()

egfr30decline_outcome <- placeholder %>% filter(lopnr %in% placeholder2$lopnr) %>% # n = 
    group_by(cohort, lopnr) %>% mutate(month = as.numeric(as.Date(datum) - index_dt)/30, # rescale to 0 months
                                       month = month - min(month)) %>% ungroup() %>% 
    group_by(cohort, lopnr) %>% mutate(max_month = max(month)) %>% ungroup() %>% 
    mutate(intercept = NA, beta = NA, solution = NA)

rm(placeholder,placeholder2)

# Make a loop, select patient one by one, run linear regression, solve linear regression, save output  
egfr30decline_outcome_c1 <- filter(egfr30decline_outcome, cohort == 1)
egfr30decline_outcome_c2 <- filter(egfr30decline_outcome, cohort == 2)

dat_c1 <- unique(egfr30decline_outcome_c1$lopnr)
dat_c2 <- unique(egfr30decline_outcome_c2$lopnr)

for(i in dat_c1){
    MyText <- as.numeric(i)
    subset <- egfr30decline_outcome_c1[egfr30decline_outcome_c1$lopnr == MyText, ]
    # linear regression
    lm <- lm(fu_egfr ~ month, data = subset)
    y = subset$threshold
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    egfr30decline_outcome_c1 <- egfr30decline_outcome_c1 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                                    beta = ifelse(lopnr == MyText, b, beta),
                                                                    solution = ifelse(lopnr == MyText, x, solution)) 
}

for(i in dat_c2){
    MyText <- as.numeric(i)
    subset <- egfr30decline_outcome_c2[egfr30decline_outcome_c2$lopnr == MyText, ]
    # linear regression
    lm <- lm(fu_egfr ~ month, data = subset)
    y = subset$threshold
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    egfr30decline_outcome_c2 <- egfr30decline_outcome_c2 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                                    beta = ifelse(lopnr == MyText, b, beta),
                                                                    solution = ifelse(lopnr == MyText, x, solution)) 
}

egfr30decline_outcome <- rbind(egfr30decline_outcome_c1, egfr30decline_outcome_c2)

# Check whether eGFR 30% decline is crossed before the last eGFR measurement 
# If slope is negative AND 30% eGFR decline is crossed between baseline and last measurement, then we have an event
# If slope is negative AND eGFR is crossed BEFORE baseline, then we also have an event.
# In this case, the time2event should be the first measurement that is below 30% eGFR decline
oc_ckd_decl <- egfr30decline_outcome %>% mutate(event_ckd_decl = ifelse((beta < 0 & solution > 0 & solution < max_month)|
                                                                            (beta < 0 & solution < 0), 1, 0)) %>% 
    filter(event_ckd_decl == 1) %>% group_by(cohort, lopnr) %>% 
    mutate(first_under30 = as.Date(min(as.Date(datum)), origin = "1970-01-01")) %>% slice(1L) %>% ungroup() %>%
    # If threshold is crossed before baseline, time2event is equal to first measurement that eGFR>30% is crossed
    mutate(time2ckd_decl = ifelse(solution > 0,solution * 30, as.numeric(first_under30 - index_dt))) %>% 
    dplyr::select(lopnr, cohort, event_ckd_decl,time2ckd_decl)

rm(lm, subset, a, b, egfr30decline_outcome_c1, egfr30decline_outcome_c2, dat_c1, dat_c2, i, MyText, x, y) 

placeholder <- sample_lab %>% dplyr::select(cohort, lopnr, index_dt, censor_dt_lab)

oc_ckd_decl <- placeholder %>% left_join(oc_ckd_decl, c("lopnr", "cohort")) %>% 
    replace_na(list(event_ckd_decl = 0)) %>% 
    mutate(time2ckd_decl = ifelse(event_ckd_decl == 0, as.numeric(censor_dt_lab - index_dt), time2ckd_decl))

# Doubling of baseline serum creatinine
placeholder <- sample_lab %>% dplyr::select(lopnr, cohort, index_dt, censor_dt, adm_creat) %>% 
    left_join(labs,"lopnr") %>% 
    filter(as.Date(datum) >= index_dt & as.Date(datum) <= censor_dt) %>% 
    arrange(cohort, lopnr, datum) %>% 
    group_by(cohort, lopnr, datum) %>% 
    mutate(creat = mean(creat)) %>% 
    slice(1L) %>% ungroup() # calculate average creatinine on one day

placeholder2 <- placeholder %>% mutate(creat_double_outcome = ifelse((creat / adm_creat) >= 2, 1, 0)) %>% 
    arrange(cohort, lopnr, desc(creat_double_outcome), datum) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    filter(!is.na(creat_double_outcome)) %>% # remove those without baseline creatinine measurement 
    filter(creat_double_outcome == 1)

creat_double_outcome <- placeholder %>% filter(lopnr %in% placeholder2$lopnr) %>% # n =
    group_by(cohort, lopnr) %>% mutate(month = as.numeric(as.Date(datum) - index_dt)/30, # rescale to 0 months
                                       month = month-min(month)) %>% ungroup() %>% 
    group_by(cohort, lopnr) %>% mutate(max_month = max(month)) %>% ungroup() %>% 
    mutate(intercept = NA, beta = NA, solution = NA,
           threshold = 2 * adm_creat)

rm(placeholder, placeholder2)

# Make a loop, select patient one by one, run linear regression, solve linear regression, save output  
creat_double_outcome_c1 <- filter(creat_double_outcome, cohort == 1)
creat_double_outcome_c2 <- filter(creat_double_outcome, cohort == 2)

dat_c1 <- unique(creat_double_outcome_c1$lopnr)
dat_c2 <- unique(creat_double_outcome_c2$lopnr)

for(i in dat_c1){
    MyText <- as.numeric(i)
    subset <- creat_double_outcome_c1[creat_double_outcome_c1$lopnr == MyText, ]
    # linear regression
    lm <- lm(creat ~ month, data = subset)
    y = subset$threshold
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    creat_double_outcome_c1 <- creat_double_outcome_c1 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                                  beta = ifelse(lopnr == MyText, b, beta),
                                                                  solution = ifelse(lopnr == MyText, x, solution)) 
}

for(i in dat_c2){
    MyText <- as.numeric(i)
    subset <- creat_double_outcome_c2[creat_double_outcome_c2$lopnr == MyText, ]
    # linear regression
    lm <- lm(creat ~ month, data = subset)
    y = subset$threshold
    a = lm[["coefficients"]][["(Intercept)"]] # Intercept
    b = lm[["coefficients"]][["month"]] # Slope
    x = (y - a) / b # Solution
    creat_double_outcome_c2 <- creat_double_outcome_c2 %>% mutate(intercept = ifelse(lopnr == MyText, a, intercept),
                                                                  beta = ifelse(lopnr == MyText, b, beta),
                                                                  solution = ifelse(lopnr == MyText, x, solution)) 
}

creat_double_outcome <- rbind(creat_double_outcome_c1, creat_double_outcome_c2)

# Check whether doubling creatinine is crossed before the last creatinine measurement 
# If slope is positive AND doubling creatinine is crossed between baseline and last measurement, then we have an event
oc_ckd_dubl <- creat_double_outcome %>% mutate(event_ckd_dubl = ifelse((beta > 0 & solution > 0 & solution < max_month), 1, 0)) %>% 
    group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    filter(event_ckd_dubl == 1) %>% 
    mutate(time2ckd_dubl = ifelse(solution > 0, solution * 30, 1)) %>% 
    dplyr::select(cohort, lopnr, event_ckd_dubl, time2ckd_dubl) 

rm(lm, subset, a, b, creat_double_outcome_c1, creat_double_outcome_c2, dat_c1, dat_c2, i, MyText, x, y, creat_double_outcome) 

placeholder <- sample_lab %>% dplyr::select(cohort, lopnr, index_dt, censor_dt_lab)

oc_ckd_dubl <- placeholder %>% left_join(oc_ckd_dubl, c("lopnr", "cohort")) %>% 
    replace_na(list(event_ckd_dubl = 0)) %>% 
    mutate(time2ckd_dubl = ifelse(event_ckd_dubl == 0, as.numeric(censor_dt_lab - index_dt), time2ckd_dubl))

oc_ckd <- sample_lab %>% dplyr::select(lopnr, cohort) %>% 
    left_join(oc_ckd_decl, c("lopnr", "cohort")) %>% 
    left_join(oc_ckd_dubl, c("lopnr", "cohort")) %>% 
    mutate(event_ckd = pmax(event_ckd_decl, event_ckd_dubl),
           time2ckd = pmin(time2ckd_decl, time2ckd_dubl)) %>% 
    dplyr::select(lopnr, cohort, event_ckd, time2ckd)

##### 6. Hyperkalemia #####
# Derive all potassium measurements
pot <- lab_values %>% filter((test == "ab_potas" | test == "potas"| test == "other_potas") & ip == 0) %>%
    dplyr::select(lopnr, datum, result) %>%
    mutate(pot_dt = as.Date(datum, origin = "1970-01-01")) %>% 
    rename(pot = result)

oc_hk <- sample_lab %>% left_join(pot, "lopnr") %>% 
    mutate(event_hk_moderate = ifelse(pot >= 5.5 & pot_dt > index_dt & pot_dt <= censor_dt_lab, 1, 0)) %>% 
    replace_na(list(event_hk_moderate = 0)) %>% 
    mutate(event_dt_moderate = as.Date(ifelse(event_hk_moderate == 1, pot_dt, censor_dt_lab), origin = "1970-01-01"),
           time2hk_moderate = as.numeric(event_dt_moderate - index_dt)) %>% 
    arrange(cohort, lopnr, desc(event_hk_moderate), event_dt_moderate) %>% 
    group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    dplyr::select(lopnr, cohort, event_hk_moderate, time2hk_moderate)

##### 7. Major adverse cardiovascular event #####
# Defined as either cardiovascular death, non-fatal stroke or non-fatal myocardial infarction
dods <- death %>% filter(diag_no == "ULORSAK") %>% 
    mutate(bdat = as.Date(dodsdat, origin = "1970-01-01")) %>%
    rename(diagnosis = diagnos)

oc_mace_death <- outcome(dods, "^I") %>% rename(event_mace_death = event, time2mace_death = time2event)
oc_mace_mi <- outcome(slv_diagnoses_long, "^I21|^I22|^I23") %>% rename(event_mace_mi = event, time2mace_mi = time2event)
oc_mace_stroke <- outcome(slv_diagnoses_long, "^I63") %>% rename(event_mace_stroke = event, time2mace_stroke = time2event)

oc_mace <- sample %>% dplyr::select(lopnr, cohort) %>% 
    left_join(oc_mace_death, c("lopnr", "cohort")) %>% 
    left_join(oc_mace_mi, c("lopnr", "cohort")) %>% 
    left_join(oc_mace_stroke, c("lopnr", "cohort")) %>% 
    mutate(event_mace = pmax(event_mace_death, event_mace_mi, event_mace_stroke),
           time2mace = pmin(time2mace_death, time2mace_mi, time2mace_stroke)) %>% 
    dplyr::select(lopnr, cohort, event_mace, time2mace)

##### 8. Composite of death and MACE #####
oc_comp <- sample %>% dplyr::select(lopnr, cohort) %>%
    left_join(oc_mort, c("lopnr", "cohort")) %>%
    left_join(oc_mace, c("lopnr", "cohort")) %>%
    mutate(event_comp = pmax(event_mort, event_mace),
           time2comp = pmin(time2death, time2mace),
           event_type = ifelse(time2comp == time2death & time2comp == time2mace, "Death",
                               ifelse(time2comp == time2mace & time2comp != time2death, "Stroke/MI",
                                      ifelse(time2comp == time2death & time2comp != time2mace, "Death", NA))),
           event_type = ifelse(event_comp == 0, "No event", event_type)) %>%
    dplyr::select(lopnr, cohort, event_comp, time2comp, event_type)

##### 8. Merging outcomes to data #####
# Also calculate time until outcome
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sample.Rdata")

sa_population6months <- sa_landmark6months %>% left_join(oc_mort, c("lopnr", "cohort")) %>% 
    left_join(oc_aki, c("lopnr", "cohort")) %>%
    left_join(oc_hf, c("lopnr", "cohort")) %>% 
    left_join(oc_eskd, c("lopnr", "cohort")) %>%
    left_join(oc_ckd, c("lopnr", "cohort")) %>% 
    left_join(oc_hk, c("lopnr", "cohort")) %>% 
    left_join(oc_mace, c("lopnr", "cohort")) %>% 
    left_join(oc_comp, c("lopnr", "cohort")) %>%
    mutate(lab_outcome = ifelse(!(is.na(event_ckd)), 1, 0),
           lab_outcome = ifelse(is.na(lab_outcome), 0, lab_outcome))

save(sa_population6months, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sa_population6months.Rdata")

timestamp()
