################################ General information ###################################
# Adverse outcomes due to RASi discontinuation after AKI hospitalization
# Roemer Janse - 11 June 2021
# Code for outcome derivation of AKI hospitalization in SCREAM2
### -------------------------------------------------------------------------------- ###

rm(list = ls())

pacman::p_load("dplyr", "tidyverse", "readr")

setwd("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes")
memory.limit(size = 60000)


##### 1. Cohort 1 merge #####
load("cohort1.Rdata")
load("covariates_c1.Rdata")

cohort1 <- cohort1 %>% left_join(covariates_c1, c("lopnr", "adm_dt", "dis_dt", "adm_egfr"))

##### 2. Cohort 2 merge #####
load("cohort2.Rdata")
load("covariates_c2.Rdata")

cohort2 <- cohort2 %>% left_join(covariates_c2, c("lopnr", "adm_dt", "dis_dt", "adm_egfr"))

##### 3. Sample merge #####
cohort <- rbind(cohort1, cohort2) %>%
    arrange(cohort, lopnr, adm_dt) %>% 
    dplyr::select(-crit_type.x) %>% 
    rename(crit_type = crit_type.y)

##### 4. Determining exposure & index date #####
setwd("~/Datasets/SCREAM2")

load("migration.Rdata")
load("death.Rdata")
load("lab_values.Rdata")

lopnrs <- cohort %>% arrange(lopnr) %>% group_by(lopnr) %>% slice(1L) %>% ungroup()
load("lmed_1.Rdata")
meds_total_1 <- lopnrs %>% left_join(lmed_1, "lopnr")
rm(lmed_1)
load("lmed_2.Rdata")
meds_total_2 <- lopnrs %>% left_join(lmed_2, "lopnr")
rm(lmed_2)
load("lmed_3.Rdata")
meds_total_3 <- lopnrs %>% left_join(lmed_3, "lopnr")
rm(lmed_3)
load("lmed_4.Rdata")
meds_total_4 <- lopnrs %>% left_join(lmed_4, "lopnr")
rm(lmed_4)
meds_total <- rbind(meds_total_1, meds_total_2, meds_total_3, meds_total_4)
rm(meds_total_1, meds_total_2, meds_total_3, meds_total_4, meds_c1, meds_c2, meds)
 
# Exposure is discontinuing (1), non-exposure is continuing (0)
# Do not sort exposure descending, as 0 indicates continuing RASi, so one 0 is enough to become non-exposed
sample <- cohort %>% mutate(index_dt = dis_dt + 90,
                            censor_dt5yr = (index_dt + (365.25 * 5))) %>%
    left_join(filter(meds_total %>% dplyr::select(lopnr, edatum, atc, antal, antnum), grepl("^C09A|^C09B|^C09C|^C09D", atc)), "lopnr") %>%
    mutate(stopping = ifelse(edatum > dis_dt & edatum <= index_dt, 0, 1)) %>% 
    # Check now if 'stopping' contains any NA: it shouldn't as all individuals are prevalent users
    arrange(cohort, lopnr, stopping) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>%
    # Censor date is first of death (from dods), migration (from demo) or administrative censoring (2019-12-31; 2018 for lab)
    left_join(migration %>% filter(hkod == "U"), "lopnr") %>%
    dplyr::select(-hdat_c) %>% 
    left_join(death %>% dplyr::select(lopnr, dodsdat), "lopnr") %>%
    mutate(censor_dt = as.Date(pmin(as.Date(hdat), as.Date(dodsdat), as.Date("2019-12-31"), censor_dt5yr, 
                                    na.rm = TRUE), origin = "1970-01-01"),
           censor_dt_lab = as.Date(pmin(as.Date(hdat), as.Date(dodsdat), as.Date("2018-12-31"), censor_dt5yr,
                                    na.rm = TRUE), origin = "1970-01-01")) %>%
    arrange(cohort, lopnr) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    filter(censor_dt > index_dt) %>% 
    dplyr::select(-dodsdat, -hdat, -hkod, -antal)

# Cohort 1: n = 10,165
# Cohort 2: n = 5,704

# Before saving we calculate eGFR at day 90 with two definitions:
# - closest to day 90 and after discharge
# - average of all measurements in outpatient care between discharge and day 90
creats <- filter(lab_values, (test == "crea" | test == "vencrea" | test == "artcrea") & ip == 0) %>% dplyr::select(lopnr, datum, result, ip)
info <- dplyr::select(sample, lopnr:cohort, index_dt)

placeholder <- info %>% left_join(creats, "lopnr") %>% filter(as.Date(datum) > dis_dt & as.Date(datum) <= index_dt)
placeholder <- placeholder %>% arrange(cohort, lopnr, datum) %>% group_by(cohort, lopnr) %>% mutate(egfr90_mean = mean(result)) %>% 
    ungroup() %>% mutate(timediff = as.numeric(index_dt - as.Date(datum))) %>% arrange(cohort, lopnr, timediff) %>% 
    group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% rename(egfr90_time = result) %>% 
    dplyr::select(lopnr, cohort, egfr90_mean, egfr90_time)

sample <- sample %>% left_join(placeholder, c("lopnr", "cohort")) %>% relocate(egfr90_mean:egfr90_time, .before = female)

load("~/Datasets/SCREAM2/demo.Rdata")
demo <- demo %>% mutate(dob = as.Date(dob, origin = "1970-01-01"))
sample <- sample %>% left_join(dplyr::select(demo, lopnr, dob), "lopnr") %>% 
    mutate(age = round((as.numeric(adm_dt - dob) / 365.25), digits = 1),
           k = ifelse(female == 1, 62, 80),
           alpha = ifelse(female == 1, -0.329, -0.411),
           egfr90_mean = ifelse(female == 1, 
                             141 * (pmin(egfr90_mean/k, 1) ^ alpha) * (pmax(egfr90_mean/k, 1) ^ (-1.209)) * (0.993 ^ age) * 1.018,
                             141 * (pmin(egfr90_mean/k, 1) ^ alpha) * (pmax(egfr90_mean/k, 1) ^ (-1.209)) * (0.993 ^ age)),
           egfr90_time = ifelse(female == 1, 
                                141 * (pmin(egfr90_time/k, 1) ^ alpha) * (pmax(egfr90_time/k, 1) ^ (-1.209)) * (0.993 ^ age) * 1.018,
                                141 * (pmin(egfr90_time/k, 1) ^ alpha) * (pmax(egfr90_time/k, 1) ^ (-1.209)) * (0.993 ^ age))) %>%
    dplyr::select(-dob, -k, -alpha) 

save(sample, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sample.Rdata")

##### 5. Sample for sensitivity analysis in which landmark is put on 6 months #####
sa_landmark6months <- cohort %>% mutate(index_dt = dis_dt + 180,
                                        censor_dt5yr = (index_dt + (365.25 * 5))) %>%
    left_join(filter(meds_total %>% dplyr::select(lopnr, edatum, atc), grepl("^C09A|^C09B|^C09C|^C09D", atc)), "lopnr") %>%
    mutate(stopping = ifelse(edatum > dis_dt & edatum <= index_dt, 0, 1)) %>% 
    # Check now if 'stopping' contains any NA: it shouldn't as all individuals are prevalent users
    arrange(cohort, lopnr, stopping) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>%
    # Censor date is first of death (from dods), migration (from demo) or administrative censoring (2019-12-31; 2018 for lab)
    left_join(migration %>% filter(hkod == "U"), "lopnr") %>%
    dplyr::select(-hdat_c) %>% 
    left_join(death %>% dplyr::select(lopnr, dodsdat), "lopnr") %>%
    mutate(censor_dt = as.Date(pmin(as.Date(hdat), as.Date(dodsdat), as.Date("2019-12-31"), censor_dt5yr,
                                    na.rm = TRUE), origin = "1970-01-01"),
           censor_dt_lab = as.Date(pmin(as.Date(hdat), as.Date(dodsdat), as.Date("2018-12-31"), censor_dt5yr,
                                    na.rm = TRUE), origin = "1970-01-01")) %>%
    arrange(cohort, lopnr) %>% group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% 
    filter(censor_dt > index_dt) %>% 
    dplyr::select(-dodsdat, -hdat, -hkod)

# Cohort 1: n = 9,409
# Cohort 2: n = 5,282

creats <- filter(lab_values, (test == "crea" | test == "vencrea" | test == "artcrea") & ip == 0) %>% dplyr::select(lopnr, datum, result, ip)
info <- dplyr::select(sa_landmark6months, lopnr:cohort, index_dt)

placeholder <- info %>% left_join(creats, "lopnr") %>% filter(as.Date(datum) > dis_dt & as.Date(datum) <= index_dt)
placeholder <- placeholder %>% arrange(cohort, lopnr, datum) %>% group_by(cohort, lopnr) %>% mutate(egfr90_mean = mean(result)) %>% 
    ungroup() %>% mutate(timediff = as.numeric(index_dt - as.Date(datum))) %>% arrange(cohort, lopnr, timediff) %>% 
    group_by(cohort, lopnr) %>% slice(1L) %>% ungroup() %>% rename(egfr90_time = result) %>% 
    dplyr::select(lopnr, cohort, egfr90_mean, egfr90_time)

sa_landmark6months <- sa_landmark6months %>% left_join(placeholder, c("lopnr", "cohort")) %>% 
    relocate(egfr90_mean:egfr90_time, .before = female)

save(sa_landmark6months, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sa_landmark6months.Rdata")




