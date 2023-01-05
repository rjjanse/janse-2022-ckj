################################ General information ###################################
# Adverse outcomes due to RASi discontinuation after AKI hospitalization
# Roemer Janse - 04 June 2021
# Code for sample derivation of AKI hospitalization in SCREAM2
### -------------------------------------------------------------------------------- ###

rm(list = ls())

pacman::p_load("dplyr", "tidyverse", "readr", "writexl")

setwd("~/Datasets/SCREAM2/")
memory.limit(size = 60000)

####### Random test sample #######
load("slv_buffer.Rdata")

set.seed(27102020)
IDs <- unique(slv_buffer$lopnr)
ID.bs <- sample(IDs, 10000, replace = FALSE) # 10000 is sample size
rows <- which(slv_buffer$lopnr == ID.bs[1])
new.id.bs <- rep(1, length(rows))

for (j in 2:length(ID.bs)){
    regels.id <- which(slv_buffer$lopnr==ID.bs[j])
    rows <- c(rows,regels.id)
    new.id.bs <- c(new.id.bs,rep(j,length(regels.id)))
}

slv_buffer <- slv_buffer[rows,]

rm(ID.bs, IDs, j, new.id.bs, regels.id, rows)

save(slv_buffer, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/random_sample.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/random_sample.Rdata")

####### General data #######
##### 1. Deriving admission and discharge dates #####
load("slv_buffer.Rdata")

slv_buffer <- slv_buffer %>% arrange(lopnr, indat, utdat) %>% group_by(lopnr, indat) %>% slice(1L) %>% ungroup()

slv_in_ut <- slv_buffer %>% dplyr::select(lopnr, indat, utdat) %>% 
    arrange(lopnr, indat, utdat) %>% 
    mutate(indat = as.Date(as.character(indat), format = "%Y%m%d"),
           utdat = as.Date(as.character(utdat), format = "%Y%m%d"))

# Remove missing data, discharge dates prior to admission date and duplicate data
adm_dis_dates <- slv_in_ut %>% filter(!(is.na(indat)) & !(is.na(utdat))) %>% # Drops 0 unique lopnrs
    filter(!(utdat < indat)) %>% # Drops 0 unique lopnrs
    arrange(lopnr, indat, utdat) %>% group_by(lopnr) %>%
    mutate(duplicate_value = ifelse(indat == lag(indat) & utdat == lag(utdat), 1, 0),
           duplicate_value = ifelse(is.na(duplicate_value), 0, duplicate_value)) %>% 
    filter(duplicate_value == 0) %>% # Drops 0 unique lopnrs
    ungroup() %>% 
    dplyr::select(-duplicate_value)

# Some admission and discharge dates do not reflect hospital admission but transfer between departments. This code adjusts for this.
# Create a loop in which each utdat that lies within 1 day of the next indat (so a 1 day of non-hospitalization is also considered a 
# transfer), is equal to that indat, or lies after that indat, is changed to the utdat corresponding to the next indat, given that the 
# first utdat is not later than the second utdat. 
repeat{
    # Determine amount of observations prior to manipulation
    x <- nrow(adm_dis_dates)
    
    # Manipulate data
    adm_dis_dates <- adm_dis_dates %>% arrange(lopnr, indat, utdat) %>% group_by(lopnr) %>% 
        mutate(utdat2 = as.Date(ifelse(as.Date(utdat) >= (lead(as.Date(indat)) - 1) & 
                                       as.Date(utdat) <= lead(as.Date(utdat)),
                                       lead(as.Date(utdat)), as.Date(utdat)), origin = "1970-01-01"),
               utdat2 = as.Date(ifelse(is.na(utdat2), as.Date(utdat), as.Date(utdat2)), origin = "1970-01-01"),
               indat = as.Date(indat, origin = "1970-01-01")) %>%
        ungroup() %>%
        dplyr::select(-utdat) %>% rename(utdat = utdat2) %>%
        arrange(lopnr, utdat, indat) %>% group_by(lopnr, utdat) %>% slice(1L) %>% ungroup()
    
    # Determine amount of observations after manipulation
    y <- nrow(adm_dis_dates)
    
    # Stop loop if no changes are made anymore in dataframe
    if(x == y){break}
}

# n = 1,419,695

save(adm_dis_dates_precheck, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/adm_dis_dates_precheck.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/adm_dis_dates_precheck.Rdata")

# Only between 2007 and 2017, only those 18 years and older, and residents of Stockholm
load("demo.Rdata")
load("migration.Rdata")

migration <- adm_dis_dates_precheck %>% left_join(migration, "lopnr") %>% filter(as.Date(hdat) < indat) %>% 
    arrange(lopnr, indat, desc(hdat)) %>% group_by(lopnr, indat) %>% slice(1L) %>% ungroup() %>% filter(hkod == "U") %>%
    dplyr::select(lopnr, hdat)

adm_dis_dates <- adm_dis_dates_precheck %>% filter(indat >= as.Date("2007-01-01") & 
                                                   utdat <= as.Date("2018-12-31")) %>% # Drops 302,926 unique lopnrs
    left_join(demo, "lopnr") %>% 
    mutate(age = round(as.numeric((indat - as.Date(dob)) / 365.25), digits = 0)) %>% 
    filter(age >= 18) %>% # Drops 130,628 unique lopnrs 
    dplyr::select(lopnr:utdat) %>% 
    left_join(migration, "lopnr") %>% 
    filter(as.Date(hdat) < indat | is.na(hdat)) %>%  # Drops 562 unique lopnrs
    dplyr::select(lopnr:utdat)

# n = 976,814

# check <- adm_dis_dates_precheck %>% filter(!(lopnr %in% adm_dis_dates$lopnr))

# Flowchart:
#       End of follow-up at:    2017             2017                              2018              2018
# - Initial sample size         n = 1,419,695                                      n = 1,419,695     
# - After date boundaries       n = 1,055,469    Drops 364,226 unique lopnrs       n = 1,116,769     Drops 302,926 unique lopnrs
# - After age                   n = 924,841      Drops 130,628 unique lopnrs       n = 977,376       Drops 139,393 unique lopnrs
# - After missing age or sex    n = 924,841      Drops 0 unique lopnrs             n = 977,376       Drops 0 unique lopnrs
# - After migration             n = 922,879      Drops 1,962 unique lopnrs         n = 976,814       Drops 562 unique lopnrs
# 264,261 lopnrs are removed due to admission dates before 2007. For 2017, 121,125 lopnrs are removed due to discharge dates after 2017.
# These two numbers may include duplicates as someone may have been admitted both before 2007 and after 2017.

save(adm_dis_dates, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/adm_dis_dates.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/adm_dis_dates.Rdata")

rm(slv_buffer, slv_in_ut, migration)

##### 2. Deriving ICD-10 code diagnosed AKI #####
load("slv_diagnoses_long.Rdata")

N17 <- slv_diagnoses_long %>% filter(grepl("^N17", diagnosis)) %>% 
    mutate(crit_type = 0,
           bdat = as.Date(bdat, origin = "1970-01-01")) %>% 
    dplyr::select(lopnr, bdat, crit_type) %>% 
    rename(diag_dt = bdat)

N17 <- N17 %>% left_join(adm_dis_dates, "lopnr") %>% 
    filter(diag_dt >= indat & diag_dt <= utdat) %>% # Drops 5,896 unique lopnrs
    arrange(lopnr, indat) %>% group_by(lopnr, indat) %>% slice(1L) %>% ungroup() %>% 
    relocate(lopnr, indat, utdat, diag_dt, crit_type)
    
####### Cohort 1 #######
# Cohort 1 consists of both KDIGO and ICD-10 AKIs
# First get all creatinine measurements into one dataset
load("lab_values.Rdata")

creatinine <- lab_values %>% filter(test == "artcrea" | test == "crea" | test == "vencrea")

## Build a creat dataset with all admission creatinine values and corresponding in hospital creatinine values
# First get all outpatient creatinine values
baseline_creat <- creatinine %>% filter(ip == 0) %>% 
    dplyr::select(lopnr, datum, result) %>% 
    mutate(datum = as.Date(datum, origin = "1970-01-01")) %>% 
    rename(creat_dt = datum, creat = result)

# Then determine all admission creats, defined as the average creatinine in the time window of 365 to 7 days prior to admission
adm_creat <- adm_dis_dates %>% left_join(baseline_creat, "lopnr") %>% 
    filter(creat_dt < (indat - 7) &
           creat_dt > (indat - 365.25)) %>% 
    arrange(lopnr, indat) %>% group_by(lopnr, indat) %>% 
    mutate(adm_creat = mean(creat)) %>% 
    slice(1L) %>% ungroup() %>% 
    dplyr::select(-creat_dt, -creat)

# Get all inpatient creatinine values
hosp_creat <- creatinine %>% filter(ip == 1) %>% 
    dplyr::select(lopnr, datum, result) %>% 
    mutate(datum = as.Date(datum, origin = "1970-01-01")) %>% 
    rename(hosp_creat_dt = datum, hosp_creat = result)

# Reduce multiple creatinine measurements per day to one per day by taking the mean
hosp_creat <- adm_dis_dates %>% left_join(hosp_creat, "lopnr") %>% 
    filter(hosp_creat_dt >= indat & 
           hosp_creat_dt <= utdat) %>% 
    arrange(lopnr, indat, hosp_creat_dt) %>% group_by(lopnr, indat, hosp_creat_dt) %>% 
    mutate(hosp_creat_mean = mean(hosp_creat)) %>% 
    slice(1L) %>% ungroup() %>%  
    mutate(hosp_creat = hosp_creat_mean) %>% 
    dplyr::select(lopnr, indat, utdat, hosp_creat, hosp_creat_dt)

# Add admission creat and hospital creat together
# We join admission creats to hospital creats because a creatinine measurement during hospitalisatoin is needed to determine KDIGO
# stages. Missing admission creats is excluded later with exclusion criteria.
creat <- hosp_creat %>% left_join(adm_creat, c("lopnr", "indat", "utdat")) %>% 
    dplyr::select(lopnr, indat, utdat, adm_creat, hosp_creat, hosp_creat_dt) %>% 
    rename(adm_dt = indat, dis_dt = utdat)

save(creat, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/creat.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/creat.Rdata")

rm(baseline_creat, adm_creat, hosp_creat)

## Add creatinine values to N17 diagnostic codes
N17 <- N17 %>% rename(adm_dt = indat, dis_dt = utdat) %>% 
    left_join(creat, c("lopnr", "adm_dt", "dis_dt")) %>% 
    arrange(lopnr, adm_dt) %>% group_by(lopnr) %>% slice(1L) %>% ungroup() %>% 
    dplyr::select(-diag_dt, -hosp_creat_dt) %>% 
    relocate(crit_type, .after = hosp_creat)

## Determining AKI using KDIGO
# # Icrit1 = first criterium of KDIGO stage I: increase of creat of 26.5 umol/L within two days during admission
# Icrit1 <- creat %>% arrange(lopnr, adm_dt, hosp_creat_dt) %>% group_by(lopnr, adm_dt) %>%
#     mutate(hosp_creat_dt_lag = lag(hosp_creat_dt),
#            hosp_creat_dt_lag2 = lag(hosp_creat_dt_lag),
#            hosp_creat_lag = lag(hosp_creat),
#            hosp_creat_lag2 = lag(hosp_creat_lag)) %>%
#     ungroup() %>%
#     mutate(Icrit1 = ifelse(hosp_creat_dt - hosp_creat_dt_lag >= 0 &
#                                  hosp_creat_dt - hosp_creat_dt_lag <= 2 &
#                                  hosp_creat - hosp_creat_lag >= 26.5, 1,
#                              ifelse(hosp_creat_dt - hosp_creat_dt_lag2 >= 0 &
#                                         hosp_creat_dt - hosp_creat_dt_lag2 <= 2 &
#                                         hosp_creat - hosp_creat_lag2 >= 26.5, 1, 0))) %>%
#     filter(Icrit1 == 1) %>%
#     arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>%
#     dplyr::select(lopnr, adm_dt, dis_dt, adm_creat, hosp_creat) %>%
#     mutate(crit_type = 1)
# 
# # IIrit2 = second criterium of KDIGO stage I: increase of >= 1.5 to 1.9 of baseline creat within 7 days after admission
# Icrit2 <- creat %>% mutate(Icrit2 = ifelse(hosp_creat >= (1.5 * adm_creat) & hosp_creat < (2.0 * adm_creat) &
#                                                    hosp_creat_dt <= (adm_dt + 7) & hosp_creat_dt >= adm_dt, 1, 0)) %>%
#     filter(Icrit2 == 1) %>%
#     arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>%
#     dplyr::select(lopnr, adm_dt, dis_dt, adm_creat, hosp_creat) %>%
#     mutate(crit_type = 1)

# IIcrit1 = first criterium of KDIGO stage II: increase of 2.0 - 2.9 times baseline creat within 7 days after admission
IIcrit1 <- creat %>% mutate(IIcrit1 = ifelse(hosp_creat >= (2.0 * adm_creat) & hosp_creat < (3.0 * adm_creat) &
                                             hosp_creat_dt >= adm_dt & hosp_creat_dt <= (adm_dt + 7), 1, 0)) %>%
    filter(IIcrit1 == 1) %>%
    arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>%
    dplyr::select(lopnr, adm_dt, dis_dt, adm_creat, hosp_creat) %>%
    mutate(crit_type = 2)

# IIIcrit1 = first criterium of KDIGO stage III: increase of creat 353.6 umol/L within two days during admission
# Take note: if the difference between two measurements is one day and the increase is not 353.6 µmol/L or more, there
# might still be an increase during the period of 48 hours. Thus the first ifelse-clause in the second mutate function
# serves to look at the increase with the measurement prior to the current one, but if this does not meet the criterium,
# the second ifelse-clause looks at the day before. If this is within 48 hours, the criterium can be still met.
# IIIcrit1 <- creat %>% arrange(lopnr, adm_dt, hosp_creat_dt) %>% group_by(lopnr, adm_dt) %>%
#     mutate(hosp_creat_dt_lag = lag(hosp_creat_dt),
#            hosp_creat_dt_lag2 = lag(hosp_creat_dt_lag),
#            hosp_creat_lag = lag(hosp_creat),
#            hosp_creat_lag2 = lag(hosp_creat_lag)) %>%
#     ungroup() %>%
#     mutate(IIIcrit1 = ifelse(hosp_creat_dt - hosp_creat_dt_lag >= 0 &
#                              hosp_creat_dt - hosp_creat_dt_lag <= 2 &
#                              hosp_creat - hosp_creat_lag >= 353.6, 1,
#                       ifelse(hosp_creat_dt - hosp_creat_dt_lag2 >= 0 &
#                              hosp_creat_dt - hosp_creat_dt_lag2 <= 2 &
#                              hosp_creat - hosp_creat_lag2 >= 353.6, 1, 0))) %>%
#     filter(IIIcrit1 == 1) %>%
#     arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>%
#     dplyr::select(lopnr, adm_dt, dis_dt, adm_creat, hosp_creat) %>%
#     mutate(crit_type = 3)

# IIIcrit1 = first criterium of KDIGO stage III: increase to creat 353.6 umol/L. This differs from the code above, but the official criterium
# is to 353.6, not a difference of 353.6. In this paper, it did not make a difference in the included population.
IIIcrit1 <- creat %>%
    mutate(IIIcrit1 = ifelse(hosp_creat >= 353.6, 1, 0)) %>%
    filter(IIIcrit1 == 1) %>%
    arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>%
    dplyr::select(lopnr, adm_dt, dis_dt, adm_creat, hosp_creat) %>%
    mutate(crit_type = 3)

# IIIcrit2 = second criterium of KDIGO stage III: increase of >= 3.0 of baseline creat within 7 days after admission
IIIcrit2 <- creat %>% mutate(IIIcrit2 = ifelse(hosp_creat >= (3.0 * adm_creat) & hosp_creat_dt <= (adm_dt + 7) &
                                               hosp_creat_dt >= adm_dt, 1, 0)) %>%
    filter(IIIcrit2 == 1) %>%
    arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>%
    dplyr::select(lopnr, adm_dt, dis_dt, adm_creat, hosp_creat) %>%
    mutate(crit_type = 3)

# IIIcrit3 = third criterum of KDIGO stage III: need for kidney replacement therapy during hospitalisation
# Here we use admission and discharge dates instead of creat, because we don't need creat. So to get the full scope
# of all patients, we want to check for kidney replacement therapy in any person
load("slv_procedurecodes.Rdata")

IIIcrit3 <- creat %>% left_join(filter(slv_procedurecodes, grepl("^DR015|^DR023", opk)), "lopnr") %>%
    mutate(IIIcrit3 = ifelse(as.Date(datum) >= adm_dt & as.Date(datum) <= dis_dt, 1, 0)) %>%
    filter(IIIcrit3 == 1) %>%
    arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>%
    dplyr::select(lopnr, adm_dt, dis_dt, adm_creat, hosp_creat) %>% 
    mutate(crit_type = 3)

## Identifying earlieast AKI
# First slice: keep one row per hospitalisation
# Second slice: keep one hospitalisation per unique person
AKIs <- rbind(IIcrit1, IIIcrit1, IIIcrit2, IIIcrit3, N17) %>%
    arrange(lopnr, adm_dt, desc(crit_type)) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>% 
    arrange(lopnr, adm_dt) %>% group_by(lopnr) %>% slice(1L) %>% ungroup()

save(AKIs, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/AKIs_cohort1.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/AKIs_cohort1.Rdata")

rm(IIcrit1, IIIcrit1, IIIcrit2, IIIcrit3, slv_procedurecodes)

# n = 41,077

## Calculate eGFR
load("demo.Rdata")

egfr <- function(){
    demo <- demo %>% mutate(dob = as.Date(dob, origin = "1970-01-01"))
    AKIs <- AKIs %>% left_join(dplyr::select(demo, lopnr, dob, female), "lopnr") %>% 
        mutate(age = round((as.numeric(adm_dt - dob) / 365.25), digits = 1),
               k = ifelse(female == 1, 62, 80),
               alpha = ifelse(female == 1, -0.329, -0.411),
               adm_egfr = ifelse(female == 1, 
                                 141 * (pmin(adm_creat/k, 1) ^ alpha) * (pmax(adm_creat/k, 1) ^ (-1.209)) * (0.993 ^ age) * 1.018,
                                 141 * (pmin(adm_creat/k, 1) ^ alpha) * (pmax(adm_creat/k, 1) ^ (-1.209)) * (0.993 ^ age))) %>%
        dplyr::select(-age, -dob, -female, -k, -alpha) %>% relocate(adm_egfr, .before = hosp_creat:crit_type)
    return(AKIs)
}

AKI <- egfr()

load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/diagnoses.Rdata")
specific_diagnoses <- dplyr::select(AKI, lopnr) %>% left_join(diagnoses, "lopnr")
load("lab_values.Rdata")

comorbidities <- function(icd_code, df){   
    comorb <- AKI %>% left_join(df, "lopnr") %>% 
        filter(grepl(icd_code, diagnosis), as.Date(bdat) <= dis_dt & !is.na(bdat)) %>%
        mutate(comorb = 1) %>% 
        arrange(lopnr, adm_dt, desc(comorb)) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% 
        ungroup() %>% 
        dplyr::select(lopnr, comorb, bdat) 
    return(comorb) 
}

# Albuminuria together with DM is defined as > 20 & < 200 mg/L & only albuminuria ≥ 200 mg/L
# https://www.farmacotherapeutischkompas.nl/vergelijken/indicatieteksten?vergelijkTeksten=chronische_nierschade 
# Here too comorbidities are included up until discharge
# - Cardiovascular disease
# - Heart failure
# - Diabetes with eGFR ≤ 60
# - Diabetes with ualb ≥ 20 | ACR ≥ 3
# - Hypertension with eGFR ≤ 60
# - Albuminuria: ualb ≥ 200 | ACR ≥ 30
indication <- function(){
    # Cardiovascular disease
    cvd_comorb <- comorbidities("^I21|^I22|^I252|^I70|^I71|^I72|^I73|^G459|^I63|^I64|^I69", specific_diagnoses) %>% 
        rename(cvd_comorb = comorb, cvd_dt = bdat) 
    
    # Heart failure
    hf_comorb <- comorbidities("^I110|^I130|^I132|^I50", specific_diagnoses) %>% rename(hf_comorb = comorb, hf_dt = bdat) 
    
    # Diabetes mellitus
    dm_comorb <- comorbidities("^E10|^E11|^E12|^E13|^E14", specific_diagnoses) %>% rename(dm_comorb = comorb, dm_dt = bdat)
    
    # Hypertension
    ht_comorb <- comorbidities("^I10|^I11|^I12|^I13|^I15", specific_diagnoses) %>% rename(ht_comorb = comorb, ht_dt = bdat)
    
    # Derive urine albuminuria
    ualb <- lab_values %>% filter(test == "ualb") %>%
        dplyr::select(lopnr, result, datum) %>%
        rename(ualb_dt = datum, ualb = result) %>% 
        mutate(ualb_dt = as.Date(ualb_dt, origin = "1970-01-01"))
    
    # Derive albuine creatinine ratios
    acr <- lab_values %>% filter(test == "uacr") %>% 
        dplyr::select(lopnr, result, datum) %>% 
        rename(acr_dt = datum, acr = result) %>% 
        mutate(acr_dt = as.Date(acr_dt, origin = "1970-01-01"))
    
    # Get indications    
    sample <- AKI %>% left_join(cvd_comorb, "lopnr") %>% left_join(hf_comorb, "lopnr") %>% left_join(dm_comorb, "lopnr") %>%
        left_join(ht_comorb, "lopnr") %>% left_join(ualb, "lopnr") %>% left_join(acr, "lopnr") %>% 
        mutate(dm_egfr = ifelse(dm_comorb == 1 & adm_egfr <= 60, 1, 0), 
               dm_ualb = ifelse(dm_comorb == 1 & ualb >= 20 & as.Date(ualb_dt) < (adm_dt - 7) & as.Date(ualb_dt) > (adm_dt - 365), 1, 0),
               dm_acr = ifelse(dm_comorb == 1 & acr >= 3 & as.Date(acr_dt) < (adm_dt - 7) & as.Date(acr_dt) > (adm_dt - 365), 1, 0),
               dm_alb = ifelse(dm_ualb == 1 | dm_acr == 1, 1, 0),
               ht_egfr = ifelse(ht_comorb == 1 & adm_egfr <= 60, 1, 0),
               ualb_c = ifelse(ualb >= 200 & as.Date(ualb_dt) < (adm_dt - 7) & as.Date(ualb_dt) > (adm_dt - 365), 1, 0),
               acr_c = ifelse(acr >= 30 & as.Date(acr_dt) < (adm_dt - 7) & as.Date(acr_dt) > (adm_dt - 365), 1, 0),
               alb = ifelse(ualb_c == 1 | acr_c == 1, 1, 0)) %>% 
        replace_na(list(cvd_comorb = 0, hf_comorb = 0, dm_comorb = 0, dm_egfr = 0, ht_egfr = 0, dm_alb = 0, alb = 0)) %>%
        arrange(lopnr, adm_dt, desc(cvd_comorb), desc(hf_comorb), desc(dm_egfr), desc(ht_egfr), desc(dm_alb), desc(alb)) %>%
        group_by(lopnr) %>% slice(1L) %>% ungroup() %>%
        filter(cvd_comorb == 1 | hf_comorb == 1 | dm_egfr == 1 | dm_alb == 1 | ht_egfr == 1 | alb == 1) %>%
        dplyr::select(lopnr:crit_type, cvd_comorb, hf_comorb, dm_egfr, dm_alb, ht_egfr, alb) 
    return(sample)
}

sample <- indication()

# n =  26,955

save(sample, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sample_cohort1.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sample_cohort1.Rdata")

rm(slv_diagnoses, ovr_diagnoses, kontakt_diagnoses, diagnoses, AKIs)

## Determining prevalent RASi users
lopnrs <- dplyr::select(sample, lopnr)

load("lmed_1.Rdata")
meds_1 <- lopnrs %>% left_join(lmed_1, "lopnr")
rm(lmed_1)

load("lmed_2.Rdata")
meds_2 <- lopnrs %>% left_join(lmed_2, "lopnr")
rm(lmed_2)

load("lmed_3.Rdata")
meds_3 <- lopnrs %>% left_join(lmed_3, "lopnr")
rm(lmed_3)

load("lmed_4.Rdata")
meds_4 <- lopnrs %>% left_join(lmed_4, "lopnr")
rm(lmed_4)

meds <- rbind(meds_1, meds_2, meds_3, meds_4) %>% 
    arrange(lopnr) 

rm(meds_1, meds_2, meds_3, meds_4)

save(meds, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/meds_c1.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/meds_c1.Rdata")

RASi_status <- function(){
    RASi_use <- sample %>% left_join(filter(meds, grepl("^C09A|^C09B|^C09C|^C09D", atc)), "lopnr") %>%
        dplyr::select(-(otyp:edatum_c)) %>%
        mutate(prev = ifelse(edatum > (adm_dt - 365.25) & edatum < adm_dt, 1, 0)) %>%
        replace_na(list(prev = 0)) %>%
        arrange(lopnr, desc(prev)) %>% group_by(lopnr) %>% slice(1L) %>% ungroup() %>%
        filter(prev == 1) %>%
        dplyr::select(lopnr:crit_type, cvd_comorb, hf_comorb, dm_egfr, dm_alb, ht_egfr, alb)
    return(RASi_use)
}

RASi_use <- RASi_status()

# n = 17,371

rm(sample)

save(RASi_use, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/RASi.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/RASi.Rdata")

## Applying exclusion criteria
load("snr_rrt.Rdata")
load("death.Rdata")

death <- death %>% arrange(lopnr, dodsdat) %>% group_by(lopnr) %>% slice(1L) %>% ungroup() %>% 
    dplyr::select(lopnr, dodsdat)

exclusion <- function(){
    RASi_use <- RASi_use %>% filter(adm_egfr >= 15 | is.na(adm_egfr))
    n1 <- length(unique(RASi_use$lopnr))
    RASi_use <- RASi_use %>% left_join(snr_rrt, "lopnr") %>% arrange(lopnr, rrt_date) %>%
        mutate(krt = ifelse(rrt == 1 & as.Date(rrt_date) < adm_dt, 1, 0)) %>%
        replace_na(list(krt = 0)) %>%
        arrange(lopnr, adm_dt, desc(krt)) %>% group_by(lopnr) %>% slice(1L) %>% ungroup() %>%
        filter(krt == 0) %>%
        dplyr::select(-rrt_date, -event_type, -krt)
    n2 <- length(unique(RASi_use$lopnr))
    t1 <- table(RASi_use$crit_type)
    RASi_use <- RASi_use %>% filter(!(is.na(adm_creat)))
    n3 <- length(unique(RASi_use$lopnr))
    t2 <- table(RASi_use$crit_type)
    RASi_use <- RASi_use %>% left_join(death, "lopnr") %>% 
        mutate(discharge_death = ifelse(dodsdat <= dis_dt, 1, 0)) %>% 
        replace_na(list(discharge_death = 0)) %>% 
        filter(discharge_death == 0)
    n4 <- length(unique(RASi_use$lopnr))
    return(list(RASi_use, n1, n2, n3, n4))
}

list <- exclusion()

RASi_use <- list[[1]]

cohort1 <- RASi_use %>% dplyr::select(lopnr, adm_dt, dis_dt, adm_creat, adm_egfr, crit_type, cvd_comorb, hf_comorb, dm_egfr, dm_alb, 
                                      ht_egfr, alb) %>% mutate(cohort = 1)

# n = 12,140

save(cohort1, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/cohort1.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/cohort1.Rdata")

rm(AKI, AKIs, creatinine, sample, specific_diagnoses, slv_diagnoses_long, meds)

####### Cohort 2 #######
## Adding eGFR to AKI's 
AKIs <- N17

AKI <- egfr()

# n = 20,363

## Getting indications
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/diagnoses.Rdata")

specific_diagnoses <- AKI %>% dplyr::select(lopnr) %>% left_join(diagnoses, "lopnr") %>% 
    mutate(bdat = as.Date(bdat, origin = "1970-01-01"))

rm(diagnoses)

sample <- indication()

save(sample, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sample_cohort2.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/sample_cohort2.Rdata")

# n = 14,239

## Extracting prevalent RASi users
lopnrs <- dplyr::select(sample, lopnr)

load("lmed_1.Rdata")
meds_1 <- lopnrs %>% left_join(lmed_1, "lopnr")
rm(lmed_1)

load("lmed_2.Rdata")
meds_2 <- lopnrs %>% left_join(lmed_2, "lopnr")
rm(lmed_2)

load("lmed_3.Rdata")
meds_3 <- lopnrs %>% left_join(lmed_3, "lopnr")
rm(lmed_3)

load("lmed_4.Rdata")
meds_4 <- lopnrs %>% left_join(lmed_4, "lopnr")
rm(lmed_4)

meds <- rbind(meds_1, meds_2, meds_3, meds_4) %>% 
    arrange(lopnr) 

rm(meds_1, meds_2, meds_3, meds_4)

save(meds, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/meds_c2.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/meds_c2.Rdata")

RASi_use <- RASi_status()

# n = 9,743

save(RASi_use, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/RASi_c2.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/RASi_c2.Rdata")

## Appying exclusion criteria
list <- exclusion()

RASi_use <- list[[1]]

cohort2 <- RASi_use %>% dplyr::select(lopnr:dis_dt, adm_creat, adm_egfr, crit_type, cvd_comorb, hf_comorb, dm_egfr, dm_alb, ht_egfr, 
                                      alb) %>% mutate(cohort = 2)

# n = 6,697

save(cohort2, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/cohort2.Rdata")
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/cohort2.Rdata")










