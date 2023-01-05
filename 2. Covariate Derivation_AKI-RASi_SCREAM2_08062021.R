################################ General information ###################################
# Adverse outcomes due to RASi discontinuation after AKI hospitalization
# Roemer Janse - 08 June 2021
# Code for covariate derivation of AKI hospitalization in SCREAM2
### -------------------------------------------------------------------------------- ###

timestamp()

rm(list = ls())

pacman::p_load("dplyr", "tidyverse")

setwd("~/Datasets/SCREAM2")
memory.limit(size = 60000)


####### Cohort 1 #######
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/cohort1.Rdata")

lopnrs <- cohort1 %>% dplyr::select(lopnr, adm_dt, dis_dt, adm_egfr)

##### 1. General information #####
load("demo.Rdata")
load("education.Rdata")

general_information <- function(){
    # Age is determined at discharge
    gi <- lopnrs %>% left_join(demo, "lopnr") %>%
        dplyr::select(lopnr, adm_dt, dis_dt, adm_egfr, dob, female) %>%
        mutate(age = as.numeric(round(((dis_dt - as.Date(dob)) / 365.25),
                                      digits = 0))) %>% 
        dplyr::select(-dob)
    
    # Age categories
    # 1 = < 50
    # 2 = 50-59
    # 3 = 60-69
    # 4 = 70-79
    # 5 = >= 80
    gi <- gi %>% mutate(age_cat = ifelse(age < 50, 1,
                                         ifelse(age >= 50 & age < 60, 2,
                                                ifelse(age >= 60 & age < 70, 3,
                                                       ifelse(age >= 70 & age < 80, 4,
                                                              ifelse(age >= 80, 5, NA))))))
    
    # eGFR categories
    gi <- gi %>% mutate(adm_egfr_cat = ifelse(adm_egfr >= 90, 1,
                                          ifelse(adm_egfr >= 60 & adm_egfr < 90, 2,
                                                 ifelse(adm_egfr >= 45 & adm_egfr < 60, 3,
                                                        ifelse(adm_egfr >= 30 & adm_egfr < 45, 4,
                                                               ifelse(adm_egfr >= 15 & adm_egfr < 30, 5,
                                                                      ifelse(adm_egfr < 15, 6, NA))))))) %>% 
        dplyr::select(-adm_egfr)
    
    # Education
    # See article by Judvigsson: https://link.springer.com/article/10.1007/s10654-019-00511-8
    # Sun2000niva.old has 7 levels, see legend of Figure 2 by Ludvigsson
    # 1 = Compulsory school, less than 9 years (did not complete compulsory education [<9 years])
    # 2 = Compulsory school, 9 years (completed compulsory education [9 years])
    # 3 = Secondary school (upper secondary [2 years])
    # 4 = Secondary school (upper secondary [3 years])
    # 5 = University (college/university < 3 years)
    # 6 = University (college/university >= 3 years)
    # 7 = University (research education)
    # Coded:
    # 1 = Compulsory school
    # 2 = Secondary school
    # 3 = University
    # 4 = missing
    education <- education %>% mutate(education = cut(sun2000niva_old, c(0, 3, 5, 10), 
                                                      labels = c("Compulsory school", "Secondary school", "University"), right = F),
                                      educ_cat = ifelse(education == "Compulsory school", 1,
                                                 ifelse(education == "Secondary school", 2,
                                                 ifelse(education == "University", 3, NA)))) %>% 
        dplyr::select(lopnr, exam_year, educ_cat)
    
    educ <- gi %>% left_join(education, "lopnr") %>% mutate(dis_year = as.numeric(format(as.Date(dis_dt, format = "%d/%m/%Y"), "%Y")),
                                                            exam_year = as.numeric(exam_year),
                                                            exam_year = ifelse(is.na(exam_year), 0, exam_year)) %>% 
        filter(exam_year <= dis_year) %>% 
        dplyr::select(lopnr, educ_cat)
    
    gi <- gi %>% left_join(educ, "lopnr") %>% replace_na(list(educ_cat = 4)) %>% dplyr::select(-dis_dt)
    
    return(gi)
}

gi_c1 <- general_information()
    
##### 2. Comorbidities #####
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/diagnoses.Rdata")

relevant_diags <- dplyr::select(lopnrs, lopnr) %>% left_join(diagnoses, "lopnr")

rm(diagnoses)

comorbidities <- function(icd_code, df){   
    comorb <- lopnrs %>% left_join(df, "lopnr") %>% 
        filter(grepl(icd_code, diagnosis), as.Date(bdat) <= dis_dt & !is.na(bdat)) %>%
        mutate(comorb = 1) %>% 
        arrange(lopnr, adm_dt, desc(comorb)) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% 
        ungroup() %>% 
        dplyr::select(lopnr, comorb, bdat) %>%
        rename(datum = bdat)
    return(comorb) 
}

comorbs <- function(){
    # Diabetes mellitus (dm)
    # Codes from Fu et al. Stopping Renin-Angiotensin System Inhibitors in Patients with Advanced CKD and Risk of Adverse Outcomes: 
    # A Nationwide Study. JASN 2021.
    comorb_dm <- comorbidities("^E10|^E11|^E12|^E13|^E14", relevant_diags) %>%
        rename(comorb_dm = comorb, datum_dm = datum)
    
    # Hypertension (ht)
    # Codes from Fu et al. Stopping Renin-Angiotensin System Inhibitors in Patients with Advanced CKD and Risk of Adverse Outcomes: 
    # A Nationwide Study. JASN 2021.
    comorb_ht <- comorbidities("^I10|^I11|^I12|^I13|^I15", relevant_diags) %>%
        rename(comorb_ht = comorb, datum_ht = datum)
    
    # Myocardial infarction (mi)
    # Codes from Fu et al. Stopping Renin-Angiotensin System Inhibitors in Patients with Advanced CKD and Risk of Adverse Outcomes: 
    # A Nationwide Study. JASN 2021.
    comorb_mi <- comorbidities("^I21|^I22|^I252", relevant_diags) %>%
        rename(comorb_mi = comorb, datum_mi = datum)
    
    # Arrythmia (arr)
    # Codes from Fu et al. Stopping Renin-Angiotensin System Inhibitors in Patients with Advanced CKD and Risk of Adverse Outcomes: 
    # A Nationwide Study. JASN 2021.
    comorb_arr <- comorbidities("^I47|^I48|^I49", relevant_diags) %>%
        rename(comorb_arr = comorb, datum_arr = datum)
    
    # Cerebrovascular disease (cd)
    # Codes from Fu et al. Comparative Effectiveness of Renin-Angiotensin System Inhibitors and Calcium Channel Blockers in Individuals 
    # With Advanced CKD: A Nationwide Observational Cohort Study. AJKD 2021.
    comorb_cd <- comorbidities("^G459|^I63|^I64|^I69", relevant_diags) %>%
        rename(comorb_cd = comorb, datum_cd = datum)
    
    # Peripheral vascular disease (pvd)
    # Codes from Fu et al. Stopping Renin-Angiotensin System Inhibitors in Patients with Advanced CKD and Risk of Adverse Outcomes: 
    # A Nationwide Study. JASN 2021.
    comorb_pvd <- comorbidities("^I70|^I71|^I72|^I73", relevant_diags) %>%
        rename(comorb_pvd = comorb, datum_pvd = datum)
    
    # Ischemic heart disease (ihd)
    # Codes from Fu et al. Stopping Renin-Angiotensin System Inhibitors in Patients with Advanced CKD and Risk of Adverse Outcomes: 
    # A Nationwide Study. JASN 2021.
    comorb_ihd <- comorbidities("^I20|^I21|^I22|^I23|^I24|^I25", relevant_diags) %>%
        rename(comorb_ihd = comorb, datum_ihd = datum)
    
    # Chronic obstructive pulmonary disease (copd)
    # Fu et al. Association of Acute Increases in Plasma Creatinine after Renin-Angiotensin Blockade with Subsequent Outcomes. CJASN 2019.
    comorb_copd <- comorbidities("^J44", relevant_diags) %>%
        rename(comorb_copd = comorb, datum_copd = datum)
    
    # Cancer (ca)
    # Pasternak et al. Use of sodium-glucose co-transporter 2 inhibitors and risk of serious renal events: Scandinavian cohort study. 
    # BMJ 2020
    comorb_ca <- comorbidities("^C0|^C1|^C2|^C3|^C40|^C41|^C43|^C45|^C46|^C46|^C47|^C48|^C49|^C5|^C6|^C7|^C8|^C9", relevant_diags) %>%
        rename(comorb_ca = comorb, datum_ca = datum)

    # Liver disease (ld)
    # On own accord from Socialstyrelsen: https://klassifikationer.socialstyrelsen.se/ICD10SE.
    comorb_ld <- comorbidities("^K70|^K71|^K72|^K73|^K74|^K75|^K76|^K77", relevant_diags) %>%
        rename(comorb_ld = comorb, datum_ld = datum)
    
    # Chronic heart failure (chf)
    # Codes from Fu et al. Comparative Effectiveness of Renin-Angiotensin System Inhibitors and Calcium Channel Blockers in Individuals 
    # With Advanced CKD: A Nationwide Observational Cohort Study. AJKD 2021.
    comorb_chf <- comorbidities("^I110|^I130|^I132|^I50", relevant_diags) %>%
        rename(comorb_chf = comorb, datum_chf = datum)
    
    # Dyslipidemia (dl)
    # Own accord from Socialstyrelsen: https://klassifikationer.socialstyrelsen.se/ICD10SE.
    comorb_dl <- comorbidities("^E78", relevant_diags) %>%
        rename(comorb_dl = comorb, datum_dl = datum)
    
    # Hypothyroidism (hyth)
    # Own accord from Socialstyrelsen: https://klassifikationer.socialstyrelsen.se/ICD10SE.
    comorb_hyth <- comorbidities("^E02|^E03|^E890", relevant_diags) %>%
        rename(comorb_hyth = comorb, datum_hyth = datum)
    
    # Extracranial hemorrhage (eh)
    # Own accord from Socialstyrelsen: https://klassifikationer.socialstyrelsen.se/ICD10SE.
    comorb_eh <- comorbidities("^I60", relevant_diags) %>%
        rename(comorb_eh = comorb, datum_eh = datum)
    
    # Valvular heart disease (vhd)
    # Own accord from Socialstyrelsen: https://klassifikationer.socialstyrelsen.se/ICD10SE.
    comorb_vhd <- comorbidities("^I34|^I35|^I36|^I37|^I38|^I39", relevant_diags) %>%
        rename(comorb_vhd = comorb, datum_vhd = datum)
    
    ### Creat comorbidities set
    # First create a list of the comorbidity objects
    comorb_objects <- list(comorb_dm, comorb_ht, comorb_mi, comorb_arr, comorb_cd, comorb_pvd, comorb_ihd, comorb_copd,
                           comorb_ca, comorb_ld, comorb_chf, comorb_dl, comorb_hyth, comorb_eh, comorb_vhd)
    
    # Then start a loop joining the list together
    comorb <- lopnrs %>% dplyr::select(lopnr)
    
    for(i in comorb_objects){
        comorb <- comorb %>% left_join(i, "lopnr")
    }
    
    # Keep the first recorded comorbidity for everyone and change NA to 0
    comorbs <- comorb %>% arrange(lopnr, desc(comorb_dm), desc(comorb_ht), desc(comorb_mi), desc(comorb_arr),
                                  desc(comorb_cd), desc(comorb_pvd), desc(comorb_ihd), desc(comorb_copd), 
                                  desc(comorb_ca), desc(comorb_ld), desc(comorb_chf), desc(comorb_dl), desc(comorb_hyth),
                                  desc(comorb_eh), desc(comorb_vhd), datum_dm, datum_ht, datum_mi, datum_arr, datum_cd, 
                                  datum_pvd, datum_ihd, datum_copd, datum_ca, datum_ld, datum_chf, datum_dl, datum_hyth,
                                  datum_eh, datum_vhd) %>%
        group_by(lopnr) %>% slice(1L) %>% ungroup() %>%
        replace_na(list(comorb_dm = 0, comorb_ht = 0, comorb_mi = 0, comorb_arr = 0, comorb_cd = 0,
                        comorb_pvd = 0, comorb_ihd = 0, comorb_copd = 0, comorb_ca = 0, comorb_ld = 0,
                        comorb_chf = 0, comorb_dl = 0, comorb_hyth = 0, comorb_eh = 0, comorb_vhd = 0))
    
    rm(comorb_dm, comorb_ht, comorb_mi, comorb_arr, comorb_cd, comorb_pvd, comorb_ihd, comorb_copd, comorb_ca, comorb_ld,
       comorb_chf, comorb_dl, comorb_hyth, comorb_eh, comorb_vhd, comorb, i, comorb_objects)
    
    return(comorbs)
}

comorbs_c1 <- comorbs()

rm(relevant_diags)

##### 3. Medication use #####
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/meds_c1.Rdata")

# For medication, we use admission date as end of ascertainment period, as medication prescribed during hospital admission
# might be initiated shorter and discontinued shortly after.
medications <- function(atc_code){
    medic <- lopnrs %>% left_join(meds, "lopnr") %>% 
        filter(grepl(atc_code, atc), edatum >= (adm_dt - 180) & edatum < adm_dt & !is.na(edatum)) %>% 
        mutate(med = 1) %>% 
        arrange(lopnr, adm_dt, desc(med)) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% 
        ungroup() %>% 
        dplyr::select(lopnr, med, edatum) %>%
        rename(datum = edatum)
    return(medic)
}

medic <- function(){
    # beta-blockers (bb)
    # Codes from Fu et al. Stopping Renin-Angiotensin System Inhibitors in Patients with Advanced CKD and Risk of Adverse Outcomes: 
    # A Nationwide Study. JASN 2021.
    medic_bb <- medications("^C07") %>% rename(med_datum_bb = datum, med_bb = med)
    
    # Calcium channel blockers (ccb)
    # Codes from Fu et al. Comparative Effectiveness of Renin-Angiotensin System Inhibitors and Calcium Channel Blockers in Individuals 
    # With Advanced CKD: A Nationwide Observational Cohort Study. AJKD 2021.
    medic_ccb <- medications("^C08") %>% rename(med_datum_ccb = datum, med_ccb = med)
    
    # Loop diuretics (ld)
    # Codes from Fu et al. Comparative Effectiveness of Renin-Angiotensin System Inhibitors and Calcium Channel Blockers in Individuals 
    # With Advanced CKD: A Nationwide Observational Cohort Study. AJKD 2021.
    medic_ld <- medications("^C03C") %>% rename(med_datum_ld = datum, med_ld = med)
    
    # Thiazide diuretics (td)
    # Codes from Fu et al. Comparative Effectiveness of Renin-Angiotensin System Inhibitors and Calcium Channel Blockers in Individuals 
    # With Advanced CKD: A Nationwide Observational Cohort Study. AJKD 2021.
    medic_td <- medications("^C03A") %>% rename(med_datum_td = datum, med_td = med)
    
    # Potassium-sparing diuretics (pd)
    # Codes from Fu et al. Comparative Effectiveness of Renin-Angiotensin System Inhibitors and Calcium Channel Blockers in Individuals 
    # With Advanced CKD: A Nationwide Observational Cohort Study. AJKD 2021.
    medic_pd <- medications("^C03D") %>% rename(med_datum_pd = datum, med_pd = med)
    
    # Nonsteroidal anti-inflammatory drugs (nsaid)
    # Fu et al. High-sensitivity C-reactive protein and the risk of chronic kidney disease progression or acute kidney injury in post-
    # myocardial infarction patients. Am Heart J 2019.
    medic_nsaid <- medications("^M01A") %>% rename(med_datum_nsaid = datum, med_nsaid = med)
    
    # Antilipids (al)
    # Codes from Fu et al. Comparative Effectiveness of Renin-Angiotensin System Inhibitors and Calcium Channel Blockers in Individuals 
    # With Advanced CKD: A Nationwide Observational Cohort Study. AJKD 2021.
    # Minus AA at the end, because that's only statins
    medic_al <- medications("^C10") %>% rename(med_datum_al = datum, med_al = med)
    
    # alfa-blockers (ab)
    # On own accord from: https://www.whocc.no/atc_ddd_index/.
    medic_ab <- medications("^C02CA|^C02LE|^G04CA|^R03AA|^R03CA") %>% rename(med_datum_ab = datum, med_ab = med)
    
    # Nitrate (nt)
    # On own accord from: https://www.whocc.no/atc_ddd_index/.
    medic_nt <- medications("^C01DA") %>% rename(med_datum_nt = datum, med_nt = med)
    
    # Aldosteron receptor antagonists (ara)
    # Fu et al. High-sensitivity C-reactive protein and the risk of chronic kidney disease progression or acute kidney injury in post-
    # myocardial infarction patients. Am Heart J 2019.
    # Minus 01 at the end (which is only spironolacton)
    medic_ara <- medications("^C03DA") %>% rename(med_datum_ara = datum, med_ara = med)
    
    # Hydralazine (hl)
    # On own accord from: https://www.whocc.no/atc_ddd_index/.
    medic_hl <- medications("^C02DB01|^C02LG01|^C02LG51|^C02DB02|^C02LG02") %>% rename(med_datum_hl = datum, med_hl = med)
    
    # Antiarrythmics (aa)
    # On own accord from: https://www.whocc.no/atc_ddd_index/.
    medic_aa <- medications("^C01B") %>% rename(med_datum_aa = datum, med_aa = med)
    
    # Digoxin (dx)
    # On own accord from: https://www.whocc.no/atc_ddd_index/.
    medic_dx <- medications("^C01AA02|^C01AA52|^C01AA05|^C01AA08") %>% rename(med_datum_dx = datum, med_dx = med)
    
    # Anticoagulants (ac)
    # On own accord from: https://www.whocc.no/atc_ddd_index/.
    medic_ac <- medications("^B01") %>% rename(med_datum_ac = datum, med_ac = med)
    
    # Diabetes therapeutics (dm)
    # Codes from Fu et al. Comparative Effectiveness of Renin-Angiotensin System Inhibitors and Calcium Channel Blockers in Individuals 
    # With Advanced CKD: A Nationwide Observational Cohort Study. AJKD 2021.
    medic_dm <- medications("^A10") %>% rename(med_datum_dm = datum, med_dm = med)
    
    # Opioids (op)
    # On own accord from: https://www.whocc.no/atc_ddd_index/.
    medic_op <- medications("^N02A") %>% rename(med_datum_op = datum, med_op = med)
    
    ### Creat medication set
    # First create a list of the medication objects
    medic_objects <- list(medic_bb, medic_ccb, medic_ld, medic_td, medic_pd, medic_nsaid, medic_al, medic_ab, medic_nt, 
                          medic_ara, medic_hl, medic_aa, medic_dx, medic_ac, medic_dm, medic_op) 
    
    # Then start a loop joining the list together
    medic <- lopnrs %>% dplyr::select(lopnr)
    
    for(i in medic_objects){
        medic <- medic %>% left_join(i, "lopnr")
    }
    
    # Keep the first recorded drug for everyone and change NA to 0
    medication <- medic %>% arrange(lopnr, desc(med_bb), desc(med_ccb), desc(med_ld), desc(med_td), desc(med_pd), 
                                    desc(med_nsaid), desc(med_al), desc(med_ab), desc(med_nt), desc(med_ara), desc(med_hl), 
                                    desc(med_aa), desc(med_dx), desc(med_ac), desc(med_dm), desc(med_op), med_datum_bb, 
                                    med_datum_ccb, med_datum_ld, med_datum_td, med_datum_pd, med_datum_nsaid, med_datum_al, 
                                    med_datum_ab, med_datum_nt, med_datum_ac, med_datum_ara, med_datum_hl, med_datum_aa, 
                                    med_datum_dx, med_datum_ac, med_datum_dm, med_datum_op) %>%
        group_by(lopnr) %>% slice(1L) %>% ungroup() %>%
        replace_na(list(med_bb = 0, med_ccb = 0, med_ld = 0, med_td = 0, med_pd = 0, med_nsaid = 0, 
                        med_al = 0, med_ab = 0, med_nt = 0, med_ara = 0, med_hl = 0, med_aa = 0, 
                        med_dx = 0, med_ac = 0, med_dm = 0, med_op = 0))
    
    rm(medic_bb, medic_ccb, medic_ld, medic_td, medic_pd, medic_nsaid, medic_al, medic_ab, medic_nt, medic_ara, medic_hl, 
       medic_aa, medic_dx, medic_ac, medic_dm, medic_op, medic_objects, medic, i)
    
    return(medication)
}

medics_c1 <- medic()

##### 4. Healthcare access in previous year #####
load("kon_diagnoses_long.Rdata")
load("ovr_diagnoses_long.Rdata")
load("slv_diagnoses_long.Rdata")

healthcare_use_all <- function(df){
    hc_use <- lopnrs %>% left_join(df, "lopnr") %>%
        filter(as.Date(bdat) > (adm_dt - 365.25) & as.Date(bdat) < adm_dt) %>%
        arrange(lopnr, bdat) %>% group_by(lopnr, bdat) %>% slice(1L) %>% ungroup() %>% 
        arrange(lopnr) %>% group_by(lopnr) %>% 
        mutate(ind = 1,
               count = sum(ind)) %>%
        slice(1L) %>% ungroup() %>%
        dplyr::select(lopnr, count) 
    return(hc_use)
}

healthcare_use_cv <- function(df){
    hc_use_cv <- lopnrs %>% left_join(filter(df, grepl("^I", diagnosis)), "lopnr") %>%
        filter(as.Date(bdat) > (adm_dt - 365.25) & as.Date(bdat) < adm_dt) %>%
        arrange(lopnr, bdat) %>% group_by(lopnr, bdat) %>% slice(1L) %>% ungroup() %>% 
        arrange(lopnr) %>% group_by(lopnr) %>% 
        mutate(ind = 1,
               count = sum(ind)) %>%
        slice(1L) %>% ungroup() %>%
        dplyr::select(lopnr, count)
    return(hc_use_cv)
}

healthcare_use <- function(){
    kon_diagnoses_long_1diag <- filter(kon_diagnoses_long, diag_no == "diag1")
    slv_diagnoses_long_1diag <- filter(slv_diagnoses_long, diag_no == "diag1")
    ovr_diagnoses_long_1diag <- filter(ovr_diagnoses_long, diag_no == "diag1")
    
    # Retrieve counts for use
    kontakt_all <- healthcare_use_all(kon_diagnoses_long_1diag) %>% rename(kontakt_all = count)
    kontakt_cv <- healthcare_use_cv(kon_diagnoses_long_1diag) %>% rename(kontakt_cv = count)
    slv_all <- healthcare_use_all(slv_diagnoses_long_1diag) %>% rename(slv_all = count)
    slv_cv <- healthcare_use_cv(slv_diagnoses_long_1diag) %>% rename(slv_cv = count)
    ovr_all <- healthcare_use_all(ovr_diagnoses_long_1diag) %>% rename(ovr_all = count)
    ovr_cv <- healthcare_use_cv(ovr_diagnoses_long_1diag) %>% rename(ovr_cv = count)
    
    # Join all together
    hcu <- dplyr::select(lopnrs, -adm_egfr) %>% left_join(kontakt_all, "lopnr") %>% left_join(kontakt_cv, "lopnr") %>% 
        left_join(slv_all, "lopnr") %>% left_join(slv_cv, "lopnr") %>%
        left_join(ovr_all, "lopnr") %>% left_join(ovr_cv, "lopnr") %>%
        replace_na(list(kontakt_all = 0, kontakt_cv = 0, slv_all = 0, slv_cv = 0, ovr_all = 0,
                        ovr_cv = 0)) %>%
        dplyr::select(-adm_dt, -dis_dt)
    
    rm(kontakt_all, kontakt_cv, slv_all, slv_cv, ovr_all, ovr_cv)
    
    ## Also as categories
    # 0 = 0 access
    # 1 = 1 access
    # 2 = 2 access
    # 3 = 3 or more access
    hcu_cat <- function(y){
        x = ifelse(y == 0, 0,
            ifelse(y == 1, 1,
            ifelse(y == 2, 2,
            ifelse(y >= 3, 3, NA))))
    }
    
    
    hcu <- hcu %>% mutate(kontakt_all_cat = hcu_cat(kontakt_all),
                          kontakt_cv_cat = hcu_cat(kontakt_cv),
                          slv_all_cat = hcu_cat(slv_all),
                          slv_cv_cat = hcu_cat(slv_cv),
                          ovr_all_cat = hcu_cat(ovr_all),
                          ovr_cv_cat = hcu_cat(ovr_cv))
    
    return(hcu)
}

hc_use_c1 <- healthcare_use()

##### 5. Reason for hospitalisation #####
# To determine cause of hospitalisation only take first diagnosis code, otherwise the total will be > 100% (especially
# with genitourinary hospitalisations)
first_diags <- slv_diagnoses_long %>% filter(diag_no == "diag1") %>%
    dplyr::select(lopnr, diagnosis, bdat) %>% 
    mutate(bdat = as.Date(bdat, origin = "1970-01-01")) %>% 
    rename(adm_dt = bdat)

hospital <- lopnrs %>% left_join(first_diags, c("lopnr", "adm_dt")) %>% 
    arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup()

hosps <- function(cause){
    hosp <- hospital %>% filter(grepl(cause, diagnosis)) %>%
        mutate(hosp = 1) %>%
        arrange(lopnr, adm_dt, desc(hosp)) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% 
        ungroup() %>%
        dplyr::select(lopnr, hosp, diagnosis)
    return(hosp)
}

hospis <- function(){
    # Cardiovascular (cv)
    hosps_cv <- hosps("^I") %>% rename(hosp_cv = hosp) %>% dplyr::select(-diagnosis)
    
    # Respiratory (rs)
    hosps_rs <- hosps("^J") %>% rename(hosp_rs = hosp) %>% dplyr::select(-diagnosis)
    
    # Gastrointestinal, endocrine and metabolic (gi)
    # Remove E875 as these are selected in the HK hospitalisation cause
    hosps_gi <- hosps("^E|^K") %>% rename(hosp_gi = hosp) %>%
        filter(diagnosis != "E875") %>% dplyr::select(-diagnosis)
    
    # Infectious disease (id)
    hosps_id <- hosps("^A|^B") %>% rename(hosp_id = hosp) %>% dplyr::select(-diagnosis)
    
    # Cancer (ca)
    hosps_ca <- hosps("^C|^D1|^D2|^D3|^D40|^D41|^D42|^D43|^D44|^D45|^D46|^D47|^D48|^D49") %>% 
        rename(hosp_ca = hosp) %>% dplyr::select(-diagnosis)
    
    # Orthopedics (op)
    hosps_op <- hosps("^M") %>% rename(hosp_op = hosp) %>% dplyr::select(-diagnosis)
    
    # Hematologic (hm)
    hosps_hm <- hosps("^D5|^D6|^D7|^D80|^D81|^D82|^D83|^D84|^D85|^D86|^D87|^D88|^D89") %>% 
        rename(hosp_hm = hosp) %>% dplyr::select(-diagnosis)
    
    # Genitourinary (gu)
    hosps_gu <- hosps("^N") %>% rename(hosp_gu = hosp) %>% dplyr::select(-diagnosis)
    
    # Injury or poisining (ip)
    hosps_ip <- hosps("^S|^T") %>% rename(hosp_ip = hosp) %>% dplyr::select(-diagnosis)
    
    # Nervous system (ns)
    hosps_ns <- hosps("^G") %>% rename(hosp_ns = hosp) %>% dplyr::select(-diagnosis)
    
    # Hyperkalemia (hk)
    hosps_hk <- hosps("^E875") %>% rename(hosp_hk = hosp) %>% dplyr::select(-diagnosis)
    
    # Other disease (od)
    hosps_od <- hosps("^R|^F|^H|^L|^O|^P|^Q|^U|^V|^W|^X|^Y|^Z") %>% rename(hosp_od = hosp) %>% dplyr::select(-diagnosis)
    
    ### Creat hospitalisations set
    # First create a list of the hospitalisation objects
    hosps_objects <- list(hosps_cv, hosps_rs, hosps_gi, hosps_id, hosps_ca, hosps_op, hosps_hm, hosps_gu, hosps_ip, 
                          hosps_ns, hosps_hk, hosps_od)
    
    # Then start a loop joining the list together
    hospi <- lopnrs %>% dplyr::select(lopnr)
    
    for(i in hosps_objects){
        hospi <- hospi %>% left_join(i, "lopnr")
    }
    
    # Change NA to 0
    hospital_causes <- hospi %>% arrange(lopnr) %>%
        replace_na(list(hosp_cv = 0, hosp_rs = 0, hosp_gi = 0, hosp_id = 0, hosp_ca = 0, 
                        hosp_op = 0, hosp_hm = 0, hosp_gu = 0, hosp_ip = 0, hosp_ns = 0, hosp_hk = 0, 
                        hosp_od = 0))
    
    rm(hosps_cv, hosps_rs, hosps_gi, hosps_id, hosps_ca, hosps_op, hosps_hm, hosps_gu, hosps_ip, hosps_ns, hosps_hk, 
       hosps_od, hosps_objects, hospi, i)
    
    return(hospital_causes)
}

hosps_c1 <- hospis()

##### 6. Procedures and conditions during hospitalisation #####
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/cohort1.Rdata")
load("slv_procedurecodes.Rdata")

proc <- function(code){
    proc <- lopnrs %>% left_join(slv_procedurecodes, "lopnr") %>%
        filter(grepl(code, opk), as.Date(datum) >= adm_dt & as.Date(datum) <= dis_dt & !is.na(datum)) %>%
        mutate(pc = 1) %>%
        arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>%
        ungroup() %>%
        dplyr::select(lopnr, pc)
    return(proc)
}

cond <- function(code){
    cond <- lopnrs %>% left_join(slv_diagnoses_long, "lopnr") %>%
        filter(grepl(code, diagnosis), as.Date(bdat) >= adm_dt & as.Date(bdat) <= dis_dt & !is.na(bdat)) %>%
        mutate(pc = 1) %>%
        arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup() %>%
        ungroup() %>%
        dplyr::select(lopnr, pc)
    return(cond)
}

procedures_conditions <- function(){
    # Sepsis (sp)
    poc_sp <- cond("^A40|^A41") %>% rename(pc_sp = pc)
    
    # Cardiac surgery (cs)
    poc_cs <- proc("^FB|^FCFB|^FD|^FE|^FF|^FG|^FH|^FJ|^FK|^FL|^FM|^FN|^FP|^FQ|^FW|^FX") %>% 
        rename(pc_cs = pc)
    
    # Cardiac catherization (cc)
    poc_cc <- proc("^TFC00|^TFC10") %>% rename(pc_cc = pc)
    
    # Abdominal aortic aneurysm (ar)
    poc_ar1 <- proc("^PDG") %>% rename(pc_ar = pc)
    
    poc_ar2 <- cond("^I713|^I714|^I715|^I716") %>% rename(pc_ar = pc)
    
    poc_ar <- rbind(poc_ar1, poc_ar2)
    
    # Pneumonia (pn)
    poc_pn <- cond("^J12|^J13|^J14|^J15|^J16|^J17|^J18") %>% rename(pc_pn = pc)
    
    # Liver failure (lf) (^K9182 doenst exist in Swedish data, only K918 which is not specific enough)
    poc_lf <- cond("^K704|^K72") %>% rename(pc_lf = pc)
    
    # Acute myocardial infarction (mi)
    poc_mi <- cond("I21") %>% rename(pc_mi = pc)
    
    # Non-cardiac surgery (ns)
    poc_ns <- proc("A|^B|^C|^D|^E|^G|^H|^J|^K|^L|^M|^N|^P|^Q|^T|^U|^X|^Y|^Z") %>% rename(pc_ns = pc)
    
    # Acute dialysis (ad)
    # poc_ad <- proc("^DR015|^DR023") %>% rename(pc_ad = pc)
    
    ### Create procedures and conditiosn set
    # First create a list of the objects
    pc_objects <- list(poc_sp, poc_cs, poc_cc, poc_ar, poc_pn, poc_lf, poc_mi, poc_ns)
    
    # Then start a loop joining the list together
    proc_cond <- lopnrs %>% dplyr::select(lopnr) 
    
    for(i in pc_objects){
        proc_cond <- proc_cond %>% left_join(i, "lopnr")
    }
    
    # Change NA to 0
    procedures_conditions <- proc_cond %>% arrange(lopnr, desc(pc_sp), desc(pc_cs), desc(pc_cc), desc(pc_ar), desc(pc_pn),
                                                   desc(pc_lf), desc(pc_mi), desc(pc_ns)) %>% group_by(lopnr) %>%
        slice(1L) %>% ungroup() %>%
        replace_na(list(pc_sp = 0, pc_cs = 0, pc_cc = 0, pc_ar = 0, pc_pn = 0, 
                        pc_lf = 0, pc_mi = 0, pc_ns = 0))
    
    rm(poc_sp, poc_cs, poc_cc, poc_ar1, poc_ar2, poc_ar, poc_pn, poc_lf, poc_mi, poc_ns, pc_objects, i)
    
    return(procedures_conditions)
}

proc_cond_c1 <- procedures_conditions()

##### 7. During index hospitalization (dih) #####
load("lab_values.Rdata")

dih <- function(cohort){
    # Duration of admission
    dur <- lopnrs %>% mutate(dur = as.numeric(dis_dt - adm_dt))
    
    # AKI as cause of admission (AKI_acoa)
    # AKI_acoa 1 = AKI was cause of admission
    # AKI_acoa 0 = AKI was not cause of admission (thus complication)
    AKIadmission <- dur %>% left_join(first_diags, c("lopnr", "adm_dt")) %>%
        mutate(AKI_acoa = ifelse(grepl("^N17", diagnosis), 1, 0))
    
    # Type of AKI
    KDIGO_stages <- AKIadmission %>% left_join(cohort, c("lopnr", "adm_dt", "dis_dt", "adm_egfr")) %>%
        dplyr::select(lopnr, adm_dt, dis_dt, dur, AKI_acoa, crit_type)
    
    # Need for acute dialysis
    AcuteDialysis <- KDIGO_stages %>% left_join(slv_procedurecodes, "lopnr") %>%
        mutate(AcutDial = ifelse(as.Date(datum) >= adm_dt & as.Date(datum) <= dis_dt &
                                     grepl("^DR015|^DR023", opk), 1, 0)) %>%
        arrange(lopnr, adm_dt, desc(AcutDial), datum) %>% group_by(lopnr, adm_dt) %>%
        slice(1L) %>% ungroup() %>%
        dplyr::select(-datum, -opk)
    
    ## Blood values of creatinine and potassium
    # Serum creatinine (peak and mean)
    creatinines_hosp <- lab_values %>% filter((test == "crea" | test == "vencrea" | test == "artcrea") & ip == 1)
    
    creats <- lopnrs %>% left_join(creatinines_hosp, "lopnr") %>%
        filter(as.Date(datum) >= adm_dt & as.Date(datum) <= dis_dt) %>%
        arrange(lopnr) %>% group_by(lopnr) %>%
        mutate(creat_mean = mean(result),
               creat_peak = max(result)) %>%
        slice(1L) %>% ungroup() %>%
        dplyr::select(-dis_dt, -(test:tid), -(Analys:vtype), -adm_egfr)
    
    # Serum potassium (peak and mean)
    potassiums_hosp <- lab_values %>% filter((test == "ab_potas" | test == "potas" | test == "other_potas") & ip == 1)
    
    potass <- lopnrs %>% left_join(potassiums_hosp, "lopnr") %>%
        filter(as.Date(datum) >= adm_dt & as.Date(datum) <= dis_dt) %>%
        arrange(lopnr) %>% group_by(lopnr) %>%
        mutate(potass_mean = mean(result),
               potass_peak = max(result)) %>%
        slice(1L) %>% ungroup() %>%
        dplyr::select(-dis_dt, -(test:tid), -(Analys:vtype), -adm_egfr)
    
    # Add together
    dih_complete <- AcuteDialysis %>% left_join(creats, c("lopnr", "adm_dt")) %>% 
        left_join(potass, c("lopnr", "adm_dt")) %>%
        dplyr::select(-adm_dt, -dis_dt, -datum_c)
    
    rm(dur, AKIadmission, AcuteDialysis, creats, potass, KDIGO_stages)
    
    return(dih_complete)
}

dih_c1 <- dih(cohort1)

##### 8. Joining everything together #####
covar <- function(gi, comorbs, medics, hc_use, hosps, proc_cond, dih){
    
    # Create a list of objects
    covar_objects <- list(gi, comorbs, medics, hc_use, hosps, proc_cond, dih)
    
    # Start a loop joining those together
    covars <- lopnrs %>% dplyr::select(lopnr, adm_egfr)
    
    for(i in covar_objects){
        covars <- covars %>% left_join(i, "lopnr")
    }
    
    covariates <- covars %>% dplyr::select(-datum_dm, -datum_ht, -datum_mi, -datum_arr, -datum_cd, -datum_pvd, -datum_ihd, 
                                    -datum_copd, -datum_ca, -datum_ld, -datum_chf, -datum_dl, -datum_hyth, -datum_eh, 
                                    -datum_vhd, -med_datum_bb, -med_datum_ccb, -med_datum_ld, -med_datum_td, 
                                    -med_datum_pd, -med_datum_nsaid, -med_datum_al, -med_datum_ab, -med_datum_nt, 
                                    -med_datum_ac, -med_datum_ara, -med_datum_hl, -med_datum_aa, -med_datum_dx,
                                    -med_datum_ac, -med_datum_dm, -med_datum_op)
    
    rm(covar_objects, covars)
    
    return(covariates)
}

covar_c1 <- covar(gi_c1, comorbs_c1, medics_c1, hc_use_c1, hosps_c1, proc_cond_c1, dih_c1) 

rm(gi_c1, comorbs_c1, medics_c1, hc_use_c1, hosps_c1, proc_cond_c1, dih_c1)

##### 9. Adding laboratory values #####
lab_value <- function(labmeasure){
    lab_val <- lopnrs %>% left_join(filter(lab_values, test == labmeasure & ip == 0), "lopnr") %>%
        filter(as.Date(datum) > (adm_dt - 180) & as.Date(datum) < (adm_dt - 7)) %>%
        arrange(lopnr) %>% group_by(lopnr) %>%
        mutate(labval = mean(result)) %>%
        slice(1L) %>% ungroup() %>%
        dplyr::select(lopnr, labval)
    return(lab_val)
}

lab <- function(covar){
    ## ACR & dipstick conversion
    # Laboratory values equal the mean of values in the six months prior to the index hospitalization
    acr_nor <- lab_value("uacr") %>% rename(acr = labval)
    
    # Dipstick for conversion
    acr_dp <- lab_value("dipstick-alb") %>% rename(dp = labval)
    
    # Adding acr's together and creating categories
    # acr categories are:
    # 1 = A1: < 3
    # 2 = A2: 3 - 29
    # 3 = A3: 30 - 69
    # 4 = A4: ≥ 70
    # 5 = missing
    acr <- lopnrs %>% left_join(acr_nor, "lopnr") %>% left_join(acr_dp, "lopnr") %>%
        dplyr::select(lopnr, acr, dp) %>%
        mutate(acr_cat = ifelse(acr < 3, 1,
                         ifelse(acr >= 3 & acr < 30, 2,
                         ifelse(acr >= 30 & acr < 70, 3,
                         ifelse(acr >= 70, 4, NA)))),
               acr_cat = ifelse(is.na(acr_cat) & dp < 1, 1,
                         ifelse(is.na(acr_cat) & dp >= 1 & dp < 2, 2,
                         ifelse(is.na(acr_cat) & dp >= 2 & dp < 4, 3,
                         ifelse(is.na(acr_cat) & dp >= 4, 4, acr_cat)))),
               acr_cat = ifelse(is.na(acr_cat), 5, acr_cat)) %>% 
        dplyr::select(-dp) 
    
    # Hemoglobin
    # hb <- lab_value("hb") %>% rename(hb = labval)
    
    # Cholesterol
    chol <- lab_value("tc") %>% rename(chol = labval)
    
    # Joining lab together
    lab <- dplyr::select(lopnrs, -adm_egfr) %>% left_join(acr, "lopnr") %>% left_join(chol, "lopnr")
    
    # Adding to other covariates
    covariates <- covar %>% left_join(lab, c("lopnr", "adm_dt"))
    
    rm(chol, lab, acr, acr_dp, acr_nor)
    
    return(covariates)
}

covariates_c1 <- lab(covar_c1) %>% relocate(adm_egfr, .after = age_cat)

save(covariates_c1, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/covariates_c1.Rdata")

rm(covar_c1, covariates_c1, cohort1, lopnrs)

timestamp()

#-------------------------------------------------------------------------------------------------------------------------#

######### Cohort 2 ######### 
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/cohort2.Rdata")

lopnrs <- cohort2 %>% dplyr::select(lopnr, adm_dt, dis_dt, adm_egfr)

gi_c2 <- general_information()

##### 2. Comorbidities #####
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/diagnoses.Rdata")

relevant_diags <- dplyr::select(lopnrs, lopnr) %>% left_join(diagnoses, "lopnr")

rm(diagnoses)

comorbs_c2 <- comorbs()

##### 3. Medication use #####
load("~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/meds_c2.Rdata")

medics_c2 <- medic()

##### 4. Healthcare access in previous year #####
hc_use_c2 <- healthcare_use()

##### 5. Reason for hospitalisation #####
first_diags <- slv_diagnoses_long %>% filter(diag_no == "diag1") %>%
    dplyr::select(lopnr, diagnosis, bdat) %>% 
    mutate(bdat = as.Date(bdat, origin = "1970-01-01")) %>% 
    rename(adm_dt = bdat)

# To determine cause of hospitalisation only take first diagnosis code, otherwise the total will be > 100% (especially
# with genitourinary hospitalisations)
hospital <- lopnrs %>% left_join(first_diags, c("lopnr", "adm_dt")) %>% 
    arrange(lopnr, adm_dt) %>% group_by(lopnr, adm_dt) %>% slice(1L) %>% ungroup()

hosps_c2 <- hospis()

##### 6. Procedures and conditions during hospitalisation #####
proc_cond_c2 <- procedures_conditions()

##### 7. During index hospitalization (dih) #####
dih_c2 <- dih(cohort2)

##### 8. Joining everything together #####
covar_c2 <- covar(gi_c2, comorbs_c2, medics_c2, hc_use_c2, hosps_c2, proc_cond_c2, dih_c2) 

rm(gi_c2, comorbs_c2, medics_c2, hc_use_c2, hosps_c2, proc_cond_c2, dih_c2)

##### 9. Adding laboratory values #####
# We add laboratory values after joining everything because for dipstick conversion we need covariate information
covariates_c2 <- lab(covar_c2)

save(covariates_c2, file = "~/Research/[] Nephrology/3. AKI_ACEi-ARB/SCREAM2/Codes/Dataframes/covariates_c2.Rdata")

timestamp()

#-------------------------------------------------------------------------------------------------------------------------#



























    
    
    
    
    
    
    