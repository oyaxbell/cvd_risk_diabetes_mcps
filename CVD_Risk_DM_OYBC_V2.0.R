# Assessment of fatal cardiovascular disease risk using data-driven diabetes subgroups and SCORE2-Diabetes in 24,943 adults in Mexico City
# Authors: Perezalonso-Espinosa, J. & Bello-Chavolla, OY.
# Last edit: 25/Feb/2026
# R version: 4.5.2 (2025-10-31)
# RStudio version: 2026.1.0.392

#### Libraries ####
require(pacman)
pacman::p_load(adjustedCurves, blandr, BlandAltmanLeh, caret, cluster, corrplot, CVrisk, datasets,
               devtools, dplyr, dummy, dcurves, Epi, factoextra, FactoMineR, finalfit,
               flextable, flexsurv, fmsb, forestmodel, fpc, ggedit,forestplot, ggimage, ggplotify,
               ggpubr, ggsankey, ggsci, globorisk, gridExtra, gtools, gtsummary, haven, keras, whoishRisk,
               lmtest, labelled, magrittr, mice, miceadds, mitools, naniar, NbClust, nephro, nhanesA,
               nortest, officer, pacman, pec, plotRCS, pROC, reticulate, RiskScorescvd,riskRegression, rms, splines,
               survival, survminer, tidytidbits, tidyr, tidyverse, timeROC, UpSetR,
               VIM, viridis, wesanderson, ggbreak,rmda,compareC,compareCstat)
#devtools::install_github("Lemonade0924/compareCstat")

header = function(df1,df2){
 n = nrow(df1)
 p = round(n/nrow(df2) * 100,1)
 h = paste0("n = ",format(n, big.mark = ",")," (",p,"%)")
 return(h)
} 
make_summary <- function(data, header_text = NULL) {
 data %>%
  set_variable_labels(.labels = var_labels) %>%
  dplyr::select(all_of(vars)) %>%
  tbl_summary(
   missing = "ifany",
   missing_text = "Missing",
   missing_stat = "{N_miss} ({p_miss}%)",
   statistic = list(all_continuous() ~ "{mean} ({sd})",
                    all_categorical() ~ "{n} ({p}%)")
  ) %>%
  modify_header(
   label = "**Variable**",
   stat_0 = if (is.null(header_text)) {
    paste0("n = ", format(nrow(data), big.mark = ","))
   } else {
    header_text})}
make_fatal_summary <- function(data, header_text = NULL) {
 data %>%
  set_variable_labels(.labels = fatal_labels) %>%
  dplyr::select(all_of(fatal_vars)) %>%
  tbl_summary(
   type = all_categorical() ~ "dichotomous",
   missing = "no",
   statistic = list(all_categorical() ~ "{n} ({p}%)")
  ) %>%
  modify_header(
   label = "**Variable**",
   stat_0 = if (is.null(header_text)) {
    paste0("n = ", format(nrow(data), big.mark = ","))
   } else {
    header_text}) %>%
  modify_table_body(
   ~ .x %>% mutate(stat_0 = if_else(label %in% c(""), "", stat_0)))}
apply_name <- function(x, model, cluster) {
 x %>% 
  mutate(model = {{model}},
         cluster = {{cluster}})
}
#### Dataset ####
#setwd("/Users/Jeronimo/Library/CloudStorage/GoogleDrive-jeroperezalonso@gmail.com/.shortcut-targets-by-id/1PnOtxLTUrx7sIMjM5Sp2ICl0Znor0vlD/Datasets/MCPS")
setwd("~/Google Drive/Mi unidad/Datasets/MCPS")

#Base MCPS
mcps_basal= read.csv("2021-004 MCPS BASELINE.csv")
#Mortality 
mcps_mort = read.csv("2022-012 MCPS MORTALITY_2022.csv")
#Laboratory
mcps_nmr = read.csv("2022-012_OmarBelloChavolla_NMR_Second_Release-Corrected/2022-012 MCPS BASE_NMR_DATA_RECALIB-corrected.csv")
#IDS
mcps_ids = read.csv("2022-012 MCPS BASELINE_IDS.csv")
#EDU_GP
mcps_edu = read.csv("2022-012 MCPS BASELINE_EDUGP.csv")
#Tags
#mcps_tags = read.csv("2022-012_OmarBelloChavolla_NMR_Second_Release-Corrected/2022-012 MCPS BASE_NMR_TAGS-corrected.csv")
#Quality control
#mcps_qc = read.csv("2022-012_OmarBelloChavolla_NMR_Second_Release-Corrected/2022-012 MCPS BASE_NMR_QC-corrrected.csv")
#Resurvey
mcps_resurvey = read.csv("2022-012 MCPS RESURVEY.csv")
#Resurvey NMR
#mcps_resurvey_nmr = read.csv("2022-012_OmarBelloChavolla_NMR_Second_Release/2022-012 MCPS RESURVEY_NMR_DATA_RECALIB.csv")

#setwd("/Users/Jeronimo/Library/CloudStorage/GoogleDrive-jeroperezalonso@gmail.com/.shortcut-targets-by-id/1RG-1Sg0NBmkz3AlzY5XL6XJ0mezx6oNg/Mexico City Cohort Study/Proyectos/Cardiovascular Risk DM")
setwd("~/Google Drive/Mi unidad/Mexico City Cohort Study/Proyectos/Cardiovascular Risk DM")

#Join datasets 
mcps = mcps_basal %>% 
 full_join(mcps_mort, by= "PATID") %>%  
 full_join(mcps_nmr, by = "PATID") %>%
 #full_join(mcps_tags, by = "PATID") %>%
 #full_join(mcps_qc, by = "PATID") %>%
 full_join(mcps_resurvey, by = "PATID") %>%
 full_join(mcps_ids, by = "PATID") %>%
 full_join(mcps_edu, by = "PATID")

rm(mcps_basal, mcps_mort, mcps_nmr, mcps_resurvey,mcps_ids,mcps_edu)

#### Filter and definitions ####
# Recodify for all of MCPS # 
mcps <- mcps %>%
 mutate(EDU_LEVEL = factor(EDU_LEVEL,
                           levels = c(1, 2, 3, 4),
                           labels = c("University/College", "High School", "Elementary", "Other"))) %>%
 mutate(BMI = (WEIGHT)/((HEIGHT/100)^2)) %>%
 mutate(SMOKEGP = ifelse(is.na(SMOKEGP),0,SMOKEGP)) %>%
 mutate(SMOKER = ifelse(SMOKEGP %in% c(3, 4, 5), 1, 0)) %>%
 mutate(Creatinine = Creatinine/88.42)  %>%
 mutate(eGFR = CKDEpi2021.creat(creatinine = Creatinine,
                                sex = MALE, age = AGE)) %>%
 mutate(Albumin = Albumin / 10) %>% #Convert from g/L to g/dl
 mutate(Glucose = Glucose * 18) %>%
 mutate(Total_C = Total_C * 38.67) %>%
 mutate(non_HDL_C = non_HDL_C * 38.67) %>%
 mutate(HDL_C = HDL_C * 38.67) %>%
 mutate(LDL_C = HDL_C * 38.67) %>%
 mutate(VLDL_C = HDL_C * 38.67) %>%
 mutate(ApoB = ApoB * 100) %>%
 mutate(Total_TG = Total_TG * 88.57) %>%
 mutate(HT_TX = ifelse(((DRUG_A1==1)|(DRUG_A2==1)|(DRUG_A3==1)|(DRUG_A4==1)|
                         (DRUG_A5==1)|(DRUG_A6==1)|(DRUG_A7==1)|(DRUG_A8==1)|
                         (DRUG_A9==1)|(DRUG_A10==1)|(DRUG_A11==1)),1,0)) %>%
 mutate(otherCV_TX = ifelse(DRUG_B1 == 1|DRUG_B2 == 1|DRUG_B3 == 1|DRUG_B4 == 1,1,0)) %>%
 mutate(thromb_TX = ifelse(DRUG_C1==1|DRUG_C2==1|DRUG_C3==1|DRUG_C4==1,1,0)) %>%
 mutate(lipid_TX = ifelse(DRUG_E1==1|DRUG_E2==1|DRUG_E3==1|DRUG_E4==1|DRUG_E5==1,1,0)) %>%
 mutate(diab_TX = ifelse(DRUG_D1==1|DRUG_D2==1|DRUG_D3==1|DRUG_D4==1,1,0)) %>%
 mutate(otherlipid_TX = ifelse(DRUG_E2==1|DRUG_E3==1|DRUG_E4==1|DRUG_E5==1,1,0)) %>%
 mutate(meanSBP = case_when(
  !is.na(SBP1) & !is.na(SBP2) & !is.na(SBP3) ~ (SBP1 + SBP2 + SBP3) / 3,
  !is.na(SBP1) & !is.na(SBP2) ~ (SBP1 + SBP2) / 2,
  !is.na(SBP1) & !is.na(SBP3) ~ (SBP1 + SBP3) / 2,
  !is.na(SBP2) & !is.na(SBP3) ~ (SBP2 + SBP3) / 2,
  !is.na(SBP1) ~ SBP1,
  !is.na(SBP2) ~ SBP2,
  !is.na(SBP3) ~ SBP3,
  TRUE ~ NA_real_
 )) %>%
 mutate(meanDBP = case_when(
  !is.na(DBP1) & !is.na(DBP2) & !is.na(DBP3) ~ (DBP1 + DBP2 + DBP3) / 3,
  !is.na(DBP1) & !is.na(DBP2) ~ (DBP1 + DBP2) / 2,
  !is.na(DBP1) & !is.na(DBP3) ~ (DBP1 + DBP3) / 2,
  !is.na(DBP2) & !is.na(DBP3) ~ (DBP2 + DBP3) / 2,
  !is.na(DBP1) ~ DBP1,
  !is.na(DBP2) ~ DBP2,
  !is.na(DBP3) ~ DBP3,
  TRUE ~ NA_real_
 )) %>%
 mutate(WHtR = WAISTC/HEIGHT) %>%
 mutate(WHR = WAISTC/HIPC) %>%
 mutate(DM_DX = BASE_DIABETES) %>% # 1: 21792, 0: 137725
 mutate(DM_noDX = ifelse(
  DM_DX == 0 & #1:8170, 0:129555
   (ifelse(is.na(DRUG_D3), FALSE, DRUG_D1 == 1) |
     ifelse(is.na(DRUG_D3), FALSE, DRUG_D2 == 1) |
     ifelse(is.na(DRUG_D3), FALSE, DRUG_D3 == 1) |
     ifelse(is.na(DRUG_D3), FALSE, DRUG_D4 == 1) |
     ifelse(is.na(BASE_HBA1C), FALSE, BASE_HBA1C >= 6.5)),1,0)) %>%
 mutate(DM = ifelse(DM_DX == 1 | DM_noDX == 1,1,0)) %>%
 ## Define cardiovascular mortality 
 # Narrow fatal CVD
 mutate(fatal_narrow_cardiac = ifelse((D001 == 1),1,0)) %>%
 mutate(fatal_narrow_cerebrovascular = ifelse((D004 == 1 | D005 == 1 | D006 ==1),1,0)) %>%
 mutate(fatal_narrow_total = ifelse((fatal_narrow_cardiac == 1 | fatal_narrow_cerebrovascular == 1),1,0)) %>%
 # Wide fatal CVD
 mutate(fatal_wide_cardiac_ischemic = ifelse((D001 ==1),1,0)) %>%
 mutate(fatal_wide_cardiac_other = ifelse(substr(ICD10_UNDERLYING, 1, 3) %in% c("I10","I11","I12","I13","I14","I15",
                                                                                "I46","I47","I48","I49","I50","I51","I52"),1,0)) %>%
 mutate(fatal_wide_cerebrovascular_stroke = ifelse((D004 == 1 | D005 == 1 | D006 ==1),1,0)) %>%
 mutate(fatal_wide_other = ifelse(substr(ICD10_UNDERLYING, 1, 3) %in% c("R96","E105","E115","E145",
                                                                        "I170","I171","I172","I173","I174",
                                                                        "I175","I176","I177","I178","I179"),1,0)) %>%
 mutate(fatal_wide_total = ifelse(fatal_wide_cardiac_ischemic == 1 | fatal_wide_cardiac_other == 1 |
                                   fatal_wide_cerebrovascular_stroke == 1 | fatal_wide_other == 1,1,0)) %>%
 mutate(mets_ir=((log((2*Glucose)+Total_TG)*BMI))/(log(HDL_C)))

### MCPS diabetes ###
mcps_dm = mcps %>% filter(DM == 1)

### MCPS nodiabetes ###
mcps_nodm = mcps %>% filter(DM == 0) %>%
 filter(BASE_HEARTATTACK == 0) %>% filter(BASE_ANGINA == 0) %>% filter(BASE_STROKE == 0) %>% #filter(BASE_PAD == 0) #Remove baseline CVD 
 filter(STATUS != "U") %>%
 filter(!((meanSBP > 270) | (meanSBP < 70))) %>%
 filter(!(BMI < 10 | BMI > 80)) %>%
 filter(!(Total_C > (20*38.67) | Total_C < (1.75*38.67) | is.na( Total_C))) %>%
 filter(!is.na(BASE_HBA1C),!is.na(WHtR),!is.na(BMI),!is.na(meanSBP),!is.na(meanDBP),!is.na(eGFR),!is.na(Total_C)) %>%
 mutate(comorbidity_count = rowSums(select(., BASE_EMPHYSEMA, BASE_ASTHMA, BASE_CKD, BASE_PEP, BASE_CIRR, 
                                           BASE_HYPERTENSION, BASE_LUNGCANCER, BASE_OTHCANCER, BASE_PROSTATECANCER, 
                                           BASE_CERVCANCER, BASE_BREASTCANCER, BASE_STOMCANCER, BASE_ORALCANCER, BASE_PAD), 
                                    na.rm = TRUE)) 

### MCPS final dataset ###
mcps_dm_fin = mcps_dm %>% 
 filter(BASE_HEARTATTACK == 0) %>% filter(BASE_ANGINA == 0) %>% filter(BASE_STROKE == 0) %>% #filter(BASE_PAD == 0) #Remove baseline CVD 
 filter(STATUS != "U") %>% # Remove those with uncertain status 
 # Recodify resurvey variables and filter non-fatal cases
 #mutate(R_HXMI= ifelse(is.na(R_HXMI),0,1)) %>% #R_HXMI: Myocardial Infarction diagnosed (1=Yes)
 #filter(!(R_HXMI == 1 & (R_HXMI_YR > YEAR_RECRUITED) & (STATUS == "A"))) %>% # Remove non-fatal MI
 #mutate(R_HXSTR= ifelse(is.na(R_HXSTR),0,1)) %>% #R_HXSTR: Stroke diagnosed (1=Yes)
 #filter(!(R_HXSTR == 1 & (R_HXSTR_YR > YEAR_RECRUITED) & (STATUS == "A"))) %>% #Remove non-fatal stroke 
 # Remove implausible values
 #filter(AGE <= 80 & AGE >= 40) %>% 
 filter(!((meanSBP > 270) | (meanSBP < 70))) %>%
 filter(!(BMI < 10 | BMI > 80)) %>%
 filter(!(Total_C > (20*38.67) | Total_C < (1.75*38.67) | is.na( Total_C))) %>%
 filter(!is.na(BASE_HBA1C),!is.na(WHtR),!is.na(BMI),!is.na(meanSBP),!is.na(meanDBP),!is.na(eGFR),!is.na(Total_C)) %>%
 # New variables for estimated age at diagnosis, number of comorbidities and medication use
 filter(!(is.na(BASE_DIABETES_DX) & BASE_DIABETES==1)) %>%
 mutate(year_diagnosis = case_when(
  BASE_DIABETES_DX == 1 ~ 1955,
  BASE_DIABETES_DX == 2 ~ 1965,
  BASE_DIABETES_DX == 3 ~ 1975,
  BASE_DIABETES_DX == 4 ~ 1985,
  BASE_DIABETES_DX == 5 ~ 1995,
  BASE_DIABETES_DX == 6 ~ 2005,
  TRUE ~ NA_real_
 )) %>%
 mutate(year_birth = YEAR_RECRUITED - AGE) %>% 
 mutate(tiempo_evol = YEAR_RECRUITED - year_diagnosis) %>%
 mutate(AGE_DX = ifelse(DM_DX == 1, AGE - tiempo_evol, AGE)) %>%
 mutate(AGE_DX =  ifelse(AGE_DX < 0 | is.na(AGE_DX), 0, AGE_DX)) %>%
 #filter(AGE_DX > 20 & AGE_DX < 100) %>%
 mutate(comorbidity_count = rowSums(select(., BASE_EMPHYSEMA, BASE_ASTHMA, BASE_CKD, BASE_PEP, BASE_CIRR, 
                                           BASE_HYPERTENSION, BASE_LUNGCANCER, BASE_OTHCANCER, BASE_PROSTATECANCER, 
                                           BASE_CERVCANCER, BASE_BREASTCANCER, BASE_STOMCANCER, BASE_ORALCANCER, BASE_PAD), 
                                    na.rm = TRUE)) 

#### Clusters ####
#setwd("/Users/Jeronimo/Library/CloudStorage/GoogleDrive-jeroperezalonso@gmail.com/.shortcut-targets-by-id/1RG-1Sg0NBmkz3AlzY5XL6XJ0mezx6oNg/Mexico City Cohort Study/Proyectos/Cardiovascular Risk DM/Diabetes Clusters/SNNN algorithm")
setwd("~/Google Drive/Mi unidad/Mexico City Cohort Study/Proyectos/Cardiovascular Risk DM/Diabetes Clusters/SNNN algorithm")
model<-keras::load_model_hdf5("Bases/mcps.h5", custom_objects = NULL, compile = FALSE)
#reticulate::use_virtualenv("r-tensorflow", required = TRUE)
mcps_cluster2<- mcps_dm_fin %>% 
 select(BMI, AGE_DX, BASE_HBA1C, meanSBP, meanDBP, WHtR)
bas<-read.csv("bas.csv")
m0<-apply(bas[,c(5:10)], 2, mean); std0<-apply(bas[,c(5:10)], 2, sd)
mcps_scale<-as.data.frame(scale(mcps_cluster2, scale=std0, center=m0))
input_array <- array(as.matrix(mcps_scale), dim = c(nrow(mcps_scale), ncol(mcps_scale)))
preds <- model$predict(input_array)
mcps_dm_fin$cluster_nolab<-max.col(preds, ties.method = "first")-1
mcps_dm_fin$cluster_nolab<-factor(mcps_dm_fin$cluster_nolab, labels = c("SIDD","SIRD","MOD","MARD"))
rm(mcps_cluster2,bas,mcps_scale,preds,m0,model,std0, input_array)
#setwd("/Users/Jeronimo/Library/CloudStorage/GoogleDrive-jeroperezalonso@gmail.com/.shortcut-targets-by-id/1RG-1Sg0NBmkz3AlzY5XL6XJ0mezx6oNg/Mexico City Cohort Study/Proyectos/Cardiovascular Risk DM")
 
setwd("~/Google Drive/Mi unidad/Mexico City Cohort Study/Proyectos/Cardiovascular Risk DM")
 
 SIDD = round(sum(mcps_dm_fin$cluster_nolab=="SIDD")/nrow(mcps_dm_fin),3)*100
 SIRD =round(sum(mcps_dm_fin$cluster_nolab=="SIRD")/nrow(mcps_dm_fin),3)*100
 MARD = round(sum(mcps_dm_fin$cluster_nolab=="MARD")/nrow(mcps_dm_fin),3)*100
 MOD = round(sum(mcps_dm_fin$cluster_nolab=="MOD")/nrow(mcps_dm_fin),3)*100

#### SCORE2-Diabetes by cluster ####
 mcps_dm_fin <- mcps_dm_fin %>%
  mutate(Risk.region = "Moderate",
         Gender = ifelse(MALE == 1, "male", "female"),
         total.chol = Total_C / 38.67,
         total.hdl = HDL_C / 38.67,
         HbA1c = (BASE_HBA1C * 10.23) - 23.5,
         SMOKER = as.integer(SMOKER),
         DM = as.integer(DM),
         AGE_DX = as.integer(AGE_DX))
 
 score2_diab <- function(x, y) {
  library(RiskScorescvd)
  x <- as.data.frame(x)
  score <- numeric(nrow(x)) 
  
  for (i in 1:nrow(x)) {
   score[i] <- SCORE2_Diabetes(
    Risk.region = x$Risk.region[i],
    Age = x$AGE[i],
    Gender = x$Gender[i],
    smoker = x$SMOKER[i],
    systolic.bp = x$meanSBP[i],
    total.chol = x$total.chol[i],
    total.hdl = x$total.hdl[i],
    diabetes = x$DM[i],
    diabetes.age = x$AGE_DX[i],
    HbA1c = x$HbA1c[i],
    eGFR = x$eGFR[i],
    classify = y
   )
  }
  
  return(score)
 }
 
 # Prepare data for SCORE2 Diabetes #
 mcps_dm_fin_SCORE2DM_RCS <- mcps_dm_fin %>%
  select(fatal_narrow_total,PERSON_YEARS)
 mcps_dm_fin_SCORE2DM_RCS$SCORE2_DM = score2_diab(mcps_dm_fin,F)
 
 # Calculate SCORE2 Diabetes
 mcps_dm_fin$SCORE2_DM = score2_diab(mcps_dm_fin,F)
 mcps_dm_fin$SCORE2_DM_class = score2_diab(mcps_dm_fin,T)
 
 mcps_dm_fin <- mcps_dm_fin %>% 
  mutate(cluster_nolab=factor(cluster_nolab)) %>%
  mutate(cluster_nolab=relevel(cluster_nolab, ref="MOD")) %>%
  mutate(tiempo_evol=AGE-AGE_DX, cluster_nolab=factor(cluster_nolab, levels=c("MOD","MARD", "SIRD", "SIDD"))) %>%
  mutate(SCORE2_DM_class=factor(SCORE2_DM_class,ordered = F)) %>%
  mutate(SCORE2_DM_class = factor(SCORE2_DM_class, levels = c("Low risk", "Moderate risk", "High risk", "Very high risk")))
 
 ## SCORE2-Diabetes categories 
 risk2<-c("Low (<5%)", "Moderate (5% to <10%)", "High (10% to <20%)", "Very high (≥20%)")
 risk3 <- c("Very high (≥20%)","High (10% to <20%)","Moderate (5% to <10%)","Low (<5%)")
 
#### Calculate SCORE2-Diabetes lp ####
 mcps_dm_fin <- mcps_dm_fin %>% 
  mutate(age_score = (AGE - 60) / 5,
         sbp_score = (meanSBP - 120) / 20,
         tcol_score = (Total_C/38.67)  - 6,
         hdl_score = ((HDL_C/38.67) - 1.3) / 0.5,
         HbA1c_mmol=10.929 * (BASE_HBA1C - 2.15),
         diab_diag = BASE_DIABETES,
         hba1c_score=(HbA1c_mmol-31)/9.34,
         cagediab=(AGE_DX-50)/5,
         lnegfr=(log(eGFR)-4.5)/0.15,
         male=MALE,
         smoke = case_when(SMOKEGP %in% c(3, 4, 5) ~ 1,
                           SMOKEGP %in% c(1, 2) ~ 0,
                           .default = NA)) %>%
  mutate(lp_score = if_else(MALE == 1,
                            0.3742*age_score + 0.6012*smoke + 0.2777*sbp_score + 0.1458*tcol_score - 0.2698*hdl_score
                            - 0.0755*smoke*age_score - 0.0255*sbp_score*age_score - 0.0281*tcol_score*age_score + 0.0426*hdl_score*age_score,
                            
                            0.4648*age_score + 0.7744*smoke + 0.3131*sbp_score + 0.1002*tcol_score - 0.2606*hdl_score
                            - 0.1088*smoke*age_score - 0.0277*sbp_score*age_score - 0.0226*tcol_score*age_score + 0.0613*hdl_score*age_score)) %>%
  mutate(lp_dm= if_else(MALE == 1,
                        0.0955*hba1c_score - 0.0591*lnegfr + 0.0058*lnegfr*lnegfr - 0.09998*cagediab 
                        - 0.0134*hba1c_score*age_score + 0.0115*lnegfr*age_score,
                        
                        0.1173*hba1c_score - 0.0640*lnegfr + 0.0062*lnegfr*lnegfr - 0.118*cagediab 
                        - 0.0196*hba1c_score*age_score + 0.0169*lnegfr*age_score),
         score2_vars= if_else(male == 1,
                              0.5368*age_score + 0.4774*smoke + 0.1322*sbp_score + 0.6457*diab_diag + 0.1102*tcol_score - 0.1087*hdl_score
                              - 0.0672*smoke*age_score - 0.0268*sbp_score*age_score - 0.0181*tcol_score*age_score + 0.0095*hdl_score*age_score-0.0983*diab_diag*age_score,
                              
                              0.6624*age_score + 0.6139*smoke + 0.1421*sbp_score + 0.8096*diab_diag + 0.1127*tcol_score - 0.1568*hdl_score
                              - 0.1122*smoke*age_score - 0.0167*sbp_score*age_score - 0.0200*tcol_score*age_score + 0.0186*hdl_score*age_score-0.1272*diab_diag*age_score),
         lp_fin=lp_score+lp_dm+score2_vars) %>% 
  mutate(score2_dm_risk_uncabil = if_else(male == 1,
                                          1 - 0.9605^exp(lp_fin),
                                          
                                          1 - 0.9776^exp(lp_fin)),
         
         #Calibration according to moderate risk region scaling factors
         score2_dm_risk = if_else(male == 1,
                                  1-exp(-exp(-0.1565 + 0.8009 * log(-log(1 - score2_dm_risk_uncabil)))),
                                  
                                  1-exp(-exp(-0.3143 + 0.7701 * log(-log(1 - score2_dm_risk_uncabil))))),
         score2_dm_risk_low = if_else(male == 1,
                                  1-exp(-exp(-0.5699 + 0.7476 * log(-log(1 - score2_dm_risk_uncabil)))),
                                      
                                  1-exp(-exp(-0.7380 + 0.7019 * log(-log(1 - score2_dm_risk_uncabil))))),
         score2_dm_risk_high = if_else(male == 1,
                                  1-exp(-exp(0.3207 + 0.9360 * log(-log(1 - score2_dm_risk_uncabil)))),
                                       
                                  1-exp(-exp(0.5710 + 0.9369 * log(-log(1 - score2_dm_risk_uncabil))))),
         score2_dm_risk_vhigh = if_else(male == 1,
                                  1-exp(-exp(0.5836 + 0.8294 * log(-log(1 - score2_dm_risk_uncabil)))),
                                        
                                  1-exp(-exp(0.9412 + 0.8329 * log(-log(1 - score2_dm_risk_uncabil))))))
 
 #10-year administrative censoring
 #CVD fatal outcomes
 mcpx_pre_cox <- survSplit(Surv(PERSON_YEARS, fatal_narrow_total) ~ ., data = mcps_dm_fin, end = "tstop", cut = c(10), episode = "cut")
 mcps_cox <- subset(mcpx_pre_cox, cut == 1)
 
 mcps_cox$score2_dm_risk <- pmin(pmax(mcps_cox$score2_dm_risk, 1e-6), 1 - 1e-6)
 mcps_cox$loglog_score2 <- log(-log(1 - mcps_cox$score2_dm_risk))
 
 mcps_cox<- mcps_cox%>% filter(!is.na(loglog_score2))
 
#### Recalibration of SCORE2-Diabetes for fatal outcomes ####
 fit_slope <- coxph(Surv(tstop, fatal_narrow_total) ~ lp_fin, data = mcps_cox)
 coef(fit_slope) 
 confint(fit_slope) # recalibration slope (beta_hat)
 # To compute predicted survival at time t0 (e.g., 10 years):
 t0 <- 10
 # baseline survival S0(t) for lp = 0
 sf0 <- survfit(fit_slope, newdata = data.frame(lp_fin = 0))
 S0_t0 <- summary(sf0, times = t0)$surv
 # final predicted risk for each subject
 mcps_cox$score2diab_recalib <- 1 - S0_t0 ^ exp(coef(fit_slope) * mcps_cox$lp_fin)
 
 mcps_dm_fin$score2diab_recalib <- 1 - S0_t0 ^ exp(coef(fit_slope) * mcps_dm_fin$lp_fin)
 mcps_cox$loglog_score2_2 <- log(-log(1 - mcps_cox$score2diab_recalib))
 
#### Discrimination of SCORE2-Diabetes ####
 
 cox_score <- cph(Surv(tstop, fatal_narrow_total) ~ loglog_score2, x = TRUE, y=TRUE, data = mcps_cox)
 
 val <- validate(cox_score, method="boot", B=300, )
 
 cox_score_mod <- coxph(Surv(tstop, fatal_narrow_total) ~ loglog_score2, x = TRUE, data = mcps_cox %>% filter(cluster_nolab=="MOD"))
 cox_score_mard <- coxph(Surv(tstop, fatal_narrow_total) ~ loglog_score2, x = TRUE, data = mcps_cox %>% filter(cluster_nolab=="MARD"))
 cox_score_sird <- coxph(Surv(tstop, fatal_narrow_total) ~ loglog_score2, x = TRUE, data = mcps_cox %>% filter(cluster_nolab=="SIRD"))
 cox_score_sidd <- coxph(Surv(tstop, fatal_narrow_total) ~ loglog_score2, x = TRUE, data = mcps_cox %>% filter(cluster_nolab=="SIDD"))
 
 cox_score2 <- coxph(Surv(tstop, fatal_narrow_total) ~ loglog_score2_2, x = TRUE, data = mcps_cox)
 cox_score_mod2 <- coxph(Surv(tstop, fatal_narrow_total) ~ loglog_score2_2, x = TRUE, data = mcps_cox %>% filter(cluster_nolab=="MOD"))
 cox_score_mard2 <- coxph(Surv(tstop, fatal_narrow_total) ~ loglog_score2_2, x = TRUE, data = mcps_cox %>% filter(cluster_nolab=="MARD"))
 cox_score_sird2 <- coxph(Surv(tstop, fatal_narrow_total) ~ loglog_score2_2, x = TRUE, data = mcps_cox %>% filter(cluster_nolab=="SIRD"))
 cox_score_sidd2 <- coxph(Surv(tstop, fatal_narrow_total) ~ loglog_score2_2, x = TRUE, data = mcps_cox %>% filter(cluster_nolab=="SIDD"))
 
 models<-list(cox_score,cox_score_mod,cox_score_mard,cox_score_sird,cox_score_sidd,
              cox_score2,cox_score_mod2,cox_score_mard2,cox_score_sird2,cox_score_sidd2)
 c_stat <- map_dbl(models, ~ round(summary(.x)$concordance[1],3))
 se <- map_dbl(models, ~ round(summary(.x)$concordance[2],3))
 c_stat_total <- tibble("cluster_nolab" = c("Overall", "MOD", "MARD", "SIRD", "SIDD",
                                            "Overall", "MOD", "MARD", "SIRD", "SIDD"),
                        "recalib"=c(rep("Original",5), rep("Recalibrated",5)),
                        "C-statistic" = c_stat,
                        "lower" = round((c_stat - 1.96*se), 3),
                        "upper" = round((c_stat + 1.96*se),3)) %>% 
  mutate("95% CI" = paste0(`C-statistic`," (", lower, "-", upper, ")")) %>% 
  select(1,2,6) 
 
 c_stat_total
 
 
#### Mean Calibration of SCORE2-Diabetes ####
 #CVD outcomes
 obs_overall <- summary(survfit(Surv(tstop, fatal_narrow_total) ~ 1,
                                data = mcps_cox),times = 10)
 
 obs_mard <- summary(survfit(Surv(tstop, fatal_narrow_total) ~ 1,
                             data = mcps_cox, subset = (cluster_nolab == "MARD")),times = 10)
 
 obs_mod <- summary(survfit(Surv(tstop, fatal_narrow_total) ~ 1,
                            data = mcps_cox, subset = (cluster_nolab == "MOD")),times = 10)
 
 obs_sird <- summary(survfit(Surv(tstop, fatal_narrow_total) ~ 1,
                             data = mcps_cox, subset = (cluster_nolab == "SIRD")),times = 10)
 
 obs_sidd <- summary(survfit(Surv(tstop, fatal_narrow_total) ~ 1,
                             data = mcps_cox, subset = (cluster_nolab == "SIDD")),times = 10)
 
 
 obs_rsk_overall <- 1 - obs_overall$surv
 obs_rsk_mard <- 1 - obs_mard$surv
 obs_rsk_mod<- 1 - obs_mod$surv
 obs_rsk_sird <- 1 - obs_sird$surv
 obs_rsk_sidd <- 1 - obs_sidd$surv
 
 mean_calib_cvd_clusters <- mcps_cox %>% 
  select(cluster_nolab, score2_dm_risk) %>% 
  group_by(cluster_nolab) %>% 
  summarise(score=mean(score2_dm_risk, na.rm=T)) %>% 
  mutate(obs_rsk = c(obs_rsk_mod,obs_rsk_mard,obs_rsk_sird,obs_rsk_sidd)) %>% 
  mutate(mean_cal = obs_rsk / score,
         lower = case_when(cluster_nolab == "MOD" ~ mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sidd$n.event)),
                           cluster_nolab == "MARD" ~ mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sird$n.event)),
                           cluster_nolab == "SIRD" ~ mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_mard$n.event)),
                           cluster_nolab == "SIDD" ~ mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_mod$n.event))),
         upper = case_when(cluster_nolab == "MOD" ~ mean_cal * exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sidd$n.event)),
                           cluster_nolab == "MARD" ~ mean_cal * exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sird$n.event)),
                           cluster_nolab == "SIRD" ~ mean_cal * exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_mard$n.event)),
                           cluster_nolab == "SIDD" ~ mean_cal * exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_mod$n.event))),
         alpha = log(mean_cal)) %>% rbind(mcps_cox %>% 
                                           select(score2_dm_risk) %>% 
                                           summarise(score=mean(score2_dm_risk, na.rm=T)) %>% 
                                           mutate(obs_rsk = c(obs_rsk_overall)) %>% 
                                           mutate(mean_cal = obs_rsk / score,
                                                  lower = mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_overall$n.event)),
                                                  upper = exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sidd$n.event)),
                                                  alpha = log(mean_cal), cluster_nolab="Overall") %>% select(7,1:6)) %>%
  mutate("95% CI" = paste0(round(mean_cal,3)," (", round(lower,3), " - ", round(upper,3), ")")) %>% select(1:3,8) %>%
  mutate(recalib=c(rep("Original",5)))
 
 mean_calib_cvd_clusters2 <- mcps_cox %>% 
  select(cluster_nolab, score2diab_recalib) %>% 
  group_by(cluster_nolab) %>% 
  summarise(score=mean(score2diab_recalib, na.rm=T)) %>% 
  mutate(obs_rsk = c(obs_rsk_mod,obs_rsk_mard,obs_rsk_sird,obs_rsk_sidd)) %>% 
  mutate(mean_cal = obs_rsk / score,
         lower = case_when(cluster_nolab == "MOD" ~ mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sidd$n.event)),
                           cluster_nolab == "MARD" ~ mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sird$n.event)),
                           cluster_nolab == "SIRD" ~ mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_mard$n.event)),
                           cluster_nolab == "SIDD" ~ mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_mod$n.event))),
         upper = case_when(cluster_nolab == "MOD" ~ mean_cal * exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sidd$n.event)),
                           cluster_nolab == "MARD" ~ mean_cal * exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sird$n.event)),
                           cluster_nolab == "SIRD" ~ mean_cal * exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_mard$n.event)),
                           cluster_nolab == "SIDD" ~ mean_cal * exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_mod$n.event))),
         alpha = log(mean_cal)) %>% rbind(mcps_cox %>% 
                                           select(score2diab_recalib) %>% 
                                           summarise(score=mean(score2diab_recalib, na.rm=T)) %>% 
                                           mutate(obs_rsk = c(obs_rsk_overall)) %>% 
                                           mutate(mean_cal = obs_rsk / score,
                                                  lower = mean_cal * exp(-qnorm(1 - 0.05 / 2) * sqrt(1 / obs_overall$n.event)),
                                                  upper = exp(+qnorm(1 - 0.05 / 2) * sqrt(1 / obs_sidd$n.event)),
                                                  alpha = log(mean_cal), cluster_nolab="Overall") %>% select(7,1:6)) %>%
  mutate("95% CI" = paste0(round(mean_cal,3)," (", round(lower,3), " - ", round(upper,3), ")")) %>% select(1:3,8) %>%
  mutate(recalib=c(rep("Recalibrated",5)))
 
 
 mean_calib_cvd_clusters<-rbind(mean_calib_cvd_clusters, mean_calib_cvd_clusters2)
 
#### Weak Calibration of SCORE2-Diabetes ####
 
 slope_fun <- function(x) {
  df <- data.frame(coef = x$coef,
                   lower = x$coef - qnorm(1 - 0.05 / 2) * sqrt(x$var),
                   upper = x$coef + qnorm(1 - 0.05 / 2) * sqrt(x$var))
  df}
 
 slope_list <- lapply(models, slope_fun)
 slope <- do.call(rbind, slope_list)
 weak_calib <- slope %>% 
  rownames_to_column(var = "models") %>% 
  mutate(cluster_nolab=c("Overall", "MOD", "MARD", "SIRD", "SIDD","Overall", "MOD", "MARD", "SIRD", "SIDD"),
         across(2:4, ~ round(.x, 4)),
         "95% CI" = paste0(coef," (", lower, " - ", upper, ")")) %>% 
  select(5,6) %>% mutate(recalib=c(rep("Original",5), rep("Recalibrated",5)))
 
#### Overall performance ####
 #Score for models with CVD outcomes  
 overall_cvd_cluster <- riskRegression::Score(models, 
                                              formula = Surv(tstop, fatal_narrow_total) ~ 1,
                                              data = mcps_cox, 
                                              conf.int = TRUE, 
                                              times = 10,
                                              cens.model = "km", 
                                              metrics = "brier",
                                              summary = "ipa") %>% 
  pluck("Brier", "score") %>% 
  select(1,3,5,6,7) %>% 
  mutate(across(2:5, ~ round(.x, 4)),
         IPA = IPA*100) %>% 
  filter(model != "Null model")
 
 ov_tot_cvd <- overall_cvd_cluster %>% 
  mutate(estimate = paste0(Brier, " (", lower, "-", upper, ")")) %>% 
  select(!contains(c("Brier", "upper", "lower"))) %>%
  mutate(cluster_nolab=c("Overall", "MOD", "MARD", "SIDD", "SIRD","Overall", "MOD", "MARD", "SIDD", "SIRD")) %>%
  select(4,3,2) %>% mutate(recalib=c(rep("Original",5), rep("Recalibrated",5)))
 
#### Moderate calibration ####
 mod_calib_cvd_df <- mcps_cox %>% 
  select(PATID,score2_dm_risk, score2diab_recalib,tstop, fatal_narrow_total) %>% 
  mutate(log_min = across(2:3, ~ log(-log(1 - .x)))) %>% 
  tidyr::unpack(log_min, names_sep = "_")
 
 vars_cvd <- list("log_min_score2_dm_risk","log_min_score2diab_recalib")
 
 calib_models_cvd_overall <- vars_cvd %>% 
  map(function(x) {
   formula <- as.formula(paste0("Surv(tstop, fatal_narrow_total) ~ rcs(", x, ",3)"))
   model <- cph(formula, x = TRUE, y = TRUE, surv = TRUE, data = mod_calib_cvd_df)
  })
 
 observed01 <- survest(calib_models_cvd_overall[[1]], times = 10, newdata = mod_calib_cvd_df)
 
 obs_score_dm_original <- data.frame(obs = 1 - c(observed01$surv),
                                     upper = 1 - c(observed01$lower),
                                     lower = 1 - c(observed01$upper),
                                     pred = mod_calib_cvd_df$score2_dm_risk)
 
 observed02 <- survest(calib_models_cvd_overall[[2]], times = 10, newdata = mod_calib_cvd_df)
 
 obs_score_dm_recalib <- data.frame(obs = 1 - c(observed02$surv),
                                    upper = 1 - c(observed02$lower),
                                    lower = 1 - c(observed02$upper),
                                    pred = mod_calib_cvd_df$score2diab_recalib)
 
 models <- list(obs_score_dm_original,obs_score_dm_recalib)
 cluster<-list("Original model", "Recalibrated")
 names <- list("Original model", "Recalibrated")
 
 obs_clust_total <-bind_rows(pmap(list(models, names, cluster), apply_name))
 
 fig_calib_total_cluster <- obs_clust_total %>% 
  ggplot(aes(x = pred, y = obs, ymin = lower, ymax = upper, color = cluster)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = "Moderate calibration of SCORE2-Diabetes",
       x = "Predicted risk",
       y = "Observed risk",
       color = "")+theme(legend.position = "top")
 
 ggsave(fig_calib_total_cluster,filename = "Figuras/SuppFigure2.jpg", 
        width = 15, height = 10,units=c("cm"),
        dpi = 300,limitsize = FALSE) 
 
 
 
 ### Per cluster ###
 
 mod_calib_cvd_df_clust <- mcps_cox %>% 
  select(cluster_nolab,PATID,score2_dm_risk, tstop, fatal_narrow_total) %>% 
  mutate(log_min = log(-log(1 - score2_dm_risk))) %>% 
  tidyr::unpack(log_min, names_sep = "_")
 
 vars_cvd <- list("log_min")
 
 calib_models_cvd_clusters <- mod_calib_cvd_df_clust %>%
  group_split(cluster_nolab) %>%
  map(function(df_cluster) {
   map(vars_cvd, function(x) {
    formula <- as.formula(
     paste0("Surv(tstop, fatal_narrow_total) ~ rcs(", x, ", 3)")
    )
    cph(formula,
        x = TRUE, y = TRUE, surv = TRUE,
        data = df_cluster)
   })
  })
 
 dd <- datadist(mod_calib_cvd_df_clust)
 options(datadist="dd")
 
 observed1 <- survest(calib_models_cvd_clusters[[1]][[1]], times = 10, newdata = mod_calib_cvd_df_clust %>% filter(cluster_nolab=="SIDD"))
 observed2 <- survest(calib_models_cvd_clusters[[2]][[1]], times = 10, newdata = mod_calib_cvd_df_clust %>% filter(cluster_nolab=="SIRD"))
 observed3 <- survest(calib_models_cvd_clusters[[3]][[1]], times = 10, newdata = mod_calib_cvd_df_clust %>% filter(cluster_nolab=="MOD"))
 observed4 <- survest(calib_models_cvd_clusters[[4]][[1]], times = 10, newdata = mod_calib_cvd_df_clust %>% filter(cluster_nolab=="MARD"))
 
 obs_score_sidd <- data.frame(obs = 1 - c(observed1$surv),
                              upper = 1 - c(observed1$lower),
                              lower = 1 - c(observed1$upper),
                              pred = mod_calib_cvd_df_clust$score2_dm_risk[mod_calib_cvd_df_clust$cluster_nolab=="SIDD"])
 
 obs_score_sird <- data.frame(obs = 1 - c(observed2$surv),
                              upper = 1 - c(observed2$lower),
                              lower = 1 - c(observed2$upper),
                              pred = mod_calib_cvd_df_clust$score2_dm_risk[mod_calib_cvd_df_clust$cluster_nolab=="SIRD"])
 
 obs_score_mod <- data.frame(obs = 1 - c(observed3$surv),
                             upper = 1 - c(observed3$lower),
                             lower = 1 - c(observed3$upper),
                             pred = mod_calib_cvd_df_clust$score2_dm_risk[mod_calib_cvd_df_clust$cluster_nolab=="MOD"])
 
 obs_score_mard <- data.frame(obs = 1 - c(observed4$surv),
                              upper = 1 - c(observed4$lower),
                              lower = 1 - c(observed4$upper),
                              pred = mod_calib_cvd_df_clust$score2_dm_risk[mod_calib_cvd_df_clust$cluster_nolab=="MARD"])
 
 models <- list(obs_score_dm_original,obs_score_sidd,obs_score_sird,obs_score_mod,obs_score_mard)
 cluster<-list("Overall", "SIDD", "SIRD", "MOD", "MARD")
 names <- list("SCORE2-Diabetes")
 
 obs_clust_total <-bind_rows(pmap(list(models, rep(names, length(models)), cluster), apply_name))
 
 fig_calib_total_cluster <- obs_clust_total %>% 
  ggplot(aes(x = pred, y = obs, ymin = lower, ymax = upper, color = cluster)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = "Moderate calibration per diabetes cluster",
       x = "Predicted risk",
       y = "Observed risk",
       color = "")+theme(legend.position = "top")
 
 ## Per cluster
 mod_calib_cvd_df_clust <- mcps_cox %>% 
  select(cluster_nolab,PATID,score2_dm_risk, score2diab_recalib, tstop, fatal_narrow_total) %>% 
  mutate(log_min = log(-log(1 - score2diab_recalib))) %>% 
  tidyr::unpack(log_min, names_sep = "_")
 
 vars_cvd <- list("log_min")
 
 calib_models_cvd_clusters <- mod_calib_cvd_df_clust %>%
  group_split(cluster_nolab) %>%
  map(function(df_cluster) {
   map(vars_cvd, function(x) {
    formula <- as.formula(
     paste0("Surv(tstop, fatal_narrow_total) ~ rcs(", x, ", 3)")
    )
    cph(formula,
        x = TRUE, y = TRUE, surv = TRUE,
        data = df_cluster)
   })
  })
 
 dd <- datadist(mod_calib_cvd_df_clust)
 options(datadist="dd")
 
 observed1 <- survest(calib_models_cvd_clusters[[1]][[1]], times = 10, newdata = mod_calib_cvd_df_clust %>% filter(cluster_nolab=="SIDD"))
 observed2 <- survest(calib_models_cvd_clusters[[2]][[1]], times = 10, newdata = mod_calib_cvd_df_clust %>% filter(cluster_nolab=="SIRD"))
 observed3 <- survest(calib_models_cvd_clusters[[3]][[1]], times = 10, newdata = mod_calib_cvd_df_clust %>% filter(cluster_nolab=="MOD"))
 observed4 <- survest(calib_models_cvd_clusters[[4]][[1]], times = 10, newdata = mod_calib_cvd_df_clust %>% filter(cluster_nolab=="MARD"))
 
 obs_score_sidd <- data.frame(obs = 1 - c(observed1$surv),
                              upper = 1 - c(observed1$lower),
                              lower = 1 - c(observed1$upper),
                              pred = mod_calib_cvd_df_clust$score2diab_recalib[mod_calib_cvd_df_clust$cluster_nolab=="SIDD"])
 
 obs_score_sird <- data.frame(obs = 1 - c(observed2$surv),
                              upper = 1 - c(observed2$lower),
                              lower = 1 - c(observed2$upper),
                              pred = mod_calib_cvd_df_clust$score2diab_recalib[mod_calib_cvd_df_clust$cluster_nolab=="SIRD"])
 
 obs_score_mod <- data.frame(obs = 1 - c(observed3$surv),
                             upper = 1 - c(observed3$lower),
                             lower = 1 - c(observed3$upper),
                             pred = mod_calib_cvd_df_clust$score2diab_recalib[mod_calib_cvd_df_clust$cluster_nolab=="MOD"])
 
 obs_score_mard <- data.frame(obs = 1 - c(observed4$surv),
                              upper = 1 - c(observed4$lower),
                              lower = 1 - c(observed4$upper),
                              pred = mod_calib_cvd_df_clust$score2diab_recalib[mod_calib_cvd_df_clust$cluster_nolab=="MARD"])
 
 models <- list(obs_score_dm_original,obs_score_dm_recalib,obs_score_sidd,obs_score_sird,obs_score_mod,obs_score_mard)
 cluster<-list("Overall original", "Overall recalibrated","SIDD recalibrated", "SIRD recalibrated", "MOD recalibrated", "MARD recalibrated")
 names <- list("SCORE2-Diabetes")
 
 obs_clust_total <-bind_rows(pmap(list(models, rep(names, length(models)), cluster), apply_name))
 
 fig_calib_total_cluster <- obs_clust_total %>% 
  ggplot(aes(x = pred, y = obs, ymin = lower, ymax = upper, color = cluster)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = 2) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = "Moderate calibration of SCORE2-Diabetes per diabetes cluster",
       x = "Predicted risk",
       y = "Observed risk",
       color = "")+theme(legend.position = "top")
 
 
#### Optimism corrected metrics ######
 
validate_score <- function(data, risk_var, time_point = 10, B = 300){
   
   # Remove NAs and ensure risk is strictly between 0 and 1 for log-log transform
   data <- data[!is.na(data[[risk_var]]), ]
   data <- data[data[[risk_var]] > 0 & data[[risk_var]] < 1, ]
   n <- nrow(data)
   
   # Convert risk -> linear predictor (the 'score')
   data$lp <- log(-log(1 - data[[risk_var]]))
   
   # Define datadist for rms
   dd <- datadist(data)
   options(datadist = "dd")
   
   fit <- cph(Surv(tstop, fatal_narrow_total) ~ lp,
              data=data, x=TRUE, y=TRUE,
              surv=TRUE, time.inc=time_point)
   
   beta_apparent <- coef(fit)["lp"]
   
   val <- validate(fit, method="boot", B=B)
   shrinkage_factor <- val["Slope", "index.corrected"]
   
   slope_corr <- beta_apparent * shrinkage_factor
   
   dxy_app  <- val["Dxy", "index.orig"]
   dxy_corr <- val["Dxy", "index.corrected"]
   c_app    <- (dxy_app + 1) / 2
   c_corr   <- (dxy_corr + 1) / 2
   
   boot_slope <- numeric(B)
   boot_c     <- numeric(B)
   
   set.seed(123)
   for(b in 1:B){
     index <- sample.int(n, replace=TRUE)
     dboot <- data[index, ]
     
     fit_b <- cph(Surv(tstop, fatal_narrow_total) ~ lp,
                  data=dboot, x=TRUE, y=TRUE)
     
     boot_slope[b] <- coef(fit_b)["lp"]
     boot_c[b]     <- (fit_b$stats["Dxy"] + 1) / 2
   }
   
   slope_ci <- quantile(boot_slope, c(0.025, 0.975), na.rm=TRUE)
   c_ci     <- quantile(boot_c, c(0.025, 0.975), na.rm=TRUE)
   
   sf <- survfit(Surv(tstop, fatal_narrow_total) ~ 1, data=data)
   obs <- summary(sf, times=time_point)
   
   obs_risk <- 1 - obs$surv
   exp_risk <- mean(data[[risk_var]])
   
   # Define OE clearly before the return
   current_OE <- obs_risk / exp_risk
   
   return(data.frame(
     OE_ratio        = current_OE,
     slope_apparent  = as.numeric(beta_apparent),
     slope_corrected = as.numeric(slope_corr),
     slope_lower95   = slope_ci[1],
     slope_upper95   = slope_ci[2],
     c_apparent      = c_app,
     c_corrected     = c_corr,
     c_lower95       = c_ci[1],
     c_upper95       = c_ci[2],
     row.names       = NULL
   ))
 }

overall_recal <- validate_score(mod_calib_cvd_df, "score2diab_recalib")
overall_recal$model <- "Recalibrated"
overall_recal$cluster <- "Overall"
clusters <- unique(mod_calib_cvd_df_clust$cluster_nolab)
 
cluster_results <- lapply(clusters, function(cl){
   
   data_sub <- subset(mod_calib_cvd_df_clust, cluster_nolab == cl)
   
   r2 <- validate_score(data_sub, "score2diab_recalib")
   r2$model <- "Recalibrated"
   r2$cluster <- cl
   
   r2
 })
 
cluster_results <- bind_rows(cluster_results)
 
final_table <- bind_rows(overall_recal,
                          cluster_results) %>%
   select(cluster, model,
          OE_ratio,
          slope_apparent,
          slope_corrected,
          c_apparent,
          c_corrected)%>%
   flextable()
 
 read_docx() %>%
   body_add_flextable(value = final_table, split = TRUE) %>%
   body_end_section_landscape() %>%
   print(target = "Tablas/Table - Optimism-corrected performance.docx")
 
 
 #### Table S1: Characteristics of Mexico City Prospective Study Cohort and Diabetes Subpopulation ####
 var_labels <- list(
  AGE = "Age (years)",
  PERSON_YEARS = "Follow-up (years)",
  EDU_LEVEL = "Educational attainment (%)",
  SMOKER = "Smoking (%)",
  BMI = "Body mass index (kg/m²)",
  meanSBP = "Systolic blood pressure (mmHg)",
  Glucose = "Glucose (mg/dl)",
  BASE_HBA1C = "Glycated HbA1c (%)",
  Total_C = "Total cholesterol (mg/dl)",
  HDL_C = "HDL cholesterol (mg/dl)",
  non_HDL_C = "Non-HDL cholesterol (mg/dl)",
  Total_TG = "Total triglycerides (mg/dl)",
  ApoB = "Apolipoprotein B (mg/dl)",
  Creatinine = "Creatinine (mg/dl)",
  eGFR = "eGFR (ml/min/1.73m²)")
 
 vars <- c("AGE","PERSON_YEARS","EDU_LEVEL","SMOKER","BMI","meanSBP","Glucose",
           "BASE_HBA1C","Total_C","HDL_C","non_HDL_C","Total_TG","ApoB",
           "Creatinine","eGFR")
 
 # All
 tab_mcps_all <- make_summary(mcps)
 
 # Men
 tab_mcps_men <- make_summary(mcps %>% filter(MALE == 1),
                              header(mcps %>% filter(MALE == 1), mcps))
 
 # Women
 tab_mcps_women <- make_summary(mcps %>% filter(MALE == 0),
                                header(mcps %>% filter(MALE == 0), mcps))
 
 # DM - All
 tab_mcps_dm_all <- make_summary(mcps_dm)
 
 # DM - Men
 tab_mcps_dm_men <- make_summary(mcps_dm %>% filter(MALE == 1),
                                 header(mcps_dm %>% filter(MALE == 1), mcps_dm))
 
 # DM - Women
 tab_mcps_dm_women <- make_summary(mcps_dm %>% filter(MALE == 0),
                                   header(mcps_dm %>% filter(MALE == 0), mcps_dm))
 
 # Merge tables
 tab_characteristics <- tbl_merge(
  list(tab_mcps_all, tab_mcps_men, tab_mcps_women,
       tab_mcps_dm_all, tab_mcps_dm_men, tab_mcps_dm_women),
  tab_spanner = c("All","Male","Female","All","Male","Female")) %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  as_flex_table() %>%
  add_header_row(
   top = TRUE,
   values = c("", 
              paste0("MCPS (n = ", format(nrow(mcps), big.mark = ","), ")"),
              paste0("MCPS with Diabetes (n = ", format(nrow(mcps_dm), big.mark = ","), ")")),
   colwidths = c(1,3,3)) %>%
  add_header_lines("Table S1. Characteristics of Mexico City Prospective Study Cohort and Diabetes Subpopulation") %>%
  bold(part = "header") %>%
  theme_vanilla() %>%
  align(j = 2:4, align = "center", part = "header") %>%
  align(j = 5:7, align = "center", part = "header") %>%
  add_footer_lines("Mean (SD) for continuous variables, n (%) for categorical variables.")
 
 # Export
 read_docx() %>%
  body_add_flextable(value = tab_characteristics, split = TRUE) %>%
  body_end_section_landscape() %>% 
  print(target = "Tablas/Table S1 - Characteristics of Mexico City Prospective Study Cohort and Diabetes Subpopulation_2.docx")
 
 rm(tab_mcps_all, tab_mcps_men, tab_mcps_women,tab_mcps_dm_all, tab_mcps_dm_men, tab_mcps_dm_women,tab_characteristics)
 
 #### Table S2: Characteristics of participants with diabetes ####
 tab_mcps_dm_all <- make_summary(mcps_dm_fin)
 
 tab_mcps_dm_men <- make_summary(mcps_dm_fin %>% filter(MALE == 1),
                                 header(mcps_dm_fin %>% filter(MALE == 1), mcps_dm_fin))
 
 tab_mcps_dm_women <- make_summary(mcps_dm_fin %>% filter(MALE == 0),
                                   header(mcps_dm_fin %>% filter(MALE == 0), mcps_dm_fin))
 
 tab_characteristics <- tbl_merge(
  list(tab_mcps_dm_all, tab_mcps_dm_men, tab_mcps_dm_women),
  tab_spanner = c("All","Male","Female")) %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  as_flex_table() %>%
  add_header_lines("Table S1.1 Characteristics of Mexico City Prospective Study Cohort with Diabetes") %>%
  bold(part = "header") %>%
  theme_vanilla() %>%
  add_footer_lines("Mean (SD) for continuous variables, n (%) for categorical variables.")
 
 read_docx() %>%
  body_add_flextable(value = tab_characteristics, split = TRUE) %>%
  body_end_section_landscape() %>% 
  print(target = "Tablas/Table S2 - Characteristics of Mexico City Prospective Study Cohort with Diabetes_1.docx")
 
 rm(tab_mcps_all, tab_mcps_men, tab_mcps_women,tab_mcps_dm_all, tab_mcps_dm_men, tab_mcps_dm_women,tab_characteristics)
 
 #### Table 1: Characteristics by clusters ####
 
 var_labels <- list(
  AGE         = "Age (years)",
  MALE        = "Male (%)",
  AGE_DX      = "Estimated age at diabetes diagnosis (years)",
  DM_noDX     = "Undiagnosed diabetes (%)",
  SMOKER      = "Smoking (%)",
  BMI         = "Body mass index (kg/m²)",
  meanSBP     = "Systolic blood pressure (mmHg)",
  meanDBP     = "Diastolic blood pressure (mmHg)",
  WHR         = "Waist-To-Hip Ratio",
  WHtR        = "Waist-To-Height Ratio",
  Glucose     = "Glucose (mg/dl)",
  BASE_HBA1C  = "Glycated HbA1c (%)",
  mets_ir     = "METS-IR",
  Creatinine  = "Creatinine (mg/dl)",
  eGFR        = "eGFR (ml/min/1.73m²)",
  Total_C     = "Total cholesterol (mg/dl)",
  HDL_C       = "HDL cholesterol (mg/dl)",
  non_HDL_C   = "Non-HDL cholesterol (mg/dl)",
  Total_TG    = "Total triglycerides (mg/dl)",
  ApoB        = "Apolipoprotein B (mg/dl)",
  DRUG_D1     = "Metformin (%)",
  DRUG_D2     = "Sulphonylurea (%)",
  DRUG_D3     = "Insulin (%)",
  DRUG_D4     = "Other antidiabetics (%)")
 
 vars <- c("AGE","MALE","AGE_DX","DM_noDX","SMOKER","BMI","meanSBP","meanDBP","WHR","WHtR",
           "Glucose","BASE_HBA1C","mets_ir","Creatinine","eGFR",
           "Total_C","HDL_C","non_HDL_C","Total_TG","ApoB",
           "DRUG_D1","DRUG_D2","DRUG_D3","DRUG_D4")
 
 tab_dm_all <- make_summary(mcps_dm_fin, header(mcps_dm_fin, mcps_dm_fin))
 tab_sidd<- make_summary(filter(mcps_dm_fin, cluster_nolab == "SIDD"),
                         header(filter(mcps_dm_fin, cluster_nolab == "SIDD"), mcps_dm_fin))
 tab_mod<- make_summary(filter(mcps_dm_fin, cluster_nolab == "MOD"),
                        header(filter(mcps_dm_fin, cluster_nolab == "MOD"), mcps_dm_fin))
 tab_mard<- make_summary(filter(mcps_dm_fin, cluster_nolab == "MARD"),
                         header(filter(mcps_dm_fin, cluster_nolab == "MARD"), mcps_dm_fin))
 tab_sird<- make_summary(filter(mcps_dm_fin, cluster_nolab == "SIRD"),
                         header(filter(mcps_dm_fin, cluster_nolab == "SIRD"), mcps_dm_fin))
 
 tab_clusters <- tbl_merge(
  list(tab_dm_all, tab_sidd, tab_mod, tab_mard, tab_sird),
  tab_spanner = c("All", "SIDD", "MOD", "MARD", "SIRD")
 ) %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  as_flex_table() %>%
  add_header_lines("Table 1. Characteristics of adult-onset diabetes subgroups in the Mexico City Prospective Study") %>%
  bold(part = "header") %>%
  theme_vanilla() %>%
  add_footer_lines(as_paragraph(
   "SIRD: severe insulin-resistant diabetes; SIDD: severe insulin-deficient diabetes; MARD: mild age-related diabetes; MOD: mild obesity-related diabetes. Clusters based on Bello-Chavolla et al. (2020) & Ahlqvist et al. (2018)"))
 
 read_docx() %>%
  body_add_flextable(value = tab_clusters, split = TRUE) %>%
  body_end_section_landscape() %>%
  print(target = "Tablas/Table 1 - Characteristics of adult-onset diabetes subgroups in MCPS.docx")
 
 rm(tab_dm_all, tab_mard, tab_mod, tab_sidd, tab_sird, var_labels, vars, tab_clusters)
 
 #### Table 2: Fatal Cardiovascular Events among diabetes subgroups in the Mexico City Prospective Study Cohort ####
 
 fatal_labels <- list(fatal_narrow_total= "Total cardiovascular events",
                      fatal_narrow_cardiac = "Ischemic cardiac events",
                      fatal_narrow_cerebrovascular = "Cerebrovascular events")
 
 fatal_vars <- c("fatal_narrow_total", "fatal_narrow_cardiac", "fatal_narrow_cerebrovascular")
 
 tab_mcps_fatal_nodiab  <- make_fatal_summary(filter(mcps_nodm, cluster_nolab == "No diabetes"), 
                                              header(filter(mcps_nodm, cluster_nolab == "No diabetes"), mcps_nodm))
 tab_mcps_fatal_predm  <- make_fatal_summary(filter(mcps_nodm, cluster_nolab == "Prediabetes"), 
                                             header(filter(mcps_nodm, cluster_nolab == "Prediabetes"), mcps_nodm))
 tab_mcps_fatal_all  <- make_fatal_summary(mcps_dm_fin, header(mcps_dm_fin, mcps_dm_fin))
 tab_mcps_fatal_sidd <- make_fatal_summary(filter(mcps_dm_fin, cluster_nolab == "SIDD"),
                                           header(filter(mcps_dm_fin, cluster_nolab == "SIDD"), mcps_dm_fin))
 tab_mcps_fatal_mod  <- make_fatal_summary(filter(mcps_dm_fin, cluster_nolab == "MOD"),
                                           header(filter(mcps_dm_fin, cluster_nolab == "MOD"), mcps_dm_fin))
 tab_mcps_fatal_mard <- make_fatal_summary(filter(mcps_dm_fin, cluster_nolab == "MARD"),
                                           header(filter(mcps_dm_fin, cluster_nolab == "MARD"), mcps_dm_fin))
 tab_mcps_fatal_sird <- make_fatal_summary(filter(mcps_dm_fin, cluster_nolab == "SIRD"),
                                           header(filter(mcps_dm_fin, cluster_nolab == "SIRD"), mcps_dm_fin))
 
 tab_fatal <- tbl_merge(
  list(tab_mcps_fatal_nodiab,
       tab_mcps_fatal_predm,
       tab_mcps_fatal_all,
       tab_mcps_fatal_sidd,
       tab_mcps_fatal_mod,
       tab_mcps_fatal_mard,
       tab_mcps_fatal_sird),
  tab_spanner = c("No diabetes","Prediabetes","All","SIDD","MOD","MARD","SIRD")
 ) %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  as_flex_table() %>%
  add_header_lines("Table 2. Fatal cardiovascular events among diabetes subgroups in the Mexico City Prospective Study Cohort") %>%
  bold(part = "header") %>%
  theme_vanilla() %>%
  footnote(
   i = 1,
   j = 1,
   value = as_paragraph(c("ICD-10: I20-I25,I60-I69.")),
   ref_symbols = c("1"))
 
 read_docx() %>%
  body_add_flextable(value = tab_fatal, split = TRUE) %>%
  body_end_section_landscape() %>%
  print(target = "Tablas/Table 2 - Fatal cardiovascular events among diabetes subgroups in the Mexico City Prospective Study Cohort_2.docx")
 
 rm(tab_fatal,tab_mcps_fatal_all,tab_mcps_fatal_sidd,tab_mcps_fatal_mod,tab_mcps_fatal_mard,tab_mcps_fatal_sird,fatal_labels,fatal_vars)
 
 #### Table S3: Fatal Cardiovascular Events in Mexico City Prospective Study Cohort and Diabetes Subpopulation ####
 prepare_fatal <- function(df) {
  df %>%
   mutate(
    DIVISION1 = NA, DIVISION2 = NA, DIVISION3 = NA, DIVISION4 = NA,
    HYPERTENSIVE_DISEASE = as.integer(substr(ICD10_UNDERLYING, 1, 3) %in% c("I10","I11","I12","I13","I14","I15")),
    CARDIAC_ARREST       = as.integer(substr(ICD10_UNDERLYING, 1, 3) %in% c("I46")),
    ARRYTHMIAS           = as.integer(substr(ICD10_UNDERLYING, 1, 3) %in% c("I47","I48","I49")),
    HEART_FAILURE        = as.integer(substr(ICD10_UNDERLYING, 1, 3) %in% c("I50")),
    OTHER_CARDIAC        = as.integer(substr(ICD10_UNDERLYING, 1, 3) %in% c("I51","I52")),
    PAD = as.integer(substr(ICD10_UNDERLYING, 1, 3) %in% 
                      c("E105","E115","E145","I170","I171","I72","I73","I74","I75","I76","I77","I78","I79"))) %>%
   set_variable_labels(
    DIVISION1 = "",
    fatal_narrow_total = "Total cardiovascular events",
    fatal_narrow_cardiac = "Ischemic cardiac events",
    fatal_narrow_cerebrovascular = "Cerebrovascular events",
    DIVISION2 = "",
    fatal_wide_total = "Total cardiovascular events",
    fatal_wide_cardiac_ischemic = "Ischemic cardiac events",
    fatal_wide_cerebrovascular_stroke = "Cerebrovascular events",
    DIVISION3 = "",
    ARRYTHMIAS = "Arrythmias",
    HEART_FAILURE = "Heart failure",
    HYPERTENSIVE_DISEASE = "Hypertensive disease",
    OTHER_CARDIAC = "Other",
    DIVISION4 = "",
    PAD = "Peripheral arterial disease"
   )
 }
 
 fatal_summary <- function(df, header_text) {
  # Define the variables of interest
  vars <- c("fatal_narrow_total","fatal_narrow_cerebrovascular","fatal_narrow_cardiac",
            "fatal_wide_total","fatal_wide_cerebrovascular_stroke","fatal_wide_cardiac_ischemic",
            "ARRYTHMIAS","HEART_FAILURE","HYPERTENSIVE_DISEASE","OTHER_CARDIAC","PAD")
  
  # Determine type per variable: dichotomous if all 0/1 (ignoring NAs), else continuous
  type_list <- lapply(vars, function(v) {
   if(all(na.omit(df[[v]]) %in% c(0,1))) "dichotomous" else "continuous"
  })
  names(type_list) <- vars
  
  # Apply tbl_summary with variable-specific types
  df %>%
   dplyr::select(DIVISION1, all_of(vars), DIVISION2, DIVISION3, DIVISION4) %>%
   tbl_summary(
    type = type_list,
    missing = "no",
    statistic = list(
     all_categorical() ~ "{n} ({p}%)",
     all_continuous()  ~ "{mean} ({sd})"
    )
   ) %>%
   modify_header(label = "**Variable**", stat_0 = header_text) %>%
   modify_table_body(~ .x %>% mutate(stat_0 = if_else(label == "", "", stat_0)))
 }
 
 tab_mcps_fatal_all      <- prepare_fatal(mcps)        %>% fatal_summary(header(mcps, mcps))
 tab_mcps_fatal_men      <- prepare_fatal(filter(mcps, MALE == 1)) %>% fatal_summary(header(filter(mcps, MALE==1), mcps))
 tab_mcps_fatal_women    <- prepare_fatal(filter(mcps, MALE == 0)) %>% fatal_summary(header(filter(mcps, MALE==0), mcps))
 tab_mcps_fatal_all_dm   <- prepare_fatal(mcps_dm_fin) %>% fatal_summary(header(mcps_dm_fin, mcps_dm_fin))
 tab_mcps_fatal_men_dm   <- prepare_fatal(filter(mcps_dm_fin, MALE == 1)) %>% fatal_summary(header(filter(mcps_dm_fin, MALE==1), mcps_dm_fin))
 tab_mcps_fatal_women_dm <- prepare_fatal(filter(mcps_dm_fin, MALE == 0)) %>% fatal_summary(header(filter(mcps_dm_fin, MALE==0), mcps_dm_fin))
 
 tab_fatal <- tbl_merge(
  list(tab_mcps_fatal_all, tab_mcps_fatal_men, tab_mcps_fatal_women,
       tab_mcps_fatal_all_dm, tab_mcps_fatal_men_dm, tab_mcps_fatal_women_dm),
  tab_spanner = c("All","Male","Female","All","Male","Female")) %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  as_flex_table() %>%
  merge_at(i = 1, j = 1:7) %>%
  mk_par(i = 1, j = 1, value = as_paragraph(as_b("Narrow definition")),part = "body") %>%
  merge_at(i = 5, j = 1:7) %>%
  mk_par(i = 5, j = 1, value = as_paragraph(as_b("Broad definition")),part = "body") %>%
  merge_at(i = 9, j = 1:7) %>%
  mk_par(i = 9, j = 1, value = as_paragraph(as_i("Other cardiac events")),part = "body") %>%
  merge_at(i = 14, j = 1:7) %>%
  mk_par(i = 14, j = 1, value = as_paragraph(as_i("Other")),part = "body") %>%
  add_header_row(top = TRUE,
                 values = c("",
                            paste0("MCPS (n = ",format(nrow(mcps), big.mark = ","),")"),
                            paste0("MCPS with Diabetes (n = ",format(nrow(mcps_dm_fin), big.mark = ","),")")),
                 colwidths = c(1,3,3)) %>%
  add_header_lines("Table S2. Fatal Cardiovascular Events in Mexico City Prospective Study Cohort and Diabetes Subpopulation") %>%
  bold(part = "header") %>%
  theme_vanilla() %>%
  align(j = 2:4, align = "center", part = "header") %>%
  align(j = 5:7, align = "center", part = "header") %>%
  align(i = c(1,5), j = 1, align = "center", part = "body") %>%
  footnote(
   i = c(1,5), j = 1,
   value = as_paragraph(c("ICD-10: I20-I25,I60-I69.",
                          "ICD-10: I20-I25,I10-I15,I46-I52,I60-I69,I170-I179, R96,E105,E115,E145.")),
   ref_symbols = c("1","2"))
 
 read_docx() %>%
  body_add_flextable(value = tab_fatal, split = TRUE) %>%
  body_end_section_landscape() %>%
  print(target = "Tablas/Table S3 - Fatal Cardiovascular Events in Mexico City Prospective Study Cohort and Diabetes Subpopulation_2.docx")
 
 rm(tab_fatal,tab_mcps_fatal_all, tab_mcps_fatal_men, tab_mcps_fatal_women,
    tab_mcps_fatal_all_dm, tab_mcps_fatal_men_dm, tab_mcps_fatal_women_dm)
 
 #### Table S4: Fatal Cardiovascular Events across diabetes subgroups ####
 # Define clusters
 clusters <- unique(mcps_dm_fin$cluster_nolab)
 cluster_tabs <- list()
 fatal_summary_safe <- function(df, header_text) {
  # Variables to summarize
  vars <- c("fatal_narrow_total","fatal_narrow_cerebrovascular","fatal_narrow_cardiac",
            "fatal_wide_total","fatal_wide_cerebrovascular_stroke","fatal_wide_cardiac_ischemic",
            "ARRYTHMIAS","HEART_FAILURE","HYPERTENSIVE_DISEASE","OTHER_CARDIAC","PAD")
  
  # Coerce each variable to factor if it has >1 unique non-NA value
  df <- df %>% mutate(across(all_of(vars), ~ {
   vals <- na.omit(.x)
   if(length(unique(vals)) > 1) factor(.x, levels = c(0,1)) else .x
  }))
  
  df %>%
   dplyr::select(DIVISION1, fatal_narrow_total, fatal_narrow_cerebrovascular, fatal_narrow_cardiac,
                 DIVISION2, fatal_wide_total, fatal_wide_cerebrovascular_stroke, fatal_wide_cardiac_ischemic,
                 DIVISION3, ARRYTHMIAS, HEART_FAILURE, HYPERTENSIVE_DISEASE, OTHER_CARDIAC,
                 DIVISION4, PAD) %>%
   tbl_summary(
    type = all_categorical() ~ "categorical",
    missing = "no",
    statistic = list(all_categorical() ~ "{n} ({p}%)",
                     all_continuous()  ~ "{mean} ({sd})")
   ) %>%
   modify_header(label = "**Variable**", stat_0 = header_text) %>%
   modify_table_body(~ .x %>% mutate(stat_0 = if_else(label == "", "", stat_0)))
 }
 
 # Loop over clusters and generate summaries
 for(clust in clusters) {
  df_cluster <- mcps_dm_fin %>% filter(cluster_nolab == clust) %>% prepare_fatal()
  header_text <- paste0("n = ", format(nrow(df_cluster), big.mark = ","), 
                        " (", round(nrow(df_cluster)/nrow(mcps_dm_fin)*100,1), "%)")
  cluster_tabs[[clust]] <- fatal_summary_safe(df_cluster, header_text)
 }
 
 # Optional: add the overall diabetes table as well
 df_all_dm <- mcps_dm_fin %>% prepare_fatal()
 header_all <- paste0("n = ", format(nrow(df_all_dm), big.mark = ","))
 cluster_tabs <- c(list("All" = fatal_summary(df_all_dm, header_all)), cluster_tabs)
 
 # Merge tables
 tab_fatal_clusters <- tbl_merge(cluster_tabs,
                                 tab_spanner = names(cluster_tabs)) %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  as_flex_table() %>%
  add_header_lines("Table: Fatal cardiovascular events by cluster in MCPS diabetes cohort") %>%
  bold(part = "header") %>%
  theme_vanilla()
 
 # Export to Word
 read_docx() %>%
  body_add_flextable(value = tab_fatal_clusters, split = TRUE) %>%
  body_end_section_landscape() %>%
  print(target = "Tablas/Table S4 - Fatal_Events_by_Cluster.docx")
 
 rm(tab_list, tab_fatal_cluster, clusters)
 
 #### Table 3: Lexis Expansion to estimate cluster-associated risk ####
 var_list<-list(
  cluster_nolab ~ "Diabetes subgroup",
  MALE ~ "Male",
  tiempo_evol ~ "Years since diagnosis",
  diab_TX ~ "Diabetes treatment",
  HT_TX ~ "Antihypertensive treatment",
  DRUG_E1 ~ "Statin treatment",
  comorbidity_count ~ "Number of comorbidities",
  EDUGP ~ "Educational attainment",
  COYOACAN ~ "Coyoacan",
  SMOKER ~ "Smoking status",
  IDS ~ "Social Development Index",
  PHYSGP ~ "Exercise",
  Total_C ~ "Total cholesterol",
  Total_TG ~ "Triglycerides")
 
 lexis_data <- mcps_dm_fin %>% filter(!is.na(YEAR_RECRUITED)) %>%mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
 lexis_data$age_risk <- lexis_data$AGE+lexis_data$PERSON_YEARS
 mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                     exit = list("period" = date_recruited + PERSON_YEARS), exit.status = fatal_narrow_total, data = lexis_data)
 mcps_fin2 <- splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
 mcps_fin2$ageout <- mcps_fin2$age + mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
 mcps_fin2$cvd_death<-ifelse(mcps_fin2$lex.Xst==T,1,0)
 mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
 mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
 mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
 
 lexis_full <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                      cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                     data=mcps_fin2,x = TRUE)
 
 table_lexiscox_overall <- tbl_regression(lexis_full, 
                                          exponentiate = TRUE,
                                          label = var_list,
                                          conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**")
 
 lexis_full_dx <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                         cluster_nolab*DM_noDX + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                        data=mcps_fin2,x = TRUE)
 
 table_lexiscox_overall_dx <- tbl_regression(lexis_full_dx, 
                                             exponentiate = TRUE,
                                             label = var_list,
                                             conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") %>% as_flex_table()
 
 read_docx() %>%
   body_add_flextable(value = table_lexiscox_overall_dx, split = TRUE) %>%
   body_end_section_landscape() %>% 
   print(target = "Tablas/Table X - Hazard Ratios for the Lexis-expanded Cox model with Undiagnosed diabetes interaction.docx")
 
 lexis_full_undx <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                           cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                          data=mcps_fin2 %>% filter(DM_noDX==1),x = TRUE)
 
 table_lexiscox_overall_undx <- tbl_regression(lexis_full_undx, 
                                               exponentiate = TRUE,
                                               label = var_list,
                                               conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                     exit = list("period" = date_recruited + PERSON_YEARS), exit.status = fatal_narrow_cardiac, data = lexis_data)
 mcps_fin2 <- splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
 mcps_fin2$ageout <- mcps_fin2$age + mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
 mcps_fin2$cvd_death<-ifelse(mcps_fin2$lex.Xst==T,1,0)
 mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
 mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
 mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
 
 lexis_full_cardiac <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                              cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 + comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                             data=mcps_fin2,x = TRUE)
 
 table_lexiscox_cardiac <- tbl_regression(lexis_full_cardiac, 
                                          exponentiate = TRUE,
                                          label = var_list,
                                          conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 lexis_full_cardiac_dx <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                                 cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                                data=mcps_fin2 %>% filter(DM_noDX==0),x = TRUE)
 
 table_lexiscox_cardiac_dx <- tbl_regression(lexis_full_cardiac_dx, 
                                             exponentiate = TRUE,
                                             label = var_list,
                                             conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 lexis_full_cardiac_undx <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                                   cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                                  data=mcps_fin2 %>% filter(DM_noDX==1),x = TRUE)
 
 table_lexiscox_cardiac_undx <- tbl_regression(lexis_full_cardiac_undx, 
                                               exponentiate = TRUE,
                                               label = var_list,
                                               conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 
 mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                     exit = list("period" = date_recruited + PERSON_YEARS), exit.status = fatal_narrow_cerebrovascular, data = lexis_data)
 mcps_fin2 <- splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
 mcps_fin2$ageout <- mcps_fin2$age + mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
 mcps_fin2$cvd_death<-ifelse(mcps_fin2$lex.Xst==T,1,0)
 mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
 mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
 mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
 
 lexis_full_cerebrovascular<- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                                     cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 + comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                                    data=mcps_fin2,x = TRUE)
 
 table_lexiscox_cerebrovascular <- tbl_regression(lexis_full_cerebrovascular, 
                                                  exponentiate = TRUE,
                                                  label = var_list,
                                                  conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 lexis_full_cerebrovascular_dx <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                                         cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                                        data=mcps_fin2 %>% filter(DM_noDX==0),x = TRUE)
 
 table_lexiscox_cerebrovascular_dx <- tbl_regression(lexis_full_cerebrovascular_dx, 
                                                     exponentiate = TRUE,
                                                     label = var_list,
                                                     conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 lexis_full_cerebrovascular_undx <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                                           cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                                          data=mcps_fin2 %>% filter(DM_noDX==1),x = TRUE)
 
 table_lexiscox_cerebrovascular_undx <- tbl_regression(lexis_full_cerebrovascular_undx, 
                                                       exponentiate = TRUE,
                                                       label = var_list,
                                                       conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 
 tab1<-tbl_merge(tbls = list(table_lexiscox_overall,table_lexiscox_cardiac,table_lexiscox_cerebrovascular),
                 tab_spanner = c("Cardiovascular death", "Ischemic cardiac death", "Cerebrovascular death")) %>%
  as_flex_table()
 
 
 read_docx() %>%
  body_add_flextable(value = tab1, split = TRUE) %>%
  body_end_section_landscape() %>% 
  print(target = "Tablas/Table 3 - Hazard Ratios for the Lexis-expanded Cox model.docx")
 
 tab1<-tbl_merge(tbls = list(table_lexiscox_overall_dx,table_lexiscox_cardiac_dx,table_lexiscox_cerebrovascular_dx),
                 tab_spanner = c("Cardiovascular death", "Ischemic cardiac death", "Cerebrovascular death")) %>%
  as_flex_table()
 
 
 read_docx() %>%
  body_add_flextable(value = tab1, split = TRUE) %>%
  body_end_section_landscape() %>% 
  print(target = "Tablas/Table 3_dx - Hazard Ratios for the Lexis-expanded Cox model.docx")
 
 tab1<-tbl_merge(tbls = list(table_lexiscox_overall_undx,table_lexiscox_cardiac_undx,table_lexiscox_cerebrovascular_undx),
                 tab_spanner = c("Cardiovascular death", "Ischemic cardiac death", "Cerebrovascular death")) %>%
  as_flex_table()
 
 
 read_docx() %>%
  body_add_flextable(value = tab1, split = TRUE) %>%
  body_end_section_landscape() %>% 
  print(target = "Tablas/Table 3_undx - Hazard Ratios for the Lexis-expanded Cox model.docx")
 
 rm(mcps_fin2, mcps_lexis)
 
 ### No diabetes ###
 
 var_list<-list(
  cluster_nolab ~ "Diabetes subgroup",
  MALE ~ "Male",
  tiempo_evol ~ "Years since diagnosis",
  diab_TX ~ "Diabetes treatment",
  HT_TX ~ "Antihypertensive treatment",
  DRUG_E1 ~ "Statin treatment",
  comorbidity_count ~ "Number of comorbidities",
  EDUGP ~ "Educational attainment",
  COYOACAN ~ "Coyoacan",
  SMOKER ~ "Smoking status",
  IDS ~ "Social Development Index",
  PHYSGP ~ "Exercise",
  Total_C ~ "Total cholesterol",
  Total_TG ~ "Triglycerides")
 
 mcps_nodm$cluster_nolab<-"No diabetes"
 mcps_nodm$tiempo_evol<-c(-1)
 mcps_nodm$predm<-between(mcps_nodm$BASE_HBA1C, left = 5.7, right = 6.4)
 mcps_nodm$cluster_nolab[mcps_nodm$predm==TRUE]<-"Prediabetes"
 mcps_fin<-mcps_dm_fin %>%
  select(cluster_nolab, BASE_HBA1C, MALE,tiempo_evol,diab_TX, HT_TX, DRUG_E1, COYOACAN, comorbidity_count, EDUGP, SMOKER, IDS, PHYSGP, Total_C, Total_TG,YEAR_RECRUITED, MONTH_RECRUITED,
         AGE, PERSON_YEARS, fatal_narrow_total,fatal_narrow_cardiac, fatal_narrow_cerebrovascular) %>%
  rbind(mcps_nodm %>%
         select(cluster_nolab, BASE_HBA1C,MALE,tiempo_evol,diab_TX, HT_TX, DRUG_E1, COYOACAN, comorbidity_count, EDUGP, SMOKER, IDS, PHYSGP, Total_C, Total_TG,YEAR_RECRUITED, MONTH_RECRUITED,
                AGE, PERSON_YEARS, fatal_narrow_total,fatal_narrow_cardiac, fatal_narrow_cerebrovascular)) %>%
  mutate(cluster_nolab=relevel(cluster_nolab, ref="No diabetes"))
 
 
 lexis_data <- mcps_fin %>% filter(!is.na(YEAR_RECRUITED)) %>%mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
 lexis_data$age_risk <- lexis_data$AGE+lexis_data$PERSON_YEARS
 mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                     exit = list("period" = date_recruited + PERSON_YEARS), exit.status = fatal_narrow_total, data = lexis_data)
 mcps_fin2 <- splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
 mcps_fin2$ageout <- mcps_fin2$age + mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
 mcps_fin2$cvd_death<-ifelse(mcps_fin2$lex.Xst==T,1,0)
 mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period) 
 
 lexis_full <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                      cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 + comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat), 
                     data=mcps_fin2,x = TRUE)
 
 table_lexiscox_overall <- tbl_regression(lexis_full, 
                                          exponentiate = TRUE,
                                          label = var_list,
                                          conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 
 coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
        cluster_nolab + strata(age.cat), data=mcps_fin2,x = TRUE) %>%tbl_regression(exponentiate = TRUE,conf.level = 0.95) %>% bold_p(t = 0.05) 
 
 
 mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                     exit = list("period" = date_recruited + PERSON_YEARS), exit.status = fatal_narrow_cardiac, data = lexis_data)
 mcps_fin2 <- splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
 mcps_fin2$ageout <- mcps_fin2$age + mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
 mcps_fin2$cvd_death<-ifelse(mcps_fin2$lex.Xst==T,1,0)
 mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
 mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
 mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
 
 lexis_full_cardiac <- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                              cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 + comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                             data=mcps_fin2,x = TRUE)
 
 table_lexiscox_cardiac <- tbl_regression(lexis_full_cardiac, 
                                          exponentiate = TRUE,
                                          label = var_list,
                                          conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                     exit = list("period" = date_recruited + PERSON_YEARS), exit.status = fatal_narrow_cerebrovascular, data = lexis_data)
 mcps_fin2 <- splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
 mcps_fin2$ageout <- mcps_fin2$age + mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
 mcps_fin2$cvd_death<-ifelse(mcps_fin2$lex.Xst==T,1,0)
 mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
 mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
 mcps_fin2<- mcps_fin2 %>% rename(
  "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
  mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
 
 lexis_full_cerebrovascular<- coxph(Surv(time_at_entry, time_at_exit, cvd_death) ~
                                     cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 + comorbidity_count+ EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG + strata(age.cat),
                                    data=mcps_fin2,x = TRUE)
 
 table_lexiscox_cerebrovascular <- tbl_regression(lexis_full_cerebrovascular, 
                                                  exponentiate = TRUE,
                                                  label = var_list,
                                                  conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 
 tab1<-tbl_merge(tbls = list(table_lexiscox_overall,table_lexiscox_cardiac,table_lexiscox_cerebrovascular),
                 tab_spanner = c("Cardiovascular death", "Ischemic cardiac death", "Cerebrovascular death")) %>%
  as_flex_table()
 
 
 read_docx() %>%
  body_add_flextable(value = tab1, split = TRUE) %>%
  body_end_section_landscape() %>% 
  print(target = "Tablas/Table 3_1 - Hazard Ratios for the Lexis-expanded Cox model.docx")
 
 rm(mcps_fin2, mcps_lexis)
 
 #### Table 4 SCORE2-Diabetes performance ####
 
 table_score2<-c_stat_total %>% left_join(mean_calib_cvd_clusters, by=c("cluster_nolab", "recalib")) %>%
  left_join(weak_calib, by=c("cluster_nolab", "recalib")) %>% left_join(ov_tot_cvd, , by=c("cluster_nolab", "recalib")) %>%
  rename(Subgroup=cluster_nolab, "c-statistic (95%CI)"=`95% CI.x`, "Observed risk"=obs_rsk, 
         "Expected risk"=score, "O/E ratio (95%CI)"=`95% CI.y`, "Calibration slope (95%CI)"=`95% CI`,
         "Brier score"=estimate, "IPA (%)"=IPA) %>%
  flextable()
 
 read_docx() %>%
  body_add_flextable(value = table_score2, split = TRUE) %>%
  body_end_section_landscape() %>%
  print(target = "Tablas/Table - SCORE2-Diabetes performance.docx")
 
 #### Figure 1: SCORE2-Diabetes Risk Categories ####
 
 risk2<-c("Low (<5%)", "Moderate (5% to <10%)", "High (10% to <20%)", "Very high (≥20%)")
 
 overall_score2<-mcps_dm_fin %>%
  filter(!is.na(score2_dm_risk)) %>%
  group_by(SCORE2_DM_class) %>% summarise(count=n(), .groups = 'drop') %>%
  mutate(SCORE2_DM_class=factor(SCORE2_DM_class, levels=c("Low risk", "Moderate risk", "High risk", "Very high risk"))) %>% 
  mutate(SCORE2_DM_class=factor(SCORE2_DM_class, labels=risk2)) %>%
  mutate(SCORE2_DM_class=factor(SCORE2_DM_class, levels=risk3)) %>%
  mutate(n = sum(count),
         per = count / n,cluster_nolab="Overall") %>% ungroup() %>%
  select(cluster_nolab, SCORE2_DM_class, count, per) 
 
 fig1 <- rbind(overall_score2, mcps_dm_fin %>%
                filter(!is.na(SCORE2_DM)) %>%
                group_by(cluster_nolab, SCORE2_DM_class) %>%
                summarise(count=n(), .groups = 'drop') %>%
                mutate(SCORE2_DM_class=factor(SCORE2_DM_class, levels=c("Low risk", "Moderate risk", "High risk", "Very high risk")),
                       SCORE2_DM_class=factor(SCORE2_DM_class, labels=risk2),
                       SCORE2_DM_class=factor(SCORE2_DM_class, levels=risk3),
                       cluster_nolab=factor(cluster_nolab, levels=c("MOD","SIRD","SIDD","MARD"))) %>%
                group_by(cluster_nolab) %>%
                mutate(n = sum(count), per = count / n) %>%
                ungroup() %>%
                select(cluster_nolab, SCORE2_DM_class, count, per)) %>%
  mutate(type=c(rep("Overall",4), rep("Diabetes subgroups",16))) %>%
  mutate(cluster_nolab = factor(as.character(cluster_nolab), levels = c("Overall","MOD","SIRD","SIDD","MARD")),
         type = factor(type, levels = c("Overall", "Diabetes subgroups"))) %>%
  ggplot(aes(x=cluster_nolab, y=per, fill=SCORE2_DM_class)) +
  geom_bar(stat="identity", color="black", width=0.5) +
  geom_text(aes(label = paste0(round(per*100,1), "%"), angle=90),
            size=2.7, col="white", position = position_stack(vjust=0.5),
            fontface="bold.italic") +
  labs(x="", y="Distribution (%)", fill="") +
  scale_fill_manual(values=c("darkred", "red", "#F0BD27", "#51B364"),
                    guide=guide_legend(reverse=TRUE)) +
  scale_y_continuous(labels=scales::percent) +
  coord_flip() +
  facet_grid(rows = vars(type),scales = "free_y",space = "free_y",switch = "y") +
  ggpubr::theme_pubclean() +
  theme(legend.position="bottom",
        plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0),
        strip.placement = "outside",
        axis.title.y = element_text(margin = margin(r = 15)),
        strip.text.y.left = element_text(angle=90, face="bold"))
 
 
 
 ggsave(fig1,filename = "Figuras/Figure1 - SCORE2-Diabetes risk categories among diabetes subgroups in MCPS.jpg", 
        width = 25, height = 15,units=c("cm"),
        dpi = 300,limitsize = FALSE) 
 
 overall_score2<-mcps_dm_fin %>%
  filter(!is.na(score2diab_recalib)) %>%
  mutate(SCORE2_DM_class_recalib=cut(score2diab_recalib, breaks=c(-Inf,0.05, 0.1, 0.20, Inf))) %>%
  group_by(SCORE2_DM_class_recalib) %>% summarise(count=n(), .groups = 'drop') %>%
  mutate(SCORE2_DM_class_recalib=factor(SCORE2_DM_class_recalib, labels =c("Low risk", "Moderate risk", "High risk", "Very high risk"))) %>% 
  mutate(SCORE2_DM_class_recalib=factor(SCORE2_DM_class_recalib, labels=risk2)) %>%
  mutate(SCORE2_DM_class_recalib=factor(SCORE2_DM_class_recalib, levels=risk3)) %>%
  mutate(n = sum(count),
         per = count / n,cluster_nolab="Overall") %>% ungroup() %>%
  select(cluster_nolab, SCORE2_DM_class_recalib, count, per) 
 
 fig1a <- rbind(overall_score2, mcps_dm_fin %>%
                 filter(!is.na(score2diab_recalib)) %>%
                 mutate(SCORE2_DM_class_recalib=cut(score2diab_recalib, breaks=c(-Inf,0.05, 0.1, 0.20, Inf))) %>%
                 
                 group_by(cluster_nolab, SCORE2_DM_class_recalib) %>%
                 summarise(count=n(), .groups = 'drop') %>%
                 mutate(SCORE2_DM_class_recalib=factor(SCORE2_DM_class_recalib, labels=c("Low risk", "Moderate risk", "High risk", "Very high risk")),
                        SCORE2_DM_class_recalib=factor(SCORE2_DM_class_recalib, labels=risk2),
                        SCORE2_DM_class_recalib=factor(SCORE2_DM_class_recalib, levels=risk3),
                        cluster_nolab=factor(cluster_nolab, levels=c("MOD","SIRD","SIDD","MARD"))) %>%
                 group_by(cluster_nolab) %>%
                 mutate(n = sum(count), per = count / n) %>%
                 ungroup() %>%
                 select(cluster_nolab, SCORE2_DM_class_recalib, count, per)) %>%
  mutate(type=c(rep("Overall",4), rep("Diabetes subgroups",16))) %>%
  mutate(cluster_nolab = factor(as.character(cluster_nolab), levels = c("Overall","MOD","SIRD","SIDD","MARD")),
         type = factor(type, levels = c("Overall", "Diabetes subgroups"))) %>%
  ggplot(aes(x=cluster_nolab, y=per, fill=SCORE2_DM_class_recalib)) +
  geom_bar(stat="identity", color="black", width=0.5) +
  geom_text(aes(label = paste0(round(per*100,1), "%"), angle=90),
            size=2.7, col="white", position = position_stack(vjust=0.5),
            fontface="bold.italic") +
  labs(x="", y="Distribution (%)", fill="") +
  scale_fill_manual(values=c("darkred", "red", "#F0BD27", "#51B364"),
                    guide=guide_legend(reverse=TRUE)) +
  scale_y_continuous(labels=scales::percent) +
  coord_flip() +
  facet_grid(rows = vars(type),scales = "free_y",space = "free_y",switch = "y") +
  ggpubr::theme_pubclean() +
  theme(legend.position="bottom",
        plot.title = element_text(size=15, face="bold", hjust=0.5, vjust=0),
        strip.placement = "outside",
        axis.title.y = element_text(margin = margin(r = 15)),
        strip.text.y.left = element_text(angle=90, face="bold"))
 
 
 
 ggsave(fig1a,filename = "Figuras/Figure1 - SCORE2-Diabetes risk categories among diabetes subgroups in MCPS_recalib.jpg", 
        width = 25, height = 15,units=c("cm"),
        dpi = 300,limitsize = FALSE) 
 
 #### Figure 2: IPW adjusted Kaplan-Meier ####
 
 treatment_model <- nnet::multinom(cluster_nolab ~ MALE + AGE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG,
                                   data = mcps_dm_fin)
 
 s_iptw <- adjustedsurv(data=mcps_dm_fin,
                        variable = "cluster_nolab",
                        ev_time = "PERSON_YEARS",
                        event = "fatal_narrow_total",
                        method = "iptw_km",
                        treatment_model = treatment_model,
                        weight_method = "glm",
                        conf_int = TRUE,
                        stabilize = TRUE)
 
 fig2a <- plot(s_iptw, risk_table = TRUE, risk_table_stratify_color = FALSE,
               risk_table_stratify = TRUE, risk_table_digits = 0, x_n_breaks = 9,
               risk_table_title_size= 11, cif=TRUE,x_breaks=c(0,5,10,15,20,24),
               gg_theme= theme_bw(), risk_table_theme = theme_pubclean() +
                theme(axis.text.y = element_text(color= c("#C77CFF","#00BFC4","#7CAE00","#F8766D"))),
               legend.position="top",
               xlab="Time in Years",additional_layers = list(scale_y_continuous(labels = scales::percent))) 
 
 fig2a$plot <- fig2a$plot +
  coord_cartesian(xlim = c(0, 24))
 
 s_iptw2 <- adjustedsurv(data=mcps_dm_fin,
                         variable = "cluster_nolab",
                         ev_time = "PERSON_YEARS",
                         event = "fatal_narrow_cardiac",
                         method = "iptw_km",
                         treatment_model = treatment_model,
                         weight_method = "glm",
                         conf_int = TRUE,
                         stabilize = TRUE)
 
 
 fig2b <- plot(s_iptw2, risk_table = TRUE, risk_table_stratify_color = FALSE,
               risk_table_stratify = TRUE, risk_table_digits = 0, x_n_breaks = 8,
               risk_table_title_size= 11, cif=TRUE,x_breaks=c(0,5,10,15,20,24),
               gg_theme= theme_bw(), risk_table_theme = theme_pubclean() +
                theme(axis.text.y = element_text(color= c("#C77CFF","#00BFC4","#7CAE00","#F8766D"))),
               legend.position="top",
               xlab="Time in Years",additional_layers = list(scale_y_continuous(labels = scales::percent)))
 
 s_iptw3 <- adjustedsurv(data=mcps_dm_fin,
                         variable = "cluster_nolab",
                         ev_time = "PERSON_YEARS",
                         event = "fatal_narrow_cerebrovascular",
                         method = "iptw_km",
                         treatment_model = treatment_model,
                         weight_method = "glm",
                         conf_int = TRUE,
                         stabilize = TRUE)
 
 fig2c <- plot(s_iptw3, risk_table = TRUE, risk_table_stratify_color = FALSE,
               risk_table_stratify = TRUE, risk_table_digits = 0, x_n_breaks = 8,
               risk_table_title_size= 11, cif=TRUE,x_breaks=c(0,5,10,15,20,24),
               gg_theme= theme_bw(), risk_table_theme = theme_pubclean() +
                theme(axis.text.y = element_text(color= c("#C77CFF","#00BFC4","#7CAE00","#F8766D"))),
               legend.position="top",
               xlab="Time in Years",additional_layers = list(scale_y_continuous(labels = scales::percent))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
 
 fig2<-ggarrange(fig2a, fig2b,fig2c, labels = "AUTO", nrow=1, ncol=3, common.legend = T)
 
 ggsave(fig2,filename = "Figuras/Figure2 -  Kaplan-Meier IPW for fatal CVD for diabetes subgroups in MCPS.jpg", 
        width = 40, 
        height = 17,
        units=c("cm"),
        dpi = 300,
        limitsize = FALSE) 
 
#### Figure 3: Cox models comparing clusters, risk factors and SCORE2-Diabetes ####
 ## Models ##
 data_model <- mcps_dm_fin %>% select("PERSON_YEARS","fatal_narrow_total","fatal_narrow_cardiac","fatal_narrow_cerebrovascular","cluster_nolab","BMI","meanSBP","meanDBP","WHR","WHtR","BASE_HBA1C","MALE",
                                      "AGE","diab_TX","comorbidity_count","score2diab_recalib","SCORE2_DM_class","COYOACAN",
                                      "SMOKER","EDUGP","IDS","tiempo_evol","PHYSGP") %>% na.omit() %>%
  mutate(cluster_nolab = relevel(cluster_nolab, ref = "MOD")) %>%
  set_variable_labels(cluster_nolab = "Cluster",
                      BMI = "BMI (kg/m²)",
                      meanSBP = "SBP (mmHg)",
                      meanDBP = "DBP (mmHg)",
                      WHR = "Waist-to-Hip Ratio",
                      WHtR = "Waist-to-Height Ratio",
                      BASE_HBA1C = "HbA1c (%)",
                      COYOACAN = "Municipality",
                      SMOKER = "Smoking",
                      EDUGP = "Educational attainment",
                      MALE = "Sex",
                      tiempo_evol = "Years since diagnosis",
                      diab_TX = "Diabetes treatment",
                      comorbidity_count = "Number of comorbidities",
                      score2diab_recalib = "SCORE2-DM",
                      SCORE2_DM_class = "SCORE2-DM Risk Category",
                      IDS = "Social Development Index",
                      PHYSGP = "Physical activity (times per week)",
                      AGE = "Age (years)") %>% as.data.frame()
 
 fit_cox <- function(formula, data) {
  coxph(formula, data = data)}
 
 fit_models_for_outcome <- function(outcome, data) {
  
  formulas <- list(
   "Risk factors" = as.formula(paste0(
    "Surv(PERSON_YEARS, ", outcome, 
    ") ~ MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP"
   )),
   "Risk factors+\nDiabetes subgroups" = as.formula(paste0(
    "Surv(PERSON_YEARS, ", outcome, 
    ") ~ cluster_nolab + MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP"
   )),
   "Risk factors+\nSCORE2-Diabetes" = as.formula(paste0(
    "Surv(PERSON_YEARS, ", outcome, 
    ") ~ score2diab_recalib + MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP"
   )),
   "Risk factors+\nSCORE2-Diabetes+\nDiabetes subgroups" = as.formula(paste0(
    "Surv(PERSON_YEARS, ", outcome, 
    ") ~ score2diab_recalib + cluster_nolab + MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP"
   ))
  )
  
  models <- map(formulas, ~fit_cox(.x, data))
  return(models)
 }
 
 ### C-statistic and AIC ###
 plot_cstat_aic <- function(models, title=NULL) {
  
  model_colors <- c(
   "Risk factors" = "#374E55FF",
   "Risk factors+\nDiabetes subgroups" = "#DF8F44FF",
   "Risk factors+\nSCORE2-Diabetes" = "#00A1D5FF",
   "Risk factors+\nSCORE2-Diabetes+\nDiabetes subgroups" = "#B24745FF")
  
  cstat_data <- map_dfr(names(models), ~{
   cc <- summary(models[[.x]])$concordance
   tibble(
    class = .x,
    c_stat = cc[1],
    lower = cc[1] - qnorm(0.975) * cc[2],
    upper = cc[1] + qnorm(0.975) * cc[2]
   )
  })
  
  fig_cstat <- ggplot(cstat_data, aes(x=class, y=c_stat, ymin=lower, ymax=upper, fill=class)) +
   geom_crossbar(width=0.3, alpha=0.85, position=position_dodge2()) +
   scale_fill_manual(values=model_colors) +
   scale_y_continuous(breaks=seq(0.65,0.78,by=0.02), limits=c(0.65,0.78)) +
   labs(x="Model", y="C-statistic (95% CI)", fill="Predictor") +
   theme_pubclean() +
   theme(legend.position="top", axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  legend <- get_legend(fig_cstat)
  fig_cstat <- fig_cstat + theme(legend.position="none")
  
  aic_data <- map_dfr(names(models), function(name) {
   tibble(
    class = name,
    AIC_value = AIC(models[[name]]))})
  
  aic_data <- aic_data %>%
   mutate(delta_aic = AIC_value - AIC_value[1]) %>%
   filter(delta_aic != 0)
  
  fig_aic <- ggplot(aic_data, aes(x=class, y=delta_aic, fill=class)) +
   geom_bar(stat="identity", col="black", width=0.5) +
   scale_fill_manual(values=model_colors) +
   labs(x="Model", y="Delta-AIC", fill="Predictor") +
   theme_pubclean() +
   theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  fig <- ggarrange(fig_cstat, fig_aic, nrow=1, ncol=2)
  
  return(list(plot = fig, legend = legend))
 }
 
 outcomes <- c("fatal_narrow_total", "fatal_narrow_cardiac", "fatal_narrow_cerebrovascular")
 models <- fit_models_for_outcome("fatal_narrow_total", data_model)

 cstat_data <- map_dfr(names(models), ~{
  cc <- summary(models[[.x]])$concordance
  tibble(
   class = .x,
   c_stat = cc[1],
   lower = cc[1] - qnorm(0.975) * cc[2],
   upper = cc[1] + qnorm(0.975) * cc[2]
  )
 })
 
 figures <- map(outcomes, function(outcome) {
  models <- fit_models_for_outcome(outcome, data_model)
  plot_cstat_aic(models)
 })
 
 fig3<-ggarrange(figures[[1]]$plot,figures[[2]]$plot,figures[[3]]$plot, labels = "AUTO", ncol=1, nrow=3)
 
 fig3 <- ggarrange(figures[[1]]$legend, fig3, nrow=2, heights=c(0.1, 1))
 
 ggsave("Figuras/Figure 3.jpg",plot=fig3, width=25, height=18, units="cm", dpi=300, limitsize=FALSE)
 
#### p-values for c-statistics ####

mod1<-coxph(Surv(PERSON_YEARS, fatal_narrow_total) ~ MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
mod2<-coxph(Surv(PERSON_YEARS, fatal_narrow_total) ~ cluster_nolab+ MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
mod3<-coxph(Surv(PERSON_YEARS, fatal_narrow_total) ~ score2diab_recalib + MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
mod4<-coxph(Surv(PERSON_YEARS, fatal_narrow_total) ~ cluster_nolab + score2diab_recalib + MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
 
compare_c_stat(mod1, mod2, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod2, mod3, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod2, mod4, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod3, mod4, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod3, mod2, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod3, mod1, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod1, mod4, data = data_model,ci_type = "perc", R = 100)
 
mod1<-coxph(Surv(PERSON_YEARS, fatal_narrow_cerebrovascular) ~ MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
mod2<-coxph(Surv(PERSON_YEARS, fatal_narrow_cerebrovascular) ~ cluster_nolab+ MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
mod3<-coxph(Surv(PERSON_YEARS, fatal_narrow_cerebrovascular) ~ score2diab_recalib + MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
mod4<-coxph(Surv(PERSON_YEARS, fatal_narrow_cerebrovascular) ~ cluster_nolab + score2diab_recalib + MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
 
compare_c_stat(mod1, mod2, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod2, mod3, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod2, mod4, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod3, mod4, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod3, mod2, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod3, mod1, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod1, mod4, data = data_model,ci_type = "perc", R = 100)
 
 
mod1<-coxph(Surv(PERSON_YEARS, fatal_narrow_cardiac) ~ MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
mod2<-coxph(Surv(PERSON_YEARS, fatal_narrow_cardiac) ~ cluster_nolab+ MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
mod3<-coxph(Surv(PERSON_YEARS, fatal_narrow_cardiac) ~ score2diab_recalib + MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
mod4<-coxph(Surv(PERSON_YEARS, fatal_narrow_cardiac) ~ cluster_nolab + score2diab_recalib + MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP, data_model)
 
compare_c_stat(mod1, mod2, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod2, mod3, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod2, mod4, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod3, mod4, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod3, mod2, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod3, mod1, data = data_model,ci_type = "perc", R = 100)
compare_c_stat(mod1, mod4, data = data_model,ci_type = "perc", R = 100)

#### Cluster stability ####
migration_dataset <- mcps_filtered_dm %>% 
                     select(PATID, AGE, AGE_DX, YEAR_RECRUITED, MONTH_RECRUITED, cluster_nolab, YEAR_RESURVEY, MONTH_RESURVEY, YEARS_BASELINE_RESURVEY,
                            R_SBP1,R_SBP2, R_SBP3, R_DBP1,R_DBP2, R_DBP3,
                            R_HEIGHT, R_WEIGHT_TANITA, R_WAISTC, R_HBA1C) %>%
 filter(!is.na(YEAR_RESURVEY)) %>% 
 mutate(R_AGE = round(AGE + YEARS_BASELINE_RESURVEY),
        R_BMI = R_WEIGHT_TANITA / (R_HEIGHT^2),
        R_WHtR = R_WAISTC/R_HEIGHT,
        R_meanSBP = case_when(
         !is.na(R_SBP1) & !is.na(R_SBP2) & !is.na(R_SBP3) ~ (R_SBP1 + R_SBP2 + R_SBP3) / 3,
         !is.na(R_SBP1) & !is.na(R_SBP2) ~ (R_SBP1 + R_SBP2) / 2,
         !is.na(R_SBP1) & !is.na(R_SBP3) ~ (R_SBP1 + R_SBP3) / 2,
         !is.na(R_SBP2) & !is.na(R_SBP3) ~ (R_SBP2 + R_SBP3) / 2,
         !is.na(R_SBP1) ~ R_SBP1,
         !is.na(R_SBP2) ~ R_SBP2,
         !is.na(R_SBP3) ~ R_SBP3,
         TRUE ~ NA_real_),
        R_meanDBP = case_when(
         !is.na(R_DBP1) & !is.na(R_DBP2) & !is.na(R_DBP3) ~ (R_DBP1 + R_DBP2 + R_DBP3) / 3,
         !is.na(R_DBP1) & !is.na(R_DBP2) ~ (R_DBP1 + R_DBP2) / 2,
         !is.na(R_DBP1) & !is.na(R_DBP3) ~ (R_DBP1 + R_DBP3) / 2,
         !is.na(R_DBP2) & !is.na(R_DBP3) ~ (R_DBP2 + R_DBP3) / 2,
         !is.na(R_DBP1) ~ R_DBP1,
         !is.na(R_DBP2) ~ R_DBP2,
         !is.na(R_DBP3) ~ R_DBP3,
         TRUE ~ NA_real_),
        baseline_date = make_date(YEAR_RECRUITED, MONTH_RECRUITED, 15),
        resurvey_date = make_date(YEAR_RESURVEY, MONTH_RESURVEY, 15),
        years_passed = as.numeric(interval(baseline_date, resurvey_date) / years(1))) %>% drop_na()

 setwd("/Users/Jeronimo/Library/CloudStorage/GoogleDrive-jeroperezalonso@gmail.com/.shortcut-targets-by-id/1RG-1Sg0NBmkz3AlzY5XL6XJ0mezx6oNg/Mexico City Cohort Study/Proyectos/Cardiovascular Risk DM/Diabetes Clusters/SNNN algorithm")
 r_model<-keras::load_model_hdf5("Bases/mcps.h5")
 #reticulate::use_virtualenv("r-tensorflow", required = TRUE)
 mcps_rcluster<- migration_dataset %>% 
  select(R_BMI, AGE_DX, R_HBA1C, R_meanSBP, R_meanDBP, R_WHtR)
 r_bas<-read.csv("bas.csv")
 r_m0<-apply(bas[,c(5:10)], 2, mean); r_std0<-apply(bas[,c(5:10)], 2, sd)
 mcps_rscale<-as.data.frame(scale(mcps_rcluster, scale=r_std0, center=r_m0))
 r_pred<-r_model %>% predict(as.matrix(mcps_rscale)) %>% k_argmax()
 migration_dataset$R_cluster<-as.vector(r_pred)
 migration_dataset$R_cluster<-factor(migration_dataset$R_cluster, levels=0:3, labels = c("SIDD","SIRD","MOD","MARD"))
 rm(mcps_rcluster,r_bas,mcps_rscale,r_pred,r_m0,r_model,r_std0)
 setwd("/Users/Jeronimo/Library/CloudStorage/GoogleDrive-jeroperezalonso@gmail.com/.shortcut-targets-by-id/1RG-1Sg0NBmkz3AlzY5XL6XJ0mezx6oNg/Mexico City Cohort Study/Proyectos/Cardiovascular Risk DM/Diabetes Clusters")

#devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)

 migration_plot <- migration_dataset %>%
  select(PATID,cluster_nolab,R_cluster)
 
 sankey_long <- migration_plot %>%
  make_long(cluster_nolab, R_cluster) %>%
 mutate(x = recode(x,"cluster_nolab" = "Baseline", "R_cluster" = "Resurvey"),
        next_x = recode(next_x,"cluster_nolab" = "Baseline", "R_cluster" = "Resurvey"))
 
 sankey_plot <- ggplot(sankey_long,aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = "black") +
  geom_sankey_label(size = 4, color = "white") +
  scale_fill_manual(values = c("SIDD" = "#C77CFF",
                               "SIRD" =  "#00BFC4",
                               "MARD" = "#7CAE00",
                               "MOD"="#F8766D")) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
   x = "",
   y = "Number of participants",
   title = "Figure SX. Diabetes subgroup stability from baseline to resurvey")
 
 baseline_summary <- migration_dataset %>%
  count(cluster_nolab, name = "n_base") %>%
  mutate(
   prop_base = round(100 * n_base / sum(n_base), 1),
   baseline = paste0(n_base, " (", prop_base, "%)")
  ) %>%
  select(diabetes_subgroup = cluster_nolab, baseline)
 
 resurvey_summary <- migration_dataset %>%
  count(R_cluster, name = "n_res", .drop = F) %>%
  mutate(
   prop_res = round(100 * n_res / sum(n_res), 1),
   resurvey = paste0(n_res, " (", prop_res, "%)")
  ) %>%
  select(diabetes_subgroup = R_cluster, resurvey)
 
 cluster_table <- left_join(baseline_summary, resurvey_summary, by = "diabetes_subgroup") %>%
  rename(
   `Diabetes subgroup` = diabetes_subgroup,
   `Baseline n (%)`    = baseline,
   `Resurvey n (%)`    = resurvey
  )
 
 table_grob <- tableGrob(
  cluster_table,
  rows = NULL,
  theme = ttheme_minimal(
   core = list(
    fg_params = list(cex = 0.95, hjust = 0.5, x = 0.5)
   ),
   colhead = list(
    fg_params = list(fontface = "bold", hjust = 0.5, x = 0.5)
   )
  )
 )
 
 
 final_figure <- grid.arrange(
  sankey_plot,
  table_grob,
  heights = c(3.5, 1.3)   # tweak if needed
 )
 
# Cluster stability overall
 stability_overall <- migration_dataset %>%
  filter(cluster_nolab == R_cluster) %>%
  summarise(
   stayed = n(),
   total = nrow(migration_dataset),
   proportion = stayed / total) %>%
  mutate(percent = round(100 * proportion, 1))
 
# Cluster stability per cluster
 cluster_stability <- migration_dataset %>%
  group_by(cluster_nolab) %>%
  summarise(
   stayed = sum(cluster_nolab == R_cluster),
   total = n(),
   proportion = stayed / total
  ) %>%
  mutate(percent = round(100 * proportion, 1))
 
 # Mean time
 mean_years_passed <- mean(migration_dataset$years_passed, na.rm = TRUE)
 sd_years_passed <- sd(migration_dataset$years_passed, na.rm = TRUE)
 time_migration <- paste0(
  round(mean_years_passed, 2),
  " ± ",
  round(sd_years_passed, 2),
  " years"
 )
 
 

#### Table SX. Characteristics of adult-onset diabetes subgroups among participants undiagnosed diabetes in the Mexico City Prospective Study ####

 mcps_dm_fin %>%
   mutate(
     DM_noDX = factor(DM_noDX,
                      levels = c(0,1),
                      labels = c("Diagnosed diabetes", "Undiagnosed diabetes")),
     cluster_nolab = factor(cluster_nolab)) %>%
   count(DM_noDX, cluster_nolab) %>%
   group_by(DM_noDX) %>%
   mutate(
     prop = n / sum(n),
     percent = scales::percent(prop, accuracy = 0.1)) %>%
   ggplot(aes(x = DM_noDX, y = prop, fill = cluster_nolab)) +
   geom_col() +
   geom_text(aes(label = percent),
             position = position_stack(vjust = 0.5),
             size = 4,
             color = "white") +
   scale_y_continuous(labels = scales::percent_format()) +
   labs(x = "",
     y = "Percentage",
     fill = "Diabetes subgroup") +
   theme_pubclean() + coord_flip()+
   theme(legend.position = "bottom")
 
 chisq.test(table(mcps_dm_fin$cluster_nolab, mcps_dm_fin$DM_noDX))
 
 mcps_dm_fin %>%
   
   ggplot(aes(x = cluster_nolab, fill = factor(DM_noDX))) +
   geom_bar(position = "dodge") +
   labs(
     x = "Diabetes subgroup",
     y = "Count",
     fill = "Undiagnosed diabetes",
     title = "Diabetes subgroups by undiagnosed diabetes status"
   ) +
   theme_minimal()
 
 var_labels <- list(
  AGE         = "Age (years)",
  MALE        = "Male (%)",
  SMOKER      = "Smoking (%)",
  BMI         = "Body mass index (kg/m²)",
  meanSBP     = "Systolic blood pressure (mmHg)",
  meanDBP     = "Diastolic blood pressure (mmHg)",
  WHR         = "Waist-To-Hip Ratio",
  WHtR        = "Waist-To-Height Ratio",
  Glucose     = "Glucose (mg/dl)",
  BASE_HBA1C  = "Glycated HbA1c (%)",
  mets_ir     = "METS-IR",
  Creatinine  = "Creatinine (mg/dl)",
  eGFR        = "eGFR (ml/min/1.73m²)",
  Total_C     = "Total cholesterol (mg/dl)",
  HDL_C       = "HDL cholesterol (mg/dl)",
  non_HDL_C   = "Non-HDL cholesterol (mg/dl)",
  Total_TG    = "Total triglycerides (mg/dl)",
  ApoB        = "Apolipoprotein B (mg/dl)",
  DRUG_D1     = "Metformin (%)")
 
 vars <- c("AGE","MALE","SMOKER","BMI","meanSBP","meanDBP","WHR","WHtR",
           "Glucose","BASE_HBA1C","mets_ir","Creatinine","eGFR",
           "Total_C","HDL_C","non_HDL_C","Total_TG","ApoB",
           "DRUG_D1")
 
 mcps_dm_nodx <- mcps_dm_fin %>%
  filter(DM_noDX==1)
 
 tab_dm_nodx <- make_summary(mcps_dm_nodx, header(mcps_dm_nodx, mcps_dm_nodx))
 
 tab_sidd_nodx <- make_summary(filter(mcps_dm_nodx, cluster_nolab == "SIDD"),
                         header(filter(mcps_dm_nodx, cluster_nolab == "SIDD"), mcps_dm_nodx))
 tab_mod_nodx <- make_summary(filter(mcps_dm_nodx, cluster_nolab == "MOD"),
                        header(filter(mcps_dm_nodx, cluster_nolab == "MOD"), mcps_dm_nodx))
 tab_mard_nodx <- make_summary(filter(mcps_dm_nodx, cluster_nolab == "MARD"),
                         header(filter(mcps_dm_nodx, cluster_nolab == "MARD"), mcps_dm_nodx))
 tab_sird_nodx<- make_summary(filter(mcps_dm_nodx, cluster_nolab == "SIRD"),
                         header(filter(mcps_dm_nodx, cluster_nolab == "SIRD"), mcps_dm_nodx))
 
 tab_clusters_nodx <- tbl_merge(
  list(tab_dm_nodx, tab_sidd_nodx, tab_mod_nodx, tab_mard_nodx, tab_sird_nodx),
  tab_spanner = c("All", "SIDD", "MOD", "MARD", "SIRD")
 ) %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  as_flex_table() %>%
  add_header_lines("Table x. Characteristics of adult-onset diabetes subgroups among participants with undiagnosed diabetes in the Mexico City Prospective Study") %>%
  bold(part = "header") %>%
  theme_vanilla() %>%
  add_footer_lines(as_paragraph(
   "SIRD: severe insulin-resistant diabetes; SIDD: severe insulin-deficient diabetes; MARD: mild age-related diabetes; MOD: mild obesity-related diabetes. Clusters based on Bello-Chavolla et al. (2020) & Ahlqvist et al. (2018)"))
 
 read_docx() %>%
  body_add_flextable(value = tab_clusters_nodx, split = TRUE) %>%
  body_end_section_landscape() #%>%
  #print(target = "Tablas/Table SX - Characteristics of adult-onset diabetes subgroups in MCPS.docx")

 
#### Figure S3. Kaplan-Meier curves for fatal cardiovascular disease events across diabetes subgroups in undiagnosed participants in the Mexico City Prospective Study ####
 treatment_model_nodx <- nnet::multinom(cluster_nolab ~ MALE + AGE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 + comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP + Total_C + Total_TG,
                                   data = mcps_dm_nodx)
 
 s_iptw_nodx <- adjustedsurv(data=mcps_dm_nodx,
                        variable = "cluster_nolab",
                        ev_time = "PERSON_YEARS",
                        event = "fatal_narrow_total",
                        method = "iptw_km",
                        treatment_model = treatment_model_nodx,
                        weight_method = "glm",
                        conf_int = TRUE,
                        stabilize = TRUE)
 
 figS3a<- plot(s_iptw_nodx, risk_table = TRUE, risk_table_stratify_color = FALSE,
               risk_table_stratify = TRUE, risk_table_digits = 0, x_n_breaks = 9,
               risk_table_title_size= 11, cif=TRUE,x_breaks=c(0,5,10,15,20,24),
               gg_theme= theme_bw(), risk_table_theme = theme_pubclean() +
                theme(axis.text.y = element_text(color= c("#C77CFF","#00BFC4","#7CAE00","#F8766D"))),
               legend.position="top",
               xlab="Time in Years",additional_layers = list(scale_y_continuous(labels = scales::percent))) 
 
 figS3a$plot <- figS3a$plot +
  coord_cartesian(xlim = c(0, 24))
 
 s_iptw2_nodx <- adjustedsurv(data=mcps_dm_nodx,
                         variable = "cluster_nolab",
                         ev_time = "PERSON_YEARS",
                         event = "fatal_narrow_cardiac",
                         method = "iptw_km",
                         treatment_model = treatment_model_nodx,
                         weight_method = "glm",
                         conf_int = TRUE,
                         stabilize = TRUE)
 
 
 figS3b <- plot(s_iptw2_nodx, risk_table = TRUE, risk_table_stratify_color = FALSE,
               risk_table_stratify = TRUE, risk_table_digits = 0, x_n_breaks = 8,
               risk_table_title_size= 11, cif=TRUE,x_breaks=c(0,5,10,15,20,24),
               gg_theme= theme_bw(), risk_table_theme = theme_pubclean() +
                theme(axis.text.y = element_text(color= c("#C77CFF","#00BFC4","#7CAE00","#F8766D"))),
               legend.position="top",
               xlab="Time in Years",additional_layers = list(scale_y_continuous(labels = scales::percent)))
 
 s_iptw3_nodx <- adjustedsurv(data=mcps_dm_nodx,
                         variable = "cluster_nolab",
                         ev_time = "PERSON_YEARS",
                         event = "fatal_narrow_cerebrovascular",
                         method = "iptw_km",
                         treatment_model = treatment_model_nodx,
                         weight_method = "glm",
                         conf_int = TRUE,
                         stabilize = TRUE)
 
 figS3c <- plot(s_iptw3_nodx, risk_table = TRUE, risk_table_stratify_color = FALSE,
               risk_table_stratify = TRUE, risk_table_digits = 0, x_n_breaks = 8,
               risk_table_title_size= 11, cif=TRUE,x_breaks=c(0,5,10,15,20,24),
               gg_theme= theme_bw(), risk_table_theme = theme_pubclean() +
                theme(axis.text.y = element_text(color= c("#C77CFF","#00BFC4","#7CAE00","#F8766D"))),
               legend.position="top",
               xlab="Time in Years",additional_layers = list(scale_y_continuous(labels = scales::percent))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
 
 figS3<-ggarrange(figS3a, figS3b,figS3c, labels = "AUTO", nrow=1, ncol=3, common.legend = T)
 
 #ggsave(figS3,filename = "Figuras/FigureS3 -  Kaplan-Meier IPW for fatal CVD for diabetes subgroups among undiagnosed particiapnts in MCPS.jpg", 
  #      width = 40, 
   #     height = 17,
    #    units=c("cm"),
     #   dpi = 300,
      #  limitsize = FALSE) 
#### Table: Fatal events by cause among diabetes subgroups in the Mexico City Prospective Study cohort ####
 mcps_mortality <- mcps_dm_fin %>% 
  select(PATID, cluster_nolab, fatal_narrow_total,fatal_narrow_cardiac,fatal_narrow_cerebrovascular,
         D019, 
         D022,
         D040,
         D044,
         D054) %>%
  mutate(all_cause = fatal_narrow_total+D019+D022+D040+D044+D054)
 
 fatal_labels <- list(all_cause = "Total events",
                      fatal_narrow_total= "Cardiovascular events",
                      fatal_narrow_cardiac = "Cardiac events",
                      fatal_narrow_cerebrovascular= "Cerebrovascular events",
                      D019 = "Renal events",
                      D022 = "Hepatobiliary events",
                      D040 = "Neoplastic events",
                      D044 = "Respiratory events",
                      D054  = "Other events")
 
 fatal_vars <- c("all_cause","fatal_narrow_total","fatal_narrow_cardiac","fatal_narrow_cerebrovascular","D019", "D022","D040","D044","D054")
 
 tab_mcps_mort_all  <- make_fatal_summary(mcps_mortality, header(mcps_mortality, mcps_mortality))
 tab_mcps_mort_sidd <- make_fatal_summary(filter(mcps_mortality, cluster_nolab == "SIDD"),
                                           header(filter(mcps_mortality, cluster_nolab == "SIDD"), mcps_mortality))
 tab_mcps_mort_mod  <- make_fatal_summary(filter(mcps_mortality, cluster_nolab == "MOD"),
                                           header(filter(mcps_mortality, cluster_nolab == "MOD"), mcps_mortality))
 tab_mcps_mort_mard <- make_fatal_summary(filter(mcps_mortality, cluster_nolab == "MARD"),
                                           header(filter(mcps_mortality, cluster_nolab == "MARD"), mcps_mortality))
 tab_mcps_mort_sird <- make_fatal_summary(filter(mcps_mortality, cluster_nolab == "SIRD"),
                                           header(filter(mcps_mortality, cluster_nolab == "SIRD"), mcps_mortality))
 
 tab_mortality <- tbl_merge(
  list(tab_mcps_mort_all,
       tab_mcps_mort_sidd,
       tab_mcps_mort_mod,
       tab_mcps_mort_mard,
       tab_mcps_mort_sird),
  tab_spanner = c("All","SIDD","MOD","MARD","SIRD")
 ) %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  as_flex_table() %>%
  add_header_lines("Fatal events by cause among diabetes subgroups in the Mexico City Prospective Study Cohort") %>%
  bold(part = "header") %>%
  theme_vanilla()
 
 rm(tab_mcps_mort_all,tab_mcps_mort_sidd,tab_mcps_mort_mod,tab_mcps_mort_mard,tab_mcps_mort_sird,fatal_labels,fatal_vars)
 
 read_docx() %>%
   body_add_flextable(value = tab_mortality, split = TRUE) %>%
   body_end_section_landscape() %>%
   print(target = "Tablas/Table - Fatal events by cause and subgroups.docx")
 
 
#### Table: Fine & Gray model ####
 var_list<-list(
  cluster_nolab ~ "Diabetes subgroup",
  MALE ~ "Male",
  tiempo_evol ~ "Years since diagnosis",
  diab_TX ~ "Diabetes treatment",
  HT_TX ~ "Antihypertensive treatment",
  DRUG_E1 ~ "Statin treatment",
  comorbidity_count ~ "Number of comorbidities",
  EDUGP ~ "Educational attainment",
  COYOACAN ~ "Coyoacan",
  SMOKER ~ "Smoking status",
  IDS ~ "Social Development Index",
  PHYSGP ~ "Exercise",
  Total_C ~ "Total cholesterol",
  Total_TG ~ "Triglycerides")
 ## Models ##
 data_fg <- mcps_dm_fin %>% 
  mutate(fatal = case_when(
  fatal_narrow_total == 1 ~ 1, 
   D019 == 1 ~ 2,
   D022 == 1 ~ 3,
   D040 == 1 ~ 4,
   D044 == 1 ~ 5,
   D054 == 1 ~ 6,
   TRUE ~ 0),
  fatal_2 = case_when(
    fatal_narrow_cerebrovascular == 1 ~ 1, 
    fatal_narrow_cardiac == 1 ~ 2,
    D019 == 1 ~ 3,
    D022 == 1 ~ 4,
    D040 == 1 ~ 5,
    D044 == 1 ~ 6,
    D054 == 1 ~ 7,
    TRUE ~ 0)) %>%
  select("PATID", "PERSON_YEARS", "cluster_nolab", "fatal","fatal_2","BMI","meanSBP","meanDBP","WHR","WHtR","BASE_HBA1C",
         "MALE","AGE","diab_TX","comorbidity_count","score2diab_recalib","SCORE2_DM_class","COYOACAN",
         "SMOKER","EDUGP","IDS","tiempo_evol","PHYSGP", "HT_TX","DRUG_E1", "Total_C", "Total_TG","YEAR_RECRUITED", "MONTH_RECRUITED") %>% na.omit() %>%
  mutate(cluster_nolab = relevel(cluster_nolab, ref = "MOD"),
         fatal = factor(fatal, levels = c(0,1,2,3,4,5,6),
          labels = c("censor",
                     "cvd_death",
                     "renal_death",
                     "hepatobiliary_death",
                     "neoplastic_death",
                     "respiratory_death",
                     "external_death")),
         fatal_2 = factor(fatal_2, levels = c(0,1,2,3,4,5,6,7),
                        labels = c("censor",
                                   "stroke_death",
                                   "cardiac_death",
                                   "renal_death",
                                   "hepatobiliary_death",
                                   "neoplastic_death",
                                   "respiratory_death",
                                   "external_death"))) %>%
  set_variable_labels(cluster_nolab = "Diabetes subgroup",
                      MALE = "Male",
                      tiempo_evol = "Years since diagnosis",
                      diab_TX = "Diabetes treatment",
                      HT_TX = "Antihypertensive treatment",
                      DRUG_E1 = "Statin treatment",
                      comorbidity_count = "Number of comorbidities",
                      EDUGP = "Educational attainment",
                      COYOACAN = "Municipality",
                      SMOKER = "Smoking",
                      IDS = "Social Development Index",
                      PHYSGP = "Physical activity (times per week)",
                      Total_C = "Total cholesterol",
                      Total_TG = "Triglycerides",
                      BMI = "BMI (kg/m²)",
                      meanSBP = "SBP (mmHg)",
                      meanDBP = "DBP (mmHg)",
                      WHR = "Waist-to-Hip Ratio",
                      WHtR = "Waist-to-Height Ratio",
                      BASE_HBA1C = "HbA1c (%)",
                      score2diab_recalib = "SCORE2-DM",
                      SCORE2_DM_class = "SCORE2-DM Risk Category",
                      AGE = "Age (years)") %>% as.data.frame()
 
 lexis_data <- data_fg %>% filter(!is.na(YEAR_RECRUITED)) %>%mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
 lexis_data$age_risk <- lexis_data$AGE+lexis_data$PERSON_YEARS
 mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                     exit = list("period" = date_recruited + PERSON_YEARS), exit.status = fatal, data = lexis_data)
 mcps_fin2 <- splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
 mcps_fin2$ageout <- mcps_fin2$age + mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
 mcps_fin2$cvd_death<-ifelse(mcps_fin2$lex.Xst==T,1,0)
 mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
 mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
 mcps_fin2<- mcps_fin2 %>% rename(
   "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
   mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
 data_fg<-mcps_fin2
 
 #CVD
 fg_cvd <- finegray(Surv(PERSON_YEARS, fatal) ~
   cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
   comorbidity_count + EDUGP + COYOACAN + SMOKER +
   IDS + PHYSGP + Total_C + Total_TG + PATID + age.cat,
  data = data_fg,
  id = PATID,
  etype = "cvd_death")
 
 fg_CVD_model <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~
   cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
   comorbidity_count + EDUGP + COYOACAN + SMOKER +
   IDS + PHYSGP + Total_C + Total_TG + cluster(PATID) + strata(age.cat),
  weight = fgwt,
  data = fg_cvd,
  robust = TRUE)
 
 table_fg_CVD <- tbl_regression(fg_CVD_model,
                                    exponentiate = TRUE,
                                    label = var_list,
                                    conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 #Renal
 fg_renal <- finegray(Surv(PERSON_YEARS, fatal) ~
                     cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
                     comorbidity_count + EDUGP + COYOACAN + SMOKER +
                     IDS + PHYSGP + Total_C + Total_TG + AGE + PATID + age.cat,
                    data = data_fg,
                    id = PATID,
                    etype = "renal_death")
 
 fg_renal_model <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~
   cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
   comorbidity_count + EDUGP + COYOACAN + SMOKER +
   IDS + PHYSGP + Total_C + Total_TG + cluster(PATID) + strata(age.cat),
  weight = fgwt,
  data = fg_renal,
  robust = TRUE)
 
 table_fg_renal <- tbl_regression(fg_renal_model,
                                exponentiate = TRUE,
                                label = var_list,
                                conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 #Hepatobiliary
 fg_hpb <- finegray(Surv(PERSON_YEARS, fatal) ~
                       cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
                       comorbidity_count + EDUGP + COYOACAN + SMOKER +
                       IDS + PHYSGP + Total_C + Total_TG + AGE + PATID + age.cat,
                      data = data_fg,
                      id = PATID,
                      etype = "hepatobiliary_death")
 
 fg_hpb_model <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~
   cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
   comorbidity_count + EDUGP + COYOACAN + SMOKER +
   IDS + PHYSGP + Total_C + Total_TG + cluster(PATID) + strata(age.cat),
  weight = fgwt,
  data = fg_hpb,
  robust = TRUE)
 
 table_fg_hpb <- tbl_regression(fg_hpb_model,
                                  exponentiate = TRUE,
                                  label = var_list,
                                  conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 #Neoplastic
 fg_neoplastic <- finegray(Surv(PERSON_YEARS, fatal) ~
                     cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
                     comorbidity_count + EDUGP + COYOACAN + SMOKER +
                     IDS + PHYSGP + Total_C + Total_TG + AGE + PATID + age.cat,
                    data = data_fg,
                    id = PATID,
                    etype = "neoplastic_death")
 
 fg_neoplastic_model <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~
   cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
   comorbidity_count + EDUGP + COYOACAN + SMOKER +
   IDS + PHYSGP + Total_C + Total_TG + cluster(PATID) + strata(age.cat),
  weight = fgwt,
  data = fg_neoplastic,
  robust = TRUE)
 
 table_fg_neoplastic <- tbl_regression(fg_neoplastic_model,
                                exponentiate = TRUE,
                                label = var_list,
                                conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 #Respiratory
 fg_respiratory <- finegray(Surv(PERSON_YEARS, fatal) ~
                            cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
                            comorbidity_count + EDUGP + COYOACAN + SMOKER +
                            IDS + PHYSGP + Total_C + Total_TG + AGE + PATID + age.cat,
                           data = data_fg,
                           id = PATID,
                           etype = "respiratory_death")
 
 fg_respiratory_model <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~
   cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
   comorbidity_count + EDUGP + COYOACAN + SMOKER +
   IDS + PHYSGP + Total_C + Total_TG + cluster(PATID) + strata(age.cat),
  weight = fgwt,
  data = fg_respiratory,
  robust = TRUE)
 
 table_fg_respiratory <- tbl_regression(fg_respiratory_model,
                                       exponentiate = TRUE,
                                       label = var_list,
                                       conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 
 
 #External
 fg_external <- finegray(Surv(PERSON_YEARS, fatal) ~
                             cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
                             comorbidity_count + EDUGP + COYOACAN + SMOKER +
                             IDS + PHYSGP + Total_C + Total_TG + AGE + PATID + age.cat,
                            data = data_fg,
                            id = PATID,
                            etype = "external_death")
 
 fg_external_model <- coxph(
  Surv(fgstart, fgstop, fgstatus) ~
   cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
   comorbidity_count + EDUGP + COYOACAN + SMOKER +
   IDS + PHYSGP + Total_C + Total_TG + cluster(PATID) + strata(age.cat),
  weight = fgwt,
  data = fg_external,
  robust = TRUE)
 
 table_fg_external <- tbl_regression(fg_external_model,
                                        exponentiate = TRUE,
                                        label = var_list,
                                        conf.level = 0.95) %>% bold_p(t = 0.05) %>%
  modify_header(label ~ "**Covariate**") 

 
 # Merge tables 
 tab_fg <- tbl_merge(tbls = list(table_fg_CVD, table_fg_renal, table_fg_hpb, table_fg_neoplastic, table_fg_respiratory, table_fg_external),
                 tab_spanner = c("Cardiovascular death", "Renal death", "Hepatobiliary death", "Neoplastic death", "Respiratory death", "Other")) %>%
  as_flex_table()
 
 read_docx() %>%
  body_add_flextable(value = tab_fg, split = TRUE) %>%
  body_end_section_landscape() %>% 
  print(target = "Tablas/Table - Subdistribution Hazard Ratios for the Fine and Gray model.docx")
 
 rm(var_list, data_fg, 
    fg_cvd, fg_CVD_model, table_fg_CVD,
    fg_renal, fg_renal_model, table_fg_renal,
    fg_hpb, fg_hpb_model, table_fg_hpb,
    fg_neoplastic, fg_neoplastic_model, table_fg_neoplastic,
    fg_respiratory, fg_respiratory_model, table_fg_respiratory,
    fg_external, fg_external_model, table_fg_external)
 

 lexis_data <- lexis_data <- data_fg %>% filter(!is.na(YEAR_RECRUITED)) %>%mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
 lexis_data$age_risk <- lexis_data$AGE+lexis_data$PERSON_YEARS
 mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                     exit = list("period" = date_recruited + PERSON_YEARS), exit.status = fatal_2, data = lexis_data)
 mcps_fin2 <- splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
 mcps_fin2$ageout <- mcps_fin2$age + mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
 mcps_fin2$cvd_death<-ifelse(mcps_fin2$lex.Xst==T,1,0)
 mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
 mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
 mcps_fin2<- mcps_fin2 %>% rename(
   "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
   mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
 data_fg<-mcps_fin2
 
 #CVD
 fg_cvd <- finegray(Surv(PERSON_YEARS, fatal_2) ~
                      cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
                      comorbidity_count + EDUGP + COYOACAN + SMOKER +
                      IDS + PHYSGP + Total_C + Total_TG + PATID + age.cat,
                    data = data_fg,
                    id = PATID,
                    etype = "cvd_death")
 
 fg_CVD_model <- coxph(
   Surv(fgstart, fgstop, fgstatus) ~
     cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
     comorbidity_count + EDUGP + COYOACAN + SMOKER +
     IDS + PHYSGP + Total_C + Total_TG + cluster(PATID) + strata(age.cat),
   weight = fgwt,
   data = fg_cvd,
   robust = TRUE)
 
 table_fg_CVD <- tbl_regression(fg_CVD_model,
                                exponentiate = TRUE,
                                label = var_list,
                                conf.level = 0.95) %>% bold_p(t = 0.05) %>%
   modify_header(label ~ "**Covariate**")  %>% filter(!is.na(YEAR_RECRUITED)) %>%mutate("date_recruited" = decimal_date(ym(paste0(YEAR_RECRUITED, MONTH_RECRUITED))))
 lexis_data$age_risk <- lexis_data$AGE+lexis_data$PERSON_YEARS
 mcps_lexis <- Lexis(entry = list("period" = date_recruited, "age" = AGE, "time_at_entry" = 0),
                     exit = list("period" = date_recruited + PERSON_YEARS), exit.status = fatal2, data = lexis_data)
 mcps_fin2 <- splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
 mcps_fin2$ageout <- mcps_fin2$age + mcps_fin2$lex.dur; mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
 mcps_fin2$cvd_death<-ifelse(mcps_fin2$lex.Xst==T,1,0)
 mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 #Deaths 35-74
 mcps_fin2$old.death<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$old.death[mcps_fin2$old.death==1 & mcps_fin2$ageout<75]<-0 #Deaths ≥75
 mcps_fin2<- mcps_fin2 %>% rename(
   "id" = lex.id, "fu_period" = lex.dur, "entry_status" = lex.Cst, "exit_status" = lex.Xst) %>%
   mutate("time_at_exit" = time_at_entry + fu_period, "age_at_exit" = age + fu_period)
 data_fg<-mcps_fin2
 
 #CVD
 cardiac_death <- finegray(Surv(PERSON_YEARS, fatal_2) ~
                      cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
                      comorbidity_count + EDUGP + COYOACAN + SMOKER +
                      IDS + PHYSGP + Total_C + Total_TG + PATID + age.cat,
                    data = data_fg,
                    id = PATID,
                    etype = "cardiac_death")
 
 fg_CVD_model <- coxph(
   Surv(fgstart, fgstop, fgstatus) ~
     cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
     comorbidity_count + EDUGP + COYOACAN + SMOKER +
     IDS + PHYSGP + Total_C + Total_TG + cluster(PATID) + strata(age.cat),
   weight = fgwt,
   data = cardiac_death,
   robust = TRUE)
 
 table_fg_CVD <- tbl_regression(fg_CVD_model,
                                exponentiate = TRUE,
                                label = var_list,
                                conf.level = 0.95) %>% bold_p(t = 0.05) %>%
   modify_header(label ~ "**Covariate**") 
 
 #Stroke
 stroke_death <- finegray(Surv(PERSON_YEARS, fatal_2) ~
                             cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
                             comorbidity_count + EDUGP + COYOACAN + SMOKER +
                             IDS + PHYSGP + Total_C + Total_TG + PATID + age.cat,
                           data = data_fg,
                           id = PATID,
                           etype = "stroke_death")
 
 fg_stroke_model <- coxph(
   Surv(fgstart, fgstop, fgstatus) ~
     cluster_nolab + MALE + tiempo_evol + diab_TX + HT_TX + DRUG_E1 +
     comorbidity_count + EDUGP + COYOACAN + SMOKER +
     IDS + PHYSGP + Total_C + Total_TG + cluster(PATID) + strata(age.cat),
   weight = fgwt,
   data = stroke_death,
   robust = TRUE)
 
 table_fg_stroke <- tbl_regression(fg_stroke_model,
                                exponentiate = TRUE,
                                label = var_list,
                                conf.level = 0.95) %>% bold_p(t = 0.05) %>%
   modify_header(label ~ "**Covariate**") 
 
 
#### Sensitivity analysis: SCORE2-Diabetes risk regions recalibration comparisons ####
 mod_calib_all_df <- mcps_cox %>% 
  select(PATID,score2_dm_risk, score2diab_recalib, score2_dm_risk_low, score2_dm_risk_high, score2_dm_risk_vhigh, tstop, fatal_narrow_total) %>% 
  mutate(across(2:6, ~ log(-log(1 - .x)), 
                .names = "log_min_{.col}"))
 
 vars_cvd <- list("log_min_score2_dm_risk","log_min_score2diab_recalib","log_min_score2_dm_risk_low","log_min_score2_dm_risk_high","log_min_score2_dm_risk_vhigh")
 
 calib_models_all_overall <- vars_cvd %>% 
  map(function(x) {
   formula <- as.formula(paste0("Surv(tstop, fatal_narrow_total) ~ rcs(", x, ",3)"))
   model <- cph(formula, x = TRUE, y = TRUE, surv = TRUE, data = mod_calib_all_df)
  })
 
 observed01 <- survest(calib_models_all_overall[[1]], times = 10, newdata = mod_calib_all_df)
 
 obs_score_dm_original <- data.frame(obs = 1 - c(observed01$surv),
                                     upper = 1 - c(observed01$lower),
                                     lower = 1 - c(observed01$upper),
                                     pred = mod_calib_all_df$score2_dm_risk)
 
 observed02 <- survest(calib_models_all_overall[[2]], times = 10, newdata = mod_calib_all_df)
 
 obs_score_dm_recalib <- data.frame(obs = 1 - c(observed02$surv),
                                    upper = 1 - c(observed02$lower),
                                    lower = 1 - c(observed02$upper),
                                    pred = mod_calib_all_df$score2diab_recalib)
 
 observed03 <- survest(calib_models_all_overall[[3]], times = 10, newdata = mod_calib_all_df)
 
 obs_score_dm_low <- data.frame(obs = 1 - c(observed03$surv),
                                    upper = 1 - c(observed03$lower),
                                    lower = 1 - c(observed03$upper),
                                    pred = mod_calib_all_df$score2_dm_risk_low)
 observed04 <- survest(calib_models_all_overall[[4]], times = 10, newdata = mod_calib_all_df)
 
 obs_score_dm_high <- data.frame(obs = 1 - c(observed04$surv),
                                    upper = 1 - c(observed04$lower),
                                    lower = 1 - c(observed04$upper),
                                    pred = mod_calib_all_df$score2_dm_risk_high)
 
 observed05 <- survest(calib_models_all_overall[[5]], times = 10, newdata = mod_calib_all_df)
 
 obs_score_dm_vhigh <- data.frame(obs = 1 - c(observed05$surv),
                                 upper = 1 - c(observed05$lower),
                                 lower = 1 - c(observed05$upper),
                                 pred = mod_calib_all_df$score2_dm_risk_vhigh)
 
 
 models <- list(obs_score_dm_low,obs_score_dm_original,obs_score_dm_high,obs_score_dm_vhigh,obs_score_dm_recalib)
 model_names <- list("Low-risk region", "Moderate-risk region", "High-risk region", "Very high-risk region", "Recalibrated")
 
 obs_total <- purrr::map2_dfr(
  models,
  model_names,
  ~ dplyr::mutate(.x, model = .y))
 
 fig_calib_all <- obs_total %>% ggplot(aes(x = pred, y = obs, ymin = lower, ymax = upper,
             color = model)) +
  geom_line() +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = 2,
              color = "red") +
  coord_cartesian(xlim = c(0, 1),
                  ylim = c(0, 1)) +
  theme_minimal() +
  labs(title = "Calibration of SCORE2-Diabetes with original risk regions and recalibration",
       x = "Predicted 10-year risk",
       y = "Observed 10-year risk",
       color = "") +
  theme(legend.position = "top",
        legend.title = element_text(face = "bold", hjust = 0.5),
        legend.title.align = 0.5)
 
 ggsave(fig_calib_all,filename = "Figuras/SFig - SCORE2-Diabetes calibration for risk regions.jpg", 
        width = 15, height = 10,units=c("cm"),
        dpi = 300,limitsize = FALSE) 
 
 rm(mod_calib_all_df,vars_cvd,calib_models_all_overall,observed01,obs_score_dm_original,
    observed02,obs_score_dm_recalib, observed03, obs_score_dm_low, observed04, obs_score_dm_high,
    observed05, obs_score_dm_vhigh, models, names, obs_total)
 
#### Number of individuals excluded ####

mcps_excl <- mcps_dm %>% #29,948
 filter(BASE_HEARTATTACK == 0) %>% filter(BASE_ANGINA == 0) %>% filter(BASE_STROKE == 0) %>% #28,538 (1,410 baseline CVD removed)
 filter(!is.na(meanSBP)) %>% #28,492 (46 removed)
 filter(!is.na(meanDBP)) %>% #29,492 (0 removed)
 filter(!(meanSBP < 70 | meanSBP > 270)) %>% #28,492
 filter(!is.na(BMI)) %>% #27,962 (530 removed)
 filter(!(BMI < 10 | BMI > 80)) %>% #27,960 (2 removed)
 filter(!is.na(Total_C)) %>% #27,064 (896 removed)
 filter(!(Total_C < (1.75*38.67) | Total_C > (20*38.67))) %>% #26,738 (326 removed)
 filter(!is.na(WHtR)) %>% #26,729 (9 removed)
 filter(!is.na(eGFR)) %>% #25,721 (1008 removed)
 filter(!is.na(6)) %>% #25,680 (41 removed)
 filter(!(is.na(BASE_DIABETES_DX) & BASE_DIABETES==1)) %>% #25,634 (46 removed)
 filter(STATUS != "U") #24,943 (691 removed)
 

# Sensitivity analysis #
 mcps_loss <- mcps_dm %>% 
  filter(BASE_HEARTATTACK == 0) %>% filter(BASE_ANGINA == 0) %>% filter(BASE_STROKE == 0) %>%
  mutate(missing = ifelse((STATUS == "A" | STATUS == "D"),0,1)) %>% 
  filter(!((meanSBP > 270) | (meanSBP < 70))) %>%
  filter(!(BMI < 10 | BMI > 80)) %>%
  filter(!(Total_C > (20*38.67) | Total_C < (1.75*38.67) | is.na( Total_C))) %>%
  filter(!is.na(BASE_HBA1C),!is.na(WHtR),!is.na(BMI),!is.na(meanSBP),!is.na(meanDBP),!is.na(eGFR),!is.na(Total_C)) %>%
  filter(!(is.na(BASE_DIABETES_DX) & BASE_DIABETES==1)) %>%
  mutate(year_diagnosis = case_when(
   BASE_DIABETES_DX == 1 ~ 1955,
   BASE_DIABETES_DX == 2 ~ 1965,
   BASE_DIABETES_DX == 3 ~ 1975,
   BASE_DIABETES_DX == 4 ~ 1985,
   BASE_DIABETES_DX == 5 ~ 1995,
   BASE_DIABETES_DX == 6 ~ 2005,
   TRUE ~ NA_real_
  )) %>%
  mutate(year_birth = YEAR_RECRUITED - AGE) %>% 
  mutate(tiempo_evol = YEAR_RECRUITED - year_diagnosis) %>%
  mutate(AGE_DX = ifelse(DM_DX == 1, AGE - tiempo_evol, AGE)) %>%
  mutate(AGE_DX =  ifelse(AGE_DX < 0 | is.na(AGE_DX), 0, AGE_DX)) %>%
  mutate(comorbidity_count = rowSums(select(., BASE_EMPHYSEMA, BASE_ASTHMA, BASE_CKD, BASE_PEP, BASE_CIRR, 
                                            BASE_HYPERTENSION, BASE_LUNGCANCER, BASE_OTHCANCER, BASE_PROSTATECANCER, 
                                            BASE_CERVCANCER, BASE_BREASTCANCER, BASE_STOMCANCER, BASE_ORALCANCER, BASE_PAD), 
                                     na.rm = TRUE)) 
 
 var_labels <- list(
  AGE         = "Age (years)",
  MALE        = "Male (%)",
  AGE_DX      = "Estimated age at diabetes diagnosis (years)",
  DM_noDX     = "Undiagnosed diabetes (%)",
  SMOKER      = "Smoking (%)",
  BMI         = "Body mass index (kg/m²)",
  meanSBP     = "Systolic blood pressure (mmHg)",
  meanDBP     = "Diastolic blood pressure (mmHg)",
  WHR         = "Waist-To-Hip Ratio",
  WHtR        = "Waist-To-Height Ratio",
  Glucose     = "Glucose (mg/dl)",
  BASE_HBA1C  = "Glycated HbA1c (%)",
  mets_ir     = "METS-IR",
  Creatinine  = "Creatinine (mg/dl)",
  eGFR        = "eGFR (ml/min/1.73m²)",
  Total_C     = "Total cholesterol (mg/dl)",
  HDL_C       = "HDL cholesterol (mg/dl)",
  non_HDL_C   = "Non-HDL cholesterol (mg/dl)",
  Total_TG    = "Total triglycerides (mg/dl)",
  ApoB        = "Apolipoprotein B (mg/dl)",
  DRUG_D1     = "Metformin (%)",
  DRUG_D2     = "Sulphonylurea (%)",
  DRUG_D3     = "Insulin (%)",
  DRUG_D4     = "Other antidiabetics (%)")
 
 vars <- c("AGE","MALE","AGE_DX","DM_noDX","SMOKER","BMI","meanSBP","meanDBP","WHR","WHtR",
           "Glucose","BASE_HBA1C","mets_ir","Creatinine","eGFR",
           "Total_C","HDL_C","non_HDL_C","Total_TG","ApoB",
           "DRUG_D1","DRUG_D2","DRUG_D3","DRUG_D4")
 
 tab_all_loss <- make_summary(filter(mcps_loss))
 tab_noloss<- make_summary(filter(mcps_loss, missing == 0),
                         header(filter(mcps_loss, missing == 0), mcps_loss))
 tab_loss<- make_summary(filter(mcps_loss, missing == 1),
                         header(filter(mcps_loss, missing == 1), mcps_loss))
 
 
 tab_final_loss <- tbl_merge(
  list(tab_all_loss,tab_noloss,tab_loss),
  tab_spanner = c("All", "Included", "Lost to follow-up")
 ) %>%
  modify_table_styling(columns = everything(), footnote = NA) %>%
  as_flex_table() %>%
  add_header_lines("Table Characteristics of participants included and lost to follow-up in the Mexico City Prospective Study") %>%
  bold(part = "header") %>%
  theme_vanilla() 
 
 read_docx() %>%
   body_add_flextable(value = tab_final_loss, split = TRUE) %>%
   body_end_section_landscape() %>%
   print(target = "Tablas/Table Characteristics of participants included and lost to follow-up in the Mexico City Prospective Study.docx")
 
 
 library(gtsummary)
 library(dplyr)
 
 tab_final_loss <-
  mcps_loss %>%
  mutate(
   missing = factor(missing,
                    levels = c(0,1),
                    labels = c("Included", "Lost to follow-up"))
  ) %>%
  select(all_of(vars), missing) %>%
  tbl_summary(
   by = missing,
   label = var_labels,
   statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n} ({p}%)"
   ),
   missing = "no"
  ) %>%
  add_p() %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels() %>%
  modify_caption("**Table. Baseline characteristics of participants included and lost to follow-up in the Mexico City Prospective Study**") %>%
  bold_p(t = 0.05) 
 
 
#### Cross validation for C-statistic ####
 
 data_model2 <- mcps_dm_fin %>% select("PERSON_YEARS","fatal_narrow_total","fatal_narrow_cardiac","fatal_narrow_cerebrovascular","cluster_nolab","BMI","meanSBP","meanDBP","WHR","WHtR","BASE_HBA1C","MALE",
                                      "AGE","diab_TX","comorbidity_count","score2diab_recalib","SCORE2_DM_class","COYOACAN",
                                      "SMOKER","EDUGP","IDS","tiempo_evol","PHYSGP") %>% na.omit() %>%
  mutate(cluster_nolab = relevel(cluster_nolab, ref = "MOD")) %>%
  set_variable_labels(cluster_nolab = "Cluster",
                      BMI = "BMI (kg/m²)",
                      meanSBP = "SBP (mmHg)",
                      meanDBP = "DBP (mmHg)",
                      WHR = "Waist-to-Hip Ratio",
                      WHtR = "Waist-to-Height Ratio",
                      BASE_HBA1C = "HbA1c (%)",
                      COYOACAN = "Municipality",
                      SMOKER = "Smoking",
                      EDUGP = "Educational attainment",
                      MALE = "Sex",
                      tiempo_evol = "Years since diagnosis",
                      diab_TX = "Diabetes treatment",
                      comorbidity_count = "Number of comorbidities",
                      score2diab_recalib = "SCORE2-DM",
                      SCORE2_DM_class = "SCORE2-DM Risk Category",
                      IDS = "Social Development Index",
                      PHYSGP = "Physical activity (times per week)",
                      AGE = "Age (years)") %>% as.data.frame()
 
 dd <- datadist(data_model2)
 options(datadist = "dd")
 
 fit_cox <- function(formula, data) {
  cph(formula, data = data,x = T, y = T, surv = T)}
 
 validate_rms_model <- function(fit, B = 300) {
  
  val <- validate(fit, method = "boot", B = B, dxy = TRUE)
  
  dxy_app  <- val["Dxy", "index.orig"]
  dxy_corr <- val["Dxy", "index.corrected"]
  optimism <- val["Dxy", "optimism"]
  
  # Convert Dxy to C-index
  c_app  <- (dxy_app  + 1) / 2
  c_corr <- (dxy_corr + 1) / 2
  opt_c  <- optimism / 2
  
  tibble(
   C_apparent  = c_app,
   C_corrected = c_corr,
   optimism = opt_c
  )
 } 
 
 fit_models_for_outcome <- function(outcome, data, B = 300) {
  
  formulas <- list(
   
   "Risk factors" =
    as.formula(paste0(
     "Surv(PERSON_YEARS, ", outcome,
     ") ~ MALE + AGE + tiempo_evol + diab_TX + comorbidity_count + 
           EDUGP + COYOACAN + SMOKER + IDS + PHYSGP"
    )),
   
   "Risk factors+\nDiabetes subgroups" =
    as.formula(paste0(
     "Surv(PERSON_YEARS, ", outcome,
     ") ~ cluster_nolab + MALE + AGE + tiempo_evol + diab_TX + 
           comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP"
    )),
   
   "Risk factors+\nSCORE2-Diabetes" =
    as.formula(paste0(
     "Surv(PERSON_YEARS, ", outcome,
     ") ~ score2diab_recalib + MALE + AGE + tiempo_evol + diab_TX + 
           comorbidity_count + EDUGP + COYOACAN + SMOKER + IDS + PHYSGP"
    )),
   
   "Risk factors+\nSCORE2-Diabetes+\nDiabetes subgroups" =
    as.formula(paste0(
     "Surv(PERSON_YEARS, ", outcome,
     ") ~ score2diab_recalib + cluster_nolab + MALE + AGE + 
           tiempo_evol + diab_TX + comorbidity_count + 
           EDUGP + COYOACAN + SMOKER + IDS + PHYSGP"
    ))
  )
  
  models <- map(formulas, ~ fit_cox(.x, data))
  
  validation <- map_dfr(names(models), function(name) {
   val <- validate_rms_model(models[[name]], B = B)
   val$class <- name
   val
  })
  
  list(models = models, validation = validation)
 }
 
 models <- fit_models_for_outcome("fatal_narrow_cardiac", data_model)
 
 ### C-statistic and AIC ###

 plot_cstat_aic <- function(model_list, B = 300) {
  
  models <- model_list$models
  validation <- model_list$validation
  
  model_colors <- c(
   "Risk factors" = "#374E55FF",
   "Risk factors+\nDiabetes subgroups" = "#DF8F44FF",
   "Risk factors+\nSCORE2-Diabetes" = "#00A1D5FF",
   "Risk factors+\nSCORE2-Diabetes+\nDiabetes subgroups" = "#B24745FF"
  )
  
  ### C-INDEX PLOT ###
  fig_cstat <- ggplot(validation,
                      aes(x = class,
                          y = C_corrected,
                          color = class)) +
   
   geom_point(size = 4) +
   
   geom_text(aes(label = sprintf("%.3f", C_corrected)),
             vjust = -1.2,
             size = 4,
             fontface = "bold",
             color = "black") +
   
   scale_color_manual(values = model_colors) +
   
   scale_y_continuous(
    limits = c(0.60, max(validation$C_corrected) + 0.05),
    expand = expansion(mult = c(0, 0.05))
   ) +
   
   labs(x = "Model",
        y = "Optimism-corrected C-index",
        color = "Model") +
   
   theme_pubclean() +
   theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
   )
  
  ### AIC ###
  aic_data <- map_dfr(names(models), function(name) {
   tibble(
    class = name,
    AIC = AIC(models[[name]])
   )
  })
  
  aic_data <- aic_data %>%
   mutate(delta_aic = AIC - min(AIC)) 
  
  fig_aic <- ggplot(aic_data,
                    aes(x = class,
                        y = delta_aic,
                        fill = class)) +
   geom_bar(stat = "identity", color = "black", width = 0.6) +
   scale_fill_manual(values = model_colors) +
   labs(x = "Model",
        y = "ΔAIC",
        fill = "Model") +
   theme_pubclean() +
   theme(legend.position = "none",
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank())
  
  fig <- ggarrange(fig_cstat, fig_aic, ncol = 2)
  
  list(plot = fig, legend = legend)
 }
 
 outcomes <- c(
  "fatal_narrow_total",
  "fatal_narrow_cardiac",
  "fatal_narrow_cerebrovascular"
 )
 
 figures <- map(outcomes, function(outcome) {
  models <- fit_models_for_outcome(outcome, data_model2, B = 300)
  plot_cstat_aic(models)
 })
 
 fig3v2 <- ggarrange(
  figures[[1]]$plot,
  figures[[2]]$plot,
  figures[[3]]$plot,
  labels = "AUTO",
  ncol = 1,
  common.legend = TRUE,
  legend = "top"
 )
 
 ggsave("Figuras/Figure 3_rev.jpg",plot=fig3v2, width=25, height=18, units="cm", dpi=300, limitsize=FALSE)
 
 