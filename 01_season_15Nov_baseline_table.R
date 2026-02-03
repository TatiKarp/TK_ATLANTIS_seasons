## This script will prepare a baseline table for the 2 seasons groups
library(dplyr)
library(lubridate)
library(data.table)
library(tableone)
library(config)


# location 
conf <- config::get()

# data
master.Table <- read.csv(file.path(conf$data_path,"Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv"))

big_master_table <-  read.csv(file.path(conf$data_path, "atlantis_patient_data.csv"), header =TRUE, na.strings=c("","NA"))
additional_data <- read.csv(file.path(conf$data_path, "Umi_dedup/Dif_expr/Asthma_groups/clinicaldata_bhr_feno_mld.csv"), header =TRUE, na.strings=c("","NA"))%>%
  dplyr::select(-c(X,FENRES, MLD_ratio))
clinical_table <- big_master_table[,c('PT','PACKNO','BMI','PHADRES', 'DUR_DIS','AGE_DIAG',
                                      'B_FEV1F','GINA', 'NUM_EX','LAMA','BIO','SYS_COR',
                                      'ICS', 'ICS_LABA', 
                                      'ICS_LABA_DDOSE_EQ', 'ICS_DDOSE_EQ', 'acq6_score','LABEOSV',
                                      'LABNEUV','LABMACV','BRONCHP', 'LYMPHOP', 'EOSP', 'MACROP', 'NEUTROP',
                                      'B_TLCPNVF', 'B_RVTLCPNVF', 'B_FEV1PNVG','B_FEV1FPNVG', 'FENRES', 'PCD',
                                      'B_R520', 'B_SCOND', 'B_SACIN', 'B_F2575', 'B_F50',
                                      'LA', 'WA','TA', 'Pi10', 'WA_TA100', 'VI_856', 'VI_950', 'lung_ratio', 'MLD_ratio')]

clinical_table <- clinical_table[!duplicated(clinical_table$PT), ] %>%
  mutate (GINA = as.factor(GINA),
          NUM_EX = as.factor(NUM_EX),
          any_ICS = if_else(ICS == "No" & ICS_LABA == "No", "No", "Yes"),
          ICS_dose_sum = if_else(is.na(ICS_DDOSE_EQ) & is.na(ICS_LABA_DDOSE_EQ), NA,
                       coalesce(ICS_DDOSE_EQ, 0) + coalesce(ICS_LABA_DDOSE_EQ, 0)))


master.Table.ATLANTIS <- master.Table %>%
  left_join(clinical_table, by = c('PT'='PT')) %>%
  left_join(additional_data, by = c('PT'='PT'))

non_normally <- c('DUR_DIS', 'AGE_DIAG','PACKNO','acq6_score', 'LABEOSV', 'BRONCHP', 'LYMPHOP',
                  'MACROP', 'NEUTROP', 'B_R520','B_SCOND', 'B_SACIN', 'VI_856', 'VI_950','FENRES', 'EOSP')
var_for_table <- c('gender','age','smoking.status', 'asthma.status','PACKNO','BMI','PHADRES', 
                   'DUR_DIS','AGE_DIAG','B_FEV1F', 'GINA', 'NUM_EX',
                   'LAMA','BIO','SYS_COR','any_ICS','ICS_dose_sum',
                   'acq6_score','LABEOSV','LABNEUV','LABMACV','BRONCHP', 
                   'LYMPHOP', 'EOSP', 'MACROP', 'NEUTROP',
                   'B_TLCPNVF', 'B_RVTLCPNVF', 'B_FEV1PNVG','B_FEV1FPNVG', 
                   'FENRES', 'PCD','bhr', 'B_R520', 'B_SCOND', 'B_SACIN', 
                   'B_F2575', 'B_F50', 'LA', 'WA','TA', 'Pi10', 'WA_TA100', 
                   'VI_856', 'VI_950', 'lung_ratio','MLD_ratio')

tabtotal <- CreateTableOne(vars = var_for_table, strata = "new_seasons" , data = master.Table.ATLANTIS)

Final_statistics <- print(tabtotal, nonnormal = non_normally, exact = "stage", formatOptions = list(big.mark = ","), quote = FALSE, noSpaces = TRUE)

write.csv(Final_statistics, file = file.path(conf$data_path, "Season/Season_new_date/Clinical_characteristics_seasons.csv"))

