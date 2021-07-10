library(readxl)
library(dplyr)
library(data.table)
library(reshape2)


drugIDs.all <- fread("AppData/allCOVID_drug_collection_v3.csv") %>% rename(Drugbank_ID = DrugBankID, DrugName = Name)
drugIDs = drugIDs.all %>%
  dplyr::select(1) %>% 
  unlist(use.names = F) %>% 
  unique()

# read net-seperation matrix
ns.m <- as.matrix(read.csv("data/network_separation_score_of_drugPair_v5.csv", row.names = 1))

# read distance-based 
zScore_D <- fread("data/distance_based_disease_proximity_v5.csv") %>% dplyr::select(-1) %>%
  as.data.frame() %>% 
  dplyr::filter(drug %in% drugIDs) %>%
  rename("DP" = "d_td", "z_DP" = "z","z_DP_pVal"="p")
zScore_FR <- fread("data/FR_based_disease_proximity_V5.csv") %>% dplyr::select(-1) %>%
  as.data.frame() %>%
  dplyr::filter(drug %in% drugIDs) %>%
  rename("FR" = "d_td", "z_FR" = "z","z_FR_pVal"="p")

res = ns.m %>% reshape2::melt() %>% rename("Net_Sep" = "value") %>%
  dplyr::filter(Var1 != Var2 & Net_Sep != -Inf) %>% rename("Drug_ID1"="Var1", "Drug_ID2"="Var2") %>%
  
  dplyr::inner_join(., drugIDs.all, by=c("Drug_ID1"="Drugbank_ID")) %>% rename("DrugName1" = "DrugName") %>%
  dplyr::inner_join(., zScore_D, by=c("Drug_ID1"="drug")) %>% rename("DP_Drug1"="DP","z_DP_Drug1" = "z_DP", "z_DP_pVal_Drug1" = "z_DP_pVal") %>%
  dplyr::inner_join(., zScore_FR, by=c("Drug_ID1"="drug")) %>% rename("FR_Drug1"="FR", "z_FR_Drug1" = "z_FR", "z_FR_pVal_Drug1" = "z_FR_pVal") %>%
  
  dplyr::inner_join(., drugIDs.all, by=c("Drug_ID2"="Drugbank_ID")) %>% rename("DrugName2" = "DrugName") %>%
  dplyr::inner_join(., zScore_D, by=c("Drug_ID2"="drug")) %>% rename("DP_Drug2"="DP","z_DP_Drug2" = "z_DP", "z_DP_pVal_Drug2" = "z_DP_pVal") %>%
  dplyr::inner_join(., zScore_FR, by=c("Drug_ID2"="drug")) %>% rename("FR_Drug2"="FR", "z_FR_Drug2" = "z_FR", "z_FR_pVal_Drug2" = "z_FR_pVal") %>%
  
  dplyr::select(c("Drug_ID1","DrugName1","DP_Drug1","z_DP_Drug1","z_DP_pVal_Drug1","FR_Drug1","z_FR_Drug1","z_FR_pVal_Drug1",
                  "Drug_ID2","DrugName2","DP_Drug2","z_DP_Drug2","z_DP_pVal_Drug2","FR_Drug2","z_FR_Drug2","z_FR_pVal_Drug2",
                  "Net_Sep"))


# --------------

res$Combination_mode_zFR <- ifelse(((res$z_FR_Drug1 < 0 & res$z_FR_Drug2 < 0) | (res$z_FR_Drug2 < 0 & res$z_FR_Drug1 < 0) & res$Net_Sep >= 0), "Complementary exposure",
                                     ifelse(((res$z_FR_Drug1 < 0 & res$z_FR_Drug2 < 0) | (res$z_FR_Drug2 < 0 & res$z_FR_Drug1 < 0) & res$Net_Sep < 0), "Overlapping/Indirect exposure",
                                            ifelse(((res$z_FR_Drug1 < 0 & res$z_FR_Drug2 >= 0) | (res$z_FR_Drug2 < 0 & res$z_FR_Drug1 >= 0) & res$Net_Sep >= 0), "Single exposure",
                                                   ifelse(((res$z_FR_Drug1 >= 0 & res$z_FR_Drug2 >= 0) | (res$z_FR_Drug2 >= 0 & res$z_FR_Drug1 >= 0) & res$Net_Sep >= 0), "Independent exposure",
                                                          ifelse(((res$z_FR_Drug1 >= 0 & res$z_FR_Drug2 >= 0) | (res$z_FR_Drug2 >= 0 & res$z_FR_Drug1 >= 0) & res$Net_Sep < 0), "Non-exposure", "NA"))))
)
res$Combination_mode_zD <- ifelse(((res$z_DP_Drug1 < 0 & res$z_DP_Drug2 < 0) | (res$z_DP_Drug2 < 0 & res$z_DP_Drug1 < 0) & res$Net_Sep >= 0), "Complementary exposure",
                                    ifelse(((res$z_DP_Drug1 < 0 & res$z_DP_Drug2 < 0) | (res$z_DP_Drug2 < 0 & res$z_DP_Drug1 < 0) & res$Net_Sep < 0), "Overlapping/Indirect exposure",
                                           ifelse(((res$z_DP_Drug1 < 0 & res$z_DP_Drug2 >= 0) | (res$z_DP_Drug2 < 0 & res$z_DP_Drug1 >= 0) & res$Net_Sep >= 0), "Single exposure",
                                                  ifelse(((res$z_DP_Drug1 >= 0 & res$z_DP_Drug2 >= 0) | (res$z_DP_Drug2 >= 0 & res$z_DP_Drug1 >= 0) & res$Net_Sep >= 0), "Independent exposure",
                                                         ifelse(((res$z_DP_Drug1 >= 0 & res$z_DP_Drug2 >= 0) | (res$z_DP_Drug2 >= 0 & res$z_DP_Drug1 >= 0) & res$Net_Sep < 0), "Non-exposure", "NA"))))
)

# res$Combination_on_nsCls_zFR <- ifelse(((res$zFR_Drug1 < 0 & res$zFR_Drug2 < 0) | (res$zFR_Drug2 < 0 & res$zFR_Drug1 < 0) & res$Net_Separation_getD >= 0), "Complementary exposure",
#                                        ifelse(((res$zFR_Drug1 < 0 & res$zFR_Drug2 < 0) | (res$zFR_Drug2 < 0 & res$zFR_Drug1 < 0) & res$Net_Separation_getD < 0), "Overlapping/Indirect exposure",
#                                               ifelse(((res$zFR_Drug1 < 0 & res$zFR_Drug2 >= 0) | (res$zFR_Drug2 < 0 & res$zFR_Drug1 >= 0) & res$Net_Separation_getD >= 0), "Single exposure",
#                                                      ifelse(((res$zFR_Drug1 >= 0 & res$zFR_Drug2 >= 0) | (res$zFR_Drug2 >= 0 & res$zFR_Drug1 >= 0) & res$Net_Separation_getD >= 0), "Independent exposure",
#                                                             ifelse(((res$zFR_Drug1 >= 0 & res$zFR_Drug2 >= 0) | (res$zFR_Drug2 >= 0 & res$zFR_Drug1 >= 0) & res$Net_Separation_getD < 0), "Non-exposure", "NA"))))
# )
# res$Combination_on_nsCls_zD <- ifelse(((res$zD_Drug1 < 0 & res$zD_Drug2 < 0) | (res$zD_Drug2 < 0 & res$zD_Drug1 < 0) & res$Net_Separation_getD >= 0), "Complementary exposure",
#                                       ifelse(((res$zD_Drug1 < 0 & res$zD_Drug2 < 0) | (res$zD_Drug2 < 0 & res$zD_Drug1 < 0) & res$Net_Separation_getD < 0), "Overlapping/Indirect exposure",
#                                              ifelse(((res$zD_Drug1 < 0 & res$zD_Drug2 >= 0) | (res$zD_Drug2 < 0 & res$zD_Drug1 >= 0) & res$Net_Separation_getD >= 0), "Single exposure",
#                                                     ifelse(((res$zD_Drug1 >= 0 & res$zD_Drug2 >= 0) | (res$zD_Drug2 >= 0 & res$zD_Drug1 >= 0) & res$Net_Separation_getD >= 0), "Independent exposure",
#                                                            ifelse(((res$zD_Drug1 >= 0 & res$zD_Drug2 >= 0) | (res$zD_Drug2 >= 0 & res$zD_Drug1 >= 0) & res$Net_Separation_getD < 0), "Non-exposure", "NA"))))
# )

fwrite(res, file="data/drug_combination_mode_v5.csv")

