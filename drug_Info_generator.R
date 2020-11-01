library(data.table)
library(dplyr)
library(stringr)

all.cv = fread(file="AppData/allCOVID_drug_collection_v3.csv") %>% 
  as.data.frame() %>% 
  dplyr::filter(!duplicated(DrugBankID))

sourceInfo = fread(file="data/allCOVID_drug_collection_v2.csv") %>% 
  as.data.frame() %>% 
  dplyr::filter(!duplicated(DrugBankID))

all.cv = dplyr::left_join(all.cv, sourceInfo%>%dplyr::select(-c("Name")), by = "DrugBankID")
all.cv[which(all.cv$DrugBankID == "DBXXXX1"),"Source"] <- "Literature"
all.cv[which(all.cv$DrugBankID == "DBXXXX1"),"Source_link"] <- "https://www.clinicaltrials.gov/ct2/show/NCT04396717"

# read drugbank meta data
db_meta_data = fread("AppData/DB_full_meta_data_5_1_6_full.csv", encoding = "UTF-8") %>% dplyr::select("DB_ID","Groups") %>% as.data.frame()
all.cv_full = dplyr::left_join(all.cv, db_meta_data, by=c("DrugBankID" = "DB_ID"))

# drug_TherapClass <- fread("AppData/TTD_drugName_TherapClass_Mapping.txt", encoding = "UTF-8") %>% dplyr::select(2,3) %>% as.data.frame()
# all.cv = dplyr::left_join(all.cv, drug_TherapClass, by = c("Name" = "DrugName"))

# have therap_class Data ----------------------
drug.info.TTD = fread("AppData/TTD_drugName_TherapClass_Indication_Mapping.txt", encoding = "UTF-8") %>% dplyr::select(2,3,4)  %>% as.data.frame()
all.cv_full = dplyr::left_join(all.cv_full, drug.info.TTD, by = c("Name" = "DrugName"))

drug.info.from.drugBank = readRDS("data/DB_full_5_1_7_slim_atc_cat.rds") %>% 
  dplyr::filter(DB_ID %in% all.cv_full$DrugBankID)

drug.info.from.drugBank$`ATC codes` <- lapply(drug.info.from.drugBank$`ATC codes`, function(x){
  aL = strsplit(x, split = "\n", fixed = T) %>% unlist(use.names = F)
  return(aL[1])
}) %>% unlist(use.names = F)
drug.info.from.drugBank$Categories <- lapply(drug.info.from.drugBank$Categories, function(x){
  aL = strsplit(x, split = "\n", fixed = T) %>% unlist(use.names = F)
  aLengths = lapply(aL, function(x) str_length(x)) %>% unlist(use.names = F)
  aDF = data.frame(cat = aL,
                   len = aLengths, 
                   stringsAsFactors = F) %>% dplyr::arrange(len)
  return(aDF[1,1])
}) %>% unlist(use.names = F)

# for drugs without any therapeutic class data, first try to use their 5-level ATC-CODE
na_idx = which(is.na(all.cv_full$Therap_Class) | all.cv_full$Therap_Class == "")
all.cv_full$Therap_Class[na_idx] = drug.info.from.drugBank[match(all.cv_full$DrugBankID[na_idx],
                                                                 drug.info.from.drugBank$DB_ID), "ATC codes"] %>%
  unlist(use.names = F)

# still some of the would be NA, use CATEGORIES data for them
na_idx = which(is.na(all.cv_full$Therap_Class) | all.cv_full$Therap_Class == "")
all.cv_full$Therap_Class[na_idx] = drug.info.from.drugBank[match(all.cv_full$DrugBankID[na_idx],
                                                                 drug.info.from.drugBank$DB_ID), "Categories"] %>%
  unlist(use.names = F)
# ----------------------


# load disease proximity score
disease.proximity.drugs = fread("data/distance_based_disease_proximity_v5.csv") %>% dplyr::select(-1) %>%
  dplyr::rename(DrugID = drug, DP = d_td, z_DP = z, z_DP_pVal = p) %>% 
  as.data.frame()

all.cv_full = dplyr::left_join(all.cv_full, disease.proximity.drugs, by=c("DrugBankID" = "DrugID"))

# load disease funcitonal relevance score
functional.relevance.drugs = fread("data/FR_based_disease_proximity_V5.csv") %>% dplyr::select(-1) %>%
  dplyr::rename(DrugID = drug, FR = d_td, z_FR = z, z_FR_pVal = p) %>% 
  as.data.frame()

all.cv_full = dplyr::left_join(all.cv_full, functional.relevance.drugs, by=c("DrugBankID" = "DrugID"))

# load disease funcitonal proximity score
functional.proximity.drugs = fread("data/Functional_proximity_score_v5.csv", header = F, skip = 1) %>% dplyr::select(-1) %>%
  dplyr::rename(DrugID = V2, Functional_Proximity = V3) %>% 
  as.data.frame()
functional.proximity.drugs$Functional_Proximity = functional.proximity.drugs$Functional_Proximity / max(functional.proximity.drugs$Functional_Proximity, na.rm = T)

  all.cv_full = dplyr::left_join(all.cv_full, functional.proximity.drugs, by=c("DrugBankID" = "DrugID"))

# may still contain duplicated drugnames with different DrugbankIDs, 
# find them using >all.cv[which(duplicated(all.cv$Name)),], and manually find unique ones using DrugBank website

# print for the supplementary file
fwrite(all.cv_full, file="data/Supplementary Files/allCOVID_drug_collection_full_v5.csv")

