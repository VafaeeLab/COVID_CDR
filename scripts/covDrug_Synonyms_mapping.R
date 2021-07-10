require(data.table)
require(dplyr)

# load drug data (IDs) from COVID-CDR
covDrugIDs = fread("AppData/allCOVID_drug_collection_full_v5.csv") %>% 
  as.data.frame() 
# load synonnym from Drugbank
syn = fread("data/drugbank vocabulary.csv") %>% as.data.frame()
# map
map = dplyr::inner_join(covDrugIDs %>% dplyr::select("DrugBankID"), 
                        syn %>% dplyr::select(c("DrugBank ID", "Common name", "Synonyms")),
                        by=c("DrugBankID"="DrugBank ID"))
# output
fwrite(map,"AppData/covDrug_Synonyms.csv")
