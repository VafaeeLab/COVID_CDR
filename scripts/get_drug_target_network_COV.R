library(data.table)
library(dplyr)

all.cv_full <- fread("AppData/allCOVID_drug_collection_full_v5.csv", encoding = "UTF-8") %>% as.data.frame()


dt <- fread("AppData/all_compiled.csv") %>% as.data.frame()

dr_dt <- base::sapply(all.cv_full$DrugBankID, function(x){
  temp <- dt[which(grepl(x, dt$`Drug IDs`, fixed = T)), 6] %>% unlist(use.names = F)
  # print(temp)
  if(length(temp) > 0)
    cbind(temp, x) %>% as.data.frame() %>% return()
}) %>% do.call("rbind",.) %>% rename(UniProtID = temp, DrugBank_ID = x)
fwrite(dr_dt, file="data/drug_target_network_V5.csv")


