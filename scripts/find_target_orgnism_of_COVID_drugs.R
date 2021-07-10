library(dplyr)
library(data.table)
all.cv = fread("data/COVID drug collection/allCOVID_drug_collection_v2.csv") %>% as.data.frame()
targets = fread("data/Drugbank/5.1.7/all.csv") %>% as.data.frame()

targets.cov.drugs = do.call(rbind, lapply(all.cv$DrugBankID, function(x){
  return(targets[which(grepl(x, targets$`Drug IDs`,fixed = T)),])
}))

targets.cov.drugs = unique(targets.cov.drugs)

# list all the organisms
print(unique(targets.cov.drugs$Species))
length(unique(targets$Species))

# testset = c("DB00303","DB00142", "DB11300")
# targets = targets[1:5, ]
# # targets$`Drug IDs`[1:5]
# 
# targets.cov.drugs = do.call(rbind, lapply(testset, function(x){
#   return(targets[which(grepl(x, targets$`Drug IDs`,fixed = T)),])
# }))
