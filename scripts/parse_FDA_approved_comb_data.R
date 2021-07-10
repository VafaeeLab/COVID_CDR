library(data.table)
library(readxl)
library(dplyr)

all.cv_full <- fread("AppData/allCOVID_drug_collection_full.csv", encoding = "UTF-8") %>% as.data.frame()

FDA_comb <- read_excel("AppData/FDA_approaved_combination_temp.xlsx") %>% dplyr::select(c(2,3)) %>% as.data.frame()
keep <- lapply(seq(1, nrow(FDA_comb)), function(i){
  x <- FDA_comb$DrugNames[i] %>% strsplit(split = " + ", fixed = T) %>% unlist()
  if(all(x %in% all.cv_full$Name))
    return(i)
  else
    return(NA)
}) %>% unlist() %>% na.omit()
FDA_comb <- FDA_comb[keep,]

FDA_comb$DrugIDs <- lapply(seq(1, nrow(FDA_comb)), function(i){
  x <- FDA_comb$DrugNames[i] %>% strsplit(split = " + ", fixed = T) %>% unlist()
  return(all.cv_full[match(x, all.cv_full$Name), "DrugBankID"] %>% paste0(.,collapse = " + "))
}) %>% unlist()

library(DT)
datatable(FDA_comb, options = list(columnDefs = list(list(
  targets = 2,
  render = JS(
    "function(data, type, row, meta) {",
    "return type === 'display' && data.length > 6 ?",
    "'<span title=\"' + data + '\">' + data.substr(0, 6) + '...</span>' : data;",
    "}")
))))
fwrite(FDA_comb[,c(3,1,2)], file = "AppData/FDA_approaved_combination.csv")
