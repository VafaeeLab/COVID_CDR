parse_COVID_Comb_DB_data <- funciton(drugbankFile, covid.drugcompFile, outfile = "Appdata/Cov_Com.drugcomDB.csv"){
  library(data.table)
  library(dplyr)
  
  cov.drugcomp = fread(covid.drugcompFile, encoding = "UTF-8") %>% dplyr::select(-c(1,2)) %>% as.data.frame()
  db = fread(drugbankFile, encoding = "UTF-8") %>% dplyr::select(c(1,2)) %>% as.data.frame()
  db$Name = db$Name %>% toupper()
  
  temp = dplyr::left_join(cov.drugcomp,db,by=c("Drug1"="Name")) %>% 
    dplyr::left_join(.,db,by=c("Drug2"="Name")) %>% dplyr::rename("DB_ID1" = "DrugBank ID.x", "DB_ID2" = "DrugBank ID.y")
  fwrite(temp %>% select(1,14,2,15,3:13), file = outfile)
}

drugbankFile = "AppData/drug links.csv"
covid.drugcompFile = "AppData/covid.drugcomp.csv"
parse_COVID_Comb_DB_data(drugbankFile, covid.drugcompFile)