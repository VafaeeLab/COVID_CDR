source("data/TTD/TTD_data_processing.R")
library(data.table)
library(dplyr)
drugID_Therapeuticlass <- TTD_file_processing(fileName = "data/TTD/P1-02-TTD_drug_download.txt",
                           valueKey = "THERCLAS")
drugID_drugName <- TTD_file_processing(fileName = "data/TTD/P1-03-TTD_crossmatching.txt",
                                       valueKey = "DRUGNAME")
drugName_Therapeuticclass <- dplyr::inner_join(drugID_drugName, drugID_Therapeuticlass, by = "ID")
colnames(drugName_Therapeuticclass) <- c("TTD_ID","DrugName","Therap_Class")

fwrite(drugName_Therapeuticclass, file = "data/TTD/TTD_drugName_TherapClass_Mapping.txt")


res <- TTD_file_processing(fileName = "data/TTD/P1-05-Drug_disease.txt",
                           valueKey = "INDICATI")
drugName_Therapeuticclass_Indication <- dplyr::left_join(drugName_Therapeuticclass, res, by = c("TTD_ID" = "ID"))
colnames(drugName_Therapeuticclass_Indication)[4] <- "Indication(s)"
fwrite(drugName_Therapeuticclass_Indication, file = "data/TTD/TTD_drugName_TherapClass_Indication_Mapping.txt")
