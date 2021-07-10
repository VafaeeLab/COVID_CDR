library(data.table)
library(dplyr)
library(tidyr)
source("initial_drug_screening.R")

# DrugBank
db.links = fread("data/Drugbank/5.1.6/drug links.csv")

db.cv = fread("data/COVID drug collection/drugbank_copyPaste.csv") 
db.cv$DrugLink = gsub("https://www.drugbank.ca/drugs/", "", db.cv$DrugLink, fixed = T)
db.cv = db.cv %>% select(c(Drug, DrugLink)) %>% as.data.frame()
db.cv = db.cv[!duplicated(t(apply(db.cv, 1, sort))),] %>% dplyr::select(c(2,1))
colnames(db.cv) = c("DrugBankID","Name")
db.cv$Source = "DrugBank"
db.cv$Source_link = "https://www.drugbank.ca/covid-19"
fwrite(db.cv, file="data/COVID drug collection/drugbank_covid.csv")

# therapeutic trackers
tt.cv = fread("data/COVID drug collection/covid_19_therapeutics_tracker.txt")
tt.cv = dplyr::inner_join(tt.cv, db.links, by = c("DrugName" = "Name")) %>% dplyr::select(c(2,1))
colnames(tt.cv) = c("DrugBankID","Name")
tt.cv$Source = "Therapeutic Tracker"
tt.cv$Source_link = "https://www.raps.org/news-and-articles/news-articles/2020/3/covid-19-therapeutics-tracker"

# WikiPedia
wiki.cv = fread("data/COVID drug collection/wikipedia_covid_drugs.txt", sep="\t")
wiki.cv = dplyr::inner_join(wiki.cv, db.links, by = c("DrugName" = "Name")) %>% select(c(2,1))
colnames(wiki.cv) = c("DrugBankID","Name")
wiki.cv$Source = "Wikipedia"
wiki.cv$Source_link = "https://en.wikipedia.org/wiki/COVID-19_drug_repurposing_research#Studies"

# Nature paper
ht.cv = initial_Drug_screeing(hop = 2) %>% na.omit()
colnames(ht.cv) = c("DrugBankID","Name")
ht.cv$Source = "Gordon et al."
ht.cv$Source_link = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2286-9/MediaObjects/41586_2020_2286_MOESM6_ESM.xlsx"

# combine all
all.cv = rbind(db.cv, tt.cv, wiki.cv, ht.cv)
all.cv = all.cv[!duplicated(t(apply(all.cv, 1, sort))),] 
all.cv = all.cv[!duplicated(all.cv$Name),]

fwrite(all.cv, file="data/COVID drug collection/allCOVID_drug_collection_v2.csv")

