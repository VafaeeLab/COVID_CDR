library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)

dcomb.db <- fread("AppData/covid.drugcomp.csv", encoding = "UTF-8") %>% dplyr::select(c("Drug1","Drug2","classification")) %>% as.data.frame() 
cov.sim.db <- fread("COV_DrugSimDB_Full_with_signif_scores.csv")
cov.sim.db$Name1 <- cov.sim.db$Name1 %>% toupper()
cov.sim.db$Name2 <- cov.sim.db$Name2 %>% toupper()
dcomb.db <- dplyr::inner_join(dcomb.db, cov.sim.db[,c("Name1","Name2","Comp_Expo")], by=c("Drug1"="Name1","Drug2"="Name2")) %>% 
  # unique()  %>% # remove duplicate row
  dplyr::filter(!is.na(Comp_Expo))
dcomb.db <- dcomb.db[!duplicated(dcomb.db[,c(1,2)]),]


# plot
data <- dcomb.db[,c(3,4)]

ggplot(data, aes(x=classification, y=Comp_Expo)) + 
  geom_boxplot()

