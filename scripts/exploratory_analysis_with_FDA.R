library(data.table)
library(dplyr)
library(ggplot2)

epsilon <- 1e-10 # A very small number


# load drugsimdb ----------------
# drugSim.db = fread("data/new_net_info_V7_pval.csv")%>% as.data.frame()
# #load FDA approved combinations
# FDA.comb = fread("data/681 FDA approved drug combination.txt") %>% as.data.frame()
# 
# # inner join (or left join)
# FDA.comb.simInfo1 = dplyr::inner_join(FDA.comb, drugSim.db, by = c("ID1" = "ID1", "ID2" = "ID2"))
# FDA.comb.simInfo2 = dplyr::inner_join(FDA.comb, drugSim.db, by = c("ID2" = "ID1", "ID1" = "ID2"))
# FDA.comb.simInfo = dplyr::bind_rows(FDA.comb.simInfo1, FDA.comb.simInfo2)
# # remove dupliate pairs
# FDA.comb.simInfo[1:2] = t(apply(FDA.comb.simInfo[1:2],1,sort))
# FDA.comb.simInfo = FDA.comb.simInfo[!duplicated(FDA.comb.simInfo[1:2]),]
# 
# FDA.comb.simInfo[,10] = 1 - FDA.comb.simInfo[,10]
# colnames(FDA.comb.simInfo)[10] = "1 - pValue"
# 
# 
# # boxplot
# boxplot(FDA.comb.simInfo[,c(3:10)])

# plot.dat = rbind(
#   cbind("Chem_Similarity", FDA.comb.simInfo[,3]),
#   cbind("Target_similarity", FDA.comb.simInfo[,4]),
#   cbind("Pathway_similarity", FDA.comb.simInfo[,5]),
#   cbind("GO_CC_Similarity", FDA.comb.simInfo[,6]),
#   cbind("GO_MF_Similarity", FDA.comb.simInfo[,7]),
#   cbind("GO_BP_Similarity", FDA.comb.simInfo[,8]),
#   cbind("Combined", FDA.comb.simInfo[,9]),
#   cbind("1 - pValue", FDA.comb.simInfo[,10])
# ) %>% as.data.frame() 
# colnames(plot.dat) = c("class","sim_val")
# plot.dat[,1] = as.character(plot.dat[,1])
# f = plot.dat[1,2]
# plot.dat[,2] = as.numeric(levels(f))[f]
# plot.dat = plot.dat %>% as_tibble()
# 
# ggplot(plot.dat, aes(x=class, y=sim_val, fill=class)) + 
#   geom_boxplot(alpha=0.3) +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette="BuPu") +
#   ylab("") + 
#   xlab("") +
#   ggtitle("All drug similarities that are approved to treat in combination by FDA") +
# theme_classic()





# load COVID-19 drugsimdb -------------------
FDA.comb = fread("data/681 FDA approved drug combination.txt") %>% as.data.frame()

drugSim.db = fread("data/cov_DrugSimDB_V5.csv")%>% as.data.frame()
drugSim.db[1:2] = t(apply(drugSim.db[1:2],1,sort))
drugSim.db = drugSim.db[!duplicated(drugSim.db[1:2]),]

# inner join (or left join)
FDA.comb.simInfo1 = dplyr::inner_join(FDA.comb, drugSim.db, by = c("ID1" = "ID1", "ID2" = "ID2"))
FDA.comb.simInfo2 = dplyr::inner_join(FDA.comb, drugSim.db, by = c("ID2" = "ID1", "ID1" = "ID2"))
FDA.comb.simInfo = dplyr::bind_rows(FDA.comb.simInfo1, FDA.comb.simInfo2)
# remove dupliate pairs
FDA.comb.simInfo[1:2] = t(apply(FDA.comb.simInfo[1:2],1,sort))
FDA.comb.simInfo = FDA.comb.simInfo[!duplicated(FDA.comb.simInfo[1:2]),]
# fwrite(FDA.comb.simInfo, "FDA_COVID_drug_combination_V5.csv")

# boxplot
FDA.comb.simInfo[,10] = 1 - FDA.comb.simInfo[,10]
colnames(FDA.comb.simInfo)[10] = "1 - pValue"


boxplot(FDA.comb.simInfo[,c(3:10)])
# plot.dat = rbind(
#   cbind("Chem_Similarity", FDA.comb.simInfo[,3]),
#   cbind("Target_similarity", FDA.comb.simInfo[,4]),
#   cbind("Pathway_similarity", FDA.comb.simInfo[,5]),
#   cbind("GO_CC_Similarity", FDA.comb.simInfo[,6]),
#   cbind("GO_MF_Similarity", FDA.comb.simInfo[,7]),
#   cbind("GO_BP_Similarity", FDA.comb.simInfo[,8]),
#   cbind("Combined", FDA.comb.simInfo[,9]),
#   cbind("1 - pValue", FDA.comb.simInfo[,10])
# ) %>% as.data.frame() 
# colnames(plot.dat) = c("class","sim_val")
# plot.dat[,1] = as.character(plot.dat[,1])
# plot.dat[,2] = as.double(plot.dat[,2])
# plot.dat = plot.dat %>% as_tibble()
# 
# ggplot(plot.dat, aes(x=class, y=sim_val, fill=class)) + 
#   geom_boxplot(alpha=0.3) +
#   theme(legend.position="none") +
#   scale_fill_brewer(palette="Dark2") +
#   ylab("") + 
#   xlab("") +
#   ggtitle("COVID-19 drug similarities that are approved to treat in combination by FDA") +
#   theme_classic()
# 



# load DDI -----
ddi = fread("data/13397 DDI from Drugbank.txt") %>% as.data.frame()
countDDI = function(ddi, aDF){
  ddi.simInfo1 = dplyr::inner_join(ddi, aDF, by = c("ID1" = "ID1", "ID2" = "ID2"))
  ddi.simInfo2 = dplyr::inner_join(ddi, aDF, by = c("ID2" = "ID1", "ID1" = "ID2"))
  ddi.simInfo = dplyr::bind_rows(ddi.simInfo1, ddi.simInfo2)
  # remove dupliate pairs
  ddi.simInfo[1:2] = t(apply(ddi.simInfo[1:2],1,sort))
  ddi.simInfo = ddi.simInfo[!duplicated(ddi.simInfo[1:2]),]
  return(nrow(ddi.simInfo))
}

# permutaion test for DDI check
obsVal = countDDI(ddi = ddi, aDF = drugSim.db)

all.cv = fread("data/COVID drug collection/allCOVID_drug_collection_v2.csv") %>% as.data.frame()
all.cv.comb = combn(all.cv[,1], 2, simplify = T)
x <- data.frame(x=all.cv[,1])
y <- data.frame(y=all.cv[,1])
all.cv.pairs = tidyr::crossing(x, y)
#remove duplicates
all.cv.pairs[1:2] = t(apply(all.cv.pairs[1:2],1,sort))
all.cv.pairs = all.cv.pairs[!duplicated(all.cv.pairs[1:2]),]
colnames(all.cv.pairs) = c("ID1","ID2")

nIter = 100000
cnt = 0
valList = c()
for(i in 1:nIter){
  keep = sample(seq(1, nrow(all.cv.pairs), by = 1), nrow(drugSim.db), replace = F)
  aDF = all.cv.pairs[keep,]
  aVal = countDDI(ddi = ddi, aDF = aDF)
  valList = c(valList, aVal)
  if(aVal <= obsVal) 
    cnt = cnt + 1
  print(i)
}
pVal = cnt/nIter
sd = sqrt(var(valList))
data <- data.frame(
     name=c("COV_simDB","Random"),
     value=c(obsVal, mean(valList)),
     sd=c(NA,sd)
 )
ggplot(data) +
  geom_bar( aes(x=name, y=value), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=name, ymin=value-sd, ymax=value+sd), width=0.4, colour="orange", alpha=0.9, size=1.3) +
  ylab("Number of drug pairs") + 
  xlab("") +
  ggtitle(paste0("Adverse drug-drug interactions (P-value: ", pVal, ")")) +
  theme_classic()
