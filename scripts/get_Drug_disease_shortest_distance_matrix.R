library(igraph)
library(data.table)
library(dplyr)
library(readxl)
library(pbapply)

all.cv_full <- fread("AppData/allCOVID_drug_collection_full_v5.csv", encoding = "UTF-8") %>% as.data.frame()

human.ppi <- fread("AppData/i2d.human.anno.ppi.Genes.csv") %>% as.data.frame() # can look for only the signaling net
ppiNet <- graph_from_data_frame(human.ppi[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

# load COVID-19-Human PPIs
cov2_human.ppi <- fread("AppData/SARS_COV2_Human_PPI_27Oct_2020.csv") %>% as.data.frame() 
colnames(cov2_human.ppi) <- c("COV","Human")
cov2_human.ppi.net <- graph_from_data_frame(cov2_human.ppi, directed = F) %>% igraph::simplify() 

full.net <- igraph::union(ppiNet, cov2_human.ppi.net)


drug.targets.all <- fread("AppData/all_compiled.csv") %>% as.data.frame() 
drug.targets <- drug.targets.all %>% dplyr::select(c(3, 12, 13))

dGenes = cov2_human.ppi$COV %>% unique() %>% as.character()

getShortestDistance <- function(drug, net, targets, disease){
  
  # t <-  targets[which(grepl(drug, targets$`Drug IDs`, fixed = T)),] %>% dplyr::filter(Species == "Humans") 
  t <-  targets[which(grepl(drug, targets$`Drug IDs`, fixed = T)),]
  keep <- which(t$`Gene Name` %in% V(net)$name)
  t <- unique(t[keep,1])
  
  if(length(t) > 0){
  dt = data.frame(from = drug,
                  to = t,
                  stringsAsFactors = F)
  
  dt.net = graph_from_data_frame(dt, directed = F) %>% igraph::simplify()
  net = igraph::union(net, dt.net)
  
  d <- disease
  keep <- which(d %in% V(net)$name)
  d <- unique(d[keep])
  
  return(shortest.paths(net, v = drug, to=d))  
  }else{
    return(NA)
  }
}


mat = do.call("rbind", pbsapply(all.cv_full$DrugBankID %>% unique(), function(aDrug){
    getShortestDistance(drug = aDrug, 
                        net = full.net, 
                        targets = drug.targets, 
                        disease = dGenes)
}))
drugNames = rownames(mat)
diseaseGeneList = colnames(mat)
mat = mat %>% as.data.frame()
colnames(mat) = diseaseGeneList
rownames(mat) = drugNames

write.csv(mat , file = "AppData/drug_SARS_COV2_shortest_distance_V5.csv", row.names = T)
