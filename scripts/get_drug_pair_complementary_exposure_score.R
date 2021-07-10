get_Complementary_Exposure_score <- function(drugA, drugB, net, targets){
  require(igraph)
  require(dplyr)
  
  tA <-  targets[which(grepl(drugA, targets$`Drug IDs`, fixed = T)),] %>% dplyr::filter(Species == "Humans")
  keep <- which(tA$`Gene Name` %in% V(net)$name)
  tA <- unique(tA[keep,1])
  
  tB <-  targets[which(grepl(drugB, targets$`Drug IDs`, fixed = T)),] %>% dplyr::filter(Species == "Humans") 
  keep <- which(tB$`Gene Name` %in% V(net)$name)
  tB <- unique(tB[keep,1])
  
  dAB <- mean(shortest.paths(net, v = tA, to=tB))
  dAA <- mean(shortest.paths(net, v = tA, to=tA))
  dBB <- mean(shortest.paths(net, v = tB, to=tB))
  sAB <- dAB - (dAA+dBB)/2
  if(is.nan(sAB))
    sAB <- NA
  return(sAB)
}

# library(igraph)
# library(readxl)
# library(dplyr)
# library(data.table)
# 
# human.ppi <- fread("AppData/i2d.human.anno.ppi.Genes.csv") %>% as.data.frame() # can look for only the signaling net
# ppiNet <- graph_from_data_frame(human.ppi[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()
# cov2_human.ppi <- read_excel("AppData/media-6.xlsx") %>% as.data.frame() %>% dplyr::select(c(1,3))
# colnames(cov2_human.ppi) <- c("COV","Human")
# filename="all"
# drug.targets.all <- fread(paste0("AppData/", filename, ".csv")) %>% as.data.frame() 
# drug.targets <- drug.targets.all %>% dplyr::select(c(3, 12, 13))
# 
# d1 = "DB00811"
# d2 = "DB00811"
# # dGenes = cov2_human.ppi$Human %>% as.character()
# # 
# get_Complementary_Exposure_score(drugA = d1, drugB = d2, net = ppiNet, targets = drug.targets)
