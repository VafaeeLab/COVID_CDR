library(igraph)
library(data.table)
library(dplyr)

getRandD <- function(a_degree, net){
  return(sample(which(degree(net) == a_degree),1, replace = F))
}
permuteTest <- function(net, t, d, d_td, N){
  r <- c()
  for (i in 1:10) {
    t_rand <- sapply(as.numeric(degree(net, v = t)), getRandD,net)
    d_rand <- sapply(as.numeric(degree(net, v = d)), getRandD,net)
    d_td_rand <- mean(shortest.paths(net, v = t_rand %>% unique(), to=d_rand %>% unique()))
    r<-c(r,d_td_rand)
  }
  r[!is.finite(r)] <- NA
  m <- mean(r, na.rm = T)
  s <- sd(r, na.rm = T)
  z <- (d_td - m)/s
  return(z)
}



getzScore <- function(drug, net, targets, disease){
  
  t <-  targets[which(grepl(drug, targets$`Drug IDs`, fixed = T)),] 
  keep <- which(t$`Gene Name` %in% V(net)$name)
  t <- unique(t[keep,1])
  
  d <- disease
  keep <- which(d %in% V(net)$name)
  d <- unique(d[keep])
  
  d_td <- mean(shortest.paths(net, v = t, to=d))
  
  z <- permuteTest(net, t, d, d_td, 2)
  p <- pnorm(-abs(z))
  
  return(cbind(z, p))
}

human.ppi <- fread("AppData/i2d.human.anno.ppi.Genes.csv") %>% as.data.frame() 
ppiNet <- graph_from_data_frame(human.ppi[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()
cov2_human.ppi <- read_excel("AppData/media-6.xlsx") %>% as.data.frame() %>% dplyr::select(c(1,3))
colnames(cov2_human.ppi) <- c("COV","Human")
# filename="all"
# drug.targets.all <- fread(paste0("AppData/", filename, ".csv")) %>% as.data.frame() 
drug_targets <- fread("AppData/all.csv") %>% as.data.frame() %>% dplyr::select(c(3, 12, 13))
dGenes = cov2_human.ppi$Human %>% as.character()

allCV.drugs <- fread("AppData/allCOVID_drug_collection_v2.csv") %>% dplyr::select(1) %>% unlist(use.names = F)
dat <- do.call("rbind", lapply(allCV.drugs, function(x){
  return(cbind(x,getzScore(drug =  x, net =  ppiNet, targets = drug_targets, disease = dGenes)))
}))
fwrite(dat, file="/srv/scratch/z3526914/COVID_drugSimDB/Data/disease_proximity_scores.csv")