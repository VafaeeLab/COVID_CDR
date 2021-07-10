# Author: A K M Azad
# This code produces z-scores (and p-values) for drug-disease proximity based on FR (functional relevance) measure
# key-param: "Net": Human-Human (I2D + more) and Human-SARS_Cov2, 
# drug-set: Drugs of interest, 
# drug-target network (includes new annotations)


library(igraph)
library(readxl)
library(dplyr)
library(data.table)
require(pbapply)

# Human-Human PPI
human.ppi <- fread("AppData/i2d.human.anno.ppi.Genes.csv", encoding = "UTF-8") %>% 
  dplyr::select(c("symbol1", "symbol2")) %>% 
  as.data.frame() 
ppiNet <- graph_from_data_frame(human.ppi, directed = F) %>% igraph::simplify()

# Human-SARS_COV2 ppi
cov2_human.ppi <- fread("AppData/SARS_COV2_Human_PPI_24Sep2020.csv", encoding = "UTF-8") %>% as.data.frame() 
colnames(cov2_human.ppi) <- c("COV","Human")
cov2_human.ppi.net <- graph_from_data_frame(cov2_human.ppi, directed = F) %>% igraph::simplify()

# extra ppi data (includes extra Human-human, 
# Human-SARS_Cov2, 
# Influenza-Human PPI (from String.DB) and 
# SARS_COV2-SARS_COV2 PPIs from literature)
ppiex.dat <- fread("./experimental data processing/new_ppi_data.csv", encoding = "UTF-8") %>% as.data.frame()
ex.ppi.net <- graph_from_data_frame(ppiex.dat, directed = F) %>% igraph::simplify()

drug_targets <- fread("AppData/all_compiled.csv") %>% as.data.frame() %>% dplyr::select(c(3, 12, 13))
dGenes = cov2_human.ppi$COV %>% unique() %>% as.character()

allCV.drugs <- fread("./experimental data processing/Drugs_for_experiment_check.csv") %>% 
  dplyr::select(1) %>% 
  unlist(use.names = F) %>% 
  unique()

full.net <- igraph::union(ppiNet, cov2_human.ppi.net, ex.ppi.net)

getRandD <- function(a_degree, net){
  return(sample(which(degree(net) == a_degree),1, replace = F))
}
permuteTest <- function(net, t, d, d_td, N){
  r <- c()
  for (i in 1:N) {
    t_rand <- sapply(as.numeric(degree(net, v = t)), getRandD,net)
    d_rand <- sapply(as.numeric(degree(net, v = d)), getRandD,net)
    d_td_rand <- get_fr(shortest.paths(net, v = t_rand %>% unique(), to=d_rand %>% unique()))
    r<-c(r,d_td_rand)
  }
  r[!is.finite(r)] <- NA
  m <- mean(r, na.rm = T)
  s <- sd(r, na.rm = T)
  z <- (d_td - m)/s
  return(z)
}


get_fr <- function(dist_matrix){
  epsilon <- 1e-3
  fr <- epsilon
  for(i in 1:length(dist_matrix)){
    if(!is.na(dist_matrix[i])){
      # sum = sum + exp(-(dist_matrix[i,j] + rank_mat[[j]]))
      fr = fr + exp(-(dist_matrix[i]))
    }
  }
  return(fr)
}

getzScore <- function(drug, net, targets, disease){
  
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
    
    d_td <- get_fr(shortest.paths(net, v = t, to=d))  
    
    z <- permuteTest(net, t, d, d_td, 1000)
    p <- pnorm(-abs(z))
    
    return(cbind(drug, d_td, z, p))
  }else{
    return(NA)
  }
}

# dat <- data.frame()

library(foreach)
library(doParallel)
library(parallel)
cl <- makeCluster(10)
registerDoParallel(cl)

start_time <- Sys.time()
# for(i in 1:length(allCV.drugs)){ 
#     dat = dat %>% rbind(.,  getzScore(allCV.drugs[i], full.net, drug_targets, dGenes))
# }
dat <- foreach(i = 1:length(allCV.drugs), .combine = 'rbind', .packages = c('igraph','dplyr'))  %dopar% {
  getzScore(allCV.drugs[i], full.net, drug_targets, dGenes)
}

fwrite(dat, file = paste0("./experimental data processing/FR_based_disease_proximity_selected_drugs.csv"), row.names = F, na = "NA")
end_time <- Sys.time()
print(end_time - start_time)
stopCluster(cl)
