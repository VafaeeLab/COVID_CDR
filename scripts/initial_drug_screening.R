library(readxl)
library(igraph)
library(data.table)
library(dplyr)
library(stringr)
source("getPPIPartners.R")

initial_Drug_screeing <- function(ppiNet = "signalink", hop = 2, ppiSupport = T){
  # load markers and PPI net
  covProt.dat <- read_excel("data/media-6.xlsx") %>% as.data.frame() %>% select(c(2,3,8))
  cov.genes <- covProt.dat$PreyGene
  length(unique(cov.genes))
  if(ppiNet == "wholePPI"){
    human.ppi <- fread("data/i2d.human.anno.ppi.Genes.csv") %>% as.data.frame() # can look for only the signaling net
    ppiNet <- graph_from_data_frame(human.ppi[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()
  }else{
    sig = fread("data/signalink_1589702515711.csv") %>% as.data.frame()
    ppiNet = graph.data.frame(sig[,c("source_name","target_name")], directed = T) %>% simplify()
  }
  # drug ---
  drug.info <- fread("data/Drugbank/5.1.6/drug links.csv") %>% as.data.frame() %>% select(c(1,2))
  # with all drugs
  filename="all"
  drug.targets <- fread(paste0("data/Drugbank/drug-target mapping/", filename, ".csv")) %>% as.data.frame() %>% select(c(2,3,6,7, 13))
  cov.genes.drug.init <- dplyr::inner_join(covProt.dat, drug.targets, by = c("PreyGene" = "Gene Name"))
  cov.drugs.init <- str_split(cov.genes.drug.init$`Drug IDs`, "; ") %>% unlist() %>% unique()
  cov.drugs.init <- drug.info[match(cov.drugs.init, drug.info$`DrugBank ID`),]
  # fwrite(cov.drugs.init, file = paste0("data/cov_init_drugs_", filename, ".csv"))
  
  # explore ppi partners of prots that not have direct Drugs; params (hop = 1 [default])
  if(ppiSupport){
    miss.prots = base::setdiff(cov.genes, cov.genes.drug.init$PreyGene)
    miss.prots.ppi.partners = getPPIpartners(ppi.net = ppiNet, geneList = miss.prots, hop=hop) %>% na.omit()
    cov.drugs.ppi.init <- NULL
    if(length(miss.prots.ppi.partners) > 0){
      miss.prots.drugs = dplyr::inner_join(data.frame(PreyGene = miss.prots.ppi.partners), 
                                           drug.targets, 
                                           by = c("PreyGene" = "Gene Name"))
      cov.drugs.ppi.init <- str_split(miss.prots.drugs$`Drug IDs`, "; ") %>% unlist() %>% unique()
      cov.drugs.ppi.init <- drug.info[match(cov.drugs.ppi.init, drug.info$`DrugBank ID`),]
      cov.drugs.init = rbind(cov.drugs.init, cov.drugs.ppi.init)
    }
    
  }
  return(cov.drugs.init)
}
# # with only active drus
# filename = "pharmacologically_active"
# drug.targets <- fread(paste0("data/Drugbank/drug-target mapping/", filename, ".csv")) %>% as.data.frame() %>% select(c(2,3,6,7, 13))
# cov.genes.drug.init <- dplyr::inner_join(covProt.dat, drug.targets, by = c("PreyGene" = "Gene Name"))
# cov.drugs.init <- str_split(cov.genes.drug.init$`Drug IDs`, "; ") %>% unlist() %>% unique()
# cov.drugs.init <- drug.info[match(cov.drugs.init, drug.info$`DrugBank ID`),]
# fwrite(cov.drugs.init, file = paste0("data/cov_init_drugs_", filename, ".csv"))
