library(data.table)
library(dplyr)
source("C:\\Users\\Azad\\OneDrive - UNSW\\Howell lab [kolling]\\Amanda_Emily\\bwl_ae_pathwayEnrichment\\Pathway_geneSetEnrichmentAnalysis.R")
pop.file.dir <- "C:\\Users\\Azad\\OneDrive - UNSW\\Howell lab [kolling]\\Amanda_Emily\\bwl_ae_pathwayEnrichment\\AppData\\Enrichr\\"

dt <- fread("AppData/all.csv") %>% as.data.frame()
all.cv <- fread("AppData/allCOVID_drug_collection_full.csv") %>% as.data.frame()
# for(i in 1:nrow(all.cv)){
mybiglist <- list()

for(i in 1:nrow(all.cv)){
  print(i)
  res <- NULL
  aDrug = all.cv$DrugBankID[i]
  t = dt[which(grepl(aDrug, dt$`Drug IDs`, fixed = T)),] %>% 
    dplyr::filter(Species == "Humans" | Species == "SARS-CoV2") %>% 
    dplyr::select(3) %>% unlist(use.names = F)
  
  # KEGG
  fileName = "KEGG_2019.txt"
  aRes1 <- Pathway_geneSetEnrichmentAnalysis(givenSet = t, pop.filepath = paste0(pop.file.dir,fileName))

    # # Reactome
  # fileName = "Reactome_2016.txt"
  # aRes2 <- Pathway_geneSetEnrichmentAnalysis(givenSet = t, pop.filepath = paste0(pop.file.dir,fileName))

  # WikiPathway
  fileName = "WikiPathways_2019.txt"
  aRes3 <- Pathway_geneSetEnrichmentAnalysis(givenSet = t, pop.filepath = paste0(pop.file.dir,fileName))

  # Biocarta
  fileName = "BioCarta_2016.txt"
  aRes4 <- Pathway_geneSetEnrichmentAnalysis(givenSet = t, pop.filepath = paste0(pop.file.dir,fileName))

  # # GO_BP
  # fileName = "GO_Biological_Process_2018.txt"
  # aRes5 <- Pathway_geneSetEnrichmentAnalysis(givenSet = t, pop.filepath = paste0(pop.file.dir,fileName))
  # 
  # # GO_CC
  # fileName = "GO_Cellular_Component_2018.txt"
  # aRes6 <- Pathway_geneSetEnrichmentAnalysis(givenSet = t, pop.filepath = paste0(pop.file.dir,fileName))
  # 
  # # GO_MF
  # fileName = "GO_Molecular_Function_2018.txt"
  # aRes7 <- Pathway_geneSetEnrichmentAnalysis(givenSet = t, pop.filepath = paste0(pop.file.dir,fileName))

  # res <- list(aRes1, aRes2, aRes3, aRes4, aRes5, aRes6, aRes7)
  res = list(aRes1, aRes3, aRes4)
  mybiglist[[aDrug]] <- res
}

saveRDS(mybiglist, file="AppData/all_Drugs_enrichments_small.rds")
