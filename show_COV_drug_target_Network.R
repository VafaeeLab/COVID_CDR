


show_COV_drug_target_Network <- function(drugs, hasSARS_CoV = T, hasInfz = T, hop = 1, mindist = 1){
  require(data.table)
  require(dplyr)
  require(ggplot2)
  require(igraph)
  require(visNetwork)
  require(readxl)
  source("get_COV2_drug_net.R")
  source("get_SARS_drug_net.R")
  source("get_Infz_drug_net.R")
  
  
  all.cv_full <- fread("AppData/allCOVID_drug_collection_full_v5.csv", encoding = "UTF-8") %>% as.data.frame()
  
  
  human.ppi <- fread("AppData/i2d.human.anno.ppi.Genes.csv", encoding = "UTF-8") %>% as.data.frame() # can look for only the signaling net
  ppiNet <- graph_from_data_frame(human.ppi[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()
  
  # add side-effect
  ddi = fread("AppData/13397 DDI from Drugbank.txt", encoding = "UTF-8") %>% as.data.frame()
  # ddi.simInfo1 = dplyr::inner_join(ddi, FDA_aprv, by = c("ID1" = "ID1", "ID2" = "ID2"))
  # ddi.simInfo2 = dplyr::inner_join(ddi, FDA_aprv, by = c("ID2" = "ID1", "ID1" = "ID2"))
  # ddi.simInfo = dplyr::bind_rows(ddi.simInfo1, ddi.simInfo2)
  # # remove dupliate pairs
  # ddi.simInfo[1:2] = t(apply(ddi.simInfo[1:2],1,sort))
  # ddi.simInfo = ddi.simInfo[!duplicated(ddi.simInfo[1:2]),]
  
  
  # load COVID-19-Human PPIs
  cov2_human.ppi <- fread("AppData/SARS_COV2_Human_PPI_27Oct_2020.csv", encoding = "UTF-8") %>% as.data.frame() 
  colnames(cov2_human.ppi) <- c("COV","Human")
  cov2_human.ppi.net <- graph_from_data_frame(cov2_human.ppi, directed = F) %>% igraph::simplify()
  
  
  # 
  
  
  # filename=
  drug.targets.all <- fread("AppData/all_compiled.csv", encoding = "UTF-8") %>% as.data.frame() 
  drug.targets <- drug.targets.all %>% dplyr::select(c(3, 12, 13))
  
  # a pair:
  to_genes = cov2_human.ppi$Human
  ppiNet.edgeDF = get.edgelist(ppiNet) %>% as.data.frame()
  
  # drugNames = drug.info[match(drugs, drug.info$`DrugBank ID`),2]
  drugNames = all.cv_full[match(drugs, all.cv_full$DrugBankID), 2]
  drug.info.TTD = all.cv_full %>% dplyr::select(c("Name", "Therap_Class", "Indication(s)"))
  
  # drug.info.TTD = fread("AppData/TTD_drugName_TherapClass_Indication_Mapping.txt", encoding = "UTF-8") %>% as.data.frame()
  drugClass = drug.info.TTD[base::match(drugNames, drug.info.TTD$Name), "Therap_Class"]
  drugIndication = drug.info.TTD[base::match(drugNames, drug.info.TTD$Name), "Indication(s)"]
  
  disease.proximity.drugs = all.cv_full %>% 
    dplyr::select(c("DrugBankID","z_DP", "z_DP_pVal")) %>% 
    dplyr::rename("DrugID" = "DrugBankID", "z_score" = "z_DP", "p_value" = "z_DP_pVal")
  
  dd.prox.score = disease.proximity.drugs[base::match(drugs, disease.proximity.drugs$DrugID),c("z_score","p_value")] 
  
  functional.proximity.drugs = all.cv_full %>% 
    dplyr::select(c("DrugBankID","Functional_Proximity"))
  
  fp.score = functional.proximity.drugs[base::match(drugs, functional.proximity.drugs$DrugBankID), c("z_FP","z_FP_pVal")] 
  
  
  
  nodes = data.frame(
    id=drugs,
    label=drugNames, 
    # image = "https://ndownloader.figshare.com/files/25144541",
    # shape = "image",
    color = rep("pink",length(drugNames)),
    title = paste0("<p>Name: <b>", drugNames, "</b></p>", 
                   "<p>Class: <b>", drugClass, "</b></p>", 
                   "<p>Indications: <b>", drugIndication, "</b></p>",
                   "<p>COVID-19 Proximity (functional): <b>", fp.score[1] %>% 
                     round(digits = 2), "(" , fp.score[2] %>% 
                     round(digits = 2), ")",  "</b></p>"),
                    "<p>COVID-19 Proximity (topological): <b>", dd.prox.score[1] %>% 
                     round(digits = 2), "(" , dd.prox.score[2] %>% 
                     round(digits = 2), ")",  "</b></p>",
    stringsAsFactors = F)
  # nodes = data.frame(
  #   id=c(d1,d2),
  #   label=c(d1,d2), image = drug_image)
  
  edges = data.frame()
  # check if DDI exists between d1 and d2
  testDF = ddi[which((ddi$ID1 %in% drugs) &
                       ddi$ID2 %in% drugs),]
  if(nrow(testDF) > 0){
    edges = edges %>% bind_rows(
      testDF
    )
  }
  # get COV2 net for all drugs
  for(i in 1:length(drugs)){
    aDrug = drugs[i]
    result <- get_COV2_drug_net(aDrug, drug.targets, to_genes, ppiNet, hop = hop, mindist = mindist, cov2_human.ppi, ppiNet.edgeDF)
    
    if(nrow(result[[1]]) > 0){
      nodes = nodes %>% bind_rows(result[[1]])
    }
    if(nrow(result[[2]]) > 0){
      edges = edges %>% bind_rows(result[[2]])
    }
  }
  
  prev_nodes = nodes[,1]
  # print(prev_nodes)
  
  # load SARS-COV-Human PPIs
  if(hasSARS_CoV){
    library(readxl)
    library(dplyr)
    SARS_uniprot <- read_excel("AppData/SARS_Uniprot.xlsx") %>% dplyr::select(1,5) %>% as.data.frame()
    
    SARS_cov_human.ppi = read_excel("AppData/SARSCOV-Human PPI_SF09072020a.xlsx") %>% 
      as.data.frame() %>% dplyr::select(c(1,3,2)) %>% na.omit()
    colnames(SARS_cov_human.ppi) <- c("Sars_COV_dgs","Sars_COV_Uniprot","Human")
    SARS_cov_human.ppi$Sars_COV_dgs <- SARS_cov_human.ppi$Sars_COV_dgs %>% toupper()
    
    # remdesivir selection -------
    SARS_cov_human.ppi <- SARS_cov_human.ppi %>% dplyr::filter(Sars_COV_dgs %in% c("SARS-COV NSP12", "SARS-COV NSP14"))
    # ---------
    
    SARS_cov_human.ppi = dplyr::inner_join(SARS_cov_human.ppi, SARS_uniprot, by=c("Sars_COV_Uniprot" = "Entry"))
    SARS_cov_human.ppi.net1 = graph_from_data_frame(SARS_cov_human.ppi[,c("Sars_COV_dgs","Gene Names")], directed = F) %>% igraph::simplify()
    SARS_cov_human.ppi.net2 = graph_from_data_frame(SARS_cov_human.ppi[,c("Sars_COV_dgs","Human")], directed = F) %>% igraph::simplify()
    SARS_cov_human.ppi.net = igraph::union(SARS_cov_human.ppi.net1, SARS_cov_human.ppi.net2) %>% igraph::simplify()
    # SARS_cov_human.ppi.net = igraph::union(SARS_cov_human.ppi.net, ppiNet)
    SARS_cov_human.ppi.net.edgeDF = get.edgelist(SARS_cov_human.ppi.net) %>% as.data.frame()
    # to_genes = SARS_cov_human.ppi$Human
    
    # drug.targets <- drug.targets.all %>% dplyr::select(c(6, 12, 13)) # uniprotID, species, drugbankIDs
    drug.targets <- drug.targets.all %>% dplyr::select(c(3, 6, 12, 13)) # GeneName, uniprotID, species, drugbankIDs
    
    # for drug1 ----------------
    for(i in 1:length(drugs)){
      aDrug = drugs[i]
      # print(aDrug)
      result <- get_SARS_drug_net(d1 =  aDrug, prev_nodes =  prev_nodes, 
                                  SARS_cov_human.ppi = SARS_cov_human.ppi, 
                                  drug.targets =  drug.targets, 
                                  drug.targets.all = drug.targets.all, 
                                  SARS_cov_human.ppi.net = SARS_cov_human.ppi.net, 
                                  ppiNet = ppiNet, hop = hop, mindist = mindist, 
                                  ppiNet.edgeDF = ppiNet.edgeDF)
      # print(result[[2]])
      
      if(nrow(result[[1]]) > 0){
        nodes = nodes %>% bind_rows(result[[1]])
      }
      if(nrow(result[[2]]) > 0){
        edges = edges %>% bind_rows(result[[2]])
      }
    }
    
  }
  
  prev_nodes = nodes[,1]
  
  # load Infz-Human PPIs
  if(hasInfz){
    library(readxl)
    library(dplyr)
    Infz_uniprot <- read_excel("AppData/Infz_Uniprot.xlsx") %>% as.data.frame()
    
    Infz_human.ppi = read_excel("AppData/Infz_HumanPPI_26_Oct.xlsx") %>% as.data.frame() 
    colnames(Infz_human.ppi) <- c("Infz_dgs","Infz_Uniprot","Human")
    Infz_human.ppi$Infz_dgs <- Infz_human.ppi$Infz_dgs %>% toupper()
    
    
    Infz_human.ppi = dplyr::inner_join(Infz_human.ppi, Infz_uniprot, by=c("Infz_Uniprot" = "UniProtID"))
    Infz_human.ppi.net1 = graph_from_data_frame(Infz_human.ppi[,c("Infz_dgs","Gene Names")], directed = F) %>% igraph::simplify()
    Infz_human.ppi.net2 = graph_from_data_frame(Infz_human.ppi[,c("Infz_dgs","Human")], directed = F) %>% igraph::simplify()
    Infz_human.ppi.net = igraph::union(Infz_human.ppi.net1, Infz_human.ppi.net2) %>% igraph::simplify()
    Infz_human.ppi.net.edgeDF = get.edgelist(Infz_human.ppi.net) %>% as.data.frame()

    # drug.targets <- drug.targets.all %>% dplyr::select(c(6, 12, 13)) # uniprotID, species, drugbankIDs
    drug.targets <- drug.targets.all %>% dplyr::select(c(3, 6, 12, 13)) # GeneName, uniprotID, species, drugbankIDs
    
    # for drug1 ----------------
    for(i in 1:length(drugs)){
      aDrug = drugs[i]
      # print(aDrug)
      result <- get_Infz_drug_net(d1 =  aDrug, prev_nodes =  prev_nodes, 
                                  Infz_human.ppi = Infz_human.ppi, 
                                  drug.targets =  drug.targets, 
                                  drug.targets.all = drug.targets.all, 
                                  Infz_human.ppi.net = Infz_human.ppi.net, 
                                  ppiNet = ppiNet, hop = hop, mindist = mindist, 
                                  ppiNet.edgeDF = ppiNet.edgeDF)
      # print(result[[2]])
      
      if(nrow(result[[1]]) > 0){
        nodes = nodes %>% bind_rows(result[[1]])
      }
      if(nrow(result[[2]]) > 0){
        edges = edges %>% bind_rows(result[[2]])
      }
    }
    
  }
  
  # distinct nodes
  nodes = nodes[!duplicated(nodes$id),]
  
  
  # distinct edges 
  edges = edges[!duplicated(edges),]
  
  # print(edges)
  
  if(nrow(edges) > 0){
    colnames(edges) = c("from","to")
    edges$id = 1:nrow(edges)
    edges$dashes = F
    # Remdesivir changes ---------------
    edges[which((edges$from == "DB14761" & edges$to == "SARS-COV NSP14") |
                  (edges$to == "DB14761" & edges$from == "SARS-COV NSP14")), 4] = T
    edges[which((edges$from == "DB14761" & edges$to == "SARS-CoV2 nsp14") |
                  (edges$to == "DB14761" & edges$from == "SARS-CoV2 nsp14")), 4] = T
    nodes[which(nodes$label == "SARS-CoV2 nsp12"),"color"] = "red"
    nodes[which(nodes$label == "SARS-CoV2 nsp14"),"color"] = "red"
    # ------
  }
  
  
  
  
  # print(edges)
  # visualize the net
  # visNetwork(nodes,edges) %>%
  #   visNodes(shapeProperties = list(useBorderWithImage = TRUE))
  
  res = NA
  return(list(nodes, edges))
}
# # 
# # # # # # # 
# # d1 = "DB14761"
# # # # d1 = "DB13609"
# d2 = "DB13609"
# # drugs = c(d1,d2)
# # hasSARS_CoV = T
# res <- show_COV_drug_target_Network(c(d2))
# # head(res[0])

# poteintial parameters of this function
# hop, mindist, ppinet (signaling/humanPPI), hasSARSCov
