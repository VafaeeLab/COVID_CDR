get_COV2_drug_net <- function(d1, 
                              drug.targets, 
                              to_genes, 
                              ppiNet, 
                              hop, 
                              mindist, 
                              cov2_human.ppi,
                              ppiNet.edgeDF){
  require(dplyr)
  require(data.table)
  source("getPPIPartners.R")
  
  
  nodes = data.frame()
  edges = data.frame()
  res = NA
  
  d1.targets = drug.targets[which(grepl(d1, drug.targets$`Drug IDs`, fixed = T)),] %>% 
    # dplyr::filter(Species == "Humans" | Species == "SARS-CoV2" | Species == "Influenza A virus")  
    dplyr::filter(Species == "Humans" | Species == "SARS-CoV2")  
  # print(d1.targets)
  if(nrow(d1.targets) > 0){
    d1.targets = d1.targets[,1]
    from_genes = d1.targets
    
    nodes = data.frame(id=d1.targets, 
                       label=d1.targets, 
                       # image = "https://ndownloader.figshare.com/files/25144634",
                       # shape = "image",
                       color = "blue",
                       title = paste0("<p>Name: <b>", d1.targets, "</b></p>"),
                       stringsAsFactors = F)
    
    
    edges = data.frame( # connect drug1 with its direct targets
      V1 = rep(d1, length(d1.targets)),
      V2 = d1.targets
    )
    temp <- getPPIpartners(ppi.net = ppiNet, geneList = from_genes, hop=hop, mindist = mindist) %>%  # mindist: 0 would include Human pray genes as direct drug-targets
      na.omit() %>% 
      intersect(to_genes) # only keep neighbours of drug targets that are Human Pray genes
    
    if(length(temp) > 0){
      # set nodes for drug1 ---------------------
      nodes = nodes %>% bind_rows(
        data.frame(id=temp, label=temp, 
                   # image = "https://ndownloader.figshare.com/files/25144628",
                   # shape = "image",
                   
                   color = "green",
                   title = paste0("<p>Name: <b>", temp, "</b></p>"),
                   stringsAsFactors = F)
      ) %>%  bind_rows(
        data.frame(id=cov2_human.ppi[which(cov2_human.ppi$Human %in% temp),1] %>% unique(),
                   label=cov2_human.ppi[which(cov2_human.ppi$Human %in% temp),1] %>% unique(), 
                   # image = "https://ndownloader.figshare.com/files/25144658",
                   # shape = "image",
                   
                   color = "red",
                   title = paste0("<p>Name: <b>", cov2_human.ppi[which(cov2_human.ppi$Human %in% temp),1] %>% unique(), "</b></p>"),
                   stringsAsFactors = F)
      )
      # set edges for drug1 ---------------------
      edges = edges %>% bind_rows( # connect drug targets with human pray genes (duplication will be there, which will be handled later)
        bind_rows(
          ppiNet.edgeDF[which(ppiNet.edgeDF$V1 %in% d1.targets & ppiNet.edgeDF$V2 %in% temp),],
          ppiNet.edgeDF[which(ppiNet.edgeDF$V2 %in% d1.targets & ppiNet.edgeDF$V1 %in% temp),]
        )
      )
      temp4 = cov2_human.ppi[which(cov2_human.ppi$Human %in% temp),]
      colnames(temp4) = c("V1","V2")
      edges = bind_rows(edges, temp4)
    }
  }
  

  return(list(nodes, edges))
}