get_Infz_drug_net <- function(d1,
                              prev_nodes,
                              Infz_human.ppi,
                              drug.targets,
                              drug.targets.all,
                              Infz_human.ppi.net,
                              ppiNet,
                              hop,
                              mindist,
                              ppiNet.edgeDF){
  require(dplyr)
  require(data.table)
  require(igraph)
  source("getPPIPartners.R")
  
  
  
  edges = data.frame()
  nodes = data.frame()
  res = NA
  
  d1.targets = drug.targets[which(grepl(d1, drug.targets$`Drug IDs`, fixed = T)),] %>% 
    dplyr::filter(Species == "Influenza A virus" | Species == "Humans") 
  
  d1.targets.Infz.gs <- d1.targets[which(d1.targets$Species == "Influenza A virus"),1] %>% unique() 
  d1.targets.Infz.dgs <- getPPIpartners(ppi.net = Infz_human.ppi.net, geneList = d1.targets.Infz.gs, hop = hop, mindist = (mindist+1)) %>% as.character()
  d1.targets.Infz.dgs.targets <- getPPIpartners(ppi.net = ppiNet, geneList = d1.targets.Infz.dgs, hop = hop, mindist = (mindist+1)) %>% 
    as.character() %>% 
    base::setdiff(.,d1.targets.Infz.gs)
  # print(d1.targets.Infz.dgs.targets) # 
  
  
  if(length(d1.targets.Infz.dgs.targets) > 0){
    
    d1.targets = d1.targets[,c(1,2,3)]    # GeneName, uniprotID, species, 
    if(length(which(d1.targets$Species == "Influenza A virus")) > 0){
      nodes = nodes %>% bind_rows(
        data.frame(id=d1.targets.Infz.gs, # gene symbol as ID
                   label=d1.targets.Infz.gs, # gene symbol as Label
                   
                   # image = "https://ndownloader.figshare.com/files/25144634",
                   # shape = "image",

                   color = "violet",
                   title = paste0("<p>Name: <b>", d1.targets.Infz.gs, "</b></p>"),
                   stringsAsFactors = F)
      ) %>% bind_rows(
      # nodes = nodes %>% bind_rows(
        data.frame(id=d1.targets.Infz.dgs, # gene symbol as ID
                   label=d1.targets.Infz.dgs, # gene symbol as Label
                   # image = "https://ndownloader.figshare.com/files/25144631",
                   # shape = "image",
                   
                   color = "green",
                   title = paste0("<p>Name: <b>", d1.targets.Infz.dgs, "</b></p>"),
                   stringsAsFactors = F)
      ) %>% bind_rows(
        data.frame(id=d1.targets.Infz.dgs.targets, # gene symbol as ID
                   label=d1.targets.Infz.dgs.targets, # gene symbol as Label
                   # image = "https://ndownloader.figshare.com/files/25144628",
                   # shape = "image",
                   
                   color = "green",
                   title = paste0("<p>Name: <b>", d1.targets.Infz.dgs.targets, "</b></p>"),
                   stringsAsFactors = F)
      )
      edges = edges %>% bind_rows(
        data.frame( # connect drug1 with its direct targets
          V1 = d1,
          V2 = d1.targets.Infz.gs
        )
      ) %>% bind_rows(
        data.frame( # connect drug1 with its direct targets
          V1 = d1.targets.Infz.gs,
          V2 = d1.targets.Infz.dgs
        )
      ) 
      
      h.ppi.portion = ppiNet.edgeDF[which(((ppiNet.edgeDF[,1] %in% d1.targets.Infz.dgs) & 
                                            (ppiNet.edgeDF[,2] %in% d1.targets.Infz.dgs.targets)) | 
                                            ((ppiNet.edgeDF[,2] %in% d1.targets.Infz.dgs) & 
                                               (ppiNet.edgeDF[,1] %in% d1.targets.Infz.dgs.targets))),]
      # print(edges)
      edges = edges %>% bind_rows(h.ppi.portion)
    }
    # temp1 <- getPPIpartners(ppi.net = Infz_human.ppi.net, geneList = d1.targets[,1], hop=1, mindist = 1)   
    # print(d1.targets.Infz.gs.targets)
    # print(getPPIpartners(ppi.net = ppiNet, geneList = d1.targets.Infz.dgs.targets, hop=hop, mindist = mindist) %>% unique())
    prev.prot.reachable.from.d1.targets.Infz.dgs.targets <- getPPIpartners(ppi.net = ppiNet, geneList = d1.targets.Infz.dgs.targets, hop=(hop+1), mindist = mindist) %>% unique() %>%
      intersect(., prev_nodes)
    # print(prev.prot.reachable.from.d1.targets.Infz.dgs.targets)
    if(length(prev.prot.reachable.from.d1.targets.Infz.dgs.targets) > 0){
      # print(prev.prot.reachable.from.d1.targets.Infz.dgs.targets)
      
      nodes = nodes  %>%  bind_rows(
        data.frame(id = prev.prot.reachable.from.d1.targets.Infz.dgs.targets,
                   label=prev.prot.reachable.from.d1.targets.Infz.dgs.targets, 
                   # image = "https://ndownloader.figshare.com/files/25144628",
                   # shape = "image",
                   
                   color = "green",
                   title = paste0("<p>Name: <b>", prev.prot.reachable.from.d1.targets.Infz.dgs.targets, "</b></p>"),
                   stringsAsFactors = F)
      )
      
      h.ppi.subnet.edgeDF <- induced.subgraph(graph = ppiNet, vids = V(ppiNet)[V(ppiNet)$name %in% c(prev.prot.reachable.from.d1.targets.Infz.dgs.targets, d1.targets.Infz.dgs.targets)]) %>%
        get.edgelist() %>% as.data.frame()
      colnames(h.ppi.subnet.edgeDF) = c("V1","V2")
      
      edges = bind_rows(edges, h.ppi.subnet.edgeDF)
      
    }
    
    edges = edges[,c(1,2)]
  }
  
  return(list(nodes, edges))
}