get_SARS_drug_net <- function(d1,
                              prev_nodes,
                              SARS_cov_human.ppi,
                              drug.targets,
                              drug.targets.all,
                              SARS_cov_human.ppi.net,
                              ppiNet,
                              hop,
                              mindist,
                              ppiNet.edgeDF){
  require(dplyr)
  require(data.table)
  require(igraph)
  source("getPPIPartners.R")
  source("get_induced_bipartite_graph.R")
  
  
  
  edges = data.frame()
  nodes = data.frame()
  res = NA
  
  d1.targets = drug.targets[which(grepl(d1, drug.targets$`Drug IDs`, fixed = T)),] %>% 
    dplyr::filter(Species == "SARS-CoV" | Species == "Humans") 
  
  d1.targets.SARS.gs <- d1.targets[which(d1.targets$Species == "SARS-CoV"),1] %>% unique() # rep
  d1.targets.SARS.dgs <- getPPIpartners(ppi.net = SARS_cov_human.ppi.net, 
                                        geneList = d1.targets.SARS.gs, 
                                        hop = hop, 
                                        mindist = (mindist+1)
  ) %>% as.character() 
  
  d1.targets.SARS.dgs.targets <- getPPIpartners(ppi.net = SARS_cov_human.ppi.net, geneList = d1.targets.SARS.dgs, hop = hop, mindist = (mindist+1)) %>% 
    as.character() %>% 
    base::setdiff(.,d1.targets.SARS.gs)
  
  print(d1.targets.SARS.gs)
  print(d1.targets.SARS.dgs)
  print(d1.targets.SARS.dgs.targets) 
  
  
  if(length(d1.targets.SARS.dgs.targets) > 0){
    
    d1.targets = d1.targets[,c(1,2,3)]    # GeneName, uniprotID, species, 
    if(length(which(d1.targets$Species == "SARS-CoV")) > 0){
      nodes = nodes %>% bind_rows(
        data.frame(id=d1.targets.SARS.dgs, # gene symbol as ID
                   label=d1.targets.SARS.dgs, # gene symbol as Label
                   # image = "https://ndownloader.figshare.com/files/25144631",
                   # shape = "image",
                   
                   color = "sienna4",
                   title = paste0("<p>Name: <b>", d1.targets.SARS.dgs, "</b></p>"),
                   stringsAsFactors = F)
      ) 
      
      edges = edges %>% bind_rows(
        data.frame( 
          V1 = d1,
          V2 = d1.targets.SARS.dgs
        )
      ) 
    }
    
    prev.prot.reachable.from.d1.targets.SARS.dgs.targets <- getPPIpartners(
      ppi.net = ppiNet, 
      geneList = d1.targets.SARS.dgs.targets, 
      hop=hop, 
      mindist = (mindist+1)
    ) %>% unique() %>%
      intersect(., prev_nodes)
    
    if(length(prev.prot.reachable.from.d1.targets.SARS.dgs.targets) > 0){
      h.ppi.subnet.edgeDF <- get_induced_bipartite_graph(ppiNet = ppiNet, 
                                                         setA = d1.targets.SARS.dgs.targets, 
                                                         setB = prev.prot.reachable.from.d1.targets.SARS.dgs.targets) %>%
        get.edgelist() %>% as.data.frame(stringsAsFactors = FALSE)
      
      colnames(h.ppi.subnet.edgeDF) = c("V1","V2")
      
      temp.filter = base::intersect(c(h.ppi.subnet.edgeDF$V1 %>% as.character(),
                                      h.ppi.subnet.edgeDF$V2 %>% as.character()), d1.targets.SARS.dgs.targets)
      
      edgeDF1 = SARS_cov_human.ppi.net %>% get.edgelist() %>% as.data.frame()
      sars.human.subnet.edgeDF = bind_rows(
        edgeDF1[which(edgeDF1$V1 %in% d1.targets.SARS.dgs & edgeDF1$V2 %in% temp.filter),],
        edgeDF1[which(edgeDF1$V2 %in% d1.targets.SARS.dgs & edgeDF1$V1 %in% temp.filter),]
      ) %>% as.data.frame(stringsAsFactors = FALSE)
      colnames(sars.human.subnet.edgeDF) = c("V1","V2")
      
      
      edges = bind_rows(edges, sars.human.subnet.edgeDF, h.ppi.subnet.edgeDF)
      
      nodelist.temp = c(
        sars.human.subnet.edgeDF$V1 %>% as.character(),
        sars.human.subnet.edgeDF$V2 %>% as.character(),
        h.ppi.subnet.edgeDF$V1 %>% as.character(),
        h.ppi.subnet.edgeDF$V2 %>% as.character()) %>% base::unique()
      
      nodes = nodes %>% bind_rows(
        data.frame(
          id = nodelist.temp,
          label = nodelist.temp,
          color = "green",
          title = paste0("<p>Name: <b>", nodelist.temp, "</b></p>"),
          stringsAsFactors = F)
      )
    }
    
    edges = edges[,c(1,2)]
  }
  return(list(nodes, edges))
}