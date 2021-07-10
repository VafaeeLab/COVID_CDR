
show_COV_drug_target_Network <- function(d1, d2, hasSARS_CoV = T){
  require(data.table)
  require(dplyr)
  require(ggplot2)
  require(igraph)
  require(visNetwork)
  require(readxl)
  source("getPPIPartners.R")
  


human.ppi <- fread("AppData/i2d.human.anno.ppi.Genes.csv") %>% as.data.frame() # can look for only the signaling net
ppiNet <- graph_from_data_frame(human.ppi[,c("symbol1", "symbol2")], directed = F) %>% igraph::simplify()

# add side-effect
ddi = fread("AppData/13397 DDI from Drugbank.txt") %>% as.data.frame()
# ddi.simInfo1 = dplyr::inner_join(ddi, FDA_aprv, by = c("ID1" = "ID1", "ID2" = "ID2"))
# ddi.simInfo2 = dplyr::inner_join(ddi, FDA_aprv, by = c("ID2" = "ID1", "ID1" = "ID2"))
# ddi.simInfo = dplyr::bind_rows(ddi.simInfo1, ddi.simInfo2)
# # remove dupliate pairs
# ddi.simInfo[1:2] = t(apply(ddi.simInfo[1:2],1,sort))
# ddi.simInfo = ddi.simInfo[!duplicated(ddi.simInfo[1:2]),]


# load COVID-19-Human PPIs
cov2_human.ppi <- read_excel("AppData/media-6.xlsx") %>% as.data.frame() %>% select(c(1,3))
colnames(cov2_human.ppi) <- c("COV","Human")
cov2_human.ppi.net <- graph_from_data_frame(cov2_human.ppi, directed = F) %>% igraph::simplify()


# 
drug.info <- fread("AppData/drug links.csv") %>% as.data.frame() %>% select(c(1,2))

filename="all"
drug.targets.all <- fread(paste0("AppData/", filename, ".csv")) %>% as.data.frame() 
drug.targets <- drug.targets.all %>% select(c(3, 12, 13))

# a pair:
to_genes = cov2_human.ppi$Human
ppiNet.edgeDF = get.edgelist(ppiNet) %>% as.data.frame()

nodes = data.frame(
  id=c(d1,d2),
  label=c(drug.info[which(drug.info$`DrugBank ID` == d1),2], 
          drug.info[which(drug.info$`DrugBank ID` == d2),2]), 
  color = c("pink","pink"))
# nodes = data.frame(
#   id=c(d1,d2),
#   label=c(d1,d2), image = drug_image)

edges = data.frame()
# check if DDI exists between d1 and d2
testDF = ddi[which((ddi$ID1 == d1 & ddi$ID2 == d2) |
                     ddi$ID2 == d1 & ddi$ID1 == d2),]
if(nrow(testDF) > 0){
  edges = edges %>% bind_rows(
    testDF
  )
}
d1.targets = drug.targets[which(grepl(d1, drug.targets$`Drug IDs`, fixed = T)),] %>% 
  dplyr::filter(Species == "Humans")  
if(nrow(d1.targets) > 0){
  d1.targets = d1.targets[,1]
  from_genes = d1.targets
  
  nodes = nodes %>% bind_rows(
    data.frame(id=d1.targets, label=d1.targets, 
               color = "blue")
  )
  edges = data.frame( # connect drug1 with its direct targets
    V1 = d1,
    V2 = d1.targets
  )
  temp <- getPPIpartners(ppi.net = ppiNet, geneList = from_genes, hop=1, mindist = 1) %>%  # mindist: 0 would include Human pray genes as direct drug-targets
    na.omit() %>% 
    intersect(to_genes) # only keep neighbours of drug targets that are Human Pray genes
  
  if(length(temp) > 0){
    # set nodes for drug1 ---------------------
    nodes = nodes %>% bind_rows(
      data.frame(id=temp, label=temp, 
                 color = "green")
    ) %>%  bind_rows(
      data.frame(id=cov2_human.ppi[which(cov2_human.ppi$Human %in% temp),1] %>% unique(),
                 label=cov2_human.ppi[which(cov2_human.ppi$Human %in% temp),1] %>% unique(), 
                 color = "red")
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


d2.targets = drug.targets[which(grepl(d2, drug.targets$`Drug IDs`, fixed = T)),] %>% 
  dplyr::filter(Species == "Humans")  
if(nrow(d2.targets) > 0){
  d2.targets = d2.targets[,1]
  from_genes = d2.targets
  nodes = nodes %>% bind_rows(
    data.frame(id=d2.targets, label=d2.targets, 
               color = "blue")
  )
  edges = edges %>% bind_rows(
    data.frame( # connect drug1 with its direct targets
      V1 = d2,
      V2 = d2.targets
    )
  )
  temp <- getPPIpartners(ppi.net = ppiNet, geneList = from_genes, hop=1, mindist = 1) %>%  # mindist: 0 would include Human pray genes as direct drug-targets
    na.omit() %>% 
    intersect(to_genes) # only keep neighbours of drug targets that are Human Pray genes
  
  # set nodes for drug2 ---------------------
  if(length(temp) > 0){
    nodes = nodes %>% bind_rows(
      data.frame(id=temp, label=temp, 
                 color = "green")
    ) %>%  bind_rows(
      data.frame(id=cov2_human.ppi[which(cov2_human.ppi$Human %in% temp),1] %>% unique(),
                 label=cov2_human.ppi[which(cov2_human.ppi$Human %in% temp),1] %>% unique(), 
                 color = "red")
    )
    # set edges for drug2 ---------------------
    edges = edges %>% bind_rows( # connect drug targets with human pray genes (duplication will be there, which will be handled later)
      bind_rows(
        ppiNet.edgeDF[which(ppiNet.edgeDF$V1 %in% d2.targets & ppiNet.edgeDF$V2 %in% temp),],
        ppiNet.edgeDF[which(ppiNet.edgeDF$V2 %in% d2.targets & ppiNet.edgeDF$V1 %in% temp),]
      )
    )
    
    temp4 = cov2_human.ppi[which(cov2_human.ppi$Human %in% temp),]
    colnames(temp4) = c("V1","V2")
    edges = bind_rows(edges, temp4)
  }
}
prev_nodes = nodes[,1]

# load SARS-COV-Human PPIs
if(hasSARS_CoV){
  library(readxl)
  library(dplyr)
  SARS_cov_human.ppi = read_excel("AppData/SARSCOV-Human PPI_SF09072020a.xlsx") %>% 
    as.data.frame() %>% select(c(3,2)) %>% na.omit()
  colnames(SARS_cov_human.ppi) <- c("Sars_COV","Human")
  SARS_cov_human.ppi.net = graph_from_data_frame(SARS_cov_human.ppi, directed = F) %>% igraph::simplify()
  # SARS_cov_human.ppi.net = igraph::union(SARS_cov_human.ppi.net, ppiNet)
  SARS_cov_human.ppi.net.edgeDF = get.edgelist(SARS_cov_human.ppi.net) %>% as.data.frame()
  # to_genes = SARS_cov_human.ppi$Human
  
  drug.targets <- drug.targets.all %>% select(c(6, 12, 13)) # uniprotID, species, drugbankIDs
  
  # for drug1 ----------------
  d1.targets = drug.targets[which(grepl(d1, drug.targets$`Drug IDs`, fixed = T)),] %>% 
    dplyr::filter(Species == "SARS-CoV" | Species == "Humans") 
  
  if(nrow(d1.targets) > 0){
    d1.targets = d1.targets[,c(1,2)]
    if(length(which(d1.targets$Species == "SARS-CoV")) > 0){
      nodes = nodes %>% bind_rows(
        data.frame(id=d1.targets[which(d1.targets$Species == "SARS-CoV"),1], 
                   label=d1.targets[which(d1.targets$Species == "SARS-CoV"),1], 
                   color = "violet")
      ) 
      edges = edges %>% bind_rows(
        data.frame( # connect drug1 with its direct targets
          V1 = d1,
          V2 = d1.targets[which(d1.targets$Species == "SARS-CoV"),1]
        )
      )
    }
  temp1 <- getPPIpartners(ppi.net = SARS_cov_human.ppi.net, geneList = d1.targets[,1], hop=1, mindist = 1)   
  if(!is.na(temp1)){
    temp2 <- getPPIpartners(ppi.net = ppiNet, geneList = temp1[,1], hop=1, mindist = 1) %>% unique() %>%
      intersect(., prev_nodes)
    if(length(temp2) > 0){
      temp = getPPIpartners(ppi.net = ppiNet, geneList = temp2, hop=1, mindist = 1) %>% intersect(temp1[,1])  %>% unique()
      if(length(temp) > 0){ 
        nodes = nodes %>% bind_rows(
          data.frame(id=temp, label=temp, 
                     color = "green")
        ) %>%  bind_rows(
          data.frame(id = temp2,
                     label=temp2, 
                     color = "green")
        )
        
        edges = edges %>% bind_rows( # connect drug targets with human pray genes (duplication will be there, which will be handled later)
          bind_rows(
            ppiNet.edgeDF[which(ppiNet.edgeDF$V1 %in% temp2 & ppiNet.edgeDF$V2 %in% temp),],
            ppiNet.edgeDF[which(ppiNet.edgeDF$V2 %in% temp2 & ppiNet.edgeDF$V1 %in% temp),]
          )
        )
        temp4 = SARS_cov_human.ppi[which(SARS_cov_human.ppi$Human %in% temp &
                                           SARS_cov_human.ppi$Sars_COV %in% d1.targets[,1]),]
        colnames(temp4) = c("V1","V2")
        edges = bind_rows(edges, temp4)
      }
    }
  }
}
  
  
  # for drug 2 --------------------
  d2.targets = drug.targets[which(grepl(d2, drug.targets$`Drug IDs`, fixed = T)),] %>% 
    dplyr::filter(Species == "SARS-CoV"  | Species == "Humans") 
  if(nrow(d2.targets) > 0){
    d2.targets = d2.targets[,c(1,2)]
    if(length(which(d2.targets$Species == "SARS-CoV")) > 0){
      nodes = nodes %>% bind_rows(
        data.frame(id=d2.targets[which(d2.targets$Species == "SARS-CoV"),1], 
                   label=d2.targets[which(d2.targets$Species == "SARS-CoV"),1], 
                   color = "violet")
      ) 
      edges = edges %>% bind_rows(
        data.frame( # connect drug1 with its direct targets
          V1 = d2,
          V2 = d2.targets[which(d2.targets$Species == "SARS-CoV"),1]
        )
      )
    }
    temp1 <- getPPIpartners(ppi.net = SARS_cov_human.ppi.net, geneList = d2.targets[,1], hop=1, mindist = 1)   # mindist: 0 would include Human pray genes as direct drug-targets
    if(!is.na(temp1)){
      temp2 <- getPPIpartners(ppi.net = ppiNet, geneList = temp1[,1], hop=1, mindist = 1) %>% unique() %>%
        intersect(., prev_nodes)
      if(length(temp2) > 0){
        temp = getPPIpartners(ppi.net = ppiNet, geneList = temp2, hop=1, mindist = 1) %>% intersect(temp1[,1])  %>% unique()
        if(length(temp) > 0){ 
          nodes = nodes %>% bind_rows(
            data.frame(id=temp, label=temp, 
                       color = "green")
          ) %>%  bind_rows(
            data.frame(id = temp2,
                       label=temp2, 
                       color = "green")
          )
          edges = edges %>% bind_rows( # connect drug targets with human pray genes (duplication will be there, which will be handled later)
            bind_rows(
              ppiNet.edgeDF[which(ppiNet.edgeDF$V1 %in% temp2 & ppiNet.edgeDF$V2 %in% temp),],
              ppiNet.edgeDF[which(ppiNet.edgeDF$V2 %in% temp2 & ppiNet.edgeDF$V1 %in% temp),]
            )
          )
          temp4 = SARS_cov_human.ppi[which(SARS_cov_human.ppi$Human %in% temp &
                                             SARS_cov_human.ppi$Sars_COV %in% d1.targets[,1]),]
          colnames(temp4) = c("V1","V2")
          edges = bind_rows(edges, temp4)
        }
      }
    }
  }
}


# distinct nodes
nodes = nodes[!duplicated(nodes$id),]

# distinct edges 
edges = edges[!duplicated(edges),]
colnames(edges) = c("from","to")
# colnames(edges) = c("source","target")


# visualize the net
# visNetwork(nodes,edges) %>%
#   visNodes(shapeProperties = list(useBorderWithImage = TRUE))

res = NA
return(list(nodes, edges))
}
# 
# # # # # 
# d1 = "DB14761"
# # d1 = "DB00175"
# d2 = "DB00641"
# res <- show_COV_drug_target_Network(d1,d2)
# # head(res[0])
