get_induced_bipartite_graph <- function(ppiNet, setA, setB){
  if(!is.null(ppiNet)){
    # step1: get induced subgraph involving setA 
    indSubnet_setA <- induced_subgraph(ppiNet, vids = setA)
    
    # step2: get induced subgraph involving setB
    indSubnet_setB <- induced_subgraph(ppiNet, vids = setB)
    
    # step3: get induced subgraph involving setA and setB
    indSubnet_bothSet <- induced_subgraph(ppiNet, vids = c(setA, setB))
    
    # step4: remove edges of indSubnet_setA from the indSubnet_bothSet
    e = get.edgelist(indSubnet_setA) %>% apply(.,1, paste, collapse="|") %>% igraph::edges()
    indSubnet_bothSet <- indSubnet_bothSet - e
    
    # step5: remove edges of indSubnet_setA from the indSubnet_bothSet
    e = get.edgelist(indSubnet_setB) %>% apply(.,1, paste, collapse="|") %>% igraph::edges()
    indSubnet_bothSet <- indSubnet_bothSet - e
    
    return(indSubnet_bothSet)
  }
}

# # ## usage
# df = data.frame(from = c('A','B','A','C','D','D','E'),
#                 to = c('B','C','E','D','E','F','F'))
# g <- graph_from_data_frame(df, directed = F)
# # V(g)$name <- c('A', 'B', 'C', 'D', 'E', 'F')
# # E(g)$name <- LETTERS[1:10]
# plot(g)
# # get_induced_bipartite_graph(g, c('A','B', 'C'), c('D', 'E', 'F')) %>% plot()
# get_induced_bipartite_graph(g, c('A', 'C'), c('D', 'E')) %>% plot()
