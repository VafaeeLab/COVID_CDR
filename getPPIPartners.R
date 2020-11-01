# Get PPI neighbours of a Protein
getNeighbours <- function(aNode, net, hop, mindist){
  #if the node doesn't exist in the net, return null
  if(sum(which(V(net)$name == aNode)) == 0)
    return(NA)
  nb <- as.character(ego(net, order=hop, nodes = aNode, mode = "all", mindist = mindist)[[1]]$name)
  return(nb)
}

getPPIpartners <- function(ppi.net, geneList, hop, mindist = 1){
  if(!is.null(geneList) | length(geneList) > 0){
    retList = sapply(geneList, FUN = function(x){
      parts = getNeighbours(x, ppi.net, hop, mindist)
        return(parts)
    }) %>% unlist(use.names = F) %>% unique() %>% return()
  }else
    return(NA)
  # print(is.null(geneList))
}
