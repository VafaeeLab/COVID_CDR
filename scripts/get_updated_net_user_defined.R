# get_updated_net_user_defined <- function(net_obj, drug_target.dt){
#   source("utils.R")
#   
#   if(length(net_obj) == 2 & nrow(drug_target.dt) != 0){
#     require(igraph)
#     require(dplyr)
#     
#     nodes = net_obj[[1]]
#     edges = net_obj[[2]]
#     
#     # work on the new network object creation (igraph) ----
#     dt.net = igraph::graph_from_data_frame(drug_target.dt[,c(1,3)], directed = F)
#     old.net = igraph::graph_from_data_frame(edges[,c(1,2)], directed = F)
#     new.net = igraph::union(dt.net, old.net)
#     # -----------------------------------------------------
#     
#     # work on the edge renewal ----------------------------
#     newDrugs = drug_target.dt[,1]
#     newTargets = drug_target.dt[,3]
#     
#     newEdges = edges %>% bind_rows(
#       data.frame(
#         from = newDrugs,
#         to = newTargets,
#         id = seq(from = nrow(edges), to = nrow(edges)+nrow(drug_target.dt), by = 1),
#         dashes = FALSE
#       )
#     )
#     # -----------------------------------------------------
#     
#     # work on the node renewal ----------------------------
#     allQuery.drugs = c(newDrugs, 
#                        nodes %>% 
#                          dplyr::filter(color = "pink") %>% 
#                          dplyr::select("color") %>% 
#                          unlist(use.names = F)
#                        )
#     # disease proximity (topological)
#     dp.t = lapply(allQuery.drugs, getzScore(., net = new.net, targets = drug_target.dt, disease = ))
#     # -----------------------------------------------------
#     
#     # enlist them and return
#   }
# }