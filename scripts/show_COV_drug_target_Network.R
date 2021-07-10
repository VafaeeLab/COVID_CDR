




#' Constructs the super-network for a drug combination combining multiple ppi
#' networks
#'
#' @param drugs
#' @param hasSARS_CoV
#' @param hasInfz
#' @param hop
#' @param mindist
#'
#' @return a list of 2 objects: nodes and edges of the super-network
#' @export
#'
#' @examples
show_COV_drug_target_Network <-
  function(drugs,
           hasSARS_CoV = T,
           hasInfz = T,
           hop,
           mindist,
           all.cv_full,
           ppiNet,
           cov2_human.ppi,
           ddi,
           drug.targets.all,
           newDrug.Names) {
    require(data.table)
    require(dplyr)
    require(ggplot2)
    require(igraph)
    require(visNetwork)
    require(readxl)
    source("get_COV2_drug_net.R")
    source("get_SARS_drug_net.R")
    source("get_Infz_drug_net.R")
    source("utils.R")
    
    
    cov2_human.ppi.net <-
      graph_from_data_frame(cov2_human.ppi, directed = F) %>% igraph::simplify()
    
    drug.targets <-
      drug.targets.all %>% dplyr::select(c(1, 3, 4)) # GeneName, species, drugbankIDs
    
    # a pair:
    to_genes = cov2_human.ppi$Human
    ppiNet.edgeDF = get.edgelist(ppiNet) %>% as.data.frame()
    
    # drugNames = drug.info[match(drugs, drug.info$`DrugBank ID`),2]
    # gathering drug information (name, class, indication, proximity scores)
    drugs.df = data.frame(DrugBankID = drugs, stringsAsFactors = F)
    temp.df = dplyr::left_join(drugs.df,
                               all.cv_full[, c(
                                 "DrugBankID",
                                 "Name",
                                 "Therap_Class",
                                 "Indication(s)",
                                 "z_DP",
                                 "z_DP_pVal",
                                 "z_FP",
                                 "z_FP_pVal"
                               )])
    # if new drugs are experimented
    print(newDrug.Names)
    if(!is.null(newDrug.Names)) {
      temp.df[match(newDrug.Names$`Drug IDs`, temp.df$DrugBankID), "Name"] = newDrug.Names$DrugNames
    }
    drugNames <-
      temp.df %>% dplyr::select("Name") %>% unlist(use.names = F)
    
    drugClass <-
      temp.df %>% dplyr::select("Therap_Class") %>% unlist(use.names = F)
    drugIndication <-
      temp.df %>% dplyr::select("Indication(s)") %>% unlist(use.names = F)
    
    dd.prox.score <-
      temp.df %>% dplyr::select(c("z_DP", "z_DP_pVal")) %>%
      dplyr::rename("z_score" = "z_DP",
                    "p_value" = "z_DP_pVal")
    fp.score <-
      temp.df %>% dplyr::select(c("DrugBankID", "z_FP", "z_FP_pVal"))
    # print(dd.prox.score[, 1])
    
    if (any(is.na(dd.prox.score$z_score))) {
      res <- get_topological_proximities_for_drugs(
        allQuery.drugs = drugs,
        ppinet = ppiNet,
        drug_target.dt = drug.targets,
        disease.genes = to_genes
      )
      dd.prox.score <- res
    }
    if (any(is.na(fp.score$z_FP))) {
      # determine which drugs (new) needs to calculate new
      newDrugs = fp.score$DrugBankID[fp.score$DrugBankID %in% drugs &
                                       is.na(fp.score$z_FP)]
      na_pos = which(fp.score$DrugBankID %in% drugs &
                       is.na(fp.score$z_FP))
      
      # need to calculate FP (functional proximities) for the new drugs
      pop.filepath = "AppData/GO_Biological_Process_2018 (copy).txt"
      num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE),
                    na.rm = TRUE)
      pathway.data = read.table(
        pop.filepath,
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE,
        fill = TRUE,
        col.names = 1:num.col,
        blank.lines.skip = TRUE,
        quote = ""
      )
      
      res <- get_functional_proximities_for_drugs(
        allQuery.drugs = newDrugs,
        ppinet = ppiNet,
        drug_target.dt = drug.targets,
        disease.genes = to_genes, 
        pathway.data = pathway.data,
        num.col = num.col
      )
      res.df = data.frame(
        DrugID = newDrugs,
        Functional_Proximity = res,
        stringsAsFactors = F
      )
      # now calculate z-scores and its p-value, for which we need to import the
      # whole population
      functional.proximity.drugs = fread(
        "AppData/Functional_proximity_score_v5.csv",
        header = F,
        skip = 1
      ) %>% dplyr::select(-1) %>%
        dplyr::rename(DrugID = V2, Functional_Proximity = V3) %>%
        as.data.frame() %>% dplyr::bind_rows(res.df)
      p <-
        pnorm(scale(functional.proximity.drugs$Functional_Proximity),
              lower.tail = F)
      # p.adj <- p.adjust(p, method = "fdr")
      functional.proximity.drugs$z_FP = scale(functional.proximity.drugs$Functional_Proximity)
      functional.proximity.drugs$z_FP_pVal = p
      
      # now replace the NAs in proper places
      fp.score[match(newDrugs, fp.score$DrugBankID), "z_FP"] <-
        functional.proximity.drugs[match(newDrugs, functional.proximity.drugs$DrugID), "z_FP"]
      fp.score[match(newDrugs, fp.score$DrugBankID), "z_FP_pVal"] <-
        functional.proximity.drugs[match(newDrugs, functional.proximity.drugs$DrugID), "z_FP_pVal"]
    }
    
    nodes = data.frame(
      id = drugs,
      label = drugNames,
      # image = "https://ndownloader.figshare.com/files/25144541",
      # shape = "image",
      color = rep("pink", length(drugNames)),
      title = paste0(
        "<p>Name: <b>",
        drugNames,
        "</b></p>",
        "<p>Class: <b>",
        drugClass,
        "</b></p>",
        "<p>Indications: <b>",
        drugIndication,
        "</b></p>",
        "<p>COVID-19 Proximity (functional): <b> z-score:",
        fp.score[, 2] %>%
          round(digits = 2),
        " (p-value:" ,
        fp.score[, 3] %>%
          round(digits = 2),
        ")",
        "</b></p>",
        "<p>COVID-19 Proximity (topological): <b>z-score:",
        dd.prox.score[, 1] %>%
          round(digits = 2),
        " (p-value:" ,
        dd.prox.score[, 2] %>%
          round(digits = 2),
        ")",
        "</b></p>"
      ),
      stringsAsFactors = F
    )
    
    edges = data.frame()
    # check if DDI exists between d1 and d2
    testDF = ddi[which((ddi$ID1 %in% drugs) &
                         ddi$ID2 %in% drugs),]
    if (nrow(testDF) > 0) {
      edges = edges %>% bind_rows(testDF)
    }
    # get COV2 net for all drugs
    # print(tail(drug.targets))
    
    for (i in 1:length(drugs)) {
      aDrug = drugs[i]
      result <-
        get_COV2_drug_net(
          aDrug,
          drug.targets,
          to_genes,
          ppiNet,
          hop = hop,
          mindist = mindist,
          cov2_human.ppi,
          ppiNet.edgeDF
        )
      # print(head(result[[2]]))
      if (nrow(result[[1]]) > 0) {
        nodes = nodes %>% bind_rows(result[[1]])
      }
      if (nrow(result[[2]]) > 0) {
        edges = edges %>% bind_rows(result[[2]])
      }
    }
    
    prev_nodes = nodes[, 1]
    
    # load SARS-COV-Human PPIs
    if (hasSARS_CoV) {
      library(readxl)
      library(dplyr)
      SARS_uniprot <-
        read_excel("AppData/SARS_Uniprot.xlsx") %>% dplyr::select(1, 5) %>% as.data.frame()
      
      SARS_cov_human.ppi = read_excel("AppData/SARSCOV-Human PPI_SF09072020a.xlsx") %>%
        as.data.frame() %>% dplyr::select(c(1, 3, 2)) %>% na.omit()
      colnames(SARS_cov_human.ppi) <-
        c("Sars_COV_dgs", "Sars_COV_Uniprot", "Human")
      SARS_cov_human.ppi$Sars_COV_dgs <-
        SARS_cov_human.ppi$Sars_COV_dgs %>% toupper()
      
      # remdesivir selection -------
      SARS_cov_human.ppi <-
        SARS_cov_human.ppi %>% dplyr::filter(Sars_COV_dgs %in% c("SARS-COV NSP12", "SARS-COV NSP14"))
      # ---------
      
      SARS_cov_human.ppi = dplyr::inner_join(SARS_cov_human.ppi,
                                             SARS_uniprot,
                                             by = c("Sars_COV_Uniprot" = "Entry"))
      SARS_cov_human.ppi.net1 = graph_from_data_frame(SARS_cov_human.ppi[, c("Sars_COV_dgs", "Gene Names")], directed = F) %>% igraph::simplify()
      SARS_cov_human.ppi.net2 = graph_from_data_frame(SARS_cov_human.ppi[, c("Sars_COV_dgs", "Human")], directed = F) %>% igraph::simplify()
      SARS_cov_human.ppi.net = igraph::union(SARS_cov_human.ppi.net1, SARS_cov_human.ppi.net2) %>% igraph::simplify()
      # SARS_cov_human.ppi.net = igraph::union(SARS_cov_human.ppi.net, ppiNet)
      SARS_cov_human.ppi.net.edgeDF = get.edgelist(SARS_cov_human.ppi.net) %>% as.data.frame()
      # to_genes = SARS_cov_human.ppi$Human
      
      # drug.targets <- drug.targets.all %>% dplyr::select(c(6, 12, 13)) # uniprotID, species, drugbankIDs
      drug.targets <-
        drug.targets.all  # GeneName, uniprotID, species, drugbankIDs
      
      # for drug1 ----------------
      for (i in 1:length(drugs)) {
        aDrug = drugs[i]
        # print(aDrug)
        result <-
          get_SARS_drug_net(
            d1 =  aDrug,
            prev_nodes =  prev_nodes,
            SARS_cov_human.ppi = SARS_cov_human.ppi,
            drug.targets =  drug.targets,
            drug.targets.all = drug.targets.all,
            SARS_cov_human.ppi.net = SARS_cov_human.ppi.net,
            ppiNet = ppiNet,
            hop = hop,
            mindist = mindist,
            ppiNet.edgeDF = ppiNet.edgeDF
          )
        # print(result[[2]])
        
        if (nrow(result[[1]]) > 0) {
          nodes = nodes %>% bind_rows(result[[1]])
        }
        if (nrow(result[[2]]) > 0) {
          edges = edges %>% bind_rows(result[[2]])
        }
      }
      
    }
    
    prev_nodes = nodes[, 1]
    
    # load Infz-Human PPIs
    if (hasInfz) {
      library(readxl)
      library(dplyr)
      Infz_uniprot <-
        read_excel("AppData/Infz_Uniprot.xlsx") %>% as.data.frame()
      
      Infz_human.ppi = read_excel("AppData/Infz_HumanPPI_26_Oct.xlsx") %>% as.data.frame()
      colnames(Infz_human.ppi) <-
        c("Infz_dgs", "Infz_Uniprot", "Human")
      Infz_human.ppi$Infz_dgs <-
        Infz_human.ppi$Infz_dgs %>% toupper()
      
      
      Infz_human.ppi = dplyr::inner_join(Infz_human.ppi,
                                         Infz_uniprot,
                                         by = c("Infz_Uniprot" = "UniProtID"))
      Infz_human.ppi.net1 = graph_from_data_frame(Infz_human.ppi[, c("Infz_dgs", "Gene Names")], directed = F) %>% igraph::simplify()
      Infz_human.ppi.net2 = graph_from_data_frame(Infz_human.ppi[, c("Infz_dgs", "Human")], directed = F) %>% igraph::simplify()
      Infz_human.ppi.net = igraph::union(Infz_human.ppi.net1, Infz_human.ppi.net2) %>% igraph::simplify()
      Infz_human.ppi.net.edgeDF = get.edgelist(Infz_human.ppi.net) %>% as.data.frame()
      
      # drug.targets <- drug.targets.all %>% dplyr::select(c(6, 12, 13)) # uniprotID, species, drugbankIDs
      drug.targets <-
        drug.targets.all  # GeneName, uniprotID, species, drugbankIDs
      
      # for drug1 ----------------
      for (i in 1:length(drugs)) {
        aDrug = drugs[i]
        # print(aDrug)
        result <-
          get_Infz_drug_net(
            d1 =  aDrug,
            prev_nodes =  prev_nodes,
            Infz_human.ppi = Infz_human.ppi,
            drug.targets =  drug.targets,
            drug.targets.all = drug.targets.all,
            Infz_human.ppi.net = Infz_human.ppi.net,
            ppiNet = ppiNet,
            hop = hop,
            mindist = mindist,
            ppiNet.edgeDF = ppiNet.edgeDF
          )
        # print(result[[2]])
        
        if (nrow(result[[1]]) > 0) {
          nodes = nodes %>% bind_rows(result[[1]])
        }
        if (nrow(result[[2]]) > 0) {
          edges = edges %>% bind_rows(result[[2]])
        }
      }
      
    }
    
    # distinct nodes
    nodes = nodes[!duplicated(nodes$id),]
    
    # distinct edges
    edges = edges[!duplicated(edges),]
    
    if (nrow(edges) > 0) {
      colnames(edges) = c("from", "to")
      edges$id = 1:nrow(edges)
      edges$dashes = F
      
      # Remdesivir changes ---------------
      edges[which((edges$from == "DB14761" &
                     edges$to == "SARS-COV NSP14") |
                    (edges$to == "DB14761" &
                       edges$from == "SARS-COV NSP14")
      ), 4] = T
      edges[which((edges$from == "DB14761" &
                     edges$to == "SARS-CoV2 nsp14") |
                    (edges$to == "DB14761" &
                       edges$from == "SARS-CoV2 nsp14")
      ), 4] = T
      nodes[which(nodes$label == "SARS-CoV2 nsp12"), "color"] = "red"
      nodes[which(nodes$label == "SARS-CoV2 nsp14"), "color"] = "red"
      # ------
    }
    
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
