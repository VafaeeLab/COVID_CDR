library(igraph)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(tidyr)
library(GOSemSim)
library(shiny)

getRandD <- function(a_degree, net) {
  return(sample(which(degree(net) == a_degree), 1, replace = F))
}
permuteTest <- function(net, t, d, d_td, N) {
  require(foreach)
  require(doParallel)
  
  r <- c()
  cl2 <- makeCluster(parallel::detectCores())
  registerDoParallel(cl2)
  r <- foreach(
    i = 1:N,
    .combine = 'c',
    .packages = c('dplyr', 'igraph'),
    .export = c('getRandD', 'net')
  ) %dopar% {
    t_rand <- sapply(as.numeric(degree(net, v = t)), getRandD, net)
    d_rand <- sapply(as.numeric(degree(net, v = d)), getRandD, net)
    d_td_rand <-
      mean(shortest.paths(net, v = t_rand %>% unique(), to = d_rand %>% unique()))
    d_td_rand
  }
  stopCluster(cl2)
  
  r[!is.finite(r)] <- NA
  m <- mean(r, na.rm = T)
  s <- sd(r, na.rm = T)
  z <- (d_td - m) / s
  return(z)
}

enrichment.test <- function(givenSet, pathway.data, num.col) {
  require(stringi)
  require(dplyr)
  require(data.table)
  require(org.Hs.eg.db)
  
  givenSet = as.character(unique(givenSet))
  
  # num.col = max(count.fields(pop.filepath, sep = "\t", blank.lines.skip = TRUE),
  #               na.rm = TRUE)
  # pathway.data = read.table(
  #   pop.filepath,
  #   sep = "\t",
  #   header = FALSE,
  #   stringsAsFactors = FALSE,
  #   fill = TRUE,
  #   col.names = 1:num.col,
  #   blank.lines.skip = TRUE,
  #   quote = ""
  # )
  pathway.data = dplyr::select(pathway.data, -2) # throw away the the column with all NA (2nd-Column)
  row.names(pathway.data) = pathway.data[, 1]
  nTerms = nrow(pathway.data)
  pop.genes = base::union(givenSet, unique(c(as.matrix(pathway.data[,-1]))))   # option-1
  
  N = length(pop.genes)
  # enrichment test using HyperGeometric test
  info.gene.overrep = data.frame(matrix(Inf, nrow = nrow(pathway.data), ncol = 8))
  K = length(givenSet)
  for (i in 1:nTerms) {
    pathway.genes.symbol = as.character(pathway.data[i, which(pathway.data[i,] != '')])[-1]
    M = length(pathway.genes.symbol)
    x.overlap.genes.symbol = intersect(pathway.genes.symbol, givenSet)
    x = length(x.overlap.genes.symbol)
    info.gene.overrep[i, 1] = pathway.data[i, 1]
    if (x > 0) {
      # info.gene.overrep[i, 2] = paste0(pathway.genes.symbol, collapse = "/")
      # info.gene.overrep[i, 3] = paste0(x.overlap.genes.symbol, collapse = "/")
      # x.overlap.genes.entrezID = c(NA)
      
      # info.gene.overrep[i, 4] = paste0(x.overlap.genes.entrezID, collapse = "/")
      # info.gene.overrep[i, 5] = x
      # info.gene.overrep[i, 6] = paste0(x, "/", K)
      
      info.gene.overrep[i, 2] = phyper(x, M, N - M, K, lower.tail = FALSE) #He wrote : overlap.genes-1
      #insert FDR val
      # print(paste0(i))
      # print(paste0(info.gene.overrep[i, 7]))
      
      info.gene.overrep[i, 3] = p.adjust(info.gene.overrep[i, 2], method = "fdr", n = nTerms)
    }
  }
  colnames(info.gene.overrep) = c(
    "Description",
    "pvalue",
    "p.adjust"
  )
  return(info.gene.overrep[which(info.gene.overrep$pvalue <= 0.05), "Description"])
}
getFPscores <- function(drug,
                        net,
                        targets,
                        disease, 
                        pathway.data, 
                        num.col) {
  
  
  # drug targets enrichment
  t <-  targets[which(grepl(drug, targets$`Drug IDs`, fixed = T)),]
  keep <- which(t$`Gene Name` %in% V(net)$name)
  t <- unique(t[keep, 1])
  t.terms = enrichment.test(givenSet = t,
                            pathway.data = pathway.data, 
                            num.col = num.col)
  
  # print(t.terms)
  if (length(t.terms) > 0) {
    t.terms <-
      t.terms %>%  strsplit(., split = "\\(|)$")[[1]][2]  # extracting the GO terms
    
    # disease genes enrichment: load from pre-computed enriched terms
    # of disease genes
    d.terms <-
      fread("AppData/covid_go_GO_BP_2020.csv", encoding = "UTF-8") %>%
      dplyr::select("Description") %>%
      unlist(use.names = F) %>%
      strsplit("\\(|)$")[[1]][2]
    
    # semantic similarity
    hsGO <- godata('org.Hs.eg.db', ont = ontTag)
    
    return(
      mgoSim(
        t.terms,
        d.terms,
        semData = hsGO,
        measure = "Wang",
        combine = "BMA"
      ))
  } else{
    return(0)
  }
}

#' calculate topology-based proximity measures of drugs with disease module.
#' Note, here the disease module considers the human (host) target genes of
#' pathogens.
#'
#' @param drug
#' @param net
#' @param targets
#' @param disease
#'
#' @return
#' @export
#'
#' @examples
getzScore <- function(drug, net, targets, disease) {
  t <-  targets[which(grepl(drug, targets$`Drug IDs`, fixed = T)), ]
  keep <- which(t$`Gene Name` %in% V(net)$name)
  t <- unique(t[keep, 1])
  
  d <- disease
  keep <- which(d %in% V(net)$name)
  d <- unique(d[keep])
  
  d_td <- mean(shortest.paths(net, v = t, to = d))
  
  z <- permuteTest(net, t, d, d_td, 2)
  p <- pnorm(z)
  
  return(cbind(z, p))
}

#' calculate topology-based proximity measures of drugs with disease module.
#' Note, here the disease module considers the human (host) target genes of
#' pathogens.
#'
#' @param allQuery.drugs
#' @param ppinet
#' @param drug_target.dt
#' @param disease.genes
#'
#' @return
#' @export
#'
#' @examples
get_topological_proximities_for_drugs <-
  function(allQuery.drugs,
           ppinet,
           drug_target.dt,
           disease.genes) {
    # notify users of the progress
    shiny::showNotification(
      id = "TP_notification",
      "Topological-Proximity calculation has been started",
      type = "message",
      duration = NULL
    )
    cl <- makeCluster(parallel::detectCores())
    registerDoParallel(cl)

    dp.tp <-
      foreach(
        i = 1:length(allQuery.drugs),
        .combine = 'rbind',
        .packages = c('dplyr', 'igraph'),
        .export = c('getzScore', 'permuteTest', 'getRandD')
      ) %dopar% {
        getzScore(
          allQuery.drugs[i],
          net = ppinet,
          targets = drug_target.dt,
          disease = disease.genes
        )
      }
    # print(dp.t)
    stopCluster(cl)
    shiny::removeNotification(id = "TP_notification")
    shiny::showNotification(
      id = "TP_notification",
      "Topological-Proximity calculation has been finished",
      type = "message",
      duration = 15
    )
    return(dp.tp)
  }

get_functional_proximities_for_drugs <-
  function(allQuery.drugs,
           ppinet,
           drug_target.dt,
           disease.genes,
           pathway.data = pathway.data,
           num.col = num.col) {
    # TODO: disease proximity (functional)
    shiny::showNotification(
      id = "FP_notification",
      "Functional-Proximity calculation has been started",
      type = "message",
      duration = NULL
    )
    cl <- makeCluster(parallel::detectCores())
    registerDoParallel(cl)
    
    dp.fn <-
      foreach(
        i = 1:length(allQuery.drugs),
        .combine = 'rbind',
        .packages = c('dplyr',
                      'igraph',
                      'GOSemSim',
                      'doParallel',
                      'foreach'),
        # .verbose = T,
        .export = c('getFPscores', 'enrichment.test', 'pathway.data', 'num.col')
      ) %dopar% {
        getFPscores(
          allQuery.drugs[i],
          net = ppinet,
          targets = drug_target.dt,
          disease = disease.genes,
          pathway.data = pathway.data,
          num.col = num.col
        )
      }
    # print(dp.fn)
    stopCluster(cl)
    shiny::removeNotification(id = "FP_notification")
    shiny::showNotification(
      "Functional-Proximity calculation has been finished",
      type = "message",
      duration = 5
    )
    return(dp.fn)
  }


get_Complementary_Exposure_score <- function(drugA, drugB, net, targets){
  require(igraph)
  require(dplyr)
  
  tA <-  targets[which(grepl(drugA, targets$`Drug IDs`, fixed = T)),] %>% dplyr::filter(Species == "Humans")
  keep <- which(tA$`Gene Name` %in% V(net)$name)
  tA <- unique(tA[keep,1])
  
  tB <-  targets[which(grepl(drugB, targets$`Drug IDs`, fixed = T)),] %>% dplyr::filter(Species == "Humans") 
  keep <- which(tB$`Gene Name` %in% V(net)$name)
  tB <- unique(tB[keep,1])
  
  dAB <- mean(shortest.paths(net, v = tA, to=tB))
  dAA <- mean(shortest.paths(net, v = tA, to=tA))
  dBB <- mean(shortest.paths(net, v = tB, to=tB))
  sAB <- dAB - (dAA+dBB)/2
  if(is.nan(sAB))
    sAB <- NA
  return(sAB)
}

#' Deprecated, use the other overloaded function
#'
#' @param net_obj
#' @param ppinet
#' @param drug_target.dt
#' @param disease.genes
#'
#' @return
#' @export
#'
#' @examples
get_updated_net_user_defined <-
  function(net_obj,
           ppinet,
           drug_target.dt,
           disease.genes) {
    if (length(net_obj) == 2 & nrow(drug_target.dt) != 0) {
      require(igraph)
      require(dplyr)
      
      nodes = net_obj[[1]]
      edges = net_obj[[2]]
      
      
      # work on the new network object creation (igraph) ----
      dt.net = igraph::graph_from_data_frame(drug_target.dt[, c(1, 3)], directed = F)
      old.net = igraph::graph_from_data_frame(edges[, c(1, 2)], directed = F)
      new.net = igraph::union(dt.net, old.net)
      # -----------------------------------------------------
      
      # work on the edge renewal ----------------------------
      newDrugs = drug_target.dt[, c(1, 2)] %>%  unique()
      newDrugIDs = newDrugs[, 1] %>% as.character()
      newDrugLabels = newDrugs[, 2] %>% as.character()
      newTargets = drug_target.dt[, 3] %>% as.character() %>% unique()
      
      newEdges = edges %>% bind_rows(data.frame(
        from = newDrugIDs,
        to = newTargets,
        id = seq(
          from = nrow(edges) + 1,
          to = nrow(edges) + nrow(drug_target.dt),
          by = 1
        ),
        dashes = FALSE
      ))
      # -----------------------------------------------------
      
      # work on the node renewal ----------------------------
      allQuery.drugs = c(
        newDrugIDs,
        nodes %>%
          dplyr::filter(color == "pink") %>%
          dplyr::select("id") %>%
          unlist(use.names = F)
      )
      # print(head(allQuery.drugs))
      
      # disease proximity (topological)
      dp.t = lapply(allQuery.drugs, function(x) {
        getzScore(x,
                  net = ppinet,
                  targets = drug_target.dt,
                  disease = disease.genes)
      }) %>% do.call("rbind", .) %>% apply(., 1, function(x) {
        paste0("z: ", x[1], " with p-value: ", x[2])
      })
      
      names(dp.t) = allQuery.drugs
      # TODO: disease proximity (functional)
      
      # update the node dataframe
      nonDrugNodes = nodes %>% dplyr::filter(color != "pink")
      drugNodes = nodes %>% dplyr::filter(color == "pink")
      newNodes = nonDrugNodes %>% bind_rows(data.frame(
        # old drug infos
        id = drugNodes$id,
        label = drugNodes$label,
        color = "pink",
        title = dp.t[drugNodes$id]
      )) %>% bind_rows(data.frame(
        # new drug infos
        id = newDrugIDs,
        label = newDrugLabels,
        color = "pink",
        title = dp.t[newDrugIDs]
      )) %>% bind_rows(data.frame(
        # new target infos
        id = newTargets,
        label = newTargets,
        color = "blue",
        title = paste0("<p>Name: <b>", newTargets, "</b></p>")
      ))
      # -----------------------------------------------------
      
      # enlist them and return
      return(list(allQuery.drugs, list(newNodes, newEdges)))
    }
  }
