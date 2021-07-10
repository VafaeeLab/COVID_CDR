mergeMatrices <- function(matList = temp){
  header = colnames(matList[[1]])
  comb = NULL
  len = length(header)
  for(i in 1:(length(matList))){
    matA = matList[[i]]
    if(!is.null(comb))
      comb = paste0(comb, "; ", matA)
    else
      comb = matA
    comb = matrix(comb, nrow = len, ncol = len)
  }
  colnames(comb) = header
  rownames(comb) = header
  return(comb)
}
getMatrix <- function(m,d1,d2){
  m = matrix(m, nrow = length(d1), ncol = length(d2))
  colnames(m) = d2
  row.names(m) = d1
  return(m)
}


library(purrr)
library(tidyr)

# load all the matrices ---------------
library(dplyr)
library(data.table)
epsilon <- 1e-10 # A very small number

# Read inputs -------------------------------------------------------------
ce.m <- as.matrix(read.csv("data/network_separation_score_of_drugPair_v5.csv", row.names = 1))
sc.m <- as.matrix(read.csv("data/chem_similarity.csv", row.names = 1))
sp.m <- as.matrix(read.csv("data/path_similarity_cov_v2.csv",  row.names = 1)) # improved: used BioCor package (for KEGG)

sp.m <- apply(sp.m, 1, function(x) replace(x, is.infinite(x),NA))

st.m <- as.matrix(read.csv("data/target_similarity.csv",  row.names = 1))
sgoCC.m <- as.matrix(read.csv("data/GO_Sim_CC_cov_combined.csv", row.names = 1))
sgoMF.m <- as.matrix(read.csv("data/GO_Sim_MF_cov_combined.csv", row.names = 1))
sgoBP.m <- as.matrix(read.csv("data/GO_Sim_BP_cov_combined.csv", row.names = 1))


diag(ce.m) <- 0; diag(sc.m) <- 0; diag(sp.m) <- 0; diag(st.m) <- 0; diag(sgoCC.m) <- 0; diag(sgoMF.m) <- 0; diag(sgoBP.m) <- 0;

# sc.m <- ifelse(is.na(sc.m), 0, sc.m) 
sc.m <- ifelse(is.na(sc.m), epsilon, sc.m) 
# epsilon instead here just be fair for some drugs
# there may not have any pathway associated
sp.m <- ifelse(is.na(sp.m), epsilon, sp.m)
st.m <- ifelse(is.na(st.m), epsilon, st.m)
# sp.m <- ifelse(is.na(sp.m), 0, sp.m)
# st.m <- ifelse(is.na(st.m), 0, st.m)
sgoCC.m <- ifelse(is.na(sgoCC.m), epsilon, sgoCC.m) 
sgoMF.m <- ifelse(is.na(sgoMF.m), epsilon, sgoMF.m) 
sgoBP.m <- ifelse(is.na(sgoBP.m), epsilon, sgoBP.m) 

m <- rowMeans(cbind(c(sc.m), c(st.m), c(sp.m), c(sgoCC.m), c(sgoMF.m), c(sgoBP.m)), na.rm = TRUE)
tmp <- m
tmp[tmp<0.05] = NA
p <- pnorm(scale(tmp),lower.tail = F)
p.adj <- p.adjust(p, method = "fdr")
rm(tmp)

# --------------------
temp = list(ce.m, sc.m, st.m, sp.m, sgoCC.m, sgoMF.m, sgoBP.m, m, p, p.adj)
temp2 <- mergeMatrices(matList = temp)
# temp3 = temp2[upper.tri(temp2, diag = T)] %>% as.vector()
temp3 = temp2 %>% as.vector()
header = colnames(temp[[1]])
all.cv = header
xx <- purrr::cross2(rownames(temp2), rownames(temp2))
xxx <- lapply(xx, function(x){
       return(paste0(x[[2]], "::", x[[1]]))
   }) %>% unlist(use.names = F)
tempMat = cbind(xxx,temp3)
tempMat = tempMat %>% 
  as.data.frame() %>% 
  tidyr::separate(1,sep="::", into=c("ID1","ID2")) %>%
  tidyr::separate(3,sep="; ", into=c("Net_sep", "Chem_sim","target_sim","Path_sim","GO_CC","GO_MF","GO_BP", "Scores", "p_value","p.adjust")) %>% 
  dplyr::filter(ID1 != ID2)

# remove duplicates and diag entries
tempMat[1:2] = t(apply(tempMat[1:2],1,sort))
tempMat = tempMat[!duplicated(tempMat[1:2]),]

drug.info <- fread("AppData/drug links.csv") %>% dplyr::select(c(1,2))
tempMat <- dplyr::inner_join(tempMat, drug.info, by=c("ID1" = "DrugBank ID"))
tempMat <- dplyr::inner_join(tempMat, drug.info, by=c("ID2" = "DrugBank ID"))
tempMat <- tempMat[,c(1,13,2,14,3:12)]
colnames(tempMat)[1:4] = c("ID1","Name1","ID2","Name2")
fwrite(tempMat, file="data/COV_DrugSimDB_Full_with_signif_scores_v5.csv")

# for testing 
# paste0(sc.m["DB04272", "DB11967"], "; ",
#        st.m["DB04272", "DB11967"], "; ",
#        sp.m["DB04272", "DB11967"], "; ",
#        sgoCC.m["DB04272", "DB11967"], "; ",
#        sgoMF.m["DB04272", "DB11967"], "; ",
#        sgoBP.m["DB04272", "DB11967"])
