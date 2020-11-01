# Approach 1

# library(UniProt.ws)
# 
# # for human # get all the human proteins GO annotations
# up <- UniProt.ws(taxId=9606)
# keys <- head(keys(up, keytype="UNIPROTKB"))
# res <- select(up, keys=c("Q04206","P59635"), columns=c("REACTOME","GO"),
#               keytype="UNIPROTKB")
# GO_Terms = str_split(unique(res$GO)[1], ";")[[1]]
# 
# up <- UniProt.ws(taxId=694009)
# keys <- head(keys(up, keytype="UNIPROTKB"))
# res <- select(up, keys=c("P59595","P0C6X7","P0C6U8","P59636","P59635","P59632","P62937","P59594","J9TC74","P59596","Q19QW4","Q19QW5","P59637","Q7TLC7","Q19QW2","Q7TFA1","P59634","P59633","Q7TFA0"), columns=c("REACTOME","GO"),
#               keytype="UNIPROTKB")
# GO_Terms = str_split(res$GO[6], ";")[[1]]
# GO_Terms

# Approach 2
require(RJSONIO)
protID = "P59632"
# json_response <- fromJSON(paste0("https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId=", protID,"&limit=200"))
json_response <- fromJSON(paste0("https://www.ebi.ac.uk/QuickGO/services/annotation/search?includeFields=goName&includeFields=taxonName&includeFields=name&includeFields=synonyms&geneProductId=", protID,"&limit=200"))
json_response <- lapply(json_response$results, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})                     
df_response <- as.data.frame(do.call("cbind", json_response)) %>% t() %>% as.data.frame()
library(DT)
datatable((df_response %>% dplyr::select(c("goName", "goId", "goAspect"))), rownames = F)
