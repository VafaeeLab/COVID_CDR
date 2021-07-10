
list2df <- function (list.object, col1 = "X1", col2 = "X2") {
  if (is.null(names(list.object))) {
    names(list.object) <- seq_along(list.object)
  }
  dat <- data.frame(x = unlist(list.object, FALSE), y = rep(names(list.object), 
                                                            sapply(list.object, length)), stringsAsFactors = FALSE, 
                    check.names = FALSE, row.names = NULL)
  colnames(dat) <- c(col1, col2)
  dat
}

TTD_file_processing <- function(fileName,
                                valueKey){

library(data.table)
library(dplyr)

# tempDat <- read.table("data/TTD/P1-02-TTD_drug_download.txt", blank.lines.skip = T, sep = "\t", header = F, stringsAsFactors = F)
tempDat <- fread(fileName, blank.lines.skip = T, stringsAsFactors = F, header = F, encoding = "UTF-8") %>% as.data.frame()
tempDat$V4 = paste0(tempDat$V2, ":::", tempDat$V3)

res <- split(tempDat$V4, tempDat$V1) %>% lapply(FUN = function(x){
  temp = strsplit(x =  x, split = ":::", fixed = T)
  
  indx = grep(pattern = valueKey, x = temp, fixed = T)
  if(length(indx) > 0){
    val = temp[indx] %>% unlist()
    if(valueKey != "INDICATI"){
      return( val[2])
    }else{
      retVal = "\n"
      for(i in 1:length(val)){
        if(i %% 2 == 0){
          indName = stri_split_regex(val[i], "\\s\\[(.*?)\\]\\s") %>% unlist()
          retVal = paste0(retVal, indName[1], "\n")
        }
      }
      return(retVal)
    }
  }else{
    return(NA)
  }
})
# res <- list2df(res)%>% dplyr::select(2,1) %>% dplyr::mutate(c("X2", "X1") <- c("ID", "Value")) %>% as.data.frame() %>% return()
res <- list2df(res)%>% dplyr::select(2,1) %>% as.data.frame()
colnames(res) <- c("ID", "Value")
return(res)

}