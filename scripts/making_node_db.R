library(data.table)
library(dplyr)
library(readxl)

# on Human-Human PPI
human.ppi <- fread("AppData/i2d.human.anno.ppi.Genes.csv") %>% as.data.frame() # can look for only the signaling net
node.db = human.ppi[,c(3,1)]
colnames(node.db) = c("GeneSymbol","UniprotID")
colnames(human.ppi)[c(4,2)] = c("GeneSymbol","UniprotID")
node.db = node.db %>% bind_rows(human.ppi[,c(4,2)])
node.db = node.db %>% unique()

# on COV2-Human PPI
cov2_genes <- read_excel("AppData/SARS_COV2_Gene_Prot_Mapping_SF_27_OCT.xlsx") %>% as.data.frame() 
colnames(cov2_genes) = c("GeneSymbol","UniprotID")
node.db = node.db %>% bind_rows(cov2_genes)

# on SARS-COV
sars_cov_genes <- read_excel("AppData/SARSCOV-Human PPI_SF09072020a.xlsx") %>% dplyr::select(c(1,3)) %>% as.data.frame()
colnames(sars_cov_genes) = c("GeneSymbol","UniprotID")
node.db = node.db %>% bind_rows(sars_cov_genes)

# on Influenza
sars_cov_genes <- fread("AppData/Influenza-Human-PPI_23_OCT.csv") %>% dplyr::select(c(1,2)) %>% as.data.frame()
colnames(sars_cov_genes) = c("GeneSymbol","UniprotID")
node.db = node.db %>% bind_rows(sars_cov_genes)

fwrite(node.db, file = "AppData/covid19_comb_therap_node_V5.db.csv")
  

