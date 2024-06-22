#############
## library ##
#############

library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)
library(org.Hs.eg.db)
library(clusterProfiler)     
library(ggalluvial)

rm(list=ls())
gc()
setwd("~/cell_cell_interaction/stMLnet_cjy/apply_in_COVID19/")

source('../code/code.R')

################
## get detail ##
################
## get cell pairs #### 

wd <- "./runscMLnet/"
files = list.files(wd)
files_tumor = files[grep(".csv",files,invert = T)]
mulNetAllList = lapply(files_tumor, function(file_tumor){
  
  readRDS(paste0(wd,file_tumor,"/scMLnet.rds"))
  
})
names(mulNetAllList) = files_tumor
mulNetAllList = mulNetAllList[!unlist(lapply(mulNetAllList, function(mulnet){nrow(mulnet$LigRec)==0}))]

df_cellpair <- names(mulNetAllList) %>% strsplit(.,"_") %>% do.call('rbind',.) %>% as.data.frame()
colnames(df_cellpair) <- c('Sender','Receiver') 
df_cellpair$CP1 <- paste(df_cellpair$Sender,df_cellpair$Receiver,sep ='_')
df_cellpair$CP2 <- paste(df_cellpair$Receiver,df_cellpair$Sender,sep ='_')

celltypes <- unique(c(df_cellpair$Sender,df_cellpair$Receiver))
cp_of_inter <- data.frame(matrix(ncol = 4,dimnames = list(c(),c('ct1','ct2','keys_of_ct1','keys_of_ct2'))))
for (i in 1:length(celltypes)) {
  
  ct1 <- celltypes[i]
  for (j in (i+1):length(celltypes)) {
    
    ct2 <- celltypes[j]
    cp1 <- paste0(ct1,'_',ct2)
    cp2 <- paste0(ct2,'_',ct1)
    if(cp1 %in% list.files('./runscMLnet/') & cp2 %in% list.files('./runscMLnet/')){
      
      cat('check in ',cp1,' and ',cp2,'\n')
      mlnet1 <- readRDS(paste0('./runscMLnet/',cp1,'/scMLnet.rds'))
      mlnet2 <- readRDS(paste0('./runscMLnet/',cp2,'/scMLnet.rds'))
      
      ct2_tgs <- unique(mlnet1$TFTar$target)
      ct1_tgs <- unique(mlnet2$TFTar$target)
      
      ct1_ligs <- unique(mlnet1$LigRec$source)
      ct2_ligs <- unique(mlnet2$LigRec$source)
      
      ct1_keys <- intersect(ct1_ligs,ct1_tgs)
      ct2_keys <- intersect(ct2_ligs,ct2_tgs)
      
      if(length(ct2_keys)>0|length(ct1_keys)>0){
        cp_of_inter <- rbind(cp_of_inter,c(ct1,ct2,length(ct1_keys),length(ct2_keys)))
      }
    }
  }
}

cp_of_inter <- na.omit(cp_of_inter)
cp_of_inter <- cp_of_inter[cp_of_inter$keys_of_ct1 != 0 & cp_of_inter$keys_of_ct2 != 0,]
rownames(cp_of_inter) <- 1:nrow(cp_of_inter)

## get genes list ####

key_of_inter <- list()
for (k in 1:nrow(cp_of_inter)) {
  
  ct1 <- cp_of_inter$ct1[k]
  ct2 <- cp_of_inter$ct2[k]
  
  mlnet1 <- readRDS(paste0('./runscMLnet/',paste(ct1,ct2,sep = '_'),'/scMLnet.rds'))
  mlnet2 <- readRDS(paste0('./runscMLnet/',paste(ct2,ct1,sep = '_'),'/scMLnet.rds'))
  
  ct2_tgs <- unique(mlnet1$TFTar$target)
  ct1_tgs <- unique(mlnet2$TFTar$target)
  
  ct1_ligs <- unique(mlnet1$LigRec$source)
  ct2_ligs <- unique(mlnet2$LigRec$source)
  
  ct1_keys <- intersect(ct1_ligs,ct1_tgs)
  ct2_keys <- intersect(ct2_ligs,ct2_tgs)
  
  key_of_inter[[paste(ct1,ct2,sep = '_')]] <- list(
    ct1_keys = ct1_keys,
    ct2_keys = ct2_keys
  )
  
}

key_of_inter$Alveoli_Macrophage
key_of_inter$Alveoli_Monocyte
key_of_inter$Macrophage_Monocyte

## transform to table ####

fbloop <- lapply(1:length(key_of_inter), function(i){
  
  ls <- key_of_inter[[i]]
  cp <- names(key_of_inter)[i]
  ct1 <- strsplit(cp,'_')[[1]][1]
  ct2 <- strsplit(cp,'_')[[1]][2]
  ct1_keys <- paste(ls$ct1_keys,collapse = ' ')
  ct2_keys <- paste(ls$ct2_keys,collapse = ' ')
  rbind(c(ct1,ct2,ct2_keys),c(ct2,ct1,ct1_keys))
  
}) %>% do.call('rbind',.) %>% as.data.frame()
colnames(fbloop) <- c('Sender','Receiver','Siganl')

cts_of_interest <- c('Alveoli','Macrophage','Monocyte')
fbloop <- fbloop[fbloop$Sender %in% cts_of_interest & fbloop$Receiver %in% cts_of_interest,]

