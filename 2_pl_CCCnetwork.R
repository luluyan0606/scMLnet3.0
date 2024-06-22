library(CellChat)
library(ggplot2)
library(ggsci)
library(tidyverse)

rm(list = ls())
gc()

setwd("/home/yll/cell_cell_interaction/apply_in_scST/seqFISH/")

load("~/cell_cell_interaction/stMLnet_cjy/apply_in_scST/giotto_seqfish_dataset/giotto_seqfish_output.rda")
###########
## color ##
###########
# celltype

celltype <- unique(df_anno$Cluster)

scales::show_col(pal_igv(palette = "default", alpha = 0.8)(15))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.8)(15)

mycolor_ct <- mycolors_nejm[1:length(celltype)]
names(mycolor_ct) <- celltype
scales::show_col(mycolor_ct)

# nodekey

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

nodekey <- c("Ligand","Receptor","TF","Target")
mycolor_key <- mycolors_locus[1:4]
names(mycolor_key) <- nodekey
scales::show_col(mycolor_key)

# nodetype

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.8)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.8)(7)

nodetype <- c("cell","Sender","Receiver")
mycolor_nt <- mycolors_locus[1:3]
names(mycolor_nt) <- nodetype
scales::show_col(mycolor_nt)

#############
## workdir ##
#############

plotdir = './visualize/'
dir.create(plotdir,recursive = T)

#################
## NetworkPlot ##
#################

wd <- "./runscMLnet/"
files = list.files(wd)
files_tumor = files[grep(".csv",files,invert = T)]
mulNetAllList = lapply(files_tumor, function(file_tumor){
  
  readRDS(paste0(wd,file_tumor,"/scMLnet.rds"))
  
})
names(mulNetAllList) = files_tumor
mulNetAllList = mulNetAllList[!unlist(lapply(mulNetAllList, function(mulnet){nrow(mulnet$LigRec)==0}))]

mulNet_tab = lapply(mulNetAllList, function(mlnet){
  
  ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
  rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
  tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
  
  res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
    merge(., tftg, by = 'TF') %>% 
    dplyr::select(Ligand, Receptor, TF, Target) %>% 
    arrange(Ligand, Receptor)
  
  res$LRpair = paste0(res$Ligand,"_",res$Receptor)
  c(length(unique(res$LRpair)),length(unique(res$Target)))
  
}) %>% do.call('rbind',.) %>% as.data.frame()

df_cellpair <- names(mulNetAllList) %>% strsplit(.,"_") %>% do.call('rbind',.) %>% as.data.frame()

mulNet_tab <- cbind(df_cellpair,mulNet_tab)
mulNet_tab <- na.omit(mulNet_tab)
colnames(mulNet_tab) <- c('cell_from','cell_to','n_LRs','n_TGs')

cts_less = unique(mulNet_tab$cell_from)[which(!unique(mulNet_tab$cell_from) %in% unique(mulNet_tab$cell_to))]
a = data.frame(cell_from = cts_less,
               cell_to = cts_less,
               n_LRs = rep(0,length(cts_less)),
               n_TGs = rep(0,length(cts_less)))

mulNet_tab <- rbind(mulNet_tab,a)
mulNet_tab$n_LRs <- as.numeric(mulNet_tab$n_LRs)
mulNet_tab$n_TGs <- as.numeric(mulNet_tab$n_TGs)

for (key in colnames(mulNet_tab)[3:4]) {
  
  tmeTab <- mulNet_tab[,c('cell_from','cell_to',key)] %>% spread(cell_to, key) %>% 
    column_to_rownames(.,var = "cell_from")
  tmeTab[is.na(tmeTab)] <- 0 
  tmeTab <- as.matrix(tmeTab)
  
  colordb <- mycolor_ct[rownames(tmeTab)]
  
  pdf(paste0(plotdir,"cellchat_networkPlot_scMLnet3.0_",key,".pdf"),height =6,width = 6)
  netVisual_circle(tmeTab, color.use = colordb,vertex.weight = rowSums(tmeTab),alpha.edge = 0.6, 
                   weight.scale = T, label.edge= F, title.name = "Number of interactions",
                   arrow.width = 1,arrow.size = 0.3,
                   text.x = 15,text.y = 1.5)
  dev.off()
  
}











