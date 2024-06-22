library(iTALK)
library(circlize)
library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)
library(ggalluvial)

rm(list=ls())
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
files = files[grep(".csv",files,invert = T)]

df_LRTG = lapply(files, function(file){
  
  mlnet = readRDS(paste0(wd,file,"/scMLnet.rds"))
  
  if (length(mlnet$LigRec)==0 ){
    message("LigRec is empty in file ", file, ", skipping...")
    return(NULL)  # Return NULL to skip further processing
    
  }else{
    ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
    rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
    tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
    
    res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
      merge(., tftg, by = 'TF') %>% 
      dplyr::select(Ligand, Receptor, TF, Target) %>% 
      arrange(Ligand, Receptor)
    
    res$LRpair = paste0(res$Ligand,"_",res$Receptor)
    
    LRTG <- lapply(unique(res$LRpair), function(lr){
      
      TG_num <- length(res$Target[res$LRpair == lr])
      df <- data.frame(source = lr %>% gsub('_.*','',.),
                       target = lr %>% gsub('.*_','',.),
                       LRpair = lr,
                       count = TG_num)
    }) %>% do.call('rbind',.)
    
    # file <- gsub('-','_',file)
    df_LigRec <- data.frame(
      source = LRTG$source,
      target = LRTG$target,
      LRpair = LRTG$LRpair,
      count = LRTG$count,
      source_group = file %>% gsub('_.*','',.),
      target_group = file %>% gsub('.*_','',.)
    )
    
  }
  
}) %>% do.call('rbind',.)


if(!is.null(df_LRTG)){
  
  # input
  res <- data.frame(ligand = df_LRTG$source,receptor = df_LRTG$target,
                    cell_from = df_LRTG$source_group,cell_to = df_LRTG$target_group, 
                    count = df_LRTG$count, comm_type = "growth factor")
  res <- res[order(res$count,decreasing=T),] 
  # plot
  if (dim(res)[1] > 20){
    # select top 20
    res <- res[1:30,]
  }
  
  min(res$count)
  
  colordb <- mycolor_ct[which(names(mycolor_ct) %in% c(unique(res$cell_from),unique(res$cell_to)))]
  scales::show_col(colordb)
  print(colordb)
  
  pdf(paste0(plotdir,"ChordPlot_v2_LRscore_scMLnet3.0_","all",".pdf"),height = 4.5,width = 4.5)
  LRPlot(res,datatype='mean count',
         cell_col=colordb,
         link.arr.lwd=res$count,
         link.arr.col="#696969", # "#808080"
         link.arr.width=0.15,
         track.height_1 = uh(1,"mm"),
         track.height_2 = uh(11,"mm"),
         annotation.height_1 = 0.015,
         annotation.height_2 = 0.01,
         text.vjust = "0.4cm")
  dev.off()
  
}



