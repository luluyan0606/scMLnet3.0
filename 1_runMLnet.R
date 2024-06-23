#############
## library ##
#############

library(Seurat)
library(tidyverse)
library(stMLnet)
library(parallel)
library(ggalluvial)

rm(list=ls())
gc()
setwd("/home/cjy/project/giotto_seqfish_dataset/")

source('../code/code.R')

###############
## get MLnet ##
###############

## load ####

load('./giotto_seqfish_output.rda')
Databases <- readRDS('../prior_knowledge/output/Databases.rds')

## data

GCMat <- df_norm
BarCluTable <- df_anno
clusters <- BarCluTable$Cluster %>% as.character() %>% unique()
clusters

Ligs_up_list <- Ligs_expr_list
ICGs_list <- lapply(ICGs_list, toupper)

str(Ligs_up_list)
str(Recs_expr_list)
str(ICGs_list)

## parameters

wd <- paste0("./runscMLnet/")
dir.create(wd,recursive = T)

## database ####

quan.cutoff = 0.98

RecTF.DB <- Databases$RecTF.DB %>% 
  .[.$score > quantile(.$score, quan.cutoff),] %>%
  dplyr::distinct(source, target)

LigRec.DB <- Databases$LigRec.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(target %in% RecTF.DB$source)

TFTG.DB <- Databases$TFTG.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(source %in% RecTF.DB$target)

## get multi-layer ####

outputDir <- getwd()
resMLnet <- runMLnet(ExprMat = ex_inputs$exprMat, AnnoMat = ex_inputs$annoMat,
                     LigClus = NULL, RecClus = 'Macrophage', Normalize = F, 
                     OutputDir = outputDir, Databases = NULL,
                     TGList=ex_inputs$tgs_of_inter, LigList=ex_inputs$ligs_of_inter, RecList=ex_inputs$recs_of_inter)




