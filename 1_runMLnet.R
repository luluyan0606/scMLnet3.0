#############
## library ##
#############

library(dplyr)
library(Seurat)
library(Giotto)  # remotes::install_github("drieslab/Giotto",  ref="v1.1.0")

# step1: Load unprecess data and proir database
# Download demo database and demo data (unprecessed_ex_inputs.rda) from [here](https://zenodo.org/records/10901775), the demo data `unprocessed_ex_inputs` contains the `exprMat`,`annoMat`, and `locaMat`, you need to calculate the `ligs_of_inter`, `recs_of_inter`, and `tgs_of_inter` according to the following code before running stMLnet.
load("/path/to/data/prior_databases.rda")
load("/path/to/data/unprocessed_ex_inputs.rda")

de_count <- ex_raw_inputs$Mat
de_cell_type <- ex_raw_inputs$annoMat

ex_inputs <- Select_Lig_Rec_TGs(ExprMat = de_count, 
                                AnnoMat = de_cell_type, 
                                Databases = ex_databases, 
                                logfc.ct = 2,
                                min.pct = 0.05, 
                                expr.ct = 0.1, 
                                pct.ct = 0.05)

# step2: run scMLnet
load('/path/to/data/prior_database.rda')
str(Databases,max.level=2)
load('/path/to/data/ex_inputs.rda')
str(ex_inputs,max.level=2)

library(Seurat)
library(tidyverse)
library(stMLnet)
library(parallel)
library(ggalluvial)

rm(list=ls())
gc()

source('path/to/code/creat_multilayer_network.R')

outputDir <- getwd()
resMLnet <- runMLnet(ExprMat = ex_inputs$exprMat, AnnoMat = ex_inputs$annoMat,
                     LigClus = NULL, RecClus = 'Macrophage', Normalize = F, 
                     OutputDir = outputDir, Databases = NULL,
                     TGList=ex_inputs$tgs_of_inter, LigList=ex_inputs$ligs_of_inter, RecList=ex_inputs$recs_of_inter)






