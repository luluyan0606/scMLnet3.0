Select_Lig_Rec_TGs <- function(ExprMat, AnnoMat, Databases, python_path,logfc.ct, min.pct, expr.ct, pct.ct){
  
  # workdir
  de_count = ExprMat
  de_cell_type = AnnoMat
  
  # Create a Seurat object 
  rownames(de_cell_type) <- de_cell_type$Barcode
  seur <- CreateSeuratObject(counts = de_count, meta.data = de_cell_type)
  
  # normalize + find HVGs + scale
  seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
  seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 1000)
  seur <- ScaleData(seur, features = rownames(seur))
  
  Idents(seur) <- seur@meta.data$Cluster
  df_DEGs <- FindAllMarkers(seur, logfc.threshold = 0.25, min.pct = 0.1)
  DEGs_list <- df_DEGs[abs(df_DEGs$avg_log2FC) >= logfc.ct & 
                       df_DEGs$p_val_adj <= 0.05 &
                       df_DEGs$pct.1 >= 0.1,]
  DEGs_list <- split(DEGs_list$gene,DEGs_list$cluster) # 按照cluster divides genes
  str(DEGs_list)

  # calculate ligand list
  ligs_in_db <- Databases$LigRec.DB$source %>% unique()
  ligs_in_db <- intersect(ligs_in_db, rownames(st_bc_A_RCTD))
  
  clusters <- seur@active.ident %>% as.character() %>% unique()
  df_markers_ligs <- lapply(clusters, function(cluster){
    
    df <- FindMarkers(seur, ident.1 = cluster, features = ligs_in_db, only.pos = T, 
                      min.pct = min.pct) 
    df$gene <- rownames(df)
    df$ident.1 <- cluster
    df
    
  }) %>% do.call('rbind',.)
  
  Ligs_up_list <- split(df_markers_ligs$gene,df_markers_ligs$ident.1)
  str(Ligs_up_list)
  
  # calculate receptor list
  data <- seur@assays$RNA@data
  BarCluTable <- seur@cell_metadata[,1:2]
  colnames(BarCluTable) <- c('Barcode','Cluster')
  
  recs_in_db <- Databases$LigRec.DB$target %>% unique()
  
  clusters <- BarCluTable$Cluster %>% as.character() %>% unique()
  
  meanExpr_of_LR <- lapply(clusters, function(cluster){
    
    cluster.ids <- BarCluTable$Barcode[BarCluTable$Cluster == cluster]
    source_mean <- rowMeans(data[,cluster.ids])
    names(source_mean) <- rownames(data)
    source_mean
    
  }) %>% do.call('cbind',.) %>% as.data.frame()
  colnames(meanExpr_of_LR) <- clusters
  
  pct_of_LR <- lapply(clusters, function(cluster){
    
    cluster.ids <- BarCluTable$Barcode[BarCluTable$Cluster == cluster]
    dat <- data[,cluster.ids]
    pct <- rowSums(dat>0)/ncol(dat)
    names(pct) <- rownames(data)
    pct
    
  }) %>% do.call('cbind',.) %>% as.data.frame()
  colnames(pct_of_LR) <- clusters
  
  Recs_expr_list <- lapply(clusters, function(cluster){
    
    recs <- rownames(data)[meanExpr_of_LR[,cluster] >= expr.ct & pct_of_LR[,cluster] >= pct.ct]
    intersect(recs, recs_in_db)
    
  })
  names(Recs_expr_list) <- clusters
  str(Recs_expr_list)
  
  ex_inputs <- list(exprMat = data, annoMat = de_cell_type,ligs_of_inter = Ligs_up_list, recs_of_inter = Recs_expr_list,
                    tgs_of_inter = ICGs_list)
  
  return(ex_inputs)
}
