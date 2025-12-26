similarity_pipe <- function(train_data, test_data,
                            feature_num, downSample_num,
                            test_layer,
                            group_label,
                            plot_width,
                            plot_height){
  
  ### a. train data ----
  train <- train_data
  Idents(train) <- "cell_type"
  
  ### b. test data ----
  test <- test_data
  Idents(test) <- "seurat_clusters"
  test <- JoinLayers(test)
  test <- NormalizeData(test)
  
  sd_train <- sort(apply(GetAssayData(train, assay = "RNA", layer = "data"), 1, sd), decreasing = T)
  sd_test <- sort(apply(GetAssayData(test, assay = "RNA", layer = test_layer), 1, sd), decreasing = T)
  
  features_used <- intersect(names(sd_train)[1:feature_num], names(sd_test)[1:feature_num])
  
  ## subset data
  subset_by_group <- function(obj, group.by = "cell_type", n = 200, seed = 2025) {
    set.seed(seed)
    
    cells_by_group <- split(colnames(obj), obj[[group.by]])
    
    selected_cells <- unlist(lapply(cells_by_group, function(cells) {
      if (length(cells) > n) {
        return(sample(cells, n))
      } else {
        return(cells)
      }
    }))
    
    obj_sub <- subset(obj, cells = selected_cells)
    
    return(obj_sub)
  }
  
  train_sub <- subset_by_group(train, group.by = "cell_type", n = downSample_num)
  test_sub <- subset_by_group(test,  group.by = "seurat_clusters", n = downSample_num)
  
  print(table(train_sub$cell_type))
  print(table(test_sub$seurat_clusters))
  
  train_sparse <- GetAssayData(train_sub, assay = "RNA", layer = "data")[features_used, ]
  test_sparse <- GetAssayData(test_sub,  assay = "RNA", layer = test_layer)[features_used, ]
  
  train_mat <- as.matrix(train_sparse)
  test_mat <- as.matrix(test_sparse)
  
  rm(train_sparse, test_sparse)
  gc()
  
  train_group <- as.character(train_sub$cell_type)
  test_group  <- as.character(test_sub$seurat_clusters)
  
  table(is.na(train_group))
  table(is.na(test_group))
  
  if (ncol(train_mat) != length(train_group)) stop("Train 维度不匹配！")
  if (ncol(test_mat)  != length(test_group))  stop("Test 维度不匹配！")
  
  ### c. calculate ----
  source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')
  
  Test_similarity <- glm.predict(
    train.data = train_mat,     
    train.group = train_group,
    genes.used = features_used,
    downsample = FALSE,
    sample.cells = 0,
    test.data = test_mat,
    test.group = test_group,
    alpha = 0.5, 
    nfolds = 5,
    seed = 123
  )
  
  heatmap_mat = Test_similarity$heatmap@matrix
  
  pdf(paste0("2.Similarity_", group_label, "_", feature_num, "_genes.pdf"), width=plot_width, height=plot_height, onefile = F)
  pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                     color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                     cluster_rows = F)
  dev.off()
  
  write.csv(heatmap_mat, paste0("2.Similarity_", group_label, "_", feature_num, "_genes.csv"), quote = FALSE, row.names = TRUE)
  save(list = c("train_mat", "test_mat", "Test_similarity", "heatmap_mat"),
       file = paste0("2.Similarity_", group_label, "_", feature_num, "_genes.RData"))
  
}







similarity_t_pipe <- function(train_data, test_data,
                              feature_num, downSample_num,
                              test_layer,
                              group_label,
                              plot_width,
                              plot_height){
  
  ### a. train data ----
  train <- train_data
  Idents(train) <- "seurat_clusters"
  train <- JoinLayers(train)
  train <- NormalizeData(train)
  
  ### b. test data ----
  test <- test_data
  Idents(test) <- "cell_type"

  sd_train <- sort(apply(GetAssayData(train, assay = "RNA", layer = "data"), 1, sd), decreasing = T)
  sd_test <- sort(apply(GetAssayData(test, assay = "RNA", layer = test_layer), 1, sd), decreasing = T)
  
  features_used <- intersect(names(sd_train)[1:feature_num], names(sd_test)[1:feature_num])
  
  ## subset data
  subset_by_group <- function(obj, group.by = "cell_type", n = 200, seed = 2025) {
    set.seed(seed)
    
    cells_by_group <- split(colnames(obj), obj[[group.by]])
    
    selected_cells <- unlist(lapply(cells_by_group, function(cells) {
      if (length(cells) > n) {
        return(sample(cells, n))
      } else {
        return(cells)
      }
    }))
    
    obj_sub <- subset(obj, cells = selected_cells)
    
    return(obj_sub)
  }
  
  train_sub <- subset_by_group(train, group.by = "seurat_clusters", n = downSample_num)
  test_sub <- subset_by_group(test,  group.by = "cell_type", n = downSample_num)
  
  print(table(train_sub$seurat_clusters))
  print(table(test_sub$cell_type))
  
  train_sparse <- GetAssayData(train_sub, assay = "RNA", layer = "data")[features_used, ]
  test_sparse <- GetAssayData(test_sub,  assay = "RNA", layer = test_layer)[features_used, ]
  
  train_mat <- as.matrix(train_sparse)
  test_mat <- as.matrix(test_sparse)
  
  rm(train_sparse, test_sparse)
  gc()
  
  train_group <- as.character(train_sub$seurat_clusters)
  test_group  <- as.character(test_sub$cell_type)
  
  table(is.na(train_group))
  table(is.na(test_group))
  
  if (ncol(train_mat) != length(train_group)) stop("Train 维度不匹配！")
  if (ncol(test_mat)  != length(test_group))  stop("Test 维度不匹配！")
  
  ### c. calculate ----
  source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')
  
  Test_similarity <- glm.predict(
    train.data = train_mat,     
    train.group = train_group,
    genes.used = features_used,
    downsample = FALSE,
    sample.cells = 0,
    test.data = test_mat,
    test.group = test_group,
    alpha = 0.5, 
    nfolds = 5,
    seed = 123
  )
  
  heatmap_mat = Test_similarity$heatmap@matrix
  
  pdf(paste0("2.Similarity_", group_label, "_", feature_num, "_genes.pdf"), width=plot_width, height=plot_height, onefile = F)
  pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                     color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                     cluster_rows = F)
  dev.off()
  
  write.csv(heatmap_mat, paste0("2.Similarity_", group_label, "_", feature_num, "_genes.csv"), quote = FALSE, row.names = TRUE)
  save(list = c("train_mat", "test_mat", "Test_similarity", "heatmap_mat"),
       file = paste0("2.Similarity_", group_label, "_", feature_num, "_genes.RData"))
  
}



similarity_scale_pipe <- function(train_data, test_data,
                                  feature_num, downSample_num,
                                  test_layer,
                                  group_label,
                                  plot_width,
                                  plot_height){
  
  ### a. train data ----
  train <- train_data
  Idents(train) <- "cell_type"
  #train <- FindVariableFeatures(train, selection.method = 'vst', nfeatures = 2000)
  
  ### b. test data ----
  test <- test_data
  Idents(test) <- "seurat_clusters"
  test <- JoinLayers(test)
  test <- NormalizeData(test)
  #test <- FindVariableFeatures(test, selection.method = 'vst', nfeatures = 2000)
  
  sd_train <- sort(apply(GetAssayData(train, assay = "RNA", layer = "data"), 1, sd), decreasing = T)
  sd_test <- sort(apply(GetAssayData(test, assay = "RNA", layer = "scale.data"), 1, sd), decreasing = T)
  
  features_used <- intersect(names(sd_train)[1:feature_num], names(sd_test)[1:feature_num])
  
  ## subset data
  subset_by_group <- function(obj, group.by = "cell_type", n = 200, seed = 2025) {
    set.seed(seed)
    
    cells_by_group <- split(colnames(obj), obj[[group.by]])
    
    selected_cells <- unlist(lapply(cells_by_group, function(cells) {
      if (length(cells) > n) {
        return(sample(cells, n))
      } else {
        return(cells)
      }
    }))
    
    obj_sub <- subset(obj, cells = selected_cells)
    
    return(obj_sub)
  }
  
  train_sub <- subset_by_group(train, group.by = "cell_type", n = downSample_num)
  test_sub <- subset_by_group(test,  group.by = "seurat_clusters", n = downSample_num)
  
  print(table(train_sub$cell_type))
  print(table(test_sub$seurat_clusters))
  
  train_sparse <- GetAssayData(train_sub, assay = "RNA", layer = "data")[features_used, ]
  test_sparse <- GetAssayData(test_sub,  assay = "RNA", layer = "scale.data")[features_used, ]
  
  train_mat <- as.matrix(train_sparse)
  test_mat <- as.matrix(test_sparse)
  
  rm(train_sparse, test_sparse)
  gc()
  
  train_group <- as.character(train_sub$cell_type)
  test_group  <- as.character(test_sub$seurat_clusters)
  
  table(is.na(train_group))
  table(is.na(test_group))
  
  if (ncol(train_mat) != length(train_group)) stop("Train 维度不匹配！")
  if (ncol(test_mat)  != length(test_group))  stop("Test 维度不匹配！")
  
  ### c. calculate ----
  source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')
  
  Test_similarity <- glm.predict(
    train.data = train_mat,     
    train.group = train_group,
    genes.used = features_used,
    downsample = FALSE,
    sample.cells = 0,
    test.data = test_mat,
    test.group = test_group,
    alpha = 0.5, 
    nfolds = 5,
    seed = 123
  )
  
  heatmap_mat = Test_similarity$heatmap@matrix
  
  pdf(paste0("2.Similarity_", group_label, "_", feature_num, "_genes.pdf"), width=plot_width, height=plot_height, onefile = F)
  pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                     color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                     cluster_rows = F)
  dev.off()
  
  write.csv(heatmap_mat, paste0("2.Similarity_", group_label, "_", feature_num, "_genes.csv"), quote = FALSE, row.names = TRUE)
  save(list = c("train_mat", "test_mat", "Test_similarity", "heatmap_mat"),
       file = paste0("2.Similarity_", group_label, "_", feature_num, "_genes.RData"))
  
}









similarity_anno_pipe <- function(train_data, test_data,
                                 feature_num, downSample_num,
                                 test_layer,
                                 group_label,
                                 plot_width,
                                 plot_height){
  
  ### a. train data ----
  train <- train_data
  Idents(train) <- "cell_type"
  
  ### b. test data ----
  test <- test_data
  Idents(test) <- "cell_type"
  test <- JoinLayers(test)
  test <- NormalizeData(test)
  sd_train <- sort(apply(GetAssayData(train, assay = "RNA", layer = "data"), 1, sd), decreasing = T)
  sd_test <- sort(apply(GetAssayData(test, assay = "RNA", layer = test_layer), 1, sd), decreasing = T)
  
  features_used <- intersect(names(sd_train)[1:feature_num], names(sd_test)[1:feature_num])
  
  ## subset data
  subset_by_group <- function(obj, group.by = "cell_type", n = 200, seed = 2025) {
    set.seed(seed)
    
    cells_by_group <- split(colnames(obj), obj[[group.by]])
    
    selected_cells <- unlist(lapply(cells_by_group, function(cells) {
      if (length(cells) > n) {
        return(sample(cells, n))
      } else {
        return(cells)
      }
    }))
    
    obj_sub <- subset(obj, cells = selected_cells)
    
    return(obj_sub)
  }
  
  train_sub <- subset_by_group(train, group.by = "cell_type", n = downSample_num)
  test_sub <- subset_by_group(test,  group.by = "cell_type", n = downSample_num)
  
  print(table(train_sub$cell_type))
  print(table(test_sub$cell_type))
  
  train_sparse <- GetAssayData(train_sub, assay = "RNA", layer = "data")[features_used, ]
  test_sparse <- GetAssayData(test_sub,  assay = "RNA", layer = test_layer)[features_used, ]
  
  train_mat <- as.matrix(train_sparse)
  test_mat <- as.matrix(test_sparse)
  
  rm(train_sparse, test_sparse)
  gc()
  
  train_group <- as.character(train_sub$cell_type)
  test_group  <- as.character(test_sub$cell_type)
  
  table(is.na(train_group))
  table(is.na(test_group))
  
  if (ncol(train_mat) != length(train_group)) stop("Train 维度不匹配！")
  if (ncol(test_mat)  != length(test_group))  stop("Test 维度不匹配！")
  
  ### c. calculate ----
  source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')
  
  Test_similarity <- glm.predict(
    train.data = train_mat,     
    train.group = train_group,
    genes.used = features_used,
    downsample = FALSE,
    sample.cells = 0,
    test.data = test_mat,
    test.group = test_group,
    alpha = 0.5, 
    nfolds = 5,
    seed = 123
  )
  
  heatmap_mat = Test_similarity$heatmap@matrix
  
  pdf(paste0("3.Similarity_", group_label, "_", feature_num, "_genes.pdf"), width=plot_width, height=plot_height, onefile = F)
  pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                     color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                     cluster_rows = F)
  dev.off()
  
  write.csv(heatmap_mat, paste0("3.Similarity_", group_label, "_", feature_num, "_genes.csv"), quote = FALSE, row.names = TRUE)
  save(list = c("train_mat", "test_mat", "Test_similarity", "heatmap_mat"),
       file = paste0("3.Similarity_", group_label, "_", feature_num, "_genes.RData"))
  
}


