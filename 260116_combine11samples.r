# Title: wuhua project analysis
# Author: Gaozy
# Time: 2026-01-16

# zygao02@comput172-general_env-R
# 0. Basic settings ----
if (dir.exists("/home/zygao02/wuhua_proj/260116/") == F) dir.create("/home/zygao02/wuhua_proj/260116/")
setwd("/home/zygao02/wuhua_proj/260116/")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggsci)
  library(scCustomize)
  library(randomcoloR)
  library(RColorBrewer)
  library(glmnet)
  library(ComplexHeatmap)
  library(readxl)
  library(harmony)
  library(stringr)
  library(data.table)
  library(circlize)
})

# ciToti1、ciToti2、ciToti3、ciToti4、E35、E45、EPSC-S4、EPSC-blastoid、TBLC-blastoid一起整合，看一下效果
ciToti1_4 <- readRDS("/home/lushi02/project/wuhua/ciToti1_4.RDS")
EPS_blastoid  <- readRDS("/home/lushi02/project/wuhua/EPS_blastoid.RDS")
TBLC_blastoid <- readRDS("/home/lushi02/project/wuhua/TBLC_blastoid.RDS")

comb <- readRDS("/home/lushi02/project/wuhua/Combined_query_data.RDS")
ciToti4 <- subset(comb, orig.ident %in% "ciToti4")
D_E35   <- subset(comb, orig.ident %in% "D_E35")
D_E45   <- subset(comb, orig.ident %in% "D_E45")
EPSC_S4 <- subset(comb, orig.ident %in% "EPSC_S4")

# merge
obj_list <- list(EPS_blastoid, TBLC_blastoid, ciToti4, D_E35, D_E45, EPSC_S4)

ids <- c("ciToti1_4", "EPS_blastoid", "TBLC_blastoid", "ciToti4", "D_E35", "D_E45", "EPSC_S4")

merged_obj <- merge(
  x = ciToti1_4,
  y = obj_list,
  add.cell.ids = ids,
  project = "Merged_Project"
)

print(merged_obj)

# QC
qc_pipeline <- function(seurat_obj, group_name){
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- JoinLayers(seurat_obj)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  seurat_obj <- UpdateSeuratObject(seurat_obj)
  
  plot_data <- seurat_obj@meta.data
  p1 <- ggplot(plot_data, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) + 
    geom_violin(trim = FALSE) + 
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle("nFeature_RNA")
  
  p2 <- ggplot(plot_data, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) + 
    geom_violin(trim = FALSE) + 
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle("nCount_RNA")
  
  p3 <- ggplot(plot_data, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) + 
    geom_violin(trim = FALSE) + 
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle("percent.mt")
  
  pdf(paste0(group_name, "_vlnplot_preQC.pdf"), width = 15, height = 3.5)
  print(p1 + p2 + p3)
  dev.off()
  
  seurat_obj_sub <- subset(seurat_obj, subset = nFeature_RNA > 200  & percent.mt < 25)
  plot_data <- seurat_obj_sub@meta.data
  p1 <- ggplot(plot_data, aes(x = orig.ident, y = nFeature_RNA, fill = orig.ident)) + 
    geom_violin(trim = FALSE) + 
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle("nFeature_RNA")
  
  p2 <- ggplot(plot_data, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) + 
    geom_violin(trim = FALSE) + 
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle("nCount_RNA")
  
  p3 <- ggplot(plot_data, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) + 
    geom_violin(trim = FALSE) + 
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle("percent.mt")
  
  pdf(paste0(group_name, "_vlnplot_postQC.pdf"), width = 15, height = 3.5)
  print(p1 + p2 + p3)
  dev.off()
  
  return(seurat_obj_sub)
}
umap_pipeline <- function(seurat_obj, group_label, plot_title_group, plot_width, plot_height){
  p1 <- DimPlot(seurat_obj, 
                reduction = "umap", 
                group.by = "orig.ident",
                pt.size = 0.5, 
                alpha = 0.6,
                label = TRUE,
                repel = TRUE,
                raster = FALSE) +
    NoAxes() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2 <- DimPlot(seurat_obj, 
                reduction = "umap", 
                group.by = "seurat_clusters",
                pt.size = 0.5, 
                alpha = 0.6,
                label = TRUE,
                repel = TRUE,
                raster = FALSE) +
    NoAxes() +
    theme(plot.title = element_text(hjust = 0.5))
  
  combined_plot <- (p1 | p2) + 
    plot_annotation(
      title = paste0("UMAP for Group ", plot_title_group),
      theme = theme(
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
      ) )
  
  ggsave(paste0("1.UMAP_", group_label, ".pdf"), plot = combined_plot, width = plot_width, height = plot_height)
  
  p3 <- DimPlot(seurat_obj, reduction = "umap", split.by = "seurat_clusters",
                pt.size = 0.5, 
                alpha = 0.6,
                label = TRUE,
                repel = TRUE,
                raster = FALSE) +
    NoAxes() +
    theme(plot.title = element_text(hjust = 0.5))
  
  p4 <- DimPlot(seurat_obj, reduction = "umap", split.by = "orig.ident",
                pt.size = 0.5, 
                alpha = 0.6,
                label = TRUE,
                repel = TRUE,
                raster = FALSE) +
    NoAxes() +
    theme(plot.title = element_text(hjust = 0.5))
  
  combined_plot_split <- (p3 / p4) + 
    plot_annotation(
      title = paste0("Split UMAP for Group ", plot_title_group),
      theme = theme(
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
      ) )
  
  ggsave(paste0("1.UMAP_", group_label, "_split.pdf"), plot = combined_plot_split, width = plot_width, height = plot_height)
  
}



merged_qc <- qc_pipeline(merged_obj, "group2_raw_merge")
merged_qc[["percent.mt"]] <- PercentageFeatureSet(merged_qc, pattern = "^mt-")

merged_qc <- NormalizeData(merged_qc) |> 
  FindVariableFeatures(nfeatures = 2000) |> 
  ScaleData(vars.to.regress = "percent.mt") |> 
  RunPCA(npcs = 50)

only_harmony <- RunHarmony(merged_qc, 
                           group.by.vars = "orig.ident", 
                           reduction = "pca", 
                           dims.use = 1:20, 
                           reduction.save = "harmony")

only_harmony <- RunUMAP(only_harmony, reduction = "harmony", dims = 1:20) |> 
  FindNeighbors(reduction = "harmony", dims = 1:20) |> 
  FindClusters(resolution = 0.5)

umap_pipeline(seurat_obj = only_harmony,
              group_label = "only_harmony", 
              plot_title_group = "7",
              plot_width = 12, plot_height = 6)

merged_qc[["percent.mt"]] <- PercentageFeatureSet(merged_qc, pattern = "^mt-")
merged_qc <- NormalizeData(merged_qc) |> FindVariableFeatures(nfeatures = 2000)

merged_qc <- CellCycleScoring(merged_qc, 
                                     s.features = str_to_title(cc.genes$s.genes), 
                                     g2m.features = str_to_title(cc.genes$g2m.genes), 
                                     set.ident = TRUE)

merged_qc_regress <- ScaleData(merged_qc, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) |> 
  RunPCA(npcs = 50)

## sample level
regress_harmony <- RunHarmony(merged_qc_regress, 
                              group.by.vars = "orig.ident", 
                              reduction = "pca", 
                              dims.use = 1:20, 
                              reduction.save = "harmony")


regress_harmony <- RunUMAP(regress_harmony, reduction = "harmony", dims = 1:20) |> 
  FindNeighbors(reduction = "harmony", dims = 1:20) |> 
  FindClusters(resolution = 0.5)

regress_harmony <- FindClusters(regress_harmony, resolution = 0.1)

umap_pipeline(seurat_obj = regress_harmony,
              group_label = "regress_harmony", 
              plot_title_group = "7",
              plot_width = 12, plot_height = 6)

only_harmony$orig.ident[only_harmony$orig.ident == "etblc"] <- "TBLC_blastoid"
regress_harmony$orig.ident[regress_harmony$orig.ident == "etblc"] <- "TBLC_blastoid"

saveRDS(only_harmony, "only_harmony.rds")

# Similarity
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(glmnet)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
options(stringsAsFactors = FALSE)

library(Hmisc)
group1_ref <- readRDS("/home/lushi02/project/wuhua/ref/group1ref_E3.5_E4.5.rds")
sub_regress_harmony <- subset(regress_harmony, orig.ident %nin% c("ciToti1", "ciToti2", "ciToti3"))

## ---------------------------
## 0) glm.predict (your version + maxit)
## ---------------------------
glm.predict <-
  function(train.data, train.group, downsample = FALSE, sample.cells = 0,
           genes.used = NA, test.data, test.group, alpha = 0.99, nfolds = 10) {
    
    require(glmnet)
    require(ComplexHeatmap)
    
    glm.fits <- list()
    glm.pred <- list()
    
    if (length(genes.used) > 1) {
      train.data <- train.data[genes.used, , drop = FALSE]
      test.data  <- test.data[genes.used,  , drop = FALSE]
      if (length(genes.used) <= 50) cat("There were less than 50 features used!\n")
    }
    
    if (sample.cells == 0 & downsample) {
      sample.cells <- max(50, min(table(train.group)))
    }
    
    if (sample.cells > 0) {
      ngroup <- length(unique(train.group))
      if (ncol(train.data) >= sample.cells * ngroup) {
        cells_used <- c()
        for (groupi in sort(unique(train.group))) {
          idx <- which(train.group == groupi)
          if (length(idx) > sample.cells) {
            cells_used <- c(cells_used, sample(idx, sample.cells))
          } else {
            cells_used <- c(cells_used, idx)
          }
        }
        train.data  <- train.data[, cells_used, drop = FALSE]
        train.group <- train.group[cells_used]
      }
    }
    
    for (groupi in sort(unique(train.group))) {
      fac <- factor(train.group == groupi)
      glm.fits[[groupi]] <-
        cv.glmnet(
          x = t(train.data),
          y = fac,
          offset = getPopulationOffset(fac),
          family = "binomial",
          intercept = FALSE,
          alpha = alpha,
          nfolds = nfolds,
          type.measure = "class",
          maxit = 1e6
        )
      
      glm.pred[[groupi]] <-
        predict(
          object = glm.fits[[groupi]],
          newx = t(test.data),
          newoffset = rep(0, ncol(test.data)),
          s = "lambda.min"
        )
    }
    
    glm.pred.df <- data.frame(do.call(cbind, glm.pred))
    colnames(glm.pred.df) <- sort(unique(train.group))
    
    glm.pred.df.prob <- (1 + exp(-glm.pred.df)) ** -1
    glm.cluster <- colnames(glm.pred.df.prob)[apply(glm.pred.df.prob, 1, which.max)]
    
    return(list(
      test.group  = test.group,
      logits      = glm.pred.df,
      probability = glm.pred.df.prob,
      cluster     = glm.cluster
    ))
  }

getPopulationOffset <- function(y) {
  if (!is.factor(y)) y <- factor(y)
  if (length(levels(y)) != 2) stop("y must be a two-level factor")
  off <- sum(y == levels(y)[2]) / length(y)
  off <- log(off / (1 - off))
  rep(off, length(y))
}

## ---------------------------
## 1) Utility: objective scoring for a ref x query similarity matrix
## ---------------------------
score_heatmap <- function(mat_plot) {
  # mat_plot: ref(rows) x query(cols)
  common <- intersect(rownames(mat_plot), colnames(mat_plot))
  if (length(common) < 2) return(NA_real_)
  
  d <- diag(mat_plot[common, common, drop = FALSE])
  
  offmax <- sapply(common, function(ct) {
    others <- setdiff(rownames(mat_plot), ct)
    if (length(others) == 0) return(NA_real_)
    max(mat_plot[others, ct, drop = TRUE], na.rm = TRUE)
  })
  
  mean(d - offmax, na.rm = TRUE)
}

## ---------------------------
## 2) Prepare train/test: filter ref types by test types
## ---------------------------
set.seed(1)  # for cv.glmnet fold randomness

train <- group1_ref
Idents(train) <- "cell_type"
train.group <- as.character(Idents(train))

test <- sub_regress_harmony
Idents(test) <- "seurat_clusters"
test.group <- as.character(Idents(test))

cat("Train cells:", ncol(train), " | Test cells:", ncol(test), "\n")
cat("Train types:", length(unique(train.group)), " | Test types:", length(unique(test.group)), "\n")

## ---------------------------
## 3) Extract matrices (NO get_expr_mat)
##     train: data
##     test : prefer scale.data if exists, else data
## ---------------------------
# Seurat v5 用 layer，v4 用 slot；这里用 tryCatch 做最小兼容，但不调用 get_expr_mat
mat_train <- tryCatch(GetAssayData(train, layer = "data"),
                      error = function(e) GetAssayData(train, slot = "data"))

mat_test <- tryCatch(GetAssayData(test, layer = "data"),
                     error = function(e) {GetAssayData(test, slot = "data")
                     })

# SD vectors
sd_train <- sort(apply(mat_train, 1, sd), decreasing = TRUE)
sd_test  <- sort(apply(mat_test,  1, sd), decreasing = TRUE)

## ---------------------------
## 4) Remove rare train types (<5)
## ---------------------------
min_cells_train <- 5
tab_train <- table(train.group)
rare_types <- names(tab_train)[tab_train < min_cells_train]
if (length(rare_types) > 0) cat("Rare train types removed:", paste(rare_types, collapse = ", "), "\n")

keep_train_cells     <- which(!train.group %in% rare_types)
mat_train_filtered   <- mat_train[, keep_train_cells, drop = FALSE]
train.group_filtered <- train.group[keep_train_cells]

features_used <- intersect(rownames(mat_train), rownames(mat_test)) %>% unique()
length(features_used)

## 3) 调用 glm.predict（注意这里 test.group 先给“每个细胞一个标签”即可）
res <- glm.predict(
  train.data   = mat_train_filtered,
  train.group  = train.group_filtered,
  downsample   = FALSE,
  sample.cells = 0,
  genes.used   = features_used,
  test.data    = mat_test,
  test.group   = test.group,  # 每个细胞一个独立 ID
  alpha        = 0.99,
  nfolds       = 10
)

## 5) 重新按 test.group 聚合，画你想要的 similarity heatmap -----
glm.predict.mean <-
  apply(res$logits, 2, function(e)
    sapply(split(e, res$test.group), mean))
glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1

col_fun <- colorRamp2(c(0, 0.5, 1), c("#e9e9e9", "white", "red"))

mat_plot <- t(glm.predict.mean.prob)

# ====== 1) 设定每个格子的边长（mm）======
cell_mm <- 16  # 你可以试 4/5/6；越大越清晰但PDF越大

hm_w <- unit(ncol(mat_plot) * cell_mm, "mm")
hm_h <- unit(nrow(mat_plot) * cell_mm, "mm")

# ====== 2) 为行名/列名预留空间（防出框）======
cn <- colnames(mat_plot)
rn <- rownames(mat_plot)
col_h <- max_text_width(cn, gp = gpar(fontsize = 14)) + unit(8, "mm")
row_w <- max_text_width(rn, gp = gpar(fontsize = 14)) + unit(8, "mm")

# ====== 3) PDF 尺寸：按格子数粗略撑大（单位：英寸）======
# 经验值：左右上下再给一些边距 + legend
pdf_w_in <- (ncol(mat_plot) * cell_mm + 140) / 25.4
pdf_h_in <- (nrow(mat_plot) * cell_mm + 120) / 25.4

pdf("glm_predict_similarity_heatmap.pdf",
    width = pdf_w_in, height = pdf_h_in, useDingbats = FALSE)

ht <- Heatmap(
  mat_plot,
  col = col_fun,
  name = "Predicted\nSimilarity",
  column_title = "query data",
  row_title = "ref data",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # ✅ 方块正方形：固定 heatmap body 的宽高比例
  width  = hm_w,
  height = hm_h,
  
  # 字体
  row_title_gp = gpar(fontsize = 16),
  column_title_gp = gpar(fontsize = 16),
  row_names_gp = gpar(fontsize = 14),
  column_names_gp = gpar(fontsize = 12), 
  
  # ✅ 防止文字出框
  row_names_max_width = row_w,
  column_names_max_height = col_h,
  
  # 列名太长就旋转
  column_names_rot = 45,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- mat_plot[i, j]
    grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = 10, col = "black"))
  }

)

draw(
  ht,
  heatmap_legend_side = "right",
  padding = unit(c(8, 35, 25, 25), "mm")  # 上 右 下 左
)

dev.off()

saveRDS(regress_harmony, "regress_harmony.rds")


regress_harmony_markers <- FindAllMarkers(regress_harmony, assay = "RNA",
                                          layer = "data", group.by = "seurat_clusters")

write.csv(regress_harmony_markers, "regress_harmony_markers.csv", quote = F, row.names = T)


regress_harmony_markers_scale <- FindAllMarkers(regress_harmony, assay = "RNA",
                                          layer = "scale.data", group.by = "seurat_clusters")

write.csv(regress_harmony_markers_scale, "regress_harmony_markers_scale.csv", quote = F, row.names = T)













