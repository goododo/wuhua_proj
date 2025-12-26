# Title: wuhua project analysis
# Author: Gaozy
# Time: 2025-12-25

# zygao02@comput172-general_env-R
# 0. Basic settings ----
if (dir.exists("/home/zygao02/wuhua_proj/251225/") == F) dir.create("/home/zygao02/wuhua_proj/251225/")
setwd("/home/zygao02/wuhua_proj/251225/")

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
})
source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

## load combined datasets & seperate ----
comb <- readRDS("/home/lushi02/project/wuhua/Combined_query_data.RDS")

ciToti4 <- subset(comb, orig.ident %in% c("ciToti4"))
ciToti7 <- subset(comb, orig.ident %in% c("ciToti7"))
ciToti8 <- subset(comb, orig.ident %in% c("ciToti8"))
ciToti10 <- subset(comb, orig.ident %in% c("ciToti10"))

BPSCEM_day8 <- subset(comb, orig.ident %in% c("BPSCEM_day8"))
D_E35 <- subset(comb, orig.ident %in% c("D_E35"))
D_E45 <- subset(comb, orig.ident %in% c("D_E45"))
D_E65 <- subset(comb, orig.ident %in% c("D_E65"))
D_E75 <- subset(comb, orig.ident %in% c("D_E75"))
D_E85 <- subset(comb, orig.ident %in% c("D_E85"))

EPSC_S4 <- subset(comb, orig.ident %in% c("EPSC_S4"))
EPSC_S7 <- subset(comb, orig.ident %in% c("EPSC_S7"))
EPSC_S8 <- subset(comb, orig.ident %in% c("EPSC_S8"))

ETiX5 <- subset(comb, orig.ident %in% c("ETiX5"))
ETiX6 <- subset(comb, orig.ident %in% c("ETiX6"))
ETiX8 <- subset(comb, orig.ident %in% c("ETiX8"))

H_E75 <- subset(comb, orig.ident %in% c("H_E75"))
H_E85 <- subset(comb, orig.ident %in% c("H_E85"))
Hanna_2022_EM_day8 <- subset(comb, orig.ident %in% c("Hanna_2022_EM_day8"))

I_E35 <- subset(comb, orig.ident %in% c("I_E35"))
iEFCEM_day6 <- subset(comb, orig.ident %in% c("iEFCEM_day6"))
iEFCEM_day8 <- subset(comb, orig.ident %in% c("iEFCEM_day8"))

TFSEM_day10 <- subset(comb, orig.ident %in% c("TFSEM_day10"))
TFSEM_day8 <- subset(comb, orig.ident %in% c("TFSEM_day8"))

W_E65 <- subset(comb, orig.ident %in% c("W_E65"))
W_E75 <- subset(comb, orig.ident %in% c("W_E75"))
W_E85 <- subset(comb, orig.ident %in% c("W_E85"))

Z_E45 <- subset(comb, orig.ident %in% c("Z_E45"))
Z_E65 <- subset(comb, orig.ident %in% c("Z_E65"))
Z_E75 <- subset(comb, orig.ident %in% c("Z_E75"))
Z_E85 <- subset(comb, orig.ident %in% c("Z_E85"))

## Harmony ----
harmony_pipe <- function(seurat_obj, reso){
  seurat_obj <- JoinLayers(seurat_obj)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "percent.mt")
  seurat_obj <- RunPCA(seurat_obj,  npcs = 50)
  seurat_obj <- RunHarmony(seurat_obj,'orig.ident')
  seurat_obj <- RunUMAP(seurat_obj,reduction = "harmony", dims=1:20)
  seurat_obj <- FindNeighbors(seurat_obj,reduction = "harmony",  dims=1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = reso)
  return(seurat_obj)
}

## Regress MT- & Cell cycle ----
regress_cellcycle_harmony_pipe <- function(seurat_obj, reso) {
  # 1) Join layers & basic QC
  seurat_obj <- JoinLayers(seurat_obj)
  
  # 计算线粒体百分比
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # 2) 归一化与高变基因识别
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # 3) 细胞周期评分
  # 使用 Seurat 内置的 S 和 G2/M 标记基因列表
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  seurat_obj <- CellCycleScoring(
    seurat_obj,
    s.features  = s.genes,
    g2m.features = g2m.genes,
    set.ident = FALSE
  )
  
  # 4) ScaleData：回归掉 mt% 与细胞周期得分
  seurat_obj <- ScaleData(
    seurat_obj,
    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
    features = rownames(seurat_obj)
  )
  
  # 5) PCA + Harmony
  seurat_obj <- RunPCA(seurat_obj, npcs = 50)
  seurat_obj <- RunHarmony(seurat_obj, 'orig.ident')
  
  # 6) UMAP & Clustering
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:20)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = reso)
  
  return(seurat_obj)
}

# 1. Group 1 ----
group1_merged <- merge(x = ciToti4, 
                       y = c(EPSC_S4, I_E35, Z_E45, D_E35, D_E45),
                       add.cell.ids = c("ciToti4", "EPSC_S4", "I_E35", "Z_E45", "D_E35", "D_E45"))

## Regress : reso = 0.06 ----
group1_cellcyc <- regress_cellcycle_harmony_pipe(group1_merged, reso = 0.04)
table(group1_cellcyc$seurat_clusters)

group1_cellcyc <- FindClusters(group1_cellcyc, resolution = 0.06)
table(group1_cellcyc$seurat_clusters)

## UMAP ----
p1 <- DimPlot(group1_cellcyc, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group1_cellcyc, 
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
    title = "UMAP for Group 1",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_group1_cellcyc.pdf", plot = combined_plot, width = 12, height = 6)

p3 <- DimPlot(group1_cellcyc, reduction = "umap", split.by = "seurat_clusters",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p4 <- DimPlot(group1_cellcyc, reduction = "umap", split.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot_split <- (p3 / p4) + 
  plot_annotation(
    title = "Split UMAP for Group 1",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_group1_cellcyc_split.pdf", plot = combined_plot_split, width = 12, height = 6)

## Similarity ----
group1_ref <- readRDS("/home/lushi02/project/wuhua/ref/group1ref_E3.5_E4.5.rds")
table(group1_ref$cell_type)

RhpcBLASctl::blas_set_num_threads(10)

similarity_pipe(train_data = group1_ref, test_data = group1_cellcyc,
                feature_num = 2000, downSample_num = 2000, 
                test_layer = "scale.data",
                plot_width = 8,
                plot_height = 2.5,
                group_label = "group1_scale")

similarity_pipe(train_data = group1_ref, test_data = group1_cellcyc,
                feature_num = 1000, downSample_num = 2000,
                test_layer = "scale.data",
                plot_width = 8,
                plot_height = 2.5,
                group_label = "group1_scale")

similarity_pipe(train_data = group1_ref, test_data = group1_cellcyc,
                feature_num = 500, downSample_num = 2000,
                test_layer = "scale.data",
                plot_width = 8,
                plot_height = 2.5,
                group_label = "group1_scale")

similarity_pipe(train_data = group1_ref, test_data = group1_cellcyc,
                feature_num = 200, downSample_num = 2000,
                test_layer = "scale.data",
                plot_width = 8,
                plot_height = 2.5,
                group_label = "group1_scale")

similarity_pipe(train_data = group1_ref, test_data = group1_cellcyc,
                feature_num = 100, downSample_num = 2000,
                test_layer = "scale.data",
                plot_width = 8,
                plot_height = 2.5,
                group_label = "group1_scale")

## Annotation according to Similarity ----
group1_cellcyc@meta.data <- group1_cellcyc@meta.data %>%
  mutate(cell_type = case_when(
    RNA_snn_res.0.06 %in% c(0, 3, 4) ~ "Trophoectoderm",
    RNA_snn_res.0.06 %in% c(1,5) ~ "Epiblast/Inner_cell_mass",
    RNA_snn_res.0.06 == 2 ~ "Primitive_endoderm",
    TRUE ~ as.character(RNA_snn_res.0.06)
  ))

group1_markers <- FindAllMarkers(group1_cellcyc, group.by = "cell_type")

write.csv(group1_markers, "group1_markers.csv", quote = FALSE, row.names = TRUE)


group1_markers <- FindAllMarkers(group1_cellcyc, group.by = "seurat_clusters")
write.csv(group1_markers, "group1_markers_seurat_clusters.csv", quote = FALSE, row.names = TRUE)








similarity_anno_pipe(train_data = group1_ref,
                     test_data = group1_cellcyc,
                     feature_num = 200,
                     downSample_num = 2000,
                     test_layer = "scale.data",
                     group_label = "group1_with_anno",
                     plot_width = 3.8,
                     plot_height = 3.5)

saveRDS(group1_cellcyc, "group1_cellcyc_with_anno.rds")

#load("3.Similarity_group1_with_anno_200_genes.RData")






# 2. Group 2 ----
group2_merged <- merge(x = ciToti7, 
                       y = c(ETiX5, Z_E65, D_E65),
                       add.cell.ids = c("ciToti7", "ETiX5", "Z_E65", "D_E65"))

## Regress ----
group2_harmony <- harmony_pipe(group2_harmony, reso = 0.04)
table(group2_harmony$seurat_clusters)

group2_harmony <- FindClusters(group2_harmony, resolution = 0.1)
table(group2_harmony$seurat_clusters)

## UMAP ----
p1 <- DimPlot(group2_harmony, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group2_harmony, 
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
    title = "UMAP for Group 2",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_group2_harmony.pdf", plot = combined_plot, width = 12, height = 6)

p3 <- DimPlot(group2_harmony, reduction = "umap", split.by = "seurat_clusters",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p4 <- DimPlot(group2_harmony, reduction = "umap", split.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot_split <- (p3 / p4) + 
  plot_annotation(
    title = "Split UMAP for Group 2",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_group2_harmony_split.pdf", plot = combined_plot_split, width = 12, height = 6)


## Similarity ----
group2_ref <- readRDS("/home/lushi02/project/wuhua/ref/group2ref_E6.5.rds")
table(group2_ref$cell_type)

RhpcBLASctl::blas_set_num_threads(10)


group2_harmony <- FindClusters(group2_harmony, resolution = 0.5)
table(group2_harmony$seurat_clusters)


similarity_pipe(train_data = group2_ref, test_data = group2_harmony,
                feature_num = 2000, downSample_num = 2000,
                group_label = "group2")


similarity_pipe(train_data = group2_ref, test_data = group2_harmony,
                feature_num = 1000, downSample_num = 2000,
                group_label = "group2")


similarity_pipe(train_data = group2_ref, test_data = group2_harmony,
                feature_num = 500, downSample_num = 2000,
                group_label = "group2")


similarity_pipe(train_data = group2_ref, test_data = group2_harmony,
                feature_num = 200, downSample_num = 2000,
                group_label = "group2")


similarity_pipe(train_data = group2_ref, test_data = group2_harmony,
                feature_num = 100, downSample_num = 2000,
                group_label = "group2")



similarity_t_pipe(train_data = group2_harmony, test_data = group2_ref,
                feature_num = 2000, downSample_num = 2000,
                group_label = "group2_t")




similarity_scale_pipe(train_data = group2_ref, test_data = group2_harmony,
                  feature_num = 2000, downSample_num = 2000,
                  group_label = "group2_scale",
                  plot_width = 6,
                  plot_height = 3)








# 3. Group 3 ----
group3_merged <- merge(x = ciToti8, 
                       y = c(EPSC_S7, ETiX6, iEFCEM_day6, TFSEM_day8, H_E75, W_E75, Z_E75, D_E75),
                       add.cell.ids = c("ciToti8", "EPSC_S7", "ETiX6", "iEFCEM_day6", "TFSEM_day8", "H_E75", "W_E75", "Z_E75", "D_E75"))

## Regress ----
group3_harmony <- harmony_pipe(group3_merged, reso = 0.2)
table(group3_harmony$seurat_clusters)

## UMAP ----
p1 <- DimPlot(group3_harmony, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group3_harmony, 
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
    title = "UMAP for Group 3",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_group3_harmony.pdf", plot = combined_plot, width = 12, height = 6)

p3 <- DimPlot(group3_harmony, reduction = "umap", split.by = "seurat_clusters",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p4 <- DimPlot(group3_harmony, reduction = "umap", split.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot_split <- (p3 / p4) + 
  plot_annotation(
    title = "Split UMAP for Group 3",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_group3_harmony_split.pdf", plot = combined_plot_split, width = 22, height = 6)


## Similarity ----
group3_ref <- readRDS("/home/lushi02/project/wuhua/ref/group3ref_E7.5.rds")
table(group3_ref$cell_type)

RhpcBLASctl::blas_set_num_threads(10)

### a. train data ----
train <- group3_ref
Idents(train) <- "cell_type"
train <- FindVariableFeatures(train, selection.method = 'vst', nfeatures = 2000)

### b. test data ----
test <- group3_harmony
Idents(test) <- "seurat_clusters"
test <- JoinLayers(test)
test <- NormalizeData(test)
test <- FindVariableFeatures(test, selection.method = 'vst', nfeatures = 2000)

sd_train <- sort(apply(GetAssayData(train, assay = "RNA", layer = "data"), 1, sd), decreasing = T)
sd_test <- sort(apply(GetAssayData(test, assay = "RNA", layer = "data"), 1, sd), decreasing = T)

features_used <- intersect(names(sd_train)[1:2000], names(sd_test)[1:2000])

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

train_sub <- subset_by_group(train, group.by = "cell_type", n = 2000)
test_sub <- subset_by_group(test,  group.by = "seurat_clusters", n = 2000)

print(table(train_sub$cell_type))
print(table(test_sub$seurat_clusters))

train_sparse <- GetAssayData(train_sub, assay = "RNA", layer = "data")[features_used, ]
test_sparse <- GetAssayData(test_sub,  assay = "RNA", layer = "data")[features_used, ]

train_mat <- as.matrix(train_sparse)
test_mat <- as.matrix(test_sparse)

### 去除 只有 1/2 个细胞的 cluster
valid_groups <- names(which(table(train_group) >= 5))
print(paste("Removing groups:", setdiff(unique(train_group), valid_groups)))

cells_to_keep <- train_group %in% valid_groups

train_mat_filtered <- train_mat[, cells_to_keep]
train_group_filtered <- train_group[cells_to_keep]

train_group_filtered <- droplevels(as.factor(train_group_filtered))

table(is.na(train_group_filtered))
table(is.na(test_group))

if (ncol(train_mat_filtered) != length(train_group_filtered)) stop("Train 维度不匹配！")
if (ncol(test_mat)  != length(test_group))  stop("Test 维度不匹配！")

### c. calculate ----
source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

Test_similarity <- glm.predict(
  train.data = train_mat_filtered,     
  train.group = train_group_filtered,
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

pdf("2.Similarity_group3.pdf", width=5, height=5.5, onefile = F)
pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                   color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                   cluster_rows = F)
dev.off()

write.csv(heatmap_mat, "2.Similarity_group3.csv", quote = FALSE, row.names = TRUE)
save(list = c("train_mat", "test_mat", "Test_similarity", "heatmap_mat"),
     file = "2.Similarity_group3.RData")
