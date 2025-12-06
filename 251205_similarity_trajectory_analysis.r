# Title: wuhua project analysis
# Author: Gaozy
# Time: 2025-11-25

# zygao02@comput172-general_env-R
# 0. Basic settings ----
if (dir.exists("/home/zygao02/wuhua_proj/251125/") == F) dir.create("/home/zygao02/wuhua_proj/251125/")
setwd("/home/zygao02/wuhua_proj/251125/")

# 1. Plot UMAP ----
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggsci)
  library(scCustomize)
  library(randomcoloR)
  library(RColorBrewer)
})

## 1) ciToti4_10 ----
seurat_obj <- readRDS("/home/lushi02/project/wuhua/ciToti4_10.RDS")

seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)

seurat_obj

p1 <- DimPlot_scCustom(seurat_obj, 
                       group.by = "seurat_clusters",
                       colors_use = "polychrome",
                       label = TRUE,
                       figure_plot = TRUE)

p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.1, 
              alpha = 0.6) +
  scale_color_d3() +
  NoAxes() +
  ggtitle("Sample Info") +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot <- (p1 | p2) + 
  plot_annotation(
    title = "ciToti4-10",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_ciToti4_10.pdf", plot = combined_plot, width = 12.5, height = 6)

## 2) mouse_atlas_qiu_li_sozen ----
seurat_obj <- readRDS("/home/lushi02/project/wuhua/mouse_atlas_qiu_li_sozen.rds")
head(seurat_obj)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

seurat_obj

p1 <- DimPlot_scCustom(seurat_obj, 
                       group.by = "seurat_clusters",
                       colors_use = "polychrome",
                       label = TRUE,
                       figure_plot = TRUE)

p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "day",
              pt.size = 0.1, 
              alpha = 0.6,
              raster = FALSE) +
  scale_color_d3(palette = "category20") +
  NoAxes() +
  ggtitle("Sample Info") +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot <- (p1 | p2) + 
  plot_annotation(
    title = "Mouse Atlas Qiu Li Sozen",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_mouse_atlas_qiu_li_sozen.pdf", plot = combined_plot, width = 12.5, height = 6)

## 3) Combined_query_data ----
seurat_obj <- readRDS("/home/lushi02/project/wuhua/Combined_query_data.RDS")
head(seurat_obj)


seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

seurat_obj

p1 <- DimPlot_scCustom(seurat_obj, 
                       group.by = "seurat_clusters",
                       colors_use = "polychrome",
                       label = TRUE,
                       figure_plot = TRUE)

p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "group",
              pt.size = 0.1, 
              alpha = 0.6,
              raster = FALSE) +
  scale_color_d3(palette = "category20") +
  NoAxes() +
  ggtitle("Sample Info") +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot <- (p1 | p2) + 
  plot_annotation(
    title = "Combined Query Data",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_Combined_query_data.pdf", plot = combined_plot, width = 12.5, height = 6)

## 4) mouse_atlas_gottgens_stelzer ----
seurat_obj <- readRDS("/home/lushi02/project/wuhua/mouse_atlas_gottgens_stelzer.rds")
head(seurat_obj)

seurat_obj <- subset(seurat_obj, subset = !is.na(cell_type))

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

seurat_obj

cell_types <- unique(seurat_obj$cell_type)
n_colors <- length(cell_types)
my_palette <- distinctColorPalette(n_colors)

p1 <- DimPlot_scCustom(seurat_obj, 
                       group.by = "seurat_clusters",
                       colors_use = "polychrome",
                       label = TRUE,
                       figure_plot = TRUE)

p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "cell_type",
              cols = my_palette,
              pt.size = 0.1, 
              alpha = 0.6,
              raster = FALSE) +
  NoAxes() +
  ggtitle("Sample Info") +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot <- (p1 | p2) + 
  plot_annotation(
    title = "Mouse Atlas Gottgens Stelzer",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_mouse_atlas_gottgens_stelzer.pdf", plot = combined_plot, width = 18, height = 6)

## 5) ciToti1_4 ----
seurat_obj <- readRDS("/home/lushi02/project/wuhua/ciToti1_4.RDS")
head(seurat_obj)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

seurat_obj

p1 <- DimPlot_scCustom(seurat_obj, 
                       group.by = "seurat_clusters",
                       colors_use = "polychrome",
                       label = TRUE,
                       figure_plot = TRUE)

p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "celltype",
              pt.size = 0.1, 
              alpha = 0.6,
              raster = FALSE) +
  scale_color_d3(palette = "category20") +
  NoAxes() +
  ggtitle("Sample Info") +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot <- (p1 | p2) + 
  plot_annotation(
    title = "ciToti1 4",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_ciToti1_4.pdf", plot = combined_plot, width = 12.5, height = 6)

# 2. Similarity ----
library(glmnet)
library(ComplexHeatmap)
source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

## 1) reference ----
ref_1 <- readRDS("/home/lushi02/project/wuhua/mouse_atlas_gottgens_stelzer.rds")
ref_1 <- subset(ref_1, subset = !is.na(cell_type))

ref_2 <- readRDS("/home/lushi02/project/wuhua/mouse_atlas_qiu_li_sozen.rds")

merged_ref <- merge(ref_1, y = ref_2, 
                    add.cell.ids = c("Ref1", "Ref2"), 
                    project = "Combined_Check")

merged_ref$Batch <- c(rep("Mouse_Atlas_Gottgens_Stelzer", ncol(ref_1)), rep("Mouse_Atlas_Qiu_Li_Sozen", ncol(ref_2)))

table(merged_ref$Batch)

merged_ref <- JoinLayers(merged_ref)
merged_ref <- NormalizeData(merged_ref)
merged_ref <- FindVariableFeatures(merged_ref, selection.method = "vst", nfeatures = 2000)
merged_ref <- ScaleData(merged_ref)
merged_ref <- RunPCA(merged_ref)
merged_ref <- RunUMAP(merged_ref, dims = 1:30)
merged_ref <- JoinLayers(merged_ref)

p1 <- DimPlot(merged_ref, reduction = "umap", group.by = "Batch") + 
  ggtitle("Colored by Batch")

merged_ref <- FindNeighbors(merged_ref, dims = 1:30)
merged_ref <- FindClusters(merged_ref, resolution = 0.5)
p2 <- DimPlot(merged_ref, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Colored by Clusters")

combined_plot <- (p1 | p2) + 
  plot_annotation(
    title = "Ref Data Batch Check",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )
ggsave("2.UMAP_Ref_Batch_Check.pdf", plot = combined_plot, width = 12.5, height = 6)

## 2) query data ----
ciToti10 <- readRDS("/home/lushi02/project/wuhua/Annotated/ciToti10.RDS")
ciToti4 <- readRDS("/home/lushi02/project/wuhua/Annotated/ciToti4.RDS")
ciToti7 <- readRDS("/home/lushi02/project/wuhua/Annotated/ciToti7.RDS")
ciToti8 <- readRDS("/home/lushi02/project/wuhua/Annotated/ciToti8.RDS")
D_E35 <- readRDS("/home/lushi02/project/wuhua/Annotated/D_E35.RDS")
D_E65 <- readRDS("/home/lushi02/project/wuhua/Annotated/D_E65.RDS")
D_E75 <- readRDS("/home/lushi02/project/wuhua/Annotated/D_E75.RDS")
D_E85 <- readRDS("/home/lushi02/project/wuhua/Annotated/D_E85.RDS")

library(Seurat)
library(Matrix)

fix_nameless_matrix <- function(obj) {
  raw_counts <- obj@assays$RNA@layers$counts
  
  gene_names <- rownames(obj) 
  cell_names <- rownames(obj@meta.data)
  
  mat_dim <- dim(raw_counts)
  message(paste("矩阵维度:", mat_dim[1], "x", mat_dim[2]))
  message(paste("基因名数量:", length(gene_names)))
  message(paste("细胞名数量:", length(cell_names)))
  
  if (length(gene_names) != mat_dim[1]) {
    stop("严重错误：基因名数量与矩阵行数不匹配！无法修复。")
  }
  if (length(cell_names) != mat_dim[2]) {
    stop("严重错误：细胞名数量与矩阵列数不匹配！无法修复。")
  }
  
  dimnames(raw_counts) <- list(gene_names, cell_names)
  
  new_obj <- CreateSeuratObject(counts = raw_counts, meta.data = obj@meta.data)
  
  message("修复成功！")
  return(new_obj)
}

ciToti4 <- fix_nameless_matrix(ciToti4)
ciToti7 <- fix_nameless_matrix(ciToti7)
D_E35 <- fix_nameless_matrix(D_E35)
D_E65 <- fix_nameless_matrix(D_E65)
D_E75 <- fix_nameless_matrix(D_E75)

merged_query <- merge(x = ciToti10, 
                      y = c(ciToti4, ciToti7, ciToti8, D_E35, D_E65, D_E75, D_E85),
                      add.cell.ids = c("ciToti10", "ciToti4", "ciToti7", "ciToti8", 
                                       "E35", "E65", "E75", "E85"))
merged_query <- JoinLayers(merged_query)

## 3) Similarity ----
RhpcBLASctl::blas_set_num_threads(10)

### a. train data ----
train <- merged_ref
train <- FindVariableFeatures(train, selection.method = 'vst', nfeatures = 2000)

### b. test data ----
test <- merged_query
test <- JoinLayers(test)
test <- NormalizeData(test)
test <- FindVariableFeatures(test, selection.method = 'vst', nfeatures = 2000)

common_genes <- intersect(rownames(train), rownames(test))
print(paste("Common genes found:", length(common_genes)))

train_sparse <- GetAssayData(train, assay = "RNA", layer = "data")[common_genes, ]
test_sparse  <- GetAssayData(test,  assay = "RNA", layer = "data")[common_genes, ]

train_mat <- as.matrix(train_sparse)
test_mat  <- as.matrix(test_sparse)

rm(train_sparse, test_sparse)
gc()

train_group <- as.character(train$cell_type_clean)
test_group  <- as.character(test$orig.ident)

table(is.na(train_group))

if (ncol(train_mat) != length(train_group)) stop("Train 维度不匹配！")
if (ncol(test_mat)  != length(test_group))  stop("Test 维度不匹配！")

### c. calculate ----
source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

Test_similarity <- glm.predict(
  train.data = train_mat,     
  train.group = train_group,
  genes.used = common_genes,
  downsample = TRUE,
  sample.cells = 200,
  test.data = test_mat,
  test.group = test_group,
  alpha = 0.5, 
  nfolds = 5,
  seed = 123
)

heatmap_mat = t(Test_similarity2$heatmap@matrix)

pdf("2.Similarity.pdf", width=18, height=5, onefile = F)
pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                   color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                   cluster_rows = F)
dev.off()


# 3. 尝试是否可以通过Similarity分析对Combined_query_data.RDS中的其他胚胎数据做annotation ？----






# 4. Trajectory time analysis ----
library(monocle3)
## 1) load data ----
ciToti1_4 <- readRDS("/home/lushi02/project/wuhua/ciToti1_4.RDS")
ciToti10 <- readRDS("/home/lushi02/project/wuhua/Annotated/ciToti10.RDS")
ciToti4 <- readRDS("/home/lushi02/project/wuhua/Annotated/ciToti4.RDS")
ciToti7 <- readRDS("/home/lushi02/project/wuhua/Annotated/ciToti7.RDS")
ciToti8 <- readRDS("/home/lushi02/project/wuhua/Annotated/ciToti8.RDS")
D_E35 <- readRDS("/home/lushi02/project/wuhua/Annotated/D_E35.RDS")
D_E65 <- readRDS("/home/lushi02/project/wuhua/Annotated/D_E65.RDS")
D_E75 <- readRDS("/home/lushi02/project/wuhua/Annotated/D_E75.RDS")
D_E85 <- readRDS("/home/lushi02/project/wuhua/Annotated/D_E85.RDS")

library(Seurat)
library(Matrix)

fix_nameless_matrix <- function(obj) {
  raw_counts <- obj@assays$RNA@layers$counts
  
  gene_names <- rownames(obj) 
  cell_names <- rownames(obj@meta.data)
  
  mat_dim <- dim(raw_counts)
  message(paste("矩阵维度:", mat_dim[1], "x", mat_dim[2]))
  message(paste("基因名数量:", length(gene_names)))
  message(paste("细胞名数量:", length(cell_names)))
  
  if (length(gene_names) != mat_dim[1]) {
    stop("基因名数量与矩阵行数不匹配！")
  }
  if (length(cell_names) != mat_dim[2]) {
    stop("细胞名数量与矩阵列数不匹配！")
  }
  
  dimnames(raw_counts) <- list(gene_names, cell_names)
  
  new_obj <- CreateSeuratObject(counts = raw_counts, meta.data = obj@meta.data)
  
  message("修复成功！")
  return(new_obj)
}

ciToti4 <- fix_nameless_matrix(ciToti4)
ciToti7 <- fix_nameless_matrix(ciToti7)
D_E35 <- fix_nameless_matrix(D_E35)
D_E65 <- fix_nameless_matrix(D_E65)
D_E75 <- fix_nameless_matrix(D_E75)

## 改 Idents
Idents(ciToti4) <- "seurat_clusters"
Idents(ciToti7) <- "seurat_clusters"
Idents(ciToti8) <- "seurat_clusters"
Idents(ciToti10) <- "seurat_clusters"
Idents(D_E35) <- "seurat_clusters"
Idents(D_E65) <- "seurat_clusters"
Idents(D_E75) <- "seurat_clusters"
Idents(D_E85) <- "seurat_clusters"

ciToti4 <- RenameIdents(ciToti4, "0" = "Intermediate", "1" = "Epiblast", "2" = "PrE", "3" = "Intermediate", 
                        "4" = "TE", "5" = "Intermediate")
ciToti7 <- RenameIdents(ciToti7,  "0" = "Epiblast", "1" = "Epiblast", "2" = "PS", "3" = "Intermediate", "4" = "PS", "5" = "Em VE",
                        "6" = "Epiblast", "7" = "Epiblast", "8" = "Epiblast", "9" = "Epiblast", "10" = "Intermediate",
                        "11" = "PS", "12" = "ExE ecto", "13" = "Primitive Node Cells", "14" = "Epiblast", "15" = "ExE ecto",
                        "16" = "Primitive Node Cells", "17" = "ExE ecto", "18" = "PGCs", "19" = "Parietal endo", "20" = "ExE VE")
ciToti8 <- RenameIdents(ciToti8, "0" = "Aminion","1"="Neuromesodermal progenitors","2"="Epiblast","3"="Mesenchyme",
                        "4"="Visceral endo","5"="ExE endo","6"="Neuromesodermal progenitors","7"="ExE ecto",
                        "8"="Haematoendothelial progenitors","9"="PS","10"="ExE endo","11"="Intermediate","12"="ExE ecto",
                        "13"="Gut progenitors","14"="Aminion","15"="Parietal endo","16"="Intermediate")
ciToti10 <- RenameIdents(ciToti10, "0" = "Mixed","1"="Aminion","2"="Neuroectoderm","3"="Allantois",
                         "4"="PS","5"="Epithelium","6"="Somites","7"="Intermediate",
                         "8"="PGCs","9"="ExE endo","10"="Erythroblast","11"="Vascular endothelium","12"="ExE endo",
                         "13"="Surface ecto","14"="Primitive erythroid cells","15"="Parietal endo","16"="Cardiacmyocyte",
                         "17"="Intermediate","18"="Gut progenitors","19"="Somites","20"="Endothelium","21"="ExE ecto",
                         "22"="Cardiac meso","23"="Haematoendothelial progenitors","24"="Notochord")
D_E35 <- RenameIdents(D_E35, "0" = "Intermediate", "1" = "Epiblast", "2" = "PrE", "3" = "Intermediate", 
                      "4" = "TE", "5" = "Intermediate", "6" = "Epiblast")
D_E65 <- RenameIdents(D_E65,  "0" = "Epiblast", "1" = "Epiblast", "2" = "PS", "3" = "Intermediate", "4" = "PS", "5" = "Em VE",
                      "6" = "Epiblast", "7" = "Epiblast", "8" = "Epiblast", "9" = "Epiblast", "10" = "Intermediate",
                      "11" = "PS", "12" = "ExE ecto", "13" = "Primitive Node Cells", "14" = "Epiblast", "15" = "ExE ecto",
                      "16" = "Primitive Node Cells", "18" = "PGCs", "19" = "Parietal endo", "20" = "ExE VE")
D_E75 <- RenameIdents(D_E75, "0" = "Aminion","1"="Neuromesodermal progenitors","2"="Epiblast","3"="Mesenchyme",
                      "4"="Visceral endo","5"="ExE endo","6"="Neuromesodermal progenitors","7"="ExE ecto",
                      "8"="Haematoendothelial progenitors","9"="PS","10"="ExE endo","11"="Intermediate","12"="ExE ecto",
                      "13"="Gut progenitors","14"="Aminion")
D_E85 <- RenameIdents(D_E85, "0" = "Mixed","1"="Aminion","2"="Neuroectoderm","3"="Allantois",
                      "4"="PS","5"="Epithelium","6"="Somites","7"="Intermediate",
                      "8"="PGCs","9"="ExE endo","10"="Erythroblast","11"="Vascular endothelium","12"="ExE endo",
                      "13"="Surface ecto","14"="Primitive erythroid cells","15"="Parietal endo","16"="Cardiacmyocyte",
                      "17"="Intermediate","18"="Gut progenitors","19"="Somites","20"="Endothelium","21"="ExE ecto",
                      "22"="Cardiac meso","23"="Haematoendothelial progenitors","24"="Notochord")

## Idents 赋值给 celltype
ciToti4$celltype <- Idents(ciToti4)
head(ciToti4)
ciToti7$celltype <- Idents(ciToti7)
head(ciToti7)
ciToti8$celltype <- Idents(ciToti8)
head(ciToti8)
ciToti10$celltype <- Idents(ciToti10)
head(ciToti10)
D_E35$celltype <- Idents(D_E35)
head(D_E35)
D_E65$celltype <- Idents(D_E65)
head(D_E65)
D_E75$celltype <- Idents(D_E75)
head(D_E75)
D_E85$celltype <- Idents(D_E85)
head(D_E85)

## 添加一列 Dataset
ciToti1_4$dataset <- "ciToti1_4"
ciToti4$dataset <- "ciToti4"
ciToti7$dataset <- "ciToti7"
ciToti8$dataset <- "ciToti8"
ciToti10$dataset <- "ciToti10"
D_E35$dataset <- "D_E35"
D_E65$dataset <- "D_E65"
D_E75$dataset <- "D_E75"
D_E85$dataset <- "D_E85"

## 2) merge data ----
merged_obj <- merge(x = ciToti1_4,
                    y = c(ciToti10, ciToti4, ciToti7, ciToti8, D_E35, D_E65, D_E75, D_E85),
                    add.cell.ids = c("ciToti1_4", "ciToti10", "ciToti4", "ciToti7", "ciToti8", 
                                     "E35", "E65", "E75", "E85"))
head(merged_obj)
table(merged_obj$celltype)
table(merged_obj$dataset)

### process
merged_obj <- JoinLayers(merged_obj)

merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)

merged_obj <- RunUMAP(merged_obj, dims = 1:30)
merged_obj <- FindNeighbors(merged_obj, dims = 1:30)

merged_obj <- FindClusters(merged_obj, resolution = 0.5)

## 3) pseudotime analysis ----
cds <- as.cell_data_set(merged_obj)
colData(cds)$celltype <- merged_obj$celltype

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

get_root_by_celltype <- function(cds, celltype_col = "celltype", start_celltype = "Totipotent"){
  
  cell_ids <- which(colData(cds)[, celltype_col] == start_celltype)
  
  if(length(cell_ids) == 0) stop("Did not find any cells with the specified start_celltype!")
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  return(root_pr_nodes)
}

root_node <- get_root_by_celltype(cds, celltype_col = "celltype", start_celltype = "Totipotent")

cds <- order_cells(cds, root_pr_nodes = root_node)

df_check <- data.frame(
  Pseudotime = pseudotime(cds),
  CellType = colData(cds)$celltype
)

df_check <- df_check[is.finite(df_check$Pseudotime), ]

## 4) plot ----
### plot cells
pdf("4.Trajectory_clusters——plot_cells.pdf", width = 15, height = 15)
plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           trajectory_graph_color = "black",
           cell_size = 0.5) +
  facet_wrap(~celltype, ncol = 5) +
  theme(legend.position = "none") +
  ggtitle("Trajectory split by Cell Type")
dev.off()

pdf("4.Trajectory_clusters——plot_cells——Pseudotime.pdf", width = 15, height = 15)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           cell_size = 0.5) +
  facet_wrap(~celltype, ncol = 5) +
  scale_color_viridis_c(option = "mako", direction = -1) +
  ggtitle("Pseudotime split by Cell Type")
dev.off()

### plot density
plot_df <- data.frame(
  Pseudotime = pseudotime(cds),
  CellType = colData(cds)$celltype
) %>% filter(is.finite(Pseudotime))

pdf("4.Trajectory_clusters——density.pdf", width = 13, height = 6)
ggplot(plot_df, aes(x = Pseudotime, fill = CellType)) +
  geom_density(alpha = 0.6, size = 0.2)+
  theme_classic() +
  scale_fill_manual(values = Seurat::DiscretePalette(length(unique(plot_df$CellType)))) + 
  labs(title = "Global Pseudotime Density Distribution",
       x = "Pseudotime", y = "Density")
dev.off()

plot_df_clean <- plot_df %>% filter(is.finite(Pseudotime))
pdf("4.Trajectory_clusters——ridge.pdf", width = 13, height = 6)
ggplot(plot_df_clean, aes(x = Pseudotime, y = CellType, fill = CellType)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.7) + 
  theme_classic() +
  scale_fill_manual(values = Seurat::DiscretePalette(length(unique(plot_df_clean$CellType)))) +
  labs(title = "Pseudotime Distribution (Ridge Plot)", x = "Pseudotime", y = "") +
  theme(legend.position = "none")
dev.off()

### heatmap
library(ComplexHeatmap)
library(circlize)
library(Seurat)

gene_subset <- VariableFeatures(merged_obj)[1:2000] 
pr_graph_test <- graph_test(subset(cds, features = gene_subset), neighbor_graph = "principal_graph", cores = 4)

top_genes <- row.names(subset(pr_graph_test, q_value < 0.05))
top_genes <- head(top_genes[order(pr_graph_test[top_genes,]$q_value)], 1000) # 取前80个

pt <- pseudotime(cds)
valid_cells <- names(pt)[is.finite(pt)]
cell_order <- valid_cells[order(pt[valid_cells])]

available_genes <- intersect(top_genes, rownames(merged_obj))
if(length(available_genes) < length(top_genes)){
  message("Warning: ", length(top_genes) - length(available_genes), " genes were missing from the matrix.")
}

mat <- GetAssayData(merged_obj, layer = "data")[available_genes, cell_order]

mat <- as.matrix(mat)
mat <- t(apply(mat, 1, scale))
colnames(mat) <- cell_order

mat[mat > 5] <- 5
mat[mat < -2.5] <- -2.5

anno_df <- colData(cds)[cell_order, c("celltype"), drop=FALSE]
unique_types <- unique(anno_df$celltype)

cell_types_vec <- merged_obj$celltype[cell_order]
cell_types_vec <- factor(cell_types_vec)

cell_types_vec <- as.character(colData(cds)[cell_order, "celltype"])
cell_types_vec[is.na(cell_types_vec) | cell_types_vec == ""] <- "Unknown"

unique_types <- sort(unique(cell_types_vec))

my_colors <- rainbow(length(unique_types))
names(my_colors) <- unique_types

print(class(my_colors))
print(any(is.na(my_colors)))
print(head(my_colors))

pt_fixed <- as.numeric(as.vector(pt[cell_order]))

ha = HeatmapAnnotation(
  Pseudotime = pt_fixed,
  CellType = cell_types_vec,
  col = list(
    Pseudotime = colorRamp2(c(0, max(pt_fixed, na.rm=TRUE)), c("white", "blue")),
    CellType = ct_colors
  ),
  simple_anno_size = unit(0.4, "cm")
)

hmcols <- colorRampPalette(c(
  "#053061", "#2166AC", "#4393C3",
  "#e0f9b5", "#ffde7d",
  "#D6604D", "#B2182B", "#67001F"
))(8)

pdf("4.Trajectory_Heatmap_Complex.pdf", width = 10, height = 12)
Heatmap(mat, 
        name = "Z-score",
        
        col = colorRamp2(c(-2, -1, 0, 1, 2, 3, 4, 5), hmcols),
        
        top_annotation = ha,
        
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_names_gp = gpar(fontsize = 8),
        
        use_raster = TRUE 
)
dev.off()

# 5. 看一下ciToti1_4.RDS中intermediate细胞的位置 ----
cds_toti4 <- cds[, colData(cds)$dataset == "ciToti1_4"]

colData(cds_toti4)$highlight <- ifelse(
  colData(cds_toti4)$celltype == "Intermediate", 
  "Intermediate", 
  "Other"
)

pdf("5.Trajectory_clusters--ciToti1_4.pdf", width = 7, height = 6)
plot_cells(cds_toti4,
           color_cells_by = "highlight",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = TRUE,
           cell_size = 0.8) +
  scale_color_manual(values = c("Intermediate" = "red", "Other" = "lightgrey")) +
  ggtitle("Location of Intermediate cells in ciToti1_4") +
  theme(legend.position = "right")
dev.off()

