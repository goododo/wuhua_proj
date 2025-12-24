# Title: wuhua project analysis
# Author: Gaozy
# Time: 2025-12-19

# zygao02@comput172-general_env-R
# 0. Basic settings ----
if (dir.exists("/home/zygao02/wuhua_proj/251219/") == F) dir.create("/home/zygao02/wuhua_proj/251219/")
setwd("/home/zygao02/wuhua_proj/251219/")

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
  library("harmony")
})

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

###
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

# 1. Group 1 ----
group1 <- merge(x = ciToti4, 
                y = c(EPSC_S4, I_E35, Z_E45, D_E35, D_E45),
                add.cell.ids = c("ciToti4", "EPSC_S4", "I_E35", "Z_E45", "D_E35", "D_E45"))

group1 <- harmony_pipe(group1, reso = 0.04)
table(group1$seurat_clusters)

p1 <- DimPlot(group1, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group1, 
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

ggsave("1.UMAP_group1.pdf", plot = combined_plot, width = 12, height = 6)

# 2. Group 2 ----
group2 <- merge(x = ciToti7, 
                y = c(ETiX5, W_E65, Z_E65, D_E65),
                add.cell.ids = c("ciToti7", "ETiX5", "W_E65", "Z_E65", "D_E65"))

group2 <- harmony_pipe(group2, reso = 0.02)
table(group2$seurat_clusters)

group2 <- FindClusters(group2, resolution = 0.04)
table(group2$seurat_clusters)

p1 <- DimPlot(group2, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group2, 
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

ggsave("1.UMAP_group2.pdf", plot = combined_plot, width = 12, height = 6)


# 3. Group 3 ----
group3 <- merge(x = ciToti8, 
                y = c(EPSC_S7, ETiX6, iEFCEM_day6, TFSEM_day8, H_E75, W_E75, Z_E75, D_E75),
                add.cell.ids = c("ciToti8", "EPSC_S7", "ETiX6", "iEFCEM_day6", "TFSEM_day8", "H_E75", "W_E75", "Z_E75", "D_E75"))

group3 <- harmony_pipe(group3, reso = 0.1)
table(group3$seurat_clusters)

group3 <- FindClusters(group3, resolution = 0.2)
table(group3$seurat_clusters)

p1 <- DimPlot(group3, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group3, 
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

ggsave("1.UMAP_group3.pdf", plot = combined_plot, width = 12, height = 6)

# 4. Group 4 ----
group4 <- merge(x = ciToti10, 
                y = c(BPSCEM_day8, EPSC_S8, ETiX8, Hanna_2022_EM_day8, iEFCEM_day8,
                      TFSEM_day10, H_E85, W_E85, Z_E85, D_E85),
                add.cell.ids = c("ciToti10", "BPSCEM_day8", "EPSC_S8", "ETiX8", "Hanna_2022_EM_day8", "iEFCEM_day8",
                                 "TFSEM_day10", "H_E85", "W_E85", "Z_E85", "D_E85"))

group4 <- harmony_pipe(group4, reso = 0.8)
table(group4$seurat_clusters)

group4 <- FindClusters(group4, resolution = 0.83)
table(group4$seurat_clusters)

p1 <- DimPlot(group4, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group4, 
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
    title = "UMAP for Group 4",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_group4.pdf", plot = combined_plot, width = 12, height = 6)

# CCA ----
cca_pipe <- function(seurat_obj_list, reso){
  seurat_obj_list <- lapply(seurat_obj_list, function(seu_obj){
    seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seu_obj <- ScaleData(seu_obj, verbose = FALSE)
    seu_obj <- RunPCA(seu_obj, verbose = FALSE)
    
    return(seu_obj)
  })
  anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, reference = c(1, 2), 
                                    reduction = "rpca", dims = 1:50)
  all.CCA <- IntegrateData(anchorset = anchors, dims = 1:50)
  #all.CCA <- JoinLayers(all.CCA)
  all.CCA <- ScaleData(all.CCA, verbose = FALSE)
  all.CCA <- RunPCA(all.CCA, verbose = FALSE)
  all.CCA <- FindNeighbors(all.CCA, dims = 1:30)
  all.CCA <- RunUMAP(all.CCA, n.neighbors=10L, dims=1:30, min.dist=0.4)
  all.CCA <- FindClusters(all.CCA, resolution = reso)
  return(all.CCA)
}

# 1. Group 1 ----
group1_list <- list(ciToti4, EPSC_S4, I_E35, Z_E45, D_E35, D_E45)
group1_cca <- cca_pipe(group1_list, reso = 0.03)
table(group1_cca$seurat_clusters)
group1_cca <- FindClusters(group1_cca, resolution = 0.04)
table(group1_cca$seurat_clusters)

p1 <- DimPlot(group1_cca, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group1_cca, 
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

ggsave("1.UMAP_group1_cca.pdf", plot = combined_plot, width = 12, height = 6)


# 2. Group 2 ----
group2_list <- list(ciToti7, ETiX5, W_E65, Z_E65, D_E65)
group2_cca <- cca_pipe(group2_list, reso = 0.05)
table(group2_cca$seurat_clusters)

group2_cca <- FindClusters(group2_cca, resolution = 0.075)
table(group2_cca$seurat_clusters)

p1 <- DimPlot(group2_cca, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group2_cca, 
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

ggsave("1.UMAP_group2_cca.pdf", plot = combined_plot, width = 12, height = 6)

# 3. Group 3 ----
group3_list <- list(ciToti8, EPSC_S7, ETiX6, iEFCEM_day6, TFSEM_day8, H_E75, W_E75, Z_E75, D_E75)
group3_cca <- cca_pipe(group3_list, reso = 0.3)
table(group3_cca$seurat_clusters)

group3_cca <- FindClusters(group3_cca, resolution = 0.173)
table(group3_cca$seurat_clusters)

p1 <- DimPlot(group3_cca, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group3_cca, 
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

ggsave("1.UMAP_group3_cca.pdf", plot = combined_plot, width = 12, height = 6)

# 4. Group 4 ----
group4_list <- list(ciToti10, BPSCEM_day8, EPSC_S8, ETiX8, Hanna_2022_EM_day8, iEFCEM_day8,
                    TFSEM_day10, H_E85, W_E85, Z_E85, D_E85)
group4_cca <- cca_pipe(group4_list, reso = 0.9)
table(group4_cca$seurat_clusters)

group4_cca <- FindClusters(group4_cca, resolution = 0.85)
table(group4_cca$seurat_clusters)

p1 <- DimPlot(group4_cca, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.5, 
              alpha = 0.6,
              label = TRUE,
              repel = TRUE,
              raster = FALSE) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(group4_cca, 
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
    title = "UMAP for Group 4",
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave("1.UMAP_group4_cca.pdf", plot = combined_plot, width = 12, height = 6)







