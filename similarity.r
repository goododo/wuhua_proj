# Title: wuhua project analysis
# Author: Gaozy
# Time: 2025-11-25

# zygao02@comput172-general_env-R
# 0. Basic settings ----
if (dir.exists("/home/zygao02/wuhua_proj/251210/") == F) dir.create("/home/zygao02/wuhua_proj/251210/")
setwd("/home/zygao02/wuhua_proj/251210/")

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
})

source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

# 1. Combined Ref v.s. Celltype devided Query ----
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

## 统一 celltype
merged_ref$cell_type <- gsub("_", " ", merged_ref$cell_type)

library(stringr)
merged_ref$cell_type <- str_to_title(merged_ref$cell_type)

table(merged_ref$cell_type)

library(dplyr)
rename_map <- c(
  "Exe Ectoderm" = "Extraembryonic Ectoderm",
  "Exe Endoderm" = "Extraembryonic Endoderm", 
  "Exe Mesoderm" = "Extraembryonic Mesoderm",
  "Exe Visceral Endoderm" = "Extraembryonic Visceral Endoderm",
  "Rostral Neurectoderm" = "Rostral Neuroectoderm",
  "Pgc" = "Primordial Germ Cells",
  "Epc Progenitors" = "EPC Progenitors", # Ectoplacental Cone
  "Tsc" = "TSC", # Trophoblast Stem Cells
  "P Tgc" = "P-TGC", # Parietal TGC
  "Spa Tgc" = "SpA-TGC", # Spiral Artery TGC
  "Spt Gly" = "SpT-Gly" # Spongiotrophoblast Glycogen
)

merged_ref$cell_type_unified <- ifelse(
  merged_ref$cell_type %in% names(rename_map), 
  rename_map[merged_ref$cell_type], 
  merged_ref$cell_type
)

merged_ref$cell_type <- merged_ref$cell_type_unified

table(merged_ref$cell_type)

Idents(merged_ref) <- merged_ref$cell_type

length(unique(Idents(merged_ref)))

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
ciToti4$cell_type <- Idents(ciToti4)
head(ciToti4)
ciToti7$cell_type <- Idents(ciToti7)
head(ciToti7)
ciToti8$cell_type <- Idents(ciToti8)
head(ciToti8)
ciToti10$cell_type <- Idents(ciToti10)
head(ciToti10)
D_E35$cell_type <- Idents(D_E35)
head(D_E35)
D_E65$cell_type <- Idents(D_E65)
head(D_E65)
D_E75$cell_type <- Idents(D_E75)
head(D_E75)
D_E85$cell_type <- Idents(D_E85)
head(D_E85)

## 添加一列 Dataset
ciToti4$dataset <- "ciToti4"
ciToti7$dataset <- "ciToti7"
ciToti8$dataset <- "ciToti8"
ciToti10$dataset <- "ciToti10"
D_E35$dataset <- "D_E35"
D_E65$dataset <- "D_E65"
D_E75$dataset <- "D_E75"
D_E85$dataset <- "D_E85"

## unify cell type
unify_map_extended <- c(
  "Aminion" = "Amnion",
  "Cardiacmyocyte" = "Cardiomyocytes",
  
  "Exe Ecto" = "Extraembryonic Ectoderm",
  "Exe Endo" = "Extraembryonic Endoderm",
  "Exe Ve" = "Extraembryonic Visceral Endoderm",
  "Em Ve" = "Embryonic Visceral Endoderm",
  "Visceral Endo" = "Visceral Endoderm",
  "Parietal Endo" = "Parietal Endoderm",
  "Surface Ecto" = "Surface Ectoderm",
  "Cardiac Meso" = "Cardiopharyngeal Mesoderm",
  
  "Ps" = "Primitive Streak",
  "Pgcs" = "Primordial Germ Cells",
  "Te" = "TE",
  "Pre" = "PrE"
)

unify_cell_names <- function(seurat_obj) {
  current_types <- seurat_obj$cell_type
  
  formatted_types <- str_to_title(current_types)
  
  final_types <- ifelse(
    formatted_types %in% names(unify_map_extended),
    unify_map_extended[formatted_types],
    formatted_types
  )
  seurat_obj$cell_type <- final_types
  
  Idents(seurat_obj) <- final_types
  
  return(seurat_obj)
}

ciToti4 <- unify_cell_names(ciToti4)
ciToti7 <- unify_cell_names(ciToti7)
ciToti8 <- unify_cell_names(ciToti8)
ciToti10 <- unify_cell_names(ciToti10)
D_E35 <- unify_cell_names(D_E35)
D_E65 <- unify_cell_names(D_E65)
D_E75 <- unify_cell_names(D_E75)
D_E85 <- unify_cell_names(D_E85)

table(ciToti4$cell_type)
table(ciToti7$cell_type)
table(ciToti8$cell_type)
table(ciToti10$cell_type)
table(D_E35$cell_type)
table(D_E65$cell_type)
table(D_E75$cell_type)
table(D_E85$cell_type)

## merge quey
merged_query <- merge(x = ciToti10, 
                      y = c(ciToti4, ciToti7, ciToti8, D_E35, D_E65, D_E75, D_E85),
                      add.cell.ids = c("ciToti10", "ciToti4", "ciToti7", "ciToti8", 
                                       "E35", "E65", "E75", "E85"))
merged_query <- JoinLayers(merged_query)

## 3) Similarity ----
RhpcBLASctl::blas_set_num_threads(10)

### a. train data ----
train <- merged_ref
Idents(train) <- "cell_type"
train <- FindVariableFeatures(train, selection.method = 'vst', nfeatures = 2000)

### b. test data ----
test <- merged_query
Idents(test) <- "cell_type"
test <- JoinLayers(test)
test <- NormalizeData(test)
test <- FindVariableFeatures(test, selection.method = 'vst', nfeatures = 2000)

common_genes <- intersect(rownames(train), rownames(test))
print(paste("Common genes found:", length(common_genes)))

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

train_sub <- subset_by_group(train, group.by = "cell_type", n = 200)
test_sub <- subset_by_group(test,  group.by = "cell_type", n = 200)

print(table(train_sub$cell_type))
print(table(test_sub$cell_type))

train_sparse <- GetAssayData(train_sub, assay = "RNA", layer = "data")[common_genes, ]
test_sparse <- GetAssayData(test_sub,  assay = "RNA", layer = "data")[common_genes, ]

train_mat <- as.matrix(train_sparse)
test_mat <- as.matrix(test_sparse)

rm(train_sparse, test_sparse)
gc()

train_group <- as.character(train_sub$cell_type)
test_group  <- as.character(test_sub$cell_type)

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
  sample.cells = 0,
  test.data = test_mat,
  test.group = test_group,
  alpha = 0.5, 
  nfolds = 5,
  seed = 123
)

heatmap_mat = Test_similarity$heatmap@matrix

pdf("1.Similarity_cell_type.pdf", width=10, height=18, onefile = F)
pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                   color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                   cluster_rows = F)
dev.off()

# 2. Combined query data ----
comb <- readRDS("/home/lushi02/project/wuhua/Combined_query_data.RDS")

this_work <- subset(comb, group %in% c("ciToti-EM", "This work"))
pub1 <- subset(comb, group %in% c("BPSC-EM"))
pub2 <- subset(comb, group %in% c("EPSC-EM"))
pub3 <- subset(comb, group %in% c("Zernicka-Goetz_2019"))
pub4 <- subset(comb, group %in% c("Zernicka-Goetz_2022", "ETiX embryoid"))
pub5 <- subset(comb, group %in% c("H_2022_EM", "Hanna_2022"))
pub6 <- subset(comb, group %in% c("Hanna_2025", "TF-SEM"))
pub7 <- subset(comb, group %in% c("Weissman_2019"))
pub8 <- subset(comb, group %in% c("Izpisua-Belmonte_2019"))
pub9 <- subset(comb, group %in% c("iEFC-EM"))

## Similarity ----
RhpcBLASctl::blas_set_num_threads(10)

### a. train data ----
train <- merged_ref
Idents(train) <- "cell_type"
train <- FindVariableFeatures(train, selection.method = 'vst', nfeatures = 2000)

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

train_sub <- subset_by_group(train, group.by = "cell_type", n = 200)

### 1) this_work ----
sub_obj <- this_work
DefaultAssay(sub_obj) <- "RNA"
sub_obj <- NormalizeData(sub_obj, verbose = FALSE)
sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sub_obj <- ScaleData(sub_obj, features = VariableFeatures(sub_obj), verbose = FALSE)
sub_obj <- RunPCA(sub_obj, npcs = 30, verbose = FALSE)

sub_obj <- FindNeighbors(sub_obj, dims = 1:20, verbose = FALSE)

sub_obj <- FindClusters(sub_obj, resolution = 0.6, verbose = FALSE)
table(sub_obj$seurat_clusters)

#### b. test data ----
test <- sub_obj
Idents(test) <- "seurat_clusters"

common_genes <- intersect(rownames(train), rownames(test))
print(paste("Common genes found:", length(common_genes)))

train_sparse <- GetAssayData(train_sub, assay = "RNA", layer = "data")[common_genes, ]
test_sparse <- GetAssayData(test,  assay = "RNA", layer = "data")[common_genes, ]

train_mat <- as.matrix(train_sparse)
test_mat <- as.matrix(test_sparse)

rm(train_sparse, test_sparse)
gc()

train_group <- as.character(train_sub$cell_type)
test_group  <- as.character(test$seurat_clusters)

table(is.na(train_group))

if (ncol(train_mat) != length(train_group)) stop("Train 维度不匹配！")
if (ncol(test_mat)  != length(test_group))  stop("Test 维度不匹配！")

#### c. calculate ----
source('/home/fengyan02/Project/YYYProject202306/BasicAnalysis/202403/20240304/script/glm.predict.R')

Test_similarity <- glm.predict(
  train.data = train_mat,     
  train.group = train_group,
  genes.used = common_genes,
  downsample = TRUE,
  sample.cells = 0,
  test.data = test_mat,
  test.group = test_group,
  alpha = 0.5, 
  nfolds = 5,
  seed = 123
)

heatmap_mat = t(Test_similarity2$heatmap@matrix)

pdf("1.Similarity_cell_type.pdf", width=18, height=5, onefile = F)
pheatmap::pheatmap(heatmap_mat, cluster_cols = F,
                   color = colorRampPalette(rev(brewer.pal(9, "Spectral")))(100), 
                   cluster_rows = F)
dev.off()








