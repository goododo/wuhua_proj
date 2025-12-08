# Title: wuhua project analysis
# Author: Gaozy
# Time: 2025-12-06

# zygao02@comput172 - general_env - R
if (dir.exists("/home/zygao02/wuhua_proj/251206/") == F) dir.create("/home/zygao02/wuhua_proj/251206/")
setwd("/home/zygao02/wuhua_proj/251206/")

# 0. Downgrade Seurat v5 ----
library(Seurat)

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

export_for_v4 <- function(obj, assay = "RNA", file) {
  obj <- JoinLayers(obj)
  a <- obj[[assay]]
  
  if ("Assay5" %in% class(a)) {
    counts <- a@layers$counts
  } else {
    counts <- GetAssayData(obj, assay = assay, slot = "counts")
  }
  
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  
  if (is.null(rownames(counts)) || any(rownames(counts) == "")) {
    rownames(counts) <- rownames(obj)
  }
  if (is.null(colnames(counts)) || any(colnames(counts) == "")) {
    colnames(counts) <- rownames(obj@meta.data)
  }
  
  export_list <- list(
    counts    = counts,
    meta.data = obj@meta.data,
    assay     = assay
  )
  
  saveRDS(export_list, file = file)
  message("已转换：", file)
  
  return(export_list)
}

## 转换
ciToti1_4_v4  <- export_for_v4(ciToti1_4, assay = "RNA", file = "ciToti1_4_v4.rds")
ciToti10_v4   <- export_for_v4(ciToti10, assay = "RNA", file = "ciToti10_v4.rds")
ciToti4_v4 <- export_for_v4(ciToti4, assay = "RNA", file = "ciToti4_v4.rds")
ciToti7_v4 <- export_for_v4(ciToti7, assay = "RNA", file = "ciToti7_v4.rds")
ciToti8_v4 <- export_for_v4(ciToti8, assay = "RNA", file = "ciToti8_v4.rds")
D_E35_v4 <- export_for_v4(D_E35, assay = "RNA", file = "D_E35_v4.rds")
D_E65_v4 <- export_for_v4(D_E65, assay = "RNA", file = "D_E65_v4.rds")
D_E75_v4 <- export_for_v4(D_E75, assay = "RNA", file = "D_E75_v4.rds")
D_E85_v4 <- export_for_v4(D_E85, assay = "RNA", file = "D_E85_v4.rds")


# zygao02@comput172 - conda activate /home/lushi02/anaconda3/envs/Han - R
# 0. Basic settings ----

.libPaths(c("/usr/lib64/R/library",
            "/usr/share/R/library",
            "/home/lushi02/.local/share/r-miniconda/envs/r-reticulate",
            "/home/lushi02/anaconda3/lib/R/library",
            "/home/lushi02/anaconda3/envs/Han/lib/R/library/ggplot2",
            "/home/lushi02/R/x86_64-redhat-linux-gnu-library/4.2"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ggsci)
  library(RColorBrewer)
  library(monocle)
  library(Matrix)
})

# 1. Trajectory analysis ----

## 1) load data ----
export <- readRDS("/home/zygao02/wuhua_proj/251206/ciToti1_4_v4.rds")
str(export) 

ciToti1_4 <- CreateSeuratObject(
  counts = export$counts,
  meta.data = export$meta.data,
  assay = export$assay
)
ciToti1_4

export <- readRDS("/home/zygao02/wuhua_proj/251206/ciToti10_v4.rds")
str(export) 

ciToti10 <- CreateSeuratObject(
  counts = export$counts,
  meta.data = export$meta.data,
  assay = export$assay
)
ciToti10

export <- readRDS("/home/zygao02/wuhua_proj/251206/ciToti4_v4.rds")
str(export) 

ciToti4 <- CreateSeuratObject(
  counts = export$counts,
  meta.data = export$meta.data,
  assay = export$assay
)
ciToti4

export <- readRDS("/home/zygao02/wuhua_proj/251206/ciToti7_v4.rds")
str(export) 

ciToti7 <- CreateSeuratObject(
  counts = export$counts,
  meta.data = export$meta.data,
  assay = export$assay
)
ciToti7

export <- readRDS("/home/zygao02/wuhua_proj/251206/ciToti8_v4.rds")
str(export) 

ciToti8 <- CreateSeuratObject(
  counts = export$counts,
  meta.data = export$meta.data,
  assay = export$assay
)
ciToti8

export <- readRDS("/home/zygao02/wuhua_proj/251206/D_E35_v4.rds")
str(export) 

D_E35 <- CreateSeuratObject(
  counts = export$counts,
  meta.data = export$meta.data,
  assay = export$assay
)
D_E35

export <- readRDS("/home/zygao02/wuhua_proj/251206/D_E65_v4.rds")
str(export) 

D_E65 <- CreateSeuratObject(
  counts = export$counts,
  meta.data = export$meta.data,
  assay = export$assay
)
D_E65

export <- readRDS("/home/zygao02/wuhua_proj/251206/D_E75_v4.rds")
str(export) 

D_E75 <- CreateSeuratObject(
  counts = export$counts,
  meta.data = export$meta.data,
  assay = export$assay
)
D_E75

export <- readRDS("/home/zygao02/wuhua_proj/251206/D_E85_v4.rds")
str(export) 

D_E85 <- CreateSeuratObject(
  counts = export$counts,
  meta.data = export$meta.data,
  assay = export$assay
)
D_E85

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

meta <- merged_obj@meta.data
exp <- merged_obj@assays$RNA@data
## downsample
set.seed(5)
cells <- rownames(meta)
downsample_cells <- as.vector(sample_n(as.data.frame(cells), 8000, replace = FALSE)$cells)

exp_data_down <- exp[,c(downsample_cells)]
meta_data_down <- meta[c(downsample_cells),]

genes <- as.data.frame(rownames(exp))
colnames(genes)<-"gene_short_name"
rownames(genes) <- genes$gene_short_name

pd <- new("AnnotatedDataFrame", data = meta_data_down)	
fd <- new("AnnotatedDataFrame", data = genes)

HSMM <- newCellDataSet(exp_data_down,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

patch_project2MST <- function() {
  ns <- asNamespace("monocle")
  f  <- get("project2MST", envir = ns)
  
  txt <- deparse(body(f))
  
  pattern <- 'if (class(projection) != "matrix") projection <- as.matrix(projection)'
  if (any(grepl("class(projection)", txt, fixed = TRUE))) {
    txt <- gsub(pattern,
                'projection <- as.matrix(projection)',
                txt, fixed = TRUE)
    body(f) <- parse(text = paste(txt, collapse = "\n"))[[1]]
    
    assignInNamespace("project2MST", f, ns = "monocle")
    message("已修补 monocle::project2MST()。")
  } else {
    message("未找到需要替换的部分，已修复。")
  }
}

patch_project2MST()

##使用clusters差异表达基因
clustering_DEG_genes <- differentialGeneTest(HSMM, fullModelFormulaStr = '~celltype', cores = 1)

HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:800]
mycds <- setOrderingFilter(HSMM, ordering_genes = HSMM_ordering_genes)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)

p1 <- plot_ordering_genes(mycds)

##使用seurat选择的高变基因
var.genes <- VariableFeatures(DC)
mycds <- setOrderingFilter(mycds, var.genes)
p2 <- plot_ordering_genes(mycds)
##使用monocle选择的高变基因
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(HSMM, disp.genes)
p3 <- plot_ordering_genes(mycds)
##结果对比
p1|p2|p3













