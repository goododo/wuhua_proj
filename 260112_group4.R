# Title: wuhua project analysis
# Author: Gaozy
# Time: 2026-01-12

# zygao02@comput172-general_env-R
# 0. Basic settings ----
if (dir.exists("/home/zygao02/wuhua_proj/260112/") == F) dir.create("/home/zygao02/wuhua_proj/260112/")
setwd("/home/zygao02/wuhua_proj/260112/")

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

## QC pipeline ----
qc_pipeline <- function(seurat_obj, group_name){
  DefaultAssay(seurat_obj) <- "RNA"
  #seurat_obj <- JoinLayers(seurat_obj)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  #seurat_obj <- UpdateSeuratObject(seurat_obj)
  
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

## UMAP pipeline ----
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
                group.by = "cell_type",
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
  
}

## load combined datasets & seperate ----
#comb <- readRDS("/home/lushi02/project/wuhua/Combined_query_data.RDS")
group4 <- readRDS("/home/lushi02/project/wuhua/group4_final.rds")

ciToti10 <- subset(group4, orig.ident %in% c("ciToti10"))
BPSCEM_day8 <- subset(group4, orig.ident %in% c("BPSCEM_day8"))
EPSC_S8 <- subset(group4, orig.ident %in% c("EPSC_S8"))
ETiX8 <- subset(group4, orig.ident %in% c("ETiX8"))
Hanna_2022_EM_day8 <- subset(group4, orig.ident %in% c("Hanna_2022_EM_day8"))
iEFCEM_day8 <- subset(group4, orig.ident %in% c("iEFCEM_day8"))
TFSEM_day10 <- subset(group4, orig.ident %in% c("TFSEM_day10"))
D_E85 <- subset(group4, orig.ident %in% c("D_E85"))

# 1. Merge group4 data ----
group4_merged <- merge(x = ciToti10, 
                       y = c(BPSCEM_day8, EPSC_S8, ETiX8, Hanna_2022_EM_day8, iEFCEM_day8, TFSEM_day10, D_E85),
                       add.cell.ids = c("ciToti10", "BPSCEM_day8", "EPSC_S8", "ETiX8", "Hanna_2022_EM_day8", "iEFCEM_day8", "TFSEM_day10", "D_E85"))

group4_merged_qc <- qc_pipeline(group4_merged, "group4_raw_merge")
group4_merged_qc[["percent.mt"]] <- PercentageFeatureSet(group4_merged_qc, pattern = "^mt-")

group4_merged_qc <- CellCycleScoring(group4_merged_qc, 
                                     s.features = str_to_title(cc.genes$s.genes), 
                                     g2m.features = str_to_title(cc.genes$g2m.genes), 
                                     set.ident = TRUE)

group4_merged_qc <- NormalizeData(group4_merged_qc) |> 
  FindVariableFeatures(nfeatures = 2000) |> 
  ScaleData(vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")) |> 
  RunPCA(npcs = 50)

# 2. Harmony ----
harmony_pipe <- function(seurat_obj, reso){
  #seurat_obj <- JoinLayers(seurat_obj)
  #seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
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

group4 <- harmony_pipe(group4_merged_qc, reso = 1.5)

umap_pipeline(seurat_obj = group4,
              group_label = "group4", 
              plot_title_group = "4",
              plot_width = 17, plot_height = 8)

p <- DimPlot(group4, reduction = "umap", split.by = "orig.ident",
             group.by = "cell_type",
             pt.size = 0.5, 
             alpha = 0.6,
             label = TRUE,
             repel = TRUE,
             raster = FALSE,
             ncol = 4) +
  NoAxes() +
  theme(plot.title = element_text(hjust = 0.5))

combined_plot_split <- (p) + 
  plot_annotation(
    title = paste0("Split UMAP for Group 4"),
    theme = theme(
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    ) )

ggsave(paste0("1.UMAP_group4_harmony_split.pdf"), plot = combined_plot_split, width = 21, height = 14)

saveRDS(group4, "group4.rds")

#library(devtools)
#install_github('theislab/kBET')
#devtools::install_github("immunogenomics/lisi")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(kBET)
  library(lisi)
  library(cluster)   # silhouette
})
options(stringsAsFactors = FALSE)

# ============================================================
# 0) Inputs (group4)
# ============================================================
obj <- group4
reduction_use <- "harmony"
ndims <- 20

stopifnot(reduction_use %in% Reductions(obj))
stopifnot(all(c("cell_type","orig.ident") %in% colnames(obj@meta.data)))

obj$cell_type  <- droplevels(factor(obj$cell_type))
obj$orig.ident <- droplevels(factor(obj$orig.ident))
Idents(obj) <- "cell_type"

pairs <- list(
  c("ciToti10",       "D_E85"),
  c("BPSCEM_day8",       "D_E85"),
  c("EPSC_S8",         "D_E85"),
  c("ETiX8",   "D_E85"),
  c("Hanna_2022_EM_day8",    "D_E85"),
  c("iEFCEM_day8",         "D_E85"),
  c("TFSEM_day10",   "D_E85")
)
pairs <- lapply(pairs, function(p) trimws(as.character(p)))

baseline <- "ciToti10"

# bootstrap settings
n_rep <- 200
k0 <- 30
perplexity <- 30

# 各指标建议 cap（每边抽多少）
n_cap_list <- list(
  kbet  = 500,
  ilisi = 500,
  asw   = 200
)

# relax settings（只对 strict 门槛卡住的项补跑）
min_k0 <- 5
min_perplexity <- 2   # 想更保守可改 5

# 忽略的 cell_type（你说 celltype=0 的不用管）
ignore_celltypes <- c("0")

# output
outdir <- "group4_metrics_3metrics"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# pdf device（没有 cairo 就回退到普通 pdf）
pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

# stage parser (D_E75 -> E7.5)
parse_stage_from_ref <- function(ref) {
  m <- stringr::str_match(ref, "^D_E(\\d+)$")
  if (is.na(m[1,2])) return(NA_character_)
  x <- m[1,2]
  if (nchar(x) >= 2) {
    paste0("E", substr(x, 1, nchar(x)-1), ".", substr(x, nchar(x), nchar(x)))
  } else {
    paste0("E", x)
  }
}

# ============================================================
# 1) Core utilities
# ============================================================
sym_downsample_cells <- function(seu, group_a, group_b, n_each, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cells_a <- WhichCells(seu, expression = orig.ident == group_a)
  cells_b <- WhichCells(seu, expression = orig.ident == group_b)
  if (length(cells_a) < n_each || length(cells_b) < n_each) return(NULL)
  take_a <- sample(cells_a, n_each)
  take_b <- sample(cells_b, n_each)
  subset(seu, cells = c(take_a, take_b))
}

calc_kbet_once <- function(seu_sub, reduction = "harmony", ndims = 20, k0 = 50) {
  emb <- Embeddings(seu_sub, reduction)[, 1:ndims, drop = FALSE]
  batch <- factor(seu_sub$orig.ident)
  keep <- !is.na(batch)
  emb <- emb[keep, , drop = FALSE]
  batch <- droplevels(batch[keep])
  if (nrow(emb) <= k0) return(NA_real_)
  res <- kBET(df = emb, batch = batch, do.pca = FALSE, k0 = k0)
  as.numeric(res$summary$kBET.observed[1])
}

calc_ilisi_once <- function(seu_sub, reduction = "harmony", ndims = 20, perplexity = 30) {
  emb <- Embeddings(seu_sub, reduction)[, 1:ndims, drop = FALSE]
  meta <- data.frame(batch = factor(seu_sub$orig.ident), row.names = colnames(seu_sub))
  keep <- !is.na(meta$batch)
  emb <- emb[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  if (nrow(emb) < 3 * perplexity) return(NA_real_)
  lisi_res <- compute_lisi(X = emb, meta_data = meta, label_colnames = "batch", perplexity = perplexity)
  mean(lisi_res[, "batch"])
}

calc_batch_asw_once <- function(seu_sub, reduction = "harmony", ndims = 20) {
  emb <- Embeddings(seu_sub, reduction)[, 1:ndims, drop = FALSE]
  batch <- factor(seu_sub$orig.ident)
  keep <- !is.na(batch)
  emb <- emb[keep, , drop = FALSE]
  batch <- droplevels(batch[keep])
  if (nrow(emb) < 10 || nlevels(batch) < 2) return(NA_real_)
  d <- dist(emb)
  sil <- cluster::silhouette(as.integer(batch), d)
  mean_sil <- mean(sil[, "sil_width"])  # [-1,1]
  (1 - mean_sil) / 2                    # [0,1], 1=best mixing
}

# ============================================================
# 2) Per-ct bootstrap runner (strict + relax)
#    returns: list(out=..., skip=...)
# ============================================================
run_metric_pair_ct_boot <- function(
    seu,
    celltype,
    group_a,
    group_b,
    metric = c("kbet", "ilisi", "asw"),
    reduction = "harmony",
    ndims = 20,
    k0 = 30,
    perplexity = 30,
    n_rep = 200,
    n_cap = 500,
    seed = 1,
    mode = c("strict","relax"),
    min_k0 = 5,
    min_perplexity = 2
) {
  metric <- match.arg(metric)
  mode   <- match.arg(mode)
  
  # 先过滤目标 celltype + 两个样本
  cells_keep <- WhichCells(
    seu,
    expression = !is.na(cell_type) & cell_type == celltype & orig.ident %in% c(group_a, group_b)
  )
  if (length(cells_keep) == 0) {
    return(list(out = NULL, skip = data.frame(
      cell_type=celltype, group_a=group_a, group_b=group_b, pair=paste(group_a,"vs",group_b),
      metric=metric, mode=mode,
      reason="no cells after filtering",
      n_a=NA_integer_, n_b=NA_integer_, n_each=NA_integer_,
      k0_used=NA_integer_, perplexity_used=NA_integer_
    )))
  }
  
  seu_ct <- subset(seu, cells = cells_keep)
  tab <- table(seu_ct$orig.ident)
  n_a <- as.integer(ifelse(group_a %in% names(tab), tab[[group_a]], 0))
  n_b <- as.integer(ifelse(group_b %in% names(tab), tab[[group_b]], 0))
  n_each <- min(n_a, n_b, n_cap)
  
  if (!is.finite(n_each) || n_each <= 1) {
    return(list(out=NULL, skip=data.frame(
      cell_type=celltype, group_a=group_a, group_b=group_b, pair=paste(group_a,"vs",group_b),
      metric=metric, mode=mode,
      reason="n_each <= 1",
      n_a=n_a, n_b=n_b, n_each=n_each,
      k0_used=NA_integer_, perplexity_used=NA_integer_
    )))
  }
  
  # --- strict/relax 参数选择 ---
  k0_used   <- k0
  perp_used <- perplexity
  
  if (mode == "relax") {
    if (metric == "kbet") {
      # 总 n=2*n_each；经验上 k0 <= (n-1)/2 = n_each-1
      k0_used <- min(k0, n_each - 1)
      k0_used <- max(min_k0, k0_used)
    }
    if (metric == "ilisi") {
      perp_max <- floor((2 * n_each) / 3)
      perp_used <- min(perplexity, perp_max)
      perp_used <- max(min_perplexity, perp_used)
    }
  }
  
  # --- thresholds ---
  if (mode == "strict") {
    if (metric == "kbet" && n_each < (k0 + 1)) {
      return(list(out=NULL, skip=data.frame(
        cell_type=celltype, group_a=group_a, group_b=group_b, pair=paste(group_a,"vs",group_b),
        metric=metric, mode=mode,
        reason=paste0("n_each < k0+1 (", n_each, "<", k0+1, ")"),
        n_a=n_a, n_b=n_b, n_each=n_each,
        k0_used=k0, perplexity_used=NA_integer_
      )))
    }
    if (metric == "ilisi" && (2 * n_each) < (3 * perplexity)) {
      return(list(out=NULL, skip=data.frame(
        cell_type=celltype, group_a=group_a, group_b=group_b, pair=paste(group_a,"vs",group_b),
        metric=metric, mode=mode,
        reason=paste0("2*n_each < 3*perplexity (", 2*n_each, "<", 3*perplexity, ")"),
        n_a=n_a, n_b=n_b, n_each=n_each,
        k0_used=NA_integer_, perplexity_used=perplexity
      )))
    }
    if (metric == "asw" && (2 * n_each) < 10) {
      return(list(out=NULL, skip=data.frame(
        cell_type=celltype, group_a=group_a, group_b=group_b, pair=paste(group_a,"vs",group_b),
        metric=metric, mode=mode,
        reason="too few cells for ASW (2*n_each<10)",
        n_a=n_a, n_b=n_b, n_each=n_each,
        k0_used=NA_integer_, perplexity_used=NA_integer_
      )))
    }
  } else {
    # relax：用调整后的参数再检查
    if (metric == "kbet" && (2 * n_each) <= k0_used) {
      return(list(out=NULL, skip=data.frame(
        cell_type=celltype, group_a=group_a, group_b=group_b, pair=paste(group_a,"vs",group_b),
        metric=metric, mode=mode,
        reason=paste0("relax still fails: 2*n_each <= k0_used (", 2*n_each, "<=", k0_used, ")"),
        n_a=n_a, n_b=n_b, n_each=n_each,
        k0_used=k0_used, perplexity_used=NA_integer_
      )))
    }
    if (metric == "ilisi" && (2 * n_each) < (3 * perp_used)) {
      return(list(out=NULL, skip=data.frame(
        cell_type=celltype, group_a=group_a, group_b=group_b, pair=paste(group_a,"vs",group_b),
        metric=metric, mode=mode,
        reason=paste0("relax still fails: 2*n_each < 3*perplexity_used (", 2*n_each, "<", 3*perp_used, ")"),
        n_a=n_a, n_b=n_b, n_each=n_each,
        k0_used=NA_integer_, perplexity_used=perp_used
      )))
    }
    if (metric == "asw" && (2 * n_each) < 10) {
      return(list(out=NULL, skip=data.frame(
        cell_type=celltype, group_a=group_a, group_b=group_b, pair=paste(group_a,"vs",group_b),
        metric=metric, mode=mode,
        reason="too few cells for ASW (2*n_each<10)",
        n_a=n_a, n_b=n_b, n_each=n_each,
        k0_used=NA_integer_, perplexity_used=NA_integer_
      )))
    }
  }
  
  # --- bootstrap ---
  vals <- numeric(0)
  for (r in seq_len(n_rep)) {
    seu_sub <- sym_downsample_cells(seu_ct, group_a, group_b, n_each, seed = seed + r)
    if (is.null(seu_sub)) next
    
    v <- switch(metric,
                kbet  = calc_kbet_once(seu_sub, reduction, ndims, k0_used),
                ilisi = calc_ilisi_once(seu_sub, reduction, ndims, perp_used),
                asw   = calc_batch_asw_once(seu_sub, reduction, ndims)
    )
    if (is.finite(v)) vals <- c(vals, v)
  }
  
  if (length(vals) < max(10, floor(n_rep * 0.3))) {
    return(list(out=NULL, skip=data.frame(
      cell_type=celltype, group_a=group_a, group_b=group_b, pair=paste(group_a,"vs",group_b),
      metric=metric, mode=mode,
      reason=paste0("too few reps used: ", length(vals), "/", n_rep),
      n_a=n_a, n_b=n_b, n_each=n_each,
      k0_used=ifelse(metric=="kbet", k0_used, NA_integer_),
      perplexity_used=ifelse(metric=="ilisi", perp_used, NA_integer_)
    )))
  }
  
  ci <- quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  
  out <- list(
    per_rep = data.frame(
      cell_type = celltype,
      group_a = group_a,
      group_b = group_b,
      pair = paste(group_a, "vs", group_b),
      metric = metric,
      mode = mode,
      rep = seq_along(vals),
      value = vals,
      n_each = n_each,
      n_a = n_a, n_b = n_b,
      k0_used = ifelse(metric == "kbet", k0_used, NA_integer_),
      perplexity_used = ifelse(metric == "ilisi", perp_used, NA_integer_)
    ),
    summary = data.frame(
      cell_type = celltype,
      group_a = group_a,
      group_b = group_b,
      pair = paste(group_a, "vs", group_b),
      metric = metric,
      mode = mode,
      mean = mean(vals),
      ci_low = ci[1],
      ci_high = ci[2],
      reps_used = length(vals),
      n_each = n_each,
      n_a = n_a, n_b = n_b,
      k0_used = ifelse(metric == "kbet", k0_used, NA_integer_),
      perplexity_used = ifelse(metric == "ilisi", perp_used, NA_integer_)
    )
  )
  
  list(out = out, skip = NULL)
}

# ============================================================
# 3.1) Global raw (你已有 run_metric_pair_global_fixed；也可继续用它)
# ============================================================
run_metric_pair_global_fixed <- function(
    seu, group_a, group_b,
    metric = c("kbet", "ilisi", "asw"),
    reduction = "harmony",
    ndims = 20,
    k0 = 30,
    perplexity = 30,
    cap_total = Inf,
    seed = 1
) {
  metric <- match.arg(metric)
  seu_sub <- subset(seu, subset = orig.ident %in% c(group_a, group_b))
  if (ncol(seu_sub) == 0) return(NULL)
  
  if (is.finite(cap_total) && ncol(seu_sub) > cap_total) {
    set.seed(seed)
    seu_sub <- subset(seu_sub, cells = sample(colnames(seu_sub), cap_total))
  }
  
  val <- switch(metric,
                kbet  = calc_kbet_once(seu_sub, reduction, ndims, k0),
                ilisi = calc_ilisi_once(seu_sub, reduction, ndims, perplexity),
                asw   = calc_batch_asw_once(seu_sub, reduction, ndims)
  )
  if (!is.finite(val)) return(NULL)
  
  data.frame(group_a=group_a, group_b=group_b, metric=metric, value=val, n_cells=ncol(seu_sub))
}

# ============================================================
# 3.2) Global symmetric bootstrap
#   - n_each_fixed: 设定每边抽多少（建议用于“跨模型公平比较”）
#   - cap_each:     每边最大抽样上限（不设 n_each_fixed 时使用）
# ============================================================
run_metric_pair_global_sym_boot <- function(
    seu, group_a, group_b,
    metric = c("kbet", "ilisi", "asw"),
    reduction = "harmony",
    ndims = 20,
    k0 = 30,
    perplexity = 30,
    n_rep = 200,
    cap_each = Inf,
    n_each_fixed = NULL,
    seed = 1
) {
  metric <- match.arg(metric)
  
  cells_a <- WhichCells(seu, expression = orig.ident == group_a)
  cells_b <- WhichCells(seu, expression = orig.ident == group_b)
  n_a <- length(cells_a); n_b <- length(cells_b)
  if (n_a == 0 || n_b == 0) return(NULL)
  
  n_each <- if (!is.null(n_each_fixed)) {
    as.integer(n_each_fixed)
  } else {
    as.integer(min(n_a, n_b, cap_each))
  }
  
  # thresholds（和你 perCT strict 对齐，避免口径变化）
  if (metric == "kbet" && n_each < (k0 + 1)) return(NULL)
  if (metric == "ilisi" && (2 * n_each) < (3 * perplexity)) return(NULL)
  if (metric == "asw"  && (2 * n_each) < 10) return(NULL)
  
  vals <- numeric(0)
  for (r in seq_len(n_rep)) {
    set.seed(seed + r)
    take_a <- sample(cells_a, n_each)
    take_b <- sample(cells_b, n_each)
    seu_sub <- subset(seu, cells = c(take_a, take_b))
    
    v <- switch(metric,
                kbet  = calc_kbet_once(seu_sub, reduction, ndims, k0),
                ilisi = calc_ilisi_once(seu_sub, reduction, ndims, perplexity),
                asw   = calc_batch_asw_once(seu_sub, reduction, ndims)
    )
    if (is.finite(v)) vals <- c(vals, v)
  }
  
  if (length(vals) < max(10, floor(n_rep * 0.3))) return(NULL)
  ci <- quantile(vals, c(0.025, 0.975), names = FALSE)
  
  list(
    summary = data.frame(
      group_a = group_a, group_b = group_b, metric = metric,
      value = mean(vals),
      ci_low = ci[1], ci_high = ci[2],
      reps_used = length(vals),
      n_each = n_each, n_a = n_a, n_b = n_b
    ),
    per_rep = data.frame(
      group_a = group_a, group_b = group_b, metric = metric,
      rep = seq_along(vals), value = vals,
      n_each = n_each, n_a = n_a, n_b = n_b
    )
  )
}

# ============================================================
# 4) Run one metric end-to-end (tables + plot)
#    - strict 主跑
#    - 对 strict 因门槛跳过的项，relax 再补跑一次（自动降 k0/perplexity）
#    - 最终结果 strict 优先；strict 没有才用 relax
# ============================================================
run_one_metric <- function(metric_use) {
  message("\n====================")
  message("Running metric: ", metric_use)
  message("====================")
  
  n_cap <- n_cap_list[[metric_use]]
  tab_ct <- table(obj$orig.ident, obj$cell_type)
  
  # ---- strict run ----
  all_rep_strict <- list()
  all_sum_strict <- list()
  all_skip_strict <- list()
  
  for (p in pairs) {
    ga <- p[1]; gb <- p[2]
    
    if (!(ga %in% rownames(tab_ct)) || !(gb %in% rownames(tab_ct))) {
      all_skip_strict[[length(all_skip_strict)+1]] <- data.frame(
        cell_type=NA_character_, group_a=ga, group_b=gb, pair=paste(ga,"vs",gb),
        metric=metric_use, mode="strict",
        reason="orig.ident not found",
        n_a=NA_integer_, n_b=NA_integer_, n_each=NA_integer_,
        k0_used=NA_integer_, perplexity_used=NA_integer_
      )
      next
    }
    
    cts_pair <- colnames(tab_ct)[ tab_ct[ga, ] > 0 & tab_ct[gb, ] > 0 ]
    cts_pair <- setdiff(cts_pair, ignore_celltypes)
    
    if (length(cts_pair) == 0) {
      all_skip_strict[[length(all_skip_strict)+1]] <- data.frame(
        cell_type=NA_character_, group_a=ga, group_b=gb, pair=paste(ga,"vs",gb),
        metric=metric_use, mode="strict",
        reason="no shared cell_type (after ignore)",
        n_a=NA_integer_, n_b=NA_integer_, n_each=NA_integer_,
        k0_used=NA_integer_, perplexity_used=NA_integer_
      )
      next
    }
    
    for (ct in cts_pair) {
      rr <- tryCatch(
        run_metric_pair_ct_boot(
          seu = obj, celltype = ct, group_a = ga, group_b = gb,
          metric = metric_use,
          reduction = reduction_use, ndims = ndims,
          k0 = k0, perplexity = perplexity,
          n_rep = n_rep, n_cap = n_cap, seed = 123,
          mode = "strict",
          min_k0 = min_k0, min_perplexity = min_perplexity
        ),
        error = function(e) list(out=NULL, skip=data.frame(
          cell_type=ct, group_a=ga, group_b=gb, pair=paste(ga,"vs",gb),
          metric=metric_use, mode="strict",
          reason=paste0("ERROR: ", conditionMessage(e)),
          n_a=NA_integer_, n_b=NA_integer_, n_each=NA_integer_,
          k0_used=NA_integer_, perplexity_used=NA_integer_
        ))
      )
      
      if (!is.null(rr$skip)) all_skip_strict[[length(all_skip_strict)+1]] <- rr$skip
      if (is.null(rr$out)) next
      
      all_rep_strict[[length(all_rep_strict)+1]] <- rr$out$per_rep
      all_sum_strict[[length(all_sum_strict)+1]] <- rr$out$summary
    }
  }
  
  res_rep_strict <- bind_rows(all_rep_strict)
  res_sum_strict <- bind_rows(all_sum_strict)
  res_skip_strict <- bind_rows(all_skip_strict)
  empty_skip_schema <- tibble::tibble(
    cell_type = character(), group_a = character(), group_b = character(),
    pair = character(), metric = character(), mode = character(), reason = character(),
    n_a = integer(), n_b = integer(), n_each = integer(),
    k0_used = integer(), perplexity_used = integer()
  )
  if (ncol(res_skip_strict) == 0) res_skip_strict <- empty_skip_schema
  
  # ---- relax补跑：只对 “strict 因门槛跳过” 的项再跑一次 ----
  need_relax <- res_skip_strict %>%
    filter(
      !is.na(cell_type),
      mode == "strict",
      (
        (metric_use == "kbet"  & grepl("^n_each < k0\\+1", reason)) |
          (metric_use == "ilisi" & grepl("^2\\*n_each < 3\\*perplexity", reason))
      )
    ) %>%
    distinct(cell_type, group_a, group_b, pair, metric)
  
  all_rep_relax <- list()
  all_sum_relax <- list()
  all_skip_relax <- list()
  
  if (nrow(need_relax) > 0) {
    for (i in seq_len(nrow(need_relax))) {
      ct <- need_relax$cell_type[i]
      ga <- need_relax$group_a[i]
      gb <- need_relax$group_b[i]
      
      rr2 <- tryCatch(
        run_metric_pair_ct_boot(
          seu = obj, celltype = ct, group_a = ga, group_b = gb,
          metric = metric_use,
          reduction = reduction_use, ndims = ndims,
          k0 = k0, perplexity = perplexity,
          n_rep = n_rep, n_cap = n_cap, seed = 999,
          mode = "relax",
          min_k0 = min_k0, min_perplexity = min_perplexity
        ),
        error = function(e) list(out=NULL, skip=data.frame(
          cell_type=ct, group_a=ga, group_b=gb, pair=paste(ga,"vs",gb),
          metric=metric_use, mode="relax",
          reason=paste0("ERROR: ", conditionMessage(e)),
          n_a=NA_integer_, n_b=NA_integer_, n_each=NA_integer_,
          k0_used=NA_integer_, perplexity_used=NA_integer_
        ))
      )
      
      if (!is.null(rr2$skip)) all_skip_relax[[length(all_skip_relax)+1]] <- rr2$skip
      if (is.null(rr2$out)) next
      
      all_rep_relax[[length(all_rep_relax)+1]] <- rr2$out$per_rep
      all_sum_relax[[length(all_sum_relax)+1]] <- rr2$out$summary
    }
  }
  
  res_rep_relax <- bind_rows(all_rep_relax)
  res_sum_relax <- bind_rows(all_sum_relax)
  res_skip_relax <- bind_rows(all_skip_relax)
  
  # ---- 合并 & 生成 FINAL（strict优先）----
  res_rep_all <- bind_rows(res_rep_strict, res_rep_relax)
  res_sum_all <- bind_rows(res_sum_strict, res_sum_relax)
  res_skip_all <- bind_rows(res_skip_strict, res_skip_relax)
  
  if (nrow(res_sum_all) == 0) {
    # 仍然要把 skip 写出去
    write.csv(res_skip_all, file.path(outdir, paste0("SKIP_log_", metric_use, ".csv")),
              row.names = FALSE, quote = FALSE)
    stop("No results for metric: ", metric_use, " (all skipped). See skip log.")
  }
  
  res_sum_final <- res_sum_all %>%
    mutate(mode_rank = ifelse(mode == "strict", 1L, 2L)) %>%
    arrange(metric, pair, cell_type, mode_rank) %>%
    group_by(metric, pair, cell_type) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(is_relaxed = (mode == "relax")) %>%
    select(-mode_rank)
  
  res_rep_final <- res_rep_all %>%
    inner_join(res_sum_final %>% select(metric, pair, cell_type, mode),
               by = c("metric","pair","cell_type","mode"))
  
  # ---- save per-ct outputs ----
  write.csv(res_sum_all,   file.path(outdir, paste0("perCT_summary_", metric_use, "_ALL.csv")),
            row.names = FALSE, quote = FALSE)
  write.csv(res_sum_final, file.path(outdir, paste0("perCT_summary_", metric_use, "_FINAL.csv")),
            row.names = FALSE, quote = FALSE)
  
  write.csv(res_rep_all,   file.path(outdir, paste0("perCT_bootReps_", metric_use, "_ALL.csv")),
            row.names = FALSE, quote = FALSE)
  write.csv(res_rep_final, file.path(outdir, paste0("perCT_bootReps_", metric_use, "_FINAL.csv")),
            row.names = FALSE, quote = FALSE)
  
  write.csv(res_skip_all,  file.path(outdir, paste0("SKIP_log_", metric_use, ".csv")),
            row.names = FALSE, quote = FALSE)
  
  # ============================================================
  # build delta (model - baseline) from FINAL bootstrap reps
  # 关键修复：delta_sum 有 model 维度，作图必须按 model facet
  # ============================================================
  res_rep2 <- res_rep_final %>%
    mutate(pair = str_squish(pair)) %>%
    separate(pair, into = c("model", "ref"), sep = "\\s+vs\\s+", remove = FALSE) %>%
    mutate(stage = parse_stage_from_ref(ref)) %>%
    filter(!is.na(stage))
  
  if (nrow(res_rep2) == 0) {
    stop("No perCT reps with parsable stage. Check ref naming like D_E75.")
  }
  
  if (!(baseline %in% unique(res_rep2$model))) {
    stop("Baseline model not found in reps: ", baseline,
         " (check pairs and baseline name)")
  }
  
  # all models from pairs (excluding baseline)
  models_all <- setdiff(sort(unique(vapply(pairs, `[`, "", 1))), baseline)
  
  # wide by model for delta
  wide <- res_rep2 %>%
    select(cell_type, stage, metric, rep, model, value) %>%
    pivot_wider(names_from = model, values_from = value)
  
  # 只对 wide 里确实存在的 model 做 delta
  models_present <- intersect(models_all, names(wide))
  
  delta_rep <- bind_rows(lapply(models_present, function(m) {
    wide %>%
      filter(!is.na(.data[[baseline]]), !is.na(.data[[m]])) %>%
      transmute(
        cell_type, stage, metric, rep,
        model = m,
        delta = .data[[m]] - .data[[baseline]]
      )
  }))
  
  delta_sum <- delta_rep %>%
    group_by(cell_type, stage, metric, model) %>%
    summarise(
      delta_mean = mean(delta),
      ci_low  = quantile(delta, 0.025),
      ci_high = quantile(delta, 0.975),
      reps_used = n(),
      .groups = "drop"
    )
  
  write.csv(delta_sum,
            file.path(outdir, paste0("delta_perCT_", metric_use, "_minus_", baseline, "_FINAL.csv")),
            row.names = FALSE, quote = FALSE)
  
  # ============================================================
  # global delta (model - baseline)
  # ============================================================
  ref_use <- unique(res_rep2$ref)
  ref_use <- ref_use[!is.na(ref_use)]
  ref_use <- ref_use[1]  # 这批基本就是 D_E75
  
  cap_global <- if (metric_use == "asw") 800 else Inf
  
  g_base <- run_metric_pair_global_fixed(
    obj, baseline, ref_use,
    metric = metric_use, reduction = reduction_use, ndims = ndims,
    k0 = k0, perplexity = perplexity,
    cap_total = cap_global, seed = 11
  )
  if (is.null(g_base)) stop("Global baseline failed: ", metric_use)
  
  res_global_delta <- bind_rows(lapply(models_all, function(m) {
    g_m <- run_metric_pair_global_fixed(
      obj, m, ref_use,
      metric = metric_use, reduction = reduction_use, ndims = ndims,
      k0 = k0, perplexity = perplexity,
      cap_total = cap_global, seed = 12
    )
    if (is.null(g_m)) return(NULL)
    data.frame(
      stage = parse_stage_from_ref(ref_use),
      model = m,
      delta_global = g_m$value - g_base$value,
      value_model = g_m$value,
      value_baseline = g_base$value,
      n_model = g_m$n_cells,
      n_baseline = g_base$n_cells
    )
  }))
  
  write.csv(res_global_delta,
            file.path(outdir, paste0("delta_global_", metric_use, "_minus_", baseline, ".csv")),
            row.names = FALSE, quote = FALSE)
  
  # ============================================================
  # plot: facet by model（否则多模型会叠在一起，看起来“全一样/全错”）
  # ============================================================
  stage_levels <- sort(unique(delta_sum$stage))
  stage_cols <- setNames(
    c("#8E8CD8", "#B10026", "#2A8C7C", "#6E6E6E", "#C57B2A")[seq_along(stage_levels)],
    stage_levels
  )
  
  delta_sum$stage <- factor(delta_sum$stage, levels = stage_levels)
  res_global_delta$stage <- factor(res_global_delta$stage, levels = stage_levels)
  
  ylab <- switch(metric_use,
                 kbet  = "ΔkBET (model - baseline)  [kBET lower=better]",
                 ilisi = "ΔiLISI (model - baseline) [higher=better]",
                 asw   = "ΔBatch-ASW (model - baseline) [higher=better]"
  )
  
  p <- ggplot(delta_sum, aes(x = cell_type, y = delta_mean, fill = stage)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5) +
    geom_hline(
      data = res_global_delta,
      aes(yintercept = delta_global, color = stage),
      linetype = "dashed",
      linewidth = 0.9,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    geom_col(width = 0.65, color = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 0.5) +
    scale_fill_manual(values = stage_cols) +
    scale_color_manual(values = stage_cols) +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = ylab, fill = NULL) +
    theme(axis.text = element_text(color = "black")) +
    facet_wrap(~model, ncol = 2, scales = "free_y")
  
  ggsave(
    filename = file.path(outdir, paste0("plot_delta_", metric_use, "_with_global_", baseline, "_FACETmodel.pdf")),
    plot = p, width = 12, height = 10, device = pdf_device
  )
  
  # 额外：把哪些条目用了 relax 单独导出，便于你核查
  relaxed_used <- res_sum_final %>% filter(is_relaxed)
  write.csv(relaxed_used,
            file.path(outdir, paste0("RELAXED_used_", metric_use, "_FINAL.csv")),
            row.names = FALSE, quote = FALSE)
  
  list(
    res_rep_all = res_rep_all,
    res_sum_all = res_sum_all,
    res_sum_final = res_sum_final,
    res_skip = res_skip_all,
    delta_sum = delta_sum,
    res_global_delta = res_global_delta,
    plot = p
  )
}

# ============================================================
# 5) Run all 3 metrics (one shot)
# ============================================================
metrics <- c("kbet","ilisi","asw")
results <- setNames(lapply(metrics, run_one_metric), metrics)

# 汇总总表（ALL/FNAL/skip）
summary_all <- bind_rows(lapply(metrics, function(m) results[[m]]$res_sum_all))
write.csv(summary_all, file.path(outdir, "perCT_summary_ALLmetrics_ALLmodes.csv"),
          row.names = FALSE, quote = FALSE)

summary_final_all <- bind_rows(lapply(metrics, function(m) results[[m]]$res_sum_final))
write.csv(summary_final_all, file.path(outdir, "perCT_summary_ALLmetrics_FINAL.csv"),
          row.names = FALSE, quote = FALSE)

skip_all <- bind_rows(lapply(metrics, function(m) results[[m]]$res_skip))
write.csv(skip_all, file.path(outdir, "SKIP_log_ALLmetrics.csv"),
          row.names = FALSE, quote = FALSE)

# 在 R 里快速检查：每个 metric 的 FINAL 是否真的不同
lapply(results, function(x) {
  head(x$res_sum_final %>% select(metric, pair, cell_type, mode, mean, k0_used, perplexity_used))
})

add_global_sym_and_replot <- function(
    metric_use,
    n_rep_global = 200,
    cap_total_global_asw = 800,    # 你之前 global ASW 用的总 cap
    cap_each_sym_kbet_ilisi = Inf, # 允许的话不截断；想控时可设 5000/8000
    sym_use_common_n_each = TRUE   # TRUE=跨模型统一 n_each（更公平）
) {
  # ---- 读 perCT delta（你已经有）----
  f_delta <- file.path(outdir, paste0("delta_perCT_", metric_use, "_minus_", baseline, "_FINAL.csv"))
  stopifnot(file.exists(f_delta))
  delta_sum <- read.csv(f_delta, stringsAsFactors = FALSE)
  
  # ref（你这里所有 pair 的 ref 都是 D_E75）
  ref_all <- unique(vapply(pairs, `[`, "", 2))
  stopifnot(length(ref_all) == 1)
  ref_use <- ref_all[1]
  stage_use <- parse_stage_from_ref(ref_use)
  
  models_all <- setdiff(sort(unique(vapply(pairs, `[`, "", 1))), baseline)
  
  # ---- 计算 global_raw（一次）----
  cap_global_raw <- if (metric_use == "asw") cap_total_global_asw else Inf
  
  g_base_raw <- run_metric_pair_global_fixed(
    obj, baseline, ref_use,
    metric = metric_use, reduction = reduction_use, ndims = ndims,
    k0 = k0, perplexity = perplexity,
    cap_total = cap_global_raw, seed = 11
  )
  if (is.null(g_base_raw)) stop("Global RAW baseline failed: ", metric_use)
  
  g_models_raw <- bind_rows(lapply(models_all, function(m) {
    g_m <- run_metric_pair_global_fixed(
      obj, m, ref_use,
      metric = metric_use, reduction = reduction_use, ndims = ndims,
      k0 = k0, perplexity = perplexity,
      cap_total = cap_global_raw, seed = 12
    )
    if (is.null(g_m)) return(NULL)
    data.frame(
      stage = stage_use, model = m,
      value_raw_model = g_m$value,
      value_raw_baseline = g_base_raw$value,
      delta_global_raw = g_m$value - g_base_raw$value,
      n_raw_model = g_m$n_cells,
      n_raw_baseline = g_base_raw$n_cells
    )
  }))
  
  # ---- 计算 global_sym（bootstrap）----
  # 关键：为了“跨模型公平比较”，建议所有模型都用同一个 n_each_fixed
  counts <- table(obj$orig.ident)
  if (sym_use_common_n_each) {
    # 每边统一抽：min( ref, baseline, 所有模型 )，这样不同模型不会因为 n_each 不同而不可比
    common_each <- min(
      as.integer(counts[[ref_use]]),
      as.integer(counts[[baseline]]),
      min(as.integer(counts[models_all])),
      if (metric_use == "asw") floor(cap_total_global_asw/2) else as.integer(cap_each_sym_kbet_ilisi)
    )
    n_each_fixed <- common_each
  } else {
    n_each_fixed <- NULL
  }
  
  cap_each_sym <- if (metric_use == "asw") floor(cap_total_global_asw/2) else cap_each_sym_kbet_ilisi
  
  g_base_sym <- run_metric_pair_global_sym_boot(
    obj, baseline, ref_use,
    metric = metric_use, reduction = reduction_use, ndims = ndims,
    k0 = k0, perplexity = perplexity,
    n_rep = n_rep_global,
    cap_each = cap_each_sym,
    n_each_fixed = n_each_fixed,
    seed = 101
  )
  if (is.null(g_base_sym)) stop("Global SYM baseline failed: ", metric_use)
  
  g_models_sym <- bind_rows(lapply(models_all, function(m) {
    g_m <- run_metric_pair_global_sym_boot(
      obj, m, ref_use,
      metric = metric_use, reduction = reduction_use, ndims = ndims,
      k0 = k0, perplexity = perplexity,
      n_rep = n_rep_global,
      cap_each = cap_each_sym,
      n_each_fixed = n_each_fixed,
      seed = 202
    )
    if (is.null(g_m)) return(NULL)
    
    data.frame(
      stage = stage_use, model = m,
      value_sym_model = g_m$summary$value,
      value_sym_baseline = g_base_sym$summary$value,
      delta_global_sym = g_m$summary$value - g_base_sym$summary$value,
      sym_ci_low  = (g_m$summary$ci_low  - g_base_sym$summary$value),
      sym_ci_high = (g_m$summary$ci_high - g_base_sym$summary$value),
      n_each_sym = g_m$summary$n_each,
      reps_used_sym = g_m$summary$reps_used
    )
  }))
  
  # 合并 raw+sym
  res_global_both <- full_join(g_models_raw, g_models_sym, by = c("stage","model"))
  
  write.csv(
    res_global_both,
    file.path(outdir, paste0("GLOBAL_both_raw_sym_", metric_use, "_minus_", baseline, ".csv")),
    row.names = FALSE, quote = FALSE
  )
  
  # ---- 重画图：两种线型区分 raw vs sym ----
  stage_levels <- sort(unique(delta_sum$stage))
  stage_cols <- setNames(
    c("#8E8CD8", "#B10026", "#2A8C7C", "#6E6E6E", "#C57B2A")[seq_along(stage_levels)],
    stage_levels
  )
  
  ylab <- switch(metric_use,
                 kbet  = "ΔkBET (model - baseline)  [kBET lower=better]",
                 ilisi = "ΔiLISI (model - baseline) [higher=better]",
                 asw   = "ΔBatch-ASW (model - baseline) [higher=better]"
  )
  
  delta_sum$stage <- factor(delta_sum$stage, levels = stage_levels)
  
  global_lines <- bind_rows(
    res_global_both %>% transmute(
      stage = factor(stage, levels = stage_levels),
      model, global_type = "global_raw", delta_global = delta_global_raw
    ),
    res_global_both %>% transmute(
      stage = factor(stage, levels = stage_levels),
      model, global_type = "global_sym", delta_global = delta_global_sym
    )
  ) %>% filter(is.finite(delta_global))
  
  p <- ggplot(delta_sum, aes(x = cell_type, y = delta_mean, fill = stage)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.5) +
    geom_hline(
      data = global_lines,
      aes(yintercept = delta_global, color = stage, linetype = global_type),
      linewidth = 0.9,
      inherit.aes = FALSE
    ) +
    geom_col(width = 0.65, color = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 0.5) +
    scale_fill_manual(values = stage_cols) +
    scale_color_manual(values = stage_cols) +
    scale_linetype_manual(values = c(global_raw = "dashed", global_sym = "dotdash")) +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = ylab, fill = NULL, linetype = NULL) +
    theme(axis.text = element_text(color = "black")) +
    facet_wrap(~model, ncol = 2, scales = "free_y")
  
  ggsave(
    filename = file.path(outdir, paste0("plot_delta_", metric_use, "_with_global_raw_sym_", baseline, ".pdf")),
    plot = p, width = 12, height = 10, device = pdf_device
  )
  
  invisible(list(global_both = res_global_both, plot = p))
}

# 对每个指标各补一遍 global_sym，并重画图
add_global_sym_and_replot("ilisi", n_rep_global = 200)
add_global_sym_and_replot("kbet",  n_rep_global = 200)
add_global_sym_and_replot("asw",   n_rep_global = 200)





