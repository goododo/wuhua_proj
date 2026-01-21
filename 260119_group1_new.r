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

group1_ref <- readRDS("/home/lushi02/project/wuhua/ref/group1ref_E3.5_E4.5.rds")
group1_new <- readRDS("/home/lushi02/project/wuhua/group1_new_noC6_noCiToti1-3.rds")

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

group1_new$cell_type <- NA
group1_new$cell_type[group1_new$seurat_clusters == 0] <- "Trophoectoderm"
group1_new$cell_type[group1_new$seurat_clusters == 1] <- "Primitive_endoderm"
group1_new$cell_type[group1_new$seurat_clusters == 2] <- "Trophoectoderm"
group1_new$cell_type[group1_new$seurat_clusters == 3] <- "Epiblast/Inner_cell_mass"
group1_new$cell_type[group1_new$seurat_clusters == 4] <- "Epiblast/Inner_cell_mass"
group1_new$cell_type[group1_new$seurat_clusters == 5] <- "Epiblast/Inner_cell_mass"

keep_ct <- sort(unique(as.character(group1_new$cell_type)))
keep_ct <- keep_ct[!is.na(keep_ct)]

group1_ref_filt <- subset(group1_ref, subset = !is.na(cell_type) & cell_type %in% keep_ct)
group1_ref_filt$cell_type <- droplevels(factor(group1_ref_filt$cell_type))

train <- group1_ref_filt
Idents(train) <- "cell_type"
train.group <- as.character(Idents(train))

test <- group1_new
Idents(test) <- "cell_type"
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

mat_test <- tryCatch(GetAssayData(test, layer = "scale.data"),
                     error = function(e) {GetAssayData(test, slot = "scale.data")
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

## ---------------------------
## 5) Small grid around topN ~ 6000
## ---------------------------
if (dir.exists("/home/zygao02/wuhua_proj/260119/group1_new_noC6_noCiToti1-3") == F) dir.create("/home/zygao02/wuhua_proj/260119/group1_new_noC6_noCiToti1-3")
outdir <- "/home/zygao02/wuhua_proj/260119/group1_new_noC6_noCiToti1-3"

# 只在 6000 左右尝试：你可以把这两个向量改得更窄/更细
Ntrain_grid <- seq(200, 5000, by = 200) 
Ntest_grid  <- seq(200, 5000, by = 200) 


# 固定其它超参（你也可改）
alpha_fix       <- 0.99
nfolds_fix      <- 10
downsample_fix  <- TRUE
sample_cells_fix <- 0  # 0 = 用 glm.predict 内部自动策略；你也可改成 300 更稳

col_fun <- colorRamp2(c(0, 0.5, 1), c("#e9e9e9", "white", "red"))

run_one <- function(Ntrain, Ntest) {
  feats <- intersect(names(sd_train)[1:min(Ntrain, length(sd_train))],
                     names(sd_test)[ 1:min(Ntest,  length(sd_test))])
  
  # 确保 features 在两边矩阵都存在
  feats <- intersect(feats, intersect(rownames(mat_train_filtered), rownames(mat_test)))
  
  if (length(feats) < 200) {
    return(list(ok = FALSE, score = NA_real_, n_used = length(feats),
                mat_plot = NULL, res = NULL))
  }
  
  res <- tryCatch(
    glm.predict(
      train.data   = mat_train_filtered,
      train.group  = train.group_filtered,
      downsample   = downsample_fix,
      sample.cells = sample_cells_fix,
      genes.used   = feats,
      test.data    = mat_test,
      test.group   = test.group,
      alpha        = alpha_fix,
      nfolds       = nfolds_fix
    ),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(list(ok = FALSE, score = NA_real_, n_used = length(feats),
                mat_plot = NULL, res = NULL))
  }
  
  glm.predict.mean <-
    apply(res$logits, 2, function(e) sapply(split(e, res$test.group), mean))
  glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
  
  mat_plot <- t(glm.predict.mean.prob)  # ref x query
  sc <- score_heatmap(mat_plot)
  
  list(ok = TRUE, score = sc, n_used = length(feats),
       mat_plot = mat_plot, res = res, feats = feats)
}

## ---------------------------
## 6) Grid run
## ---------------------------
rank_tbl <- list()
best_obj <- NULL
best_score <- -Inf

k <- 0
for (Ntr in Ntrain_grid) {
  for (Nte in Ntest_grid) {
    k <- k + 1
    cat(sprintf("Run %02d | Ntrain=%d, Ntest=%d ... ", k, Ntr, Nte))
    
    tmp <- run_one(Ntr, Nte)
    cat(sprintf("n_used=%d, score=%s\n", tmp$n_used, format(round(tmp$score, 4), nsmall = 4)))
    
    rank_tbl[[k]] <- data.frame(
      Ntrain = Ntr,
      Ntest  = Nte,
      n_used = tmp$n_used,
      score  = tmp$score,
      ok     = tmp$ok,
      stringsAsFactors = FALSE
    )
    
    if (isTRUE(tmp$ok) && !is.na(tmp$score) && tmp$score > best_score) {
      best_score <- tmp$score
      best_obj <- tmp
      best_Ntrain <- Ntr
      best_Ntest  <- Nte
    }
  }
}

rank_df <- bind_rows(rank_tbl) %>% arrange(desc(score))

if (dir.exists("/home/zygao02/wuhua_proj/260119/group1_new_noC6_noCiToti1-3") == F) dir.create("/home/zygao02/wuhua_proj/260119/group1_new_noC6_noCiToti1-3")
outdir <- "/home/zygao02/wuhua_proj/260119/group1_new_noC6_noCiToti1-3"

write.csv(rank_df, file.path(outdir, "ranked_results_topN_around200-5000.csv"),
          row.names = FALSE, quote = FALSE)

cat("\nSaved ranking to:", file.path(outdir, "ranked_results_topN_around200-5000.csv"), "\n")

## ---------------------------
## 7) Save BEST heatmap + csv + RData
## ---------------------------
if (!is.null(best_obj)) {
  cat(sprintf("\nBEST: Ntrain=%d, Ntest=%d | n_used=%d | score=%.4f\n",
              best_Ntrain, best_Ntest, best_obj$n_used, best_score))
  
  best_mat <- best_obj$mat_plot
  
  pdf(file.path(outdir, "BEST_heatmap200-5000.pdf"),
      width = 15, height = 8, useDingbats = FALSE)
  
  ht <- Heatmap(
    best_mat,
    col = col_fun,
    name = "Predicted\nSimilarity",
    column_title = paste0("query data (cell_type) | Ntest=", best_Ntest),
    row_title = paste0("ref data (cell_type2) | Ntrain=", best_Ntrain),
    show_row_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title_gp = gpar(fontsize = 16),
    column_title_gp = gpar(fontsize = 16),
    row_names_gp = gpar(fontsize = 12),
    column_names_gp = gpar(fontsize = 12),
    cell_fun = function(j, i, x, y, width, height, fill) {
      v <- best_mat[i, j]
      grid.text(sprintf("%.2f", v), x, y, gp = gpar(fontsize = 6))
    }
  )
  
  draw(ht, heatmap_legend_side = "right",
       padding = unit(c(5, 25, 5, 5), "mm"))
  dev.off()
  
  write.csv(best_mat, file.path(outdir, "BEST_heatmap200-5000.csv"),
            row.names = TRUE, quote = FALSE)
  
  save(
    list = c("group1_ref_filt", "train", "test",
             "mat_train_filtered", "mat_test",
             "sd_train", "sd_test",
             "best_Ntrain", "best_Ntest", "best_score",
             "best_obj", "rank_df"),
    file = file.path(outdir, "BEST_run200-5000.RData")
  )
  
  cat("Saved BEST heatmap/csv/RData to:", outdir, "\n")
} else {
  cat("\nNo successful run found (all failed or too few features).\n")
}

mat_plot <- as.matrix(read.csv(file.path(outdir, "BEST_heatmap200-5000.csv"),
                               row.names = 1, check.names = FALSE))

col_fun <- colorRamp2(c(0, 0.4, 1), c("#e9e9e9", "white", "red"))

# ====== 1) 设定每个格子的边长（mm）======
cell_mm <- 8  # 你可以试 4/5/6；越大越清晰但PDF越大

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

pdf(file.path(outdir, "BEST_heatmap_square200-5000.pdf"),
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
  column_names_rot = 45
)

draw(
  ht,
  heatmap_legend_side = "right",
  padding = unit(c(8, 35, 25, 25), "mm")  # 上 右 下 左
)

dev.off()


load(file.path(outdir, "BEST_run5000-20000.RData"))

# best_obj$res 里应当有 logits / probability（每个细胞一行）
# 保险起见：把“行顺序”对齐到 test 的细胞（= mat_test 的列名）
stopifnot(exists("best_obj"), exists("test"), exists("mat_test"))

logits <- best_obj$res$logits
prob   <- best_obj$res$probability

# 若 logits/prob 没有 cellname 行名，则强制赋值为 mat_test 的列名（= test cells）
if (is.null(rownames(logits)) || all(rownames(logits) == as.character(seq_len(nrow(logits))))) {
  rownames(logits) <- colnames(mat_test)
}
if (is.null(rownames(prob)) || all(rownames(prob) == as.character(seq_len(nrow(prob))))) {
  rownames(prob) <- colnames(mat_test)
}

# meta：从 test 取 orig.ident 和 cell_type
meta <- test@meta.data %>%
  dplyr::select(orig.ident, cell_type) %>%
  tibble::rownames_to_column("cell")

# 对齐到 logits 的细胞集合
meta <- meta %>% filter(cell %in% rownames(logits))
logits <- logits[meta$cell, , drop = FALSE]

col_fun <- colorRamp2(c(0, 0.5, 1), c("#e9e9e9", "white", "red"))

# 画正方形 heatmap 的函数（复用你之前逻辑）
draw_square_heatmap <- function(mat_plot, pdf_file, cell_mm = 5) {
  hm_w <- unit(ncol(mat_plot) * cell_mm, "mm")
  hm_h <- unit(nrow(mat_plot) * cell_mm, "mm")
  
  cn <- colnames(mat_plot); rn <- rownames(mat_plot)
  col_h <- max_text_width(cn, gp = gpar(fontsize = 14)) + unit(8, "mm")
  row_w <- max_text_width(rn, gp = gpar(fontsize = 14)) + unit(8, "mm")
  
  pdf_w_in <- (ncol(mat_plot) * cell_mm + 160) / 25.4
  pdf_h_in <- (nrow(mat_plot) * cell_mm + 140) / 25.4
  
  pdf(pdf_file, width = pdf_w_in, height = pdf_h_in, useDingbats = FALSE)
  ht <- Heatmap(
    mat_plot,
    col = col_fun,
    name = "Predicted\nSimilarity",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    width  = hm_w,
    height = hm_h,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 12),
    row_names_max_width = row_w,
    column_names_max_height = col_h,
    column_names_rot = 45
  )
  draw(ht, heatmap_legend_side = "right",
       padding = unit(c(8, 35, 25, 25), "mm"))
  dev.off()
}

# ====== 核心：按 orig.ident 分组重新聚合 ======
samples <- sort(unique(meta$orig.ident))
dir.create(file.path(outdir, "by_orig_ident_fromBEST"), showWarnings = FALSE)

for (s in samples) {
  meta_s <- meta %>% filter(orig.ident == s, !is.na(cell_type))
  if (nrow(meta_s) < 50) next
  
  logits_s <- logits[meta_s$cell, , drop = FALSE]
  test_group_s <- meta_s$cell_type
  
  # 1) 先按 cell_type 对每个 ref-type 的 logits 求均值
  #    得到： query(cell_type) × ref(cell_type2)
  glm.predict.mean <- apply(logits_s, 2, function(e) sapply(split(e, test_group_s), mean))
  
  # 2) 转 prob
  glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
  
  # 3) 转成 ref×query（跟你原来一致）
  mat_plot_s <- t(glm.predict.mean.prob)
  
  # 保存
  write.csv(mat_plot_s,
            file.path(outdir, "by_orig_ident_fromBEST", paste0("heatmap_", s, ".csv")),
            quote = FALSE)
  
  draw_square_heatmap(
    mat_plot_s,
    pdf_file = file.path(outdir, "by_orig_ident_fromBEST", paste0("heatmap_", s, "_square.pdf")),
    cell_mm = 5
  )
  
  cat("[OK]", s, "cells=", nrow(meta_s),
      " | ref=", nrow(mat_plot_s), " query=", ncol(mat_plot_s), "\n")
}


# =========================
# A) One-page BIG heatmap (Cell-style)
#   - row font smaller
#   - move whole plot LEFT a bit (by reducing left padding)
# =========================

indir  <- file.path(outdir, "by_orig_ident_fromBEST")

load(file.path(outdir, "BEST_run5000-20000.RData"))
ref_ct   <- rownames(best_obj$mat_plot)    # 17
query_ct <- colnames(best_obj$mat_plot)    # 17

samples <- c("ciToti4","D_E35","D_E45","EPS_blastoid","EPSC_S4",
             "TBLC_blastoid")

# ---- Cell-style muted palette ----
sample_cols <- c(
  "ciToti4"       = "#B24745",
  "D_E35"         = "#6E6E6E",
  "D_E45"       = "#6A5AA8",
  "EPS_blastoid"         = "#3B74A7",
  "EPSC_S4"         = "#2A8C7C",
  "TBLC_blastoid"   = "#C57B2A"
)
miss <- setdiff(samples, names(sample_cols))
if (length(miss) > 0) sample_cols[miss] <- "#999999"
sample_cols <- sample_cols[samples]

# ---- Heatmap color ----
col_fun <- colorRamp2(
  c(0, 0.5, 1),
  c("#F2F2F2", "#FFFFFF", "#B2182B")
)

# ---- Read per-sample matrices: by_orig_ident_fromBEST/heatmap_<sample>.csv ----
read_and_fix_17x17 <- function(s) {
  f <- file.path(indir, paste0("heatmap_", s, ".csv"))
  if (!file.exists(f)) stop("Missing file: ", f)
  mat <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))
  
  mat_full <- matrix(NA_real_,
                     nrow = length(ref_ct),
                     ncol = length(query_ct),
                     dimnames = list(ref_ct, query_ct))
  
  common_r <- intersect(rownames(mat), ref_ct)
  common_c <- intersect(colnames(mat), query_ct)
  if (length(common_r) > 0 && length(common_c) > 0) {
    mat_full[common_r, common_c] <- mat[common_r, common_c, drop = FALSE]
  }
  mat_full
}

mats <- lapply(samples, read_and_fix_17x17)

big_mat <- do.call(cbind, mats)
colnames(big_mat) <- unlist(lapply(samples, function(s) paste0(s, "|", query_ct)))

col_split <- factor(rep(samples, each = length(query_ct)), levels = samples)

top_ha <- HeatmapAnnotation(
  orig.ident = col_split,
  col = list(orig.ident = sample_cols),
  show_annotation_name = FALSE,
  annotation_height = unit(4, "mm"),
  annotation_legend_param = list(
    title = "orig.ident",
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

# ---- Size ----
cell_mm <- 8
hm_w <- unit(ncol(big_mat) * cell_mm, "mm")
hm_h <- unit(nrow(big_mat) * cell_mm, "mm")

pdf_file <- file.path(indir, "ALL_6samples_BIG_onepage_CellStyle.pdf")
pdf(pdf_file,
    width  = (ncol(big_mat) * cell_mm + 125) / 25.4,
    height = (nrow(big_mat) * cell_mm + 115) / 25.4,
    useDingbats = FALSE)

ht <- Heatmap(
  big_mat,
  col = col_fun,
  name = "Predicted\nSimilarity",
  top_annotation = top_ha,
  
  column_split = col_split,
  column_gap = unit(2, "mm"),
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),   # ✅ 1) 行字体更小（10 -> 8）
  row_title_gp = gpar(fontsize = 13, fontface = "bold"),
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  
  width  = hm_w,
  height = hm_h,
  
  na_col = "#F0F0F0",
  
  column_title = "Query: 6 samples × 3 cell types (columns split by orig.ident)",
  row_title    = "Reference cell types (3)"
)

draw(
  ht,
  heatmap_legend_side = "right",
  padding = unit(c(6, 25, 6, 6), "mm")  # ✅ 2) 整体向左：左边距 10mm -> 6mm
)

dev.off()

cat("Saved one-page BIG heatmap to:\n  ", pdf_file, "\n", sep = "")
cat("Dim(big_mat) = ", nrow(big_mat), " x ", ncol(big_mat), "\n", sep = "")



