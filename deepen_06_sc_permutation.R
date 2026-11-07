#!/usr/bin/env Rscript
# Deepen-06: 单细胞侧的置换检验
# 报告给出"髓系 0.386 高于 T 细胞 0.068"与"肌成纤维 -> SPP1+ 权重 0.014"，
# 但没有零模型。这里做标签置换，给经验 p 与效应量 CI。
suppressPackageStartupMessages({
  library(data.table); library(survival)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
set.seed(20260924)

# 用与主报告共识表同一口径的子集对象（15003 细胞：髓系 1639、T 5660 ...）
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("02_results/scrna/gse132465_scored_subset.rds")
meta <- NULL
if (inherits(obj, "Seurat")) {
  meta <- as.data.table(obj@meta.data)
} else if (is.data.frame(obj)) {
  meta <- as.data.table(obj)
} else if (is.list(obj)) {
  cand <- obj[vapply(obj, is.data.frame, logical(1))]
  if (length(cand)) meta <- as.data.table(cand[[which.max(vapply(cand, nrow, integer(1)))]])
}
if (is.null(meta)) stop("cannot extract metadata from gse132465_scored.rds")
log_msg("meta: ", nrow(meta), " cells; cols: ", paste(head(names(meta), 25), collapse = ", "))

ct_col <- names(meta)[grepl("celltype|cell_type|CellType|major", names(meta), ignore.case = TRUE)][1]
cons_col <- names(meta)[grepl("consensus", names(meta), ignore.case = TRUE)][1]
snai1_col <- names(meta)[grepl("^SNAI1$|SNAI1$", names(meta), ignore.case = TRUE)][1]
log_msg("using ct=", ct_col, " consensus=", cons_col, " snai1=", snai1_col)

meta <- meta[!is.na(get(ct_col)) & !is.na(get(cons_col))]
ct <- as.character(meta[[ct_col]]); cons <- as.numeric(meta[[cons_col]])
log_msg("cell types: ", paste(sort(unique(ct)), collapse = " | "))

perm_diff <- function(ct, vals, a, b, B = 2000) {
  obs <- mean(vals[ct == a]) - mean(vals[ct == b])
  null <- vapply(seq_len(B), function(i) {
    sh <- sample(ct)
    mean(vals[sh == a]) - mean(vals[sh == b])
  }, numeric(1))
  p_two <- (1 + sum(abs(null) >= abs(obs))) / (1 + B)
  list(obs = obs, null_mean = mean(null), null_sd = sd(null),
       z = (obs - mean(null)) / sd(null), p = p_two,
       ci = quantile(null, c(0.025, 0.975)))
}

res <- list()
pairs <- list(c("Myeloids", "T cells"), c("Myeloids", "Stromal cells"), c("Stromal cells", "T cells"))
for (pr in pairs) {
  if (!(pr[1] %in% ct && pr[2] %in% ct)) next
  r <- perm_diff(ct, cons, pr[1], pr[2])
  res[[length(res) + 1]] <- data.table(statistic = "13-gene program consensus",
                                       contrast = paste0(pr[1], " - ", pr[2]),
                                       observed = r$obs, null_mean = r$null_mean, null_sd = r$null_sd,
                                       z = r$z, p_permutation = r$p,
                                       null_2.5 = r$ci[1], null_97.5 = r$ci[2],
                                       n_a = sum(ct == pr[1]), n_b = sum(ct == pr[2]))
}
snai1_vec <- NULL
if ("SNAI1" %in% rownames(obj)) {
  snai1_vec <- tryCatch(as.numeric(GetAssayData(obj, layer = "data")["SNAI1", ]),
                        error = function(e) as.numeric(GetAssayData(obj, assay = "RNA", slot = "data")["SNAI1", ]))
}
if (!is.na(snai1_col)) {
  v <- as.numeric(meta[[snai1_col]])
  ok <- !is.na(v); ct2 <- ct[ok]; v <- v[ok]
  for (pr in pairs) {
    if (!(pr[1] %in% ct2 && pr[2] %in% ct2)) next
    r <- perm_diff(ct2, v, pr[1], pr[2])
    res[[length(res) + 1]] <- data.table(statistic = "SNAI1 gene mean",
                                         contrast = paste0(pr[1], " - ", pr[2]),
                                         observed = r$obs, null_mean = r$null_mean, null_sd = r$null_sd,
                                         z = r$z, p_permutation = r$p,
                                         null_2.5 = r$ci[1], null_97.5 = r$ci[2],
                                         n_a = sum(ct2 == pr[1]), n_b = sum(ct2 == pr[2]))
  }
}
if (!is.null(snai1_vec)) {
  v <- snai1_vec; ok <- !is.na(v); ct2 <- ct[ok]; v <- v[ok]
  for (pr in pairs) {
    if (!(pr[1] %in% ct2 && pr[2] %in% ct2)) next
    r <- perm_diff(ct2, v, pr[1], pr[2])
    res[[length(res) + 1]] <- data.table(statistic = "SNAI1 normalized expression",
                                         contrast = paste0(pr[1], " - ", pr[2]),
                                         observed = r$obs, null_mean = r$null_mean, null_sd = r$null_sd,
                                         z = r$z, p_permutation = r$p,
                                         null_2.5 = r$ci[1], null_97.5 = r$ci[2],
                                         n_a = sum(ct2 == pr[1]), n_b = sum(ct2 == pr[2]))
  }
}
tab <- rbindlist(res)
fwrite(tab, "02_results/deepen/scrna_permutation.csv")
print(tab)

# ---- 亚群比较：SPP1+ vs 促炎髓系（Wilcoxon + Cliff's delta + 置换） ----
sub_obj <- try(readRDS("02_results/scrna/gse132465_scored_subset.rds"), silent = TRUE)
if (!inherits(sub_obj, "try-error")) {
  m2 <- if (inherits(sub_obj, "Seurat")) as.data.table(sub_obj@meta.data) else as.data.table(sub_obj)
  st_col <- names(m2)[grepl("subtype|subtype_label|cluster", names(m2), ignore.case = TRUE)][1]
  if (!is.na(st_col)) {
    st <- as.character(m2[[st_col]])
    v <- as.numeric(m2[[cons_col]])
    a <- "SPP1+"; b <- "Myeloid_Inflammatory"
    cand_a <- unique(st)[grepl("SPP1", unique(st), ignore.case = TRUE)]
    cand_b <- unique(st)[grepl("Inflamm", unique(st), ignore.case = TRUE)]
    if (length(cand_a) && length(cand_b)) {
      a <- cand_a[1]; b <- cand_b[1]
      va <- v[st == a & !is.na(v)]; vb <- v[st == b & !is.na(v)]
      w <- wilcox.test(va, vb)
      cliff <- mean(outer(va, vb, ">")) - mean(outer(va, vb, "<"))
      fwrite(data.table(contrast = paste0(a, " vs ", b), n_a = length(va), n_b = length(vb),
                        median_a = median(va), median_b = median(vb),
                        wilcox_p = w$p.value, cliffs_delta = cliff),
             "02_results/deepen/scrna_subtype_contrast.csv")
      log_msg("subtype contrast written: ", a, " vs ", b)
    }
  }
}
log_msg("DONE deepen_06")
