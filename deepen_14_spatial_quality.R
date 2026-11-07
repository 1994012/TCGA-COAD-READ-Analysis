#!/usr/bin/env Rscript
# Deepen-14: 空间异例切片 P2/P3 的质量与构成复核 —— spot 数、UMI 深度、基因数、
# 上皮（肿瘤）含量代理指标，与同队列其余切片对比，给出异质性的具体解释。
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(ggplot2)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
dir.create("03_figures", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
theme_set(theme_classic(base_family = "sans"))

read_spots <- function(proc_dir) {
  mdir <- file.path(proc_dir, "filtered_feature_bc_matrix")
  bc <- readLines(gzfile(file.path(mdir, "barcodes.tsv.gz")))
  feat <- read.delim(gzfile(file.path(mdir, "features.tsv.gz")), header = FALSE)
  m <- readMM(gzfile(file.path(mdir, "matrix.mtx.gz")))
  m <- as(as(m, "dgTMatrix"), "dgCMatrix")  # 该数据 mtx 已是 genes x spots
  dimnames(m) <- list(as.character(feat$V2), bc)
  m
}

# 与 step10b 相同的 QC 阈值：nFeature > 200 & nCount > 500
# 上皮/肿瘤含量代理：EPCAM + KRT 家族平均表达（log 归一化后），spot 得分高于
# 全样本均值 + 1 SD 记为肿瘤高表达 spot。
epi_genes <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CEACAM5", "VIL1")

rows <- list()
for (s in c("Ajou_Visium_P1", "Ajou_Visium_P2", "Ajou_Visium_P3", "Ajou_Visium_P4")) {
  d <- file.path("00_data/spatial/gse226997_processed", s)
  m <- read_spots(d)
  cnt <- Matrix::colSums(m); feat_n <- Matrix::colSums(m > 0)
  keep <- feat_n > 200 & cnt > 500
  mk <- rownames(m) %in% epi_genes
  lib <- cnt
  mk_log <- t(log1p(t(m[mk, , drop = FALSE]) / lib * 1e4))
  epi_score <- as.numeric(Matrix::colMeans(mk_log))
  thr <- mean(epi_score) + sd(epi_score)
  rows[[s]] <- data.table(
    sample = s, cohort = "GSE226997",
    spots_raw = ncol(m), spots_qc = sum(keep),
    qc_pass_pct = 100 * sum(keep) / ncol(m),
    median_umi = median(cnt), median_genes = median(feat_n),
    median_umi_qc = median(cnt[keep]),
    epi_score_mean = mean(epi_score), epi_score_sd = sd(epi_score),
    tumor_spot_pct = 100 * mean(epi_score > thr))
  log_msg(s, ": raw=", ncol(m), " qc=", sum(keep),
          " medUMI=", round(median(cnt)), " tumor%=", round(100*mean(epi_score > thr), 1))
}
# GSE267401 作对照（如目录存在）
g2 <- c("CTC17M", "CTC17P", "CTC21M", "CTC21P")
for (s in g2) {
  d <- file.path("00_data/spatial/gse267401_processed", s)
  if (!dir.exists(file.path(d, "filtered_feature_bc_matrix"))) next
  m <- read_spots(d)
  cnt <- Matrix::colSums(m); feat_n <- Matrix::colSums(m > 0)
  keep <- feat_n > 200 & cnt > 500
  mk <- rownames(m) %in% epi_genes
  lib <- cnt
  mk_log <- t(log1p(t(m[mk, , drop = FALSE]) / lib * 1e4))
  epi_score <- as.numeric(Matrix::colMeans(mk_log))
  thr <- mean(epi_score) + sd(epi_score)
  rows[[s]] <- data.table(
    sample = s, cohort = "GSE267401",
    spots_raw = ncol(m), spots_qc = sum(keep),
    qc_pass_pct = 100 * sum(keep) / ncol(m),
    median_umi = median(cnt), median_genes = median(feat_n),
    median_umi_qc = median(cnt[keep]),
    epi_score_mean = mean(epi_score), epi_score_sd = sd(epi_score),
    tumor_spot_pct = 100 * mean(epi_score > thr))
  log_msg(s, ": raw=", ncol(m), " qc=", sum(keep),
          " medUMI=", round(median(cnt)), " tumor%=", round(100*mean(epi_score > thr), 1))
}
out <- rbindlist(rows)
fwrite(out, "02_results/deepen/spatial_quality_check.csv")

# ---- 图：三联条形（QC 通过率、中位 UMI、肿瘤 spot 占比）----
out[, slabel := gsub("Ajou_Visium_", "", sample)]
out[, slabel := factor(slabel, levels = out[cohort == "GSE226997", slabel])]
p1 <- ggplot(out[cohort == "GSE226997"], aes(slabel, qc_pass_pct, fill = slabel)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f", qc_pass_pct)), vjust = -0.35, size = 3.2) +
  scale_fill_manual(values = c(P1 = "#4393C3", P2 = "#D6604D", P3 = "#D6604D", P4 = "#4393C3")) +
  labs(x = NULL, y = "QC 通过率（%）", title = "spot 质量")
p2 <- ggplot(out[cohort == "GSE226997"], aes(slabel, median_umi, fill = slabel)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = comma <- sprintf("%d", round(median_umi))), vjust = -0.35, size = 3.2) +
  scale_fill_manual(values = c(P1 = "#4393C3", P2 = "#D6604D", P3 = "#D6604D", P4 = "#4393C3")) +
  labs(x = NULL, y = "中位 UMI / spot", title = "测序深度")
p3 <- ggplot(out[cohort == "GSE226997"], aes(slabel, tumor_spot_pct, fill = slabel)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", tumor_spot_pct)), vjust = -0.35, size = 3.2) +
  scale_fill_manual(values = c(P1 = "#4393C3", P2 = "#D6604D", P3 = "#D6604D", P4 = "#4393C3")) +
  labs(x = NULL, y = "肿瘤高表达 spot 占比（%）", title = "肿瘤区域构成")
if (requireNamespace("patchwork", quietly = TRUE)) {
  pp <- (p1 | p2 | p3) + patchwork::plot_annotation(
    title = "GSE226997 逐切片质量与构成复核（红：异例切片 P2 / P3）",
    theme = theme(plot.title = element_text(size = 12, face = "bold")))
  ggsave("03_figures/FigD6_spatial_quality.png", pp, width = 9.6, height = 3.6, dpi = 300)
} else {
  cowplot::plot_grid(p1, p2, p3, nrow = 1)
  ggsave("03_figures/FigD6_spatial_quality.png",
         cowplot::plot_grid(p1, p2, p3, nrow = 1), width = 9.6, height = 3.6, dpi = 300)
}
log_msg("done")
