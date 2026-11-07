#!/usr/bin/env Rscript
# Deepen-16: 修复落盘签名口径 —— 将 spot_signature_caf_myeloid.csv（GSE267401，
# 18288 spots）的 sig 列从 17 基因扰动签名改为 13 基因程序打分。
# 原文件备份为 spot_signature_caf_myeloid_17gene.csv；caf/mye 列不变。
suppressPackageStartupMessages(library(data.table))
setwd("/root/projects/crc_caf_myeloid")
old <- fread("02_results/spatial/spot_signature_caf_myeloid.csv")
new <- fread("02_results/deepen/spatial_scores_both_signatures.csv")
ng <- new[cohort == "GSE267401"]
stopifnot(nrow(old) == nrow(ng))
# GSE267401 四张切片的顺序与行序校验：逐 sample 比较 caf/mye/sig17 序列
for (s in unique(old$sample)) {
  a <- old[sample == s]; b <- ng[slice == s]
  stopifnot(nrow(a) == nrow(b))
  stopifnot(max(abs(a$caf - b$caf)) < 1e-8, max(abs(a$mye - b$mye)) < 1e-8,
            max(abs(a$sig - b$sig17)) < 1e-8)
  cat(s, "aligned,", nrow(a), "spots, max diff =",
      max(abs(a$sig - b$sig17)), "\n")
}
if (!file.exists("02_results/spatial/spot_signature_caf_myeloid_17gene.csv")) {
  fwrite(old, "02_results/spatial/spot_signature_caf_myeloid_17gene.csv")
}
for (s in unique(old$sample)) {
  old[sample == s, sig := ng[slice == s, sig13]]
}
fwrite(old, "02_results/spatial/spot_signature_caf_myeloid.csv")
cat("fixed. sig column now = 13-gene program score.\n")
