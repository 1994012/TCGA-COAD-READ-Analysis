#!/usr/bin/env Rscript
# Deepen-19: 空间侧 SPP1 共享基因敏感性分析
#
# 背景：髓系评分标志基因含 SPP1(CD68/LYZ/C1QA/SPP1)，而 SPP1 同时是 13 基因签名成员。
#       sig13 与 mye4 共用 SPP1，rho(sig,mye) 可能部分来源于共享基因而非生物学共定位。
# 本脚本重算髓系评分（剔除 SPP1），检验该结论是否依赖共享基因。
#
# 口径与 deepen_08_spatial_13gene.R 完全一致：SCT data layer，标志基因取算术均值（非 beta 加权）。
# 敏感性变体：
#   sig13  x mye4  原口径（含 SPP1 共享，基线）
#   sig13  x mye3  髓系剔除 SPP1  → 直接检验共享基因贡献（主分析）
#   sig12  x mye3  签名与髓系均剔除 SPP1 → 完全无共享基因（最严格）
# 另附髓系标志基因逐一剔除（LOO）与逐基因诊断。
suppressPackageStartupMessages({
  library(data.table); library(Seurat); library(Matrix)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
set.seed(20260924)

beta <- fread("02_results/model/model_genes_beta.csv")$gene
sig12 <- setdiff(beta, "SPP1")
caf_markers  <- c("COL1A1", "FAP", "PDGFRB", "LUM")
mye_markers  <- c("SPP1", "CD68", "LYZ", "C1QA")
mye_markers_noSPP1 <- setdiff(mye_markers, "SPP1")

LAB_A <- "A 原始口径（签名13 x 髓系4，含SPP1共享）"
LAB_B <- "B 髓系剔除SPP1（签名13 x 髓系3）"
LAB_C <- "C 签名与髓系均剔除SPP1（12 x 3）"

log_msg("sig13 : ", paste(beta, collapse = ","))
log_msg("sig12 : ", paste(sig12, collapse = ","))
log_msg("mye4  : ", paste(mye_markers, collapse = ","))
log_msg("mye3  : ", paste(mye_markers_noSPP1, collapse = ","))
stopifnot("SPP1" %in% beta, "SPP1" %in% mye_markers)

mean_score <- function(exp, genes) {  sg <- intersect(genes, rownames(exp))
  if (length(sg) < 1) return(rep(NA_real_, ncol(exp)))
  as.numeric(Matrix::colMeans(exp[sg, , drop = FALSE]))
}
get_exp <- function(so) {
  a <- if ("SCT" %in% names(so@assays)) "SCT" else DefaultAssay(so)
  lyr <- tryCatch(Layers(so[[a]]), error = function(e) "data")
  lyr <- if (all(c("data", "counts", "scale.data") %in% lyr)) "data" else lyr[1]
  GetAssayData(so, assay = a, layer = lyr)
}
pspear <- function(x, y, z) {
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  cor(lm(rx ~ rz)$residuals, lm(ry ~ rz)$residuals)
}

## ---- 1. 汇总 8 张切片的 spot 级分数 ----
collect <- function(objs, cohort) {
  rbindlist(lapply(seq_along(objs), function(i) {
    so <- objs[[i]]; exp <- get_exp(so)
    nm <- if (!is.null(so$sample)) as.character(so$sample)[1] else names(objs)[i]
    d <- data.table(cohort = cohort, slice = nm,
                    sig13 = mean_score(exp, beta),
                    sig12 = mean_score(exp, sig12),
                    caf   = mean_score(exp, caf_markers),
                    mye4  = mean_score(exp, mye_markers),
                    mye3  = mean_score(exp, mye_markers_noSPP1))
    for (g in c("SPP1", mye_markers_noSPP1)) d[[g]] <- as.numeric(exp[g, ])
    d
  }))
}

objs1 <- readRDS("02_results/spatial/spatial_objects.rds")
if (!is.list(objs1)) objs1 <- list(objs1)
objs2 <- readRDS("02_results/spatial/gse226997_objects.rds")
if (!is.list(objs2)) objs2 <- list(objs2)
log_msg("GSE267401 objects: ", length(objs1), "  GSE226997 objects: ", length(objs2))

pool <- rbindlist(list(collect(objs1, "GSE267401"), collect(objs2, "GSE226997")))
fwrite(pool, "02_results/deepen/spatial_scores_spp1.csv")
log_msg("pooled spots: ", nrow(pool))

## ---- 1b. 回源核验：用原方法重算的 sig13/caf/mye4 应与 deepen_08 落盘值逐 spot 一致 ----
# 两个表来自同一对象遍历顺序，故以 (cohort, slice, spot_id) 对齐；避免按组 merge 造成笛卡尔积。
pool[, spot_id := seq_len(.N), by = .(cohort, slice)]
prev <- fread("02_results/deepen/spatial_scores_both_signatures.csv")
prev[, spot_id := seq_len(.N), by = .(cohort, slice)]
stopifnot(nrow(prev) == nrow(pool))
m <- merge(pool, prev[, .(cohort, slice, spot_id, sig13_old = sig13, caf_old = caf, mye_old = mye)],
           by = c("cohort", "slice", "spot_id"))
stopifnot(nrow(m) == nrow(pool))
recon <- m[, .(n = .N,
               max_abs_d_sig13 = max(abs(sig13 - sig13_old)),
               max_abs_d_caf   = max(abs(caf   - caf_old)),
               max_abs_d_mye   = max(abs(mye4  - mye_old))), by = .(cohort, slice)]
fwrite(recon, "02_results/deepen/spatial_spp1_reconciliation.csv")
log_msg("回源核验 max|d| (应≈0): sig13=", max(recon$max_abs_d_sig13),
        "  caf=", max(recon$max_abs_d_caf), "  mye=", max(recon$max_abs_d_mye))
stopifnot(max(recon$max_abs_d_sig13) < 1e-9,
          max(recon$max_abs_d_caf)   < 1e-9,
          max(recon$max_abs_d_mye)   < 1e-9)

## ---- 2. 逐切片 × 定义 的主表 ----
analyse <- function(d, sigcol, myecol, label, signame, myename) {
  s <- d[complete.cases(d[, c(sigcol, myecol, "caf"), with = FALSE])]
  s <- s[, .(slice, sig = get(sigcol), caf, mye = get(myecol), cohort)]
  rbindlist(lapply(unique(s$slice), function(sl) {
    s1 <- s[slice == sl]
    if (nrow(s1) < 50) return(NULL)
    r_sm <- suppressWarnings(cor(s1$sig, s1$mye, method = "spearman"))
    r_sc <- suppressWarnings(cor(s1$sig, s1$caf, method = "spearman"))
    r_mc <- suppressWarnings(cor(s1$caf, s1$mye, method = "spearman"))
    B <- 1000
    bs <- replicate(B, {
      i <- sample(seq_len(nrow(s1)), replace = TRUE)
      c(cor(s1$sig[i], s1$mye[i], method = "spearman"),
        cor(s1$sig[i], s1$caf[i], method = "spearman"))
    })
    d_obs <- r_sm - r_sc
    d_bs <- bs[1, ] - bs[2, ]
    null <- replicate(2000, {
      sw <- sample(c(TRUE, FALSE), nrow(s1), replace = TRUE)
      mye2 <- s1$mye; caf2 <- s1$caf
      mye2[sw] <- s1$caf[sw]; caf2[sw] <- s1$mye[sw]
      cor(s1$sig, mye2, method = "spearman") - cor(s1$sig, caf2, method = "spearman")
    })
    p_perm <- (1 + sum(abs(null) >= abs(d_obs))) / (1 + 2000)
    q <- quantile(s1$sig, 0.75); hi <- s1$sig >= q
    data.table(definition = label, cohort = s1$cohort[1], slice = sl,
               signature_genes = signame, myeloid_genes = myename,
               shared_genes = length(intersect(
                 if (signame == "13") beta else sig12,
                 if (myename == "4") mye_markers else mye_markers_noSPP1)),
               n_spots = nrow(s1),
               rho_sig_mye = r_sm, rho_sig_mye_lo = quantile(bs[1, ], 0.025),
               rho_sig_mye_hi = quantile(bs[1, ], 0.975),
               rho_sig_caf = r_sc, rho_sig_caf_lo = quantile(bs[2, ], 0.025),
               rho_sig_caf_hi = quantile(bs[2, ], 0.975),
               rho_caf_mye = r_mc,
               diff_mye_minus_caf = d_obs,
               diff_lo = quantile(d_bs, 0.025), diff_hi = quantile(d_bs, 0.975),
               p_paired_permutation = p_perm,
               partial_sig_mye_given_caf = pspear(s1$sig, s1$mye, s1$caf),
               partial_sig_caf_given_mye = pspear(s1$sig, s1$caf, s1$mye),
               hotspot_prop_mye = mean(s1$mye[hi] >= quantile(s1$mye, 0.75)),
               hotspot_prop_caf = mean(s1$caf[hi] >= quantile(s1$caf, 0.75)))
  }))
}

res <- rbindlist(list(
  analyse(pool, "sig13", "mye4", LAB_A, "13", "4"),
  analyse(pool, "sig13", "mye3", LAB_B, "13", "3"),
  analyse(pool, "sig12", "mye3", LAB_C, "12", "3")
))
fwrite(res, "02_results/deepen/spatial_spp1_sensitivity.csv")
log_msg("--- 主表（rho_sig_mye / rho_sig_caf / diff / partial） ---")
print(res[, .(definition, cohort, slice, n_spots, rho_sig_mye, rho_sig_caf,
              diff_mye_minus_caf, diff_lo, diff_hi, p_paired_permutation,
              partial_sig_mye_given_caf)])

## ---- 3. 合并表：同一切片下 SPP1 造成的相关增量 ----
A <- res[definition == LAB_A, .(cohort, slice, n_spots,
                               rho_mye_A = rho_sig_mye, rho_caf_A = rho_sig_caf,
                               partial_A = partial_sig_mye_given_caf,
                               partial_caf_A = partial_sig_caf_given_mye)]
B <- res[definition == LAB_B, .(cohort, slice,
                               rho_mye_B = rho_sig_mye,
                               partial_B = partial_sig_mye_given_caf)]
C <- res[definition == LAB_C, .(cohort, slice,
                               rho_mye_C = rho_sig_mye, rho_caf_C = rho_sig_caf,
                               partial_C = partial_sig_mye_given_caf)]
w <- merge(merge(A, B, by = c("cohort", "slice")), C, by = c("cohort", "slice"))
stopifnot(nrow(w) == 8)
w[, rho_sig_caf := rho_caf_A]        # CAF 评分与髓系定义无关，三种定义下相同
w[, delta_rho_spp1 := rho_mye_A - rho_mye_B]          # 共享 SPP1 贡献的相关增量
w[, retained_frac := rho_mye_B / rho_mye_A]          # 剔除 SPP1 后保留的比例
w[, myeloid_dominant_A := rho_mye_A > rho_caf_A]
w[, myeloid_dominant_B := rho_mye_B > rho_caf_A]
w[, myeloid_dominant_C := rho_mye_C > rho_caf_C]
w[, partial_pos_A := partial_A > 0]
w[, partial_pos_B := partial_B > 0]
w[, partial_pos_C := partial_C > 0]
setorder(w, cohort, slice)
fwrite(w, "02_results/deepen/spatial_spp1_comparison.csv")
log_msg("--- 共享 SPP1 的相关增量 ---")
print(w[, .(cohort, slice, n_spots, rho_mye_A, rho_mye_B, rho_mye_C, rho_sig_caf,
            delta_rho_spp1, retained_frac, partial_A, partial_B, partial_C)])
log_msg("髓系占优切片数 A/B/C: ", sum(w$myeloid_dominant_A), "/",
        sum(w$myeloid_dominant_B), "/", sum(w$myeloid_dominant_C))
log_msg("偏相关为正切片数 A/B/C: ", sum(w$partial_pos_A), "/",
        sum(w$partial_pos_B), "/", sum(w$partial_pos_C))

## ---- 4. 队列合并（pooled） ----
pooled <- rbindlist(lapply(c("GSE267401", "GSE226997"), function(ch) {
  s <- pool[cohort == ch & complete.cases(sig13, sig12, caf, mye4, mye3)]
  rbindlist(list(
    data.table(cohort = ch, definition = "A 原始口径（签名13 x 髓系4，含SPP1共享）",
               signature_genes = "13", myeloid_genes = "4", n = nrow(s),
               pooled_rho_sig_mye = cor(s$sig13, s$mye4, method = "spearman"),
               pooled_rho_sig_caf = cor(s$sig13, s$caf,  method = "spearman"),
               pooled_rho_caf_mye = cor(s$caf,   s$mye4, method = "spearman")),
    data.table(cohort = ch, definition = "B 髓系剔除SPP1（签名13 x 髓系3）",
               signature_genes = "13", myeloid_genes = "3", n = nrow(s),
               pooled_rho_sig_mye = cor(s$sig13, s$mye3, method = "spearman"),
               pooled_rho_sig_caf = cor(s$sig13, s$caf,  method = "spearman"),
               pooled_rho_caf_mye = cor(s$caf,   s$mye3, method = "spearman")),
    data.table(cohort = ch, definition = "C 签名与髓系均剔除SPP1（12 x 3）",
               signature_genes = "12", myeloid_genes = "3", n = nrow(s),
               pooled_rho_sig_mye = cor(s$sig12, s$mye3, method = "spearman"),
               pooled_rho_sig_caf = cor(s$sig12, s$caf,  method = "spearman"),
               pooled_rho_caf_mye = cor(s$caf,   s$mye3, method = "spearman"))
  ))
}))
fwrite(pooled, "02_results/deepen/spatial_spp1_pooled.csv")
print(pooled)

## ---- 5. 逐基因诊断：共享基因到底贡献了多少 ----
diag <- rbindlist(lapply(unique(pool$slice), function(sl) {
  s <- pool[slice == sl]
  if (nrow(s) < 50) return(NULL)
  data.table(cohort = s$cohort[1], slice = sl, n_spots = nrow(s),
    rho_sig13_SPP1 = cor(s$sig13, s$SPP1, method = "spearman"),
    rho_sig12_SPP1 = cor(s$sig12, s$SPP1, method = "spearman"),
    rho_caf_SPP1   = cor(s$caf,   s$SPP1, method = "spearman"),
    rho_sig13_CD68 = cor(s$sig13, s$CD68, method = "spearman"),
    rho_sig13_LYZ  = cor(s$sig13, s$LYZ,  method = "spearman"),
    rho_sig13_C1QA = cor(s$sig13, s$C1QA, method = "spearman"),
    rho_sig12_CD68 = cor(s$sig12, s$CD68, method = "spearman"),
    rho_sig12_LYZ  = cor(s$sig12, s$LYZ,  method = "spearman"),
    rho_sig12_C1QA = cor(s$sig12, s$C1QA, method = "spearman"),
    # SPP1 对签名评分方差的解释比例（单基因线性 R^2）
    spp1_r2_sig13 = cor(s$sig13, s$SPP1, method = "spearman")^2,
    n_spp1_pos = sum(s$SPP1 > 0))
}))
fwrite(diag, "02_results/deepen/spatial_spp1_gene_diagnostics.csv")
log_msg("--- 逐基因诊断 ---")
print(diag)

## ---- 6. 髓系标志基因逐一剔除（LOO）：检验结论是否只依赖 SPP1 ----
loo <- rbindlist(lapply(unique(pool$slice), function(sl) {
  s <- pool[slice == sl]
  if (nrow(s) < 50) return(NULL)
  rbindlist(lapply(c(NA, mye_markers), function(drop) {
    gs <- if (is.na(drop)) mye_markers else setdiff(mye_markers, drop)
    sc <- as.numeric(Matrix::rowMeans(s[, ..gs, with = FALSE]))
    data.table(cohort = s$cohort[1], slice = sl, dropped = if (is.na(drop)) "none (4 genes)" else drop,
               myeloid_genes = paste(gs, collapse = "/"),
               rho_sig13_mye = cor(s$sig13, sc, method = "spearman"),
               rho_sig12_mye = cor(s$sig12, sc, method = "spearman"),
               rho_caf_mye   = cor(s$caf,   sc, method = "spearman"),
               partial_sig13_mye_given_caf = pspear(s$sig13, sc, s$caf))
  }))
}))
fwrite(loo, "02_results/deepen/spatial_mye_markers_loo.csv")
log_msg("--- 髓系标志基因 LOO ---")
print(dcast(loo, cohort + slice ~ dropped, value.var = "rho_sig13_mye"))

log_msg("DONE deepen_19")
