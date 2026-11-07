#!/usr/bin/env Rscript
# Deepen-08: 空间分数按报告口径（13 基因回退）重算，再与落盘版（17 基因扰动签名）对照
# 背景：02_results/spatial/step10_run_meta.csv 记录的最终运行 signature_source=perturbation_signature(17 genes)，
#       与报告正文"空转使用 13 基因回退"不一致。这里两套都算，避免把 17 基因基质签名的结论读成 13 基因程序。
suppressPackageStartupMessages({
  library(data.table); library(Seurat)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
set.seed(20260924)

beta <- fread("02_results/model/model_genes_beta.csv")$gene
pert <- readLines("02_results/perturbation/perturbation_signature.txt")
pert <- pert[nchar(pert) > 0]
caf_markers <- c("COL1A1", "FAP", "PDGFRB", "LUM")
mye_markers <- c("SPP1", "CD68", "LYZ", "C1QA")
log_msg("13-gene: ", paste(beta, collapse = ","))
log_msg("17-gene perturbation signature: ", paste(pert, collapse = ","))

mean_score <- function(exp, genes) {
  sg <- intersect(genes, rownames(exp))
  if (length(sg) < 1) return(rep(NA_real_, ncol(exp)))
  as.numeric(Matrix::colMeans(exp[sg, , drop = FALSE]))
}
get_exp <- function(so) {
  a <- if ("SCT" %in% names(so@assays)) "SCT" else DefaultAssay(so)
  lyr <- tryCatch(Layers(so[[a]]), error = function(e) "data")
  lyr <- if (all(c("data", "counts", "scale.data") %in% lyr)) "data" else lyr[1]
  GetAssayData(so, assay = a, layer = lyr)
}

objs <- readRDS("02_results/spatial/spatial_objects.rds")
if (!is.list(objs)) objs <- list(objs)
log_msg("objects: ", length(objs))

pool <- rbindlist(lapply(seq_along(objs), function(i) {
  so <- objs[[i]]
  exp <- get_exp(so)
  nm <- if (!is.null(so$sample)) as.character(so$sample)[1] else names(objs)[i]
  data.table(cohort = "GSE267401", slice = nm,
             sig13 = mean_score(exp, beta),
             sig17 = mean_score(exp, pert),
             caf = mean_score(exp, caf_markers),
             mye = mean_score(exp, mye_markers))
}))
objs2 <- try(readRDS("02_results/spatial/gse226997_objects.rds"), silent = TRUE)
if (!inherits(objs2, "try-error")) {
  if (!is.list(objs2)) objs2 <- list(objs2)
  pool2 <- rbindlist(lapply(seq_along(objs2), function(i) {
    so <- objs2[[i]]
    exp <- get_exp(so)
    nm <- if (!is.null(so$sample)) as.character(so$sample)[1] else names(objs2)[i]
    data.table(cohort = "GSE226997", slice = nm,
               sig13 = mean_score(exp, beta),
               sig17 = mean_score(exp, pert),
               caf = mean_score(exp, caf_markers),
               mye = mean_score(exp, mye_markers))
  }))
  pool <- rbind(pool, pool2)
}
fwrite(pool, "02_results/deepen/spatial_scores_both_signatures.csv")
log_msg("pooled spots: ", nrow(pool))

analyse <- function(d, coh, sigcol, label) {
  s <- d[cohort == coh & complete.cases(d[, c("sig13", "sig17", "caf", "mye"), with = FALSE])]
  s <- s[, .(slice, sig = get(sigcol), caf, mye)]
  out <- rbindlist(lapply(unique(s$slice), function(sl) {
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
    d_obs <- r_sm - r_sc; d_bs <- bs[1, ] - bs[2, ]
    null <- replicate(2000, {
      sw <- sample(c(TRUE, FALSE), nrow(s1), replace = TRUE)
      mye2 <- s1$mye; caf2 <- s1$caf
      mye2[sw] <- s1$caf[sw]; caf2[sw] <- s1$mye[sw]
      cor(s1$sig, mye2, method = "spearman") - cor(s1$sig, caf2, method = "spearman")
    })
    p_perm <- (1 + sum(abs(null) >= abs(d_obs))) / (1 + 2000)
    pspear <- function(x, y, z) {
      rx <- rank(x); ry <- rank(y); rz <- rank(z)
      cor(lm(rx ~ rz)$residuals, lm(ry ~ rz)$residuals)
    }
    q <- quantile(s1$sig, 0.75); hi <- s1$sig >= q
    data.table(signature = label, cohort = coh, slice = sl, n_spots = nrow(s1),
               rho_sig_mye = r_sm, rho_sig_mye_lo = quantile(bs[1, ], 0.025),
               rho_sig_mye_hi = quantile(bs[1, ], 0.975),
               rho_sig_caf = r_sc, rho_sig_caf_lo = quantile(bs[2, ], 0.025),
               rho_sig_caf_hi = quantile(bs[2, ], 0.975),
               rho_caf_mye = r_mc, diff_mye_minus_caf = d_obs,
               diff_lo = quantile(d_bs, 0.025), diff_hi = quantile(d_bs, 0.975),
               p_paired_permutation = p_perm,
               partial_sig_mye_given_caf = pspear(s1$sig, s1$mye, s1$caf),
               partial_sig_caf_given_mye = pspear(s1$sig, s1$caf, s1$mye),
               hotspot_prop_mye = mean(s1$mye[hi] >= quantile(s1$mye, 0.75)),
               hotspot_prop_caf = mean(s1$caf[hi] >= quantile(s1$caf, 0.75)))
  }))
  out
}

res <- rbindlist(list(
  analyse(pool, "GSE267401", "sig13", "13-gene program (report claim)"),
  analyse(pool, "GSE267401", "sig17", "17-gene perturbation signature (persisted)"),
  analyse(pool, "GSE226997", "sig13", "13-gene program (report claim)"),
  analyse(pool, "GSE226997", "sig17", "17-gene perturbation signature (persisted)")
))
fwrite(res, "02_results/deepen/spatial_permutation_both_signatures.csv")
print(res[, .(signature, cohort, slice, n_spots, rho_sig_mye, rho_sig_caf, rho_caf_mye,
              partial_sig_mye_given_caf, partial_sig_caf_given_mye,
              hotspot_prop_mye, hotspot_prop_caf, p_paired_permutation)])

# pooled 对照
pooled <- rbindlist(lapply(c("GSE267401", "GSE226997"), function(ch) {
  s <- pool[cohort == ch]
  rbindlist(lapply(c("sig13", "sig17"), function(sc) {
    x <- s[[sc]]
    data.table(cohort = ch, signature = sc, n = nrow(s),
               pooled_rho_sig_mye = cor(x, s$mye, method = "spearman", use = "complete.obs"),
               pooled_rho_sig_caf = cor(x, s$caf, method = "spearman", use = "complete.obs"),
               pooled_rho_caf_mye = cor(s$caf, s$mye, method = "spearman", use = "complete.obs"))
  }))
}))
fwrite(pooled, "02_results/deepen/spatial_pooled_both_signatures.csv")
print(pooled)
log_msg("DONE deepen_08")
