#!/usr/bin/env Rscript
# Deepen-07: 空间侧的零模型与偏相关
# 报告说"未做置换零模型，不作为旁分泌证据"。这里补三项：
#  (1) 每切片程序-髓系 / 程序-CAF 的 Spearman 及其 bootstrap CI
#  (2) 配对置换检验：同一 spot 上交换 CAF 与髓系分数，检验"程序更靠近谁"
#  (3) 偏 Spearman：ρ(程序,髓系 | CAF) 与 ρ(程序,CAF | 髓系)
suppressPackageStartupMessages({ library(data.table) })
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
set.seed(20260924)

pspear <- function(x, y, z) {
  # 秩次回归残差法：ρ(x,y | z)
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  ex <- lm(rx ~ rz)$residuals; ey <- lm(ry ~ rz)$residuals
  cor(ex, ey, method = "pearson")
}

analyse <- function(file, tag) {
  d <- fread(file)
  need <- c("sample", "sig", "caf", "mye")
  if (!all(need %in% names(d))) {
    log_msg("skip ", tag, ": columns ", paste(names(d), collapse = ",")); return(NULL)
  }
  out <- rbindlist(lapply(unique(d$sample), function(s) {
    s1 <- d[sample == s & complete.cases(d[, ..need])]
    if (nrow(s1) < 50) return(NULL)
    r_sm <- suppressWarnings(cor(s1$sig, s1$mye, method = "spearman"))
    r_sc <- suppressWarnings(cor(s1$sig, s1$caf, method = "spearman"))
    r_mc <- suppressWarnings(cor(s1$mye, s1$caf, method = "spearman"))
    # bootstrap CI
    B <- 1000
    bs <- replicate(B, {
      i <- sample(seq_len(nrow(s1)), replace = TRUE)
      c(cor(s1$sig[i], s1$mye[i], method = "spearman"),
        cor(s1$sig[i], s1$caf[i], method = "spearman"))
    })
    d_obs <- r_sm - r_sc
    d_bs <- bs[1, ] - bs[2, ]
    # 配对置换：交换 mye / caf
    Bp <- 2000
    null <- replicate(Bp, {
      sw <- sample(c(TRUE, FALSE), nrow(s1), replace = TRUE)
      mye2 <- s1$mye; caf2 <- s1$caf
      mye2[sw] <- s1$caf[sw]; caf2[sw] <- s1$mye[sw]
      cor(s1$sig, mye2, method = "spearman") - cor(s1$sig, caf2, method = "spearman")
    })
    p_perm <- (1 + sum(abs(null) >= abs(d_obs))) / (1 + Bp)
    # hotspot 重叠（切片内 75% 分位）
    q <- quantile(s1$sig, 0.75)
    hi <- s1$sig >= q
    qm <- quantile(s1$mye, 0.75); qc <- quantile(s1$caf, 0.75)
    prop_m <- mean(s1$mye[hi] >= qm); prop_c <- mean(s1$caf[hi] >= qc)
    data.table(cohort = tag, slice = s, n_spots = nrow(s1),
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
               hotspot_prop_mye = prop_m, hotspot_prop_caf = prop_c)
  }))
  out
}

res <- rbindlist(Filter(Negate(is.null), list(
  analyse("02_results/spatial/spot_signature_caf_myeloid.csv", "GSE267401"),
  analyse("02_results/spatial/gse226997_spot_signature_caf_myeloid.csv", "GSE226997")
)))
fwrite(res, "02_results/deepen/spatial_permutation_partial.csv")
print(res)

# ---- 图：分切片 ρ(程序,髓系) vs ρ(程序,CAF) ----
if (nrow(res)) {
  png("03_figures/deepen_fig_spatial_rho.png", width = 2200, height = 1300, res = 220)
  d <- res
  xpos <- seq_len(nrow(d))
  plot(xpos, d$rho_sig_mye, ylim = c(-0.3, 0.9), pch = 19, col = "#c0392b",
       xlab = "", ylab = "Spearman rho (13-gene program vs marker)", xaxt = "n",
       main = "Per-slice program correlation: myeloid vs CAF")
  arrows(xpos, d$rho_sig_mye_lo, xpos, d$rho_sig_mye_hi, code = 3, angle = 90, length = 0.05, col = "#c0392b")
  points(xpos + 0.15, d$rho_sig_caf, pch = 17, col = "#2471a3")
  arrows(xpos + 0.15, d$rho_sig_caf_lo, xpos + 0.15, d$rho_sig_caf_hi, code = 3, angle = 90, length = 0.05, col = "#2471a3")
  abline(h = 0, lty = 2)
  axis(1, at = xpos + 0.075, labels = paste0(d$cohort, "\n", d$slice), cex.axis = 0.7, las = 2)
  legend("topright", legend = c("program ~ myeloid", "program ~ CAF"),
         pch = c(19, 17), col = c("#c0392b", "#2471a3"), bty = "n")
  dev.off()
  log_msg("wrote deepen_fig_spatial_rho.png")
}
log_msg("DONE deepen_07")
