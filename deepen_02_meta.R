#!/usr/bin/env Rscript
# Deepen-02: 多队列连续 Cox 的随机效应 meta（TCGA + GSE39582 + GSE17536 + GSE17537）
# 目的：把"外部只过 KM、判别掉到 0.53~0.58"改成可合并的效应量，并给异质性。
suppressPackageStartupMessages({
  library(data.table); library(survival); library(timeROC); library(metafor)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

cohort_stat <- function(tag, f, time_col = "os_time", event_col = "os_event",
                        marker_col = "riskScore", filter_fun = NULL,
                        time_div = 12) {   # TCGA 原始单位为日(÷365.25)，GEO 为月(÷12)
  d <- fread(f)
  if (!is.null(filter_fun)) d <- filter_fun(d)
  t <- d[[time_col]]; ev <- d[[event_col]]; mk <- d[[marker_col]]
  ok <- !is.na(t) & !is.na(ev) & !is.na(mk)
  t <- t[ok]; ev <- ev[ok]; mk <- mk[ok]
  cf <- coxph(Surv(t, ev) ~ mk); s <- summary(cf)
  grp <- ifelse(mk >= median(mk), "High", "Low")
  km <- survdiff(Surv(t, ev) ~ grp)
  p_km <- pchisq(km$chisq, df = 1, lower.tail = FALSE)
  aucs <- rep(NA_real_, 3)
  try({
    tmax <- max(t[ev == 1]) / time_div
    tt <- c(1, 3, 5); tt <- tt[tt < tmax]
    tr <- timeROC(T = t / time_div, delta = ev, marker = mk, cause = 1, times = tt)
    aucs[1:length(tt)] <- as.numeric(tr$AUC)
  }, silent = TRUE)
  data.table(cohort = tag, n = length(t), events = sum(ev),
             logHR = s$coefficients[1, 1], se = s$coefficients[1, 3],
             HR = s$coefficients[1, 2], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4],
             p_cont = s$coefficients[1, 5], km_p = p_km,
             auc_1y = aucs[1], auc_3y = aucs[2], auc_5y = aucs[3])
}

tumor_only <- function(d) {
  # 剔除含 "Non Tumoral" 的样本（GSE39582 的 17 例非肿瘤）
  cols <- names(d)[vapply(d, function(x) any(grepl("Non Tumoral", as.character(x), fixed = TRUE)), logical(1))]
  if (length(cols)) {
    hit <- Reduce(`|`, lapply(cols, function(cc) grepl("Non Tumoral", as.character(d[[cc]]), fixed = TRUE)))
    log_msg("GSE39582 drop Non Tumoral: ", sum(hit), " via cols: ", paste(cols, collapse = ","))
    d <- d[!hit]
  }
  d
}

tab <- rbindlist(list(
  cohort_stat("TCGA-COADREAD (train)", "02_results/model/tcga_riskscore.csv", time_div = 365.25),
  cohort_stat("GSE39582 (tumor)", "02_results/validation/gse39582_riskscore.csv", filter_fun = tumor_only),
  cohort_stat("GSE17536", "02_results/validation/gse17536_riskscore.csv"),
  cohort_stat("GSE17537", "02_results/deepen/gse17537_riskscore.csv")
))
fwrite(tab, "02_results/deepen/meta_cohort_table.csv")
print(tab)

# ---- 随机效应 meta：全部队列 + 仅外部队列 ----
meta_fit <- function(sub, tag) {
  m <- rma(yi = logHR, sei = se, data = sub, method = "REML", slab = cohort)
  pred <- predict(m)
  out <- data.table(analysis = tag, k = m$k,
                    hr = exp(m$b[1]), lo = exp(m$ci.lb), hi = exp(m$ci.ub),
                    p = m$pval, I2 = m$I2, tau2 = m$tau2, Q_p = m$QEp,
                    pred_lo = exp(pred$pi.lb), pred_hi = exp(pred$pi.ub))
  return(list(fit = m, out = out))
}
allm <- meta_fit(tab, "all cohorts (train + external)")
ext  <- meta_fit(tab[cohort != "TCGA-COADREAD (train)"], "external only")
res <- rbindlist(list(allm$out, ext$out))
fwrite(res, "02_results/deepen/meta_random_effects.csv")
print(res)

# ---- 逐队列留一 ----
loo <- rbindlist(lapply(seq_len(nrow(tab)), function(i) {
  m <- rma(yi = logHR, sei = se, data = tab[-i], method = "REML")
  data.table(left_out = tab$cohort[i], hr = exp(m$b[1]), lo = exp(m$ci.lb),
             hi = exp(m$ci.ub), p = m$pval, I2 = m$I2, tau2 = m$tau2)
}))
fwrite(loo, "02_results/deepen/meta_leave_one_out.csv")
print(loo)

# ---- 森林图数据 + 图 ----
png("03_figures/deepen_fig_meta_forest.png", width = 2200, height = 1300, res = 220)
par(mar = c(5, 16, 3, 3))
forest(allm$fit, xlab = "Hazard ratio (continuous 13-gene score, per unit)",
       slab = tab$cohort, atransf = exp, at = log(c(0.5, 1, 2, 4)),
       cex = 0.85, psize = 1.1)
text(-2.6, nrow(tab) + 2.4, "Cohort", pos = 4, cex = 0.9, font = 2)
text(log(4.6), nrow(tab) + 2.4, "HR [95% CI]", pos = 2, cex = 0.9, font = 2)
dev.off()
log_msg("wrote deepen_fig_meta_forest.png")
log_msg("DONE deepen_02")
