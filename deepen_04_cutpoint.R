#!/usr/bin/env Rscript
# Deepen-04: 切点敏感性 —— 中位 / 三分位(T1 vs T3) / maximally selected rank statistic
# 最优切点只作探索性诊断，不改主方法的中位切分口径。
suppressPackageStartupMessages({
  library(data.table); library(survival); library(timeROC); library(maxstat)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

km_p <- function(t, ev, grp) {
  ok <- !is.na(t) & !is.na(ev) & !is.na(grp)
  s <- survdiff(Surv(t[ok], ev[ok]) ~ grp[ok])
  pchisq(s$chisq, df = 1, lower.tail = FALSE)
}
cox_hr <- function(t, ev, x) {
  ok <- !is.na(t) & !is.na(ev) & !is.na(x)
  cf <- coxph(Surv(t[ok], ev[ok]) ~ x[ok]); s <- summary(cf)
  c(s$coefficients[1, 2], s$conf.int[1, 3], s$conf.int[1, 4], s$coefficients[1, 5])
}

cut_variants <- function(tag, t, ev, mk) {
  t_y <- t / 12
  # 1) 中位
  g_med <- ifelse(mk >= median(mk), "High", "Low")
  n_med <- table(g_med)
  # 2) 三分位 T1 vs T3
  q <- quantile(mk, c(1/3, 2/3), na.rm = TRUE)
  keep <- mk <= q[1] | mk >= q[2]
  g_t13 <- rep(NA_character_, length(mk)); g_t13[mk <= q[1]] <- "Low(T1)"; g_t13[mk >= q[2]] <- "High(T3)"
  # 3) 最优切点（探索性）：20%-80% 网格 + maxstat 校正 p
  cand <- quantile(mk, seq(0.2, 0.8, by = 0.01), na.rm = TRUE)
  ps <- vapply(cand, function(cc) km_p(t_y, ev, ifelse(mk >= cc, "High", "Low")), numeric(1))
  cut_opt <- cand[which.min(ps)]; p_opt <- min(ps)
  g_opt <- ifelse(mk >= cut_opt, "High", "Low")
  hr_opt <- cox_hr(t_y, ev, as.numeric(g_opt == "High"))
  best_p <- NA_real_
  try({
    dfm <- data.frame(t_y = t_y, ev = ev, mk = mk)
    ms <- maxstat.test(Surv(t_y, ev) ~ mk, data = dfm, smethod = "LogRank",
                       pmethod = "condMC", B = 9999, minprop = 0.2, maxprop = 0.8)
    best_p <- ms$p.value
  }, silent = TRUE)
  hr_med <- cox_hr(t_y, ev, as.numeric(g_med == "High"))
  hr_t13 <- cox_hr(t_y[keep], ev[keep], as.numeric(g_t13[keep] == "High(T3)"))
  rbindlist(list(
    data.table(cohort = tag, cut = "median", n = length(mk), events = sum(ev),
               n_high = as.integer(n_med[["High"]]), n_low = as.integer(n_med[["Low"]]),
               km_p = km_p(t_y, ev, g_med), hr = hr_med[1], lo = hr_med[2], hi = hr_med[3], p = hr_med[4],
               cutpoint = median(mk), p_maxstat_adjusted = NA_real_),
    data.table(cohort = tag, cut = "tertile T1 vs T3", n = sum(keep), events = sum(ev[keep]),
               n_high = sum(g_t13[keep] == "High(T3)"), n_low = sum(g_t13[keep] == "Low(T1)"),
               km_p = km_p(t_y[keep], ev[keep], g_t13[keep]), hr = hr_t13[1], lo = hr_t13[2],
               hi = hr_t13[3], p = hr_t13[4], cutpoint = q[2], p_maxstat_adjusted = NA_real_),
    data.table(cohort = tag, cut = "maxstat (exploratory)", n = length(mk), events = sum(ev),
               n_high = ifelse(is.na(cut_opt), NA_integer_, sum(mk >= cut_opt)),
               n_low = ifelse(is.na(cut_opt), NA_integer_, sum(mk < cut_opt)),
               km_p = p_opt, hr = hr_opt[1], lo = hr_opt[2], hi = hr_opt[3], p = hr_opt[4],
               cutpoint = cut_opt, p_maxstat_adjusted = best_p)
  ))
}

d39582 <- fread("02_results/validation/gse39582_riskscore.csv")
cols <- names(d39582)[vapply(d39582, function(x) any(grepl("Non Tumoral", as.character(x), fixed = TRUE)), logical(1))]
if (length(cols)) {
  hit <- Reduce(`|`, lapply(cols, function(cc) grepl("Non Tumoral", as.character(d39582[[cc]]), fixed = TRUE)))
  d39582 <- d39582[!hit]
}
d17536 <- fread("02_results/validation/gse17536_riskscore.csv")
d17537 <- fread("02_results/deepen/gse17537_riskscore.csv")
dtcga  <- fread("02_results/model/tcga_riskscore.csv")

res <- rbindlist(list(
  cut_variants("TCGA (train)", dtcga$os_time / 30.4375, dtcga$os_event, dtcga$riskScore),
  cut_variants("GSE39582 (tumor)", d39582$os_time, d39582$os_event, d39582$riskScore),
  cut_variants("GSE17536", d17536$os_time, d17536$os_event, d17536$riskScore),
  cut_variants("GSE17537", d17537$os_time, d17537$os_event, d17537$riskScore)
))
fwrite(res, "02_results/deepen/cutpoint_sensitivity.csv")
print(res)
log_msg("DONE deepen_04")
