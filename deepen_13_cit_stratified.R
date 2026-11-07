#!/usr/bin/env Rscript
# Deepen-13: GSE39582 CIT 分层内部验证 —— 在每个 CIT 亚型内部检验风险评分是否仍有效，
# 并做 riskScore x CIT 交互检验，回应多因素分析中 CIT 共线（HR 衰减 56%）的问题。
suppressPackageStartupMessages({
  library(data.table); library(survival); library(ggplot2)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
dir.create("03_figures", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
theme_set(theme_classic(base_family = "sans"))

d3 <- fread("02_results/validation/gse39582_riskscore.csv")
cols <- names(d3)[vapply(d3, function(x) any(grepl("Non Tumoral", as.character(x), fixed = TRUE)), logical(1))]
if (length(cols)) {
  hit <- Reduce(`|`, lapply(cols, function(cc) grepl("Non Tumoral", as.character(d3[[cc]]), fixed = TRUE)))
  d3 <- d3[!hit]
}
d3 <- d3[!is.na(os_time) & !is.na(os_event) & !is.na(riskScore) & cit.molecularsubtype != "N/A"]
d3[, cit_f := factor(cit.molecularsubtype)]
d3[, groupH := factor(group, levels = c("Low", "High"))]
log_msg("n=", nrow(d3), " subtypes: ", paste(names(table(d3$cit_f)), collapse = ","))

# ---- 每亚型内部：连续 riskScore 的 Cox ----
res <- lapply(sort(levels(d3$cit_f)), function(s) {
  dd <- d3[cit_f == s]
  if (nrow(dd) < 20 || sum(dd$os_event) < 5) {
    return(data.table(subtype = s, n = nrow(dd), events = sum(dd$os_event),
                      hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_,
                      hr_high = NA_real_, lo_high = NA_real_, hi_high = NA_real_,
                      p_high = NA_real_, p_logrank = NA_real_))
  }
  cf <- summary(coxph(Surv(os_time, os_event) ~ riskScore, data = dd))
  cf2 <- try(summary(coxph(Surv(os_time, os_event) ~ groupH, data = dd)), silent = TRUE)
  lr <- try(survdiff(Surv(os_time, os_event) ~ groupH, data = dd), silent = TRUE)
  p_lr <- if (inherits(lr, "try-error")) NA_real_ else pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  data.table(subtype = s, n = nrow(dd), events = sum(dd$os_event),
             hr = cf$coefficients[1, 2], lo = cf$conf.int[1, 3], hi = cf$conf.int[1, 4],
             p = cf$coefficients[1, 5],
             hr_high = if (inherits(cf2, "try-error")) NA_real_ else cf2$coefficients[1, 2],
             lo_high = if (inherits(cf2, "try-error")) NA_real_ else cf2$conf.int[1, 3],
             hi_high = if (inherits(cf2, "try-error")) NA_real_ else cf2$conf.int[1, 4],
             p_high = if (inherits(cf2, "try-error")) NA_real_ else cf2$coefficients[1, 5],
             p_logrank = p_lr)
})
res <- rbindlist(res)

# ---- 交互检验：riskScore x CIT（因子）----
fit0 <- coxph(Surv(os_time, os_event) ~ riskScore + cit_f, data = d3)
fit1 <- coxph(Surv(os_time, os_event) ~ riskScore * cit_f, data = d3)
lrt <- anova(fit0, fit1, test = "LRT")
# 每单位 SD 的量级，便于横向比较
sd_rs <- sd(d3$riskScore)
res[, hr_sd := hr^sd_rs]

out <- data.table(table_name = "GSE39582 CIT stratified",
                  interaction_lrt_p = lrt$`Pr(>|Chi|)`[2],
                  riskScore_sd = sd_rs)
fwrite(res, "02_results/deepen/cit_stratified.csv")
fwrite(out, "02_results/deepen/cit_stratified_meta.csv")
log_msg("interaction LRT p = ", format(lrt$`Pr(>|Chi|)`[2], digits = 3))

# ---- 森林图 ----
lab <- res[, .(subtype, n, events, hr, lo, hi, p)]
lab[, lab := sprintf("%s (n=%d, %d ev)", subtype, n, events)]
overall <- data.table(lab = "GSE39582 全队列", hr = 1.20, lo = 1.05, hi = 1.37, p = 0.006)
pf <- rbindlist(list(lab[, .(lab, hr, lo, hi, p)], overall))
pf[, sig := p < 0.05]
pf[, lab := factor(lab, levels = rev(pf$lab))]
p <- ggplot(pf, aes(x = hr, y = lab)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi, color = sig), size = 0.55, linewidth = 0.7) +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey55"), guide = "none") +
  scale_x_log10() +
  labs(x = "HR（每单位风险评分，对数轴）", y = NULL,
       title = "GSE39582：按 CIT 分子亚型分层的风险评分效应") +
  coord_cartesian(xlim = c(0.4, 3.5))
ggsave("03_figures/FigD5_cit_stratified.png", p, width = 7.2, height = 4.2, dpi = 300)
log_msg("done")
