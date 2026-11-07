#!/usr/bin/env Rscript
# Deepen-03: 训练队列判别度的乐观校正（bootstrap B=1000）
# 报告只给了 C-index 的 200 次 bootstrap 分位区间；这里补 1000 次的乐观校正 C-index 与 1/3/5 年 AUC。
suppressPackageStartupMessages({
  library(data.table); library(survival); library(timeROC)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
set.seed(20260924)

beta_tab <- fread("02_results/model/model_genes_beta.csv")
genes <- beta_tab$gene
clin <- fread("02_results/model/tcga_riskscore.csv")
clin <- clin[!is.na(os_time) & !is.na(os_event) & !is.na(riskScore)]

expr <- readRDS("02_results/deg/tcga_voom_expr.rds")
if (is.list(expr) && !is.null(expr$E)) expr <- expr$E
expr <- as.matrix(expr)
log_msg("expr: ", nrow(expr), " x ", ncol(expr))
miss <- setdiff(genes, rownames(expr))
if (length(miss)) stop("genes missing in expr: ", paste(miss, collapse = ", "))

# 对齐：风险表样本 x 表达矩阵列（表达矩阵列为 TCGA-xx-xxxx-01B 全长 barcode）
colnames(expr) <- substr(colnames(expr), 1, 12)
clin[, case_id := substr(as.character(case_id), 1, 12)]
common <- intersect(clin$case_id, colnames(expr))
log_msg("overlap samples: ", length(common))
e13 <- t(expr[genes, common, drop = FALSE])
d <- clin[match(common, case_id)]
d <- d[!is.na(os_time) & !is.na(os_event)]
e13 <- e13[rownames(e13) %in% d$case_id, , drop = FALSE]
d <- d[match(rownames(e13), case_id)]
stopifnot(nrow(d) == nrow(e13))
log_msg("analysis n=", nrow(d), " events=", sum(d$os_event))

t_mon <- d$os_time / 30.4375        # 报告口径：月
t_yr  <- d$os_time / 365.25
ev <- d$os_event

# lp 越大风险越高 => 生存时间越短，需要 reverse=TRUE 才是常规方向的 C 指数
c_index <- function(lp, time, event) {
  concordancefit(Surv(time, event), lp, reverse = TRUE)$concordance
}
auc_t <- function(marker, time, event, times = c(1, 3, 5)) {
  # 注意：timeROC 在 times 只给单个值时返回的 AUC 长度为 2（首位为 NA），
  # 因此这里一次性传入全部时点，再按 times 对齐取值。
  mx <- max(time[event == 1], na.rm = TRUE)
  out <- rep(NA_real_, length(times))
  tt <- times[is.finite(mx) & times < mx]
  if (!length(tt)) return(out)
  r <- timeROC(T = time, delta = event, marker = marker, cause = 1, times = tt, iid = FALSE)
  vals <- as.numeric(r$AUC); tn <- as.numeric(r$times)
  for (k in seq_along(tt)) {
    j <- which(abs(tn - tt[k]) < 1e-8)[1]
    i <- which(abs(times - tt[k]) < 1e-8)[1]
    if (!is.na(j) && !is.na(i)) out[i] <- vals[j]
  }
  out
}

fit_lp <- function(X, time, event) {
  df <- data.frame(as.data.frame(X)); df$.t <- time; df$.e <- event
  cf <- coxph(Surv(.t, .e) ~ ., data = df)
  list(lp = as.numeric(predict(cf, newdata = df, type = "lp")), cf = cf)
}

# ---- apparent ----
app <- fit_lp(e13, t_yr, ev)
app_C <- c_index(app$lp, t_yr, ev)
app_A <- auc_t(app$lp, t_yr, ev)
log_msg(sprintf("apparent C=%.3f  AUC1/3/5=%.3f/%.3f/%.3f", app_C, app_A[1], app_A[2], app_A[3]))

B <- 1000
opt <- matrix(NA_real_, nrow = B, ncol = 4)   # C, AUC1, AUC3, AUC5
for (b in seq_len(B)) {
  idx <- sample(seq_len(nrow(e13)), replace = TRUE)
  Xb <- e13[idx, , drop = FALSE]; tb <- t_yr[idx]; eb <- ev[idx]
  if (sum(eb) < 10) next
  fb <- try(fit_lp(Xb, tb, eb), silent = TRUE)
  if (inherits(fb, "try-error")) next
  lp_boot <- fb$lp
  dfb <- data.frame(as.data.frame(e13)); colnames(dfb) <- colnames(e13)
  lp_orig <- try(as.numeric(predict(fb$cf, newdata = dfb, type = "lp")), silent = TRUE)
  if (inherits(lp_orig, "try-error")) next
  cb <- c_index(lp_boot, tb, eb); co <- c_index(lp_orig, t_yr, ev)
  ab <- auc_t(lp_boot, tb, eb); ao <- auc_t(lp_orig, t_yr, ev)
  opt[b, ] <- c(cb - co, ab - ao)
  if (b %% 100 == 0) log_msg("bootstrap ", b, "/", B)
}
opt <- opt[rowSums(is.na(opt)) < 4, , drop = FALSE]
log_msg("usable bootstrap replicates: ", nrow(opt))
m_opt <- colMeans(opt, na.rm = TRUE)
corr <- c(app_C - m_opt[1], app_A - m_opt[2:4])

out <- data.table(
  metric = c("Harrell C", "AUC 1y", "AUC 3y", "AUC 5y"),
  apparent = c(app_C, app_A),
  optimism = m_opt,
  optimism_corrected = corr,
  opt_2.5pct = apply(opt, 2, quantile, 0.025, na.rm = TRUE),
  opt_97.5pct = apply(opt, 2, quantile, 0.975, na.rm = TRUE),
  n_replicates_used = colSums(!is.na(opt))
)
fwrite(out, "02_results/deepen/tcga_optimism_bootstrap.csv")
print(out)

png("03_figures/deepen_fig_optimism.png", width = 2000, height = 1200, res = 220)
op <- par(mfrow = c(1, 2))
hist(opt[, 1], breaks = 40, main = "Optimism: Harrell C", xlab = "bootstrap optimism", col = "#8fb3d9")
abline(v = m_opt[1], col = "red", lwd = 2)
hist(opt[, 4], breaks = 40, main = "Optimism: 5-year AUC", xlab = "bootstrap optimism", col = "#f2b880")
abline(v = m_opt[4], col = "red", lwd = 2)
par(op)
dev.off()
log_msg("DONE deepen_03")
