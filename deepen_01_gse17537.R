#!/usr/bin/env Rscript
# Deepen-01: 第三个外部队列 GSE17537 (GPL570, VMC 55 例, GSE17538 的另一个中心)
# 只投训炼得到的 13 基因线性预测子，不重估系数。系列矩阵已为 log2 尺度，不再跑 CEL/RMA。
suppressPackageStartupMessages({
  library(data.table); library(hgu133plus2.db); library(AnnotationDbi)
  library(survival); library(timeROC)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

beta_tab <- fread("02_results/model/model_genes_beta.csv")
beta <- setNames(beta_tab$beta, beta_tab$gene)
genes <- names(beta)
log_msg("model genes: ", paste(genes, collapse = ", "))

sm <- "00_data/geo_bulk/series/GSE17537_series_matrix.txt.gz"
lines <- readLines(gzfile(sm))
b <- grep("^!series_matrix_table_begin", lines)
e <- grep("^!series_matrix_table_end", lines)
stopifnot(length(b) == 1, length(e) == 1)
tbl <- read.delim(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"),
                  check.names = FALSE, quote = "\"")
probe <- tbl[[1]]
mat <- as.matrix(tbl[, -1, drop = FALSE])
rownames(mat) <- probe
log_msg("series matrix: ", nrow(mat), " probes x ", ncol(mat), " samples; value range ",
        sprintf("%.2f ~ %.2f", min(mat), max(mat)))

# ---- 临床：逐样本解析 key: value ----
acc <- lines[grep("^!Sample_geo_accession", lines)]
gsm <- trimws(gsub("\"", "", strsplit(acc, "\t")[[1]][-1]))
ch  <- lines[grep("^!Sample_characteristics_ch1", lines)]
clin_list <- lapply(seq_along(gsm), function(i) list())
field_map <- list()
for (ln in ch) {
  vals <- trimws(gsub("\"", "", strsplit(ln, "\t")[[1]][-1]))
  for (i in seq_along(vals)) {
    v <- vals[i]
    if (is.na(v) || v == "" || !grepl(":", v)) next
    k <- trimws(sub(":.*$", "", v))
    val <- trimws(sub("^[^:]*:", "", v))
    clin_list[[i]][[k]] <- val
    field_map[[k]] <- TRUE
  }
}
log_msg("clinical keys: ", paste(names(field_map), collapse = " | "))
getf <- function(i, key) {
  v <- clin_list[[i]][[key]]
  if (is.null(v)) NA_character_ else v
}
clin <- data.table(
  gsm = gsm,
  age = as.numeric(sapply(seq_along(gsm), getf, "age")),
  gender = sapply(seq_along(gsm), getf, "gender"),
  ajcc_stage = sapply(seq_along(gsm), getf, "ajcc_stage"),
  os_raw = sapply(seq_along(gsm), getf, "overall_event (death from any cause)"),
  os_time = as.numeric(sapply(seq_along(gsm), getf, "overall survival follow-up time")),
  dfs_raw = sapply(seq_along(gsm), getf, "dfs_event (disease free survival; cancer recurrence)"),
  dfs_time = as.numeric(sapply(seq_along(gsm), getf, "dfs_time"))
)
clin[, os_event := ifelse(tolower(os_raw) == "death", 1, ifelse(tolower(os_raw) == "no death", 0, NA))]
clin[, dfs_event := ifelse(grepl("recurrence", tolower(dfs_raw), fixed = TRUE) & !grepl("no", tolower(dfs_raw)), 1,
                           ifelse(tolower(dfs_raw) == "no recurrence", 0, NA))]
log_msg("OS events: ", sum(clin$os_event == 1, na.rm = TRUE), " / n=", sum(!is.na(clin$os_event)))
log_msg("DFS events: ", sum(clin$dfs_event == 1, na.rm = TRUE), " / n=", sum(!is.na(clin$dfs_event)))

# ---- 探针 -> 基因符号（max-mean） ----
sym <- AnnotationDbi::mapIds(hgu133plus2.db, rownames(mat), "SYMBOL", "PROBEID")
keep <- !is.na(sym)
mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
o <- order(rowMeans(mat), decreasing = TRUE)
mat <- mat[o, , drop = FALSE]; sym <- sym[o]
mat <- mat[!duplicated(sym), , drop = FALSE]
rownames(mat) <- sym[!duplicated(sym)]
log_msg("symbol matrix: ", nrow(mat), " genes x ", ncol(mat), " samples")

miss <- setdiff(genes, rownames(mat))
if (length(miss)) stop("missing model genes: ", paste(miss, collapse = ", "))
e13 <- mat[genes, , drop = FALSE]

# 两种尺度：(a) 原始 log2（与主法口径一致）(b) 队列内 z-score（平台尺度对照）
risk_raw <- as.numeric(t(e13) %*% beta)
z <- t(scale(t(e13)))                     # 每基因在队列内标准化
risk_z <- as.numeric(t(z) %*% beta)

dt <- copy(clin)
dt[, riskScore := risk_raw]
dt[, riskScore_z := risk_z]
dt[, group := ifelse(riskScore >= median(riskScore), "High", "Low")]
fwrite(dt, "02_results/deepen/gse17537_riskscore.csv")
log_msg("wrote gse17537_riskscore.csv; n=", nrow(dt))

# ---- 统计量 ----
stat_row <- function(tag, time, event, marker) {
  ok <- !is.na(time) & !is.na(event) & !is.na(marker)
  t <- time[ok]; ev <- event[ok]; mk <- marker[ok]
  grp <- ifelse(mk >= median(mk), "High", "Low")
  km <- survdiff(Surv(t, ev) ~ grp)
  p_km <- pchisq(km$chisq, df = 1, lower.tail = FALSE)
  cf <- coxph(Surv(t, ev) ~ mk)
  s <- summary(cf)
  hr <- s$coefficients[1, 2]; lo <- s$conf.int[1, 3]; hi <- s$conf.int[1, 4]
  p_cont <- s$coefficients[1, 5]
  # 校准斜率：把转移过来的线性预测子（riskScore）直接进 Cox，其系数即 van Houwelingen 斜率
  slope <- s$coefficients[1, 1]
  slope_ci <- c(s$coefficients[1, 1] - 1.959964 * s$coefficients[1, 3],
                s$coefficients[1, 1] + 1.959964 * s$coefficients[1, 3])
  aucs <- NA_real_
  try({
    tr <- timeROC(T = t / 12, delta = ev, marker = mk, cause = 1,
                  times = c(1, 3, 5), iid = TRUE)
    aucs <- tr$AUC
  }, silent = TRUE)
  data.table(config = tag, n = length(t), events = sum(ev),
             km_p = p_km, hr_cont = hr, hr_low = lo, hr_high = hi, p_cont = p_cont,
             calib_slope = slope, calib_low = slope_ci[1], calib_high = slope_ci[2],
             auc_1y = aucs[1], auc_3y = aucs[2], auc_5y = aucs[3])
}
res <- rbindlist(list(
  stat_row("GSE17537 OS (raw log2)", dt$os_time, dt$os_event, dt$riskScore),
  stat_row("GSE17537 OS (cohort z-score)", dt$os_time, dt$os_event, dt$riskScore_z),
  stat_row("GSE17537 DFS (companion)", dt$dfs_time, dt$dfs_event, dt$riskScore)
))
fwrite(res, "02_results/deepen/gse17537_stats.csv")
print(res)
log_msg("DONE deepen_01")
