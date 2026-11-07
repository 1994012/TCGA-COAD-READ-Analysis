#!/usr/bin/env Rscript
# Deepen-15: 风险评分 vs CIT 分型的预测能力对比（Harrell C、3 年 time-dependent AUC、
# 3 年 NRI/IDI），回应局限四的补充分析清单第 (2)(3) 条。
suppressPackageStartupMessages({
  library(data.table); library(survival); library(timeROC)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

d <- fread("02_results/validation/gse39582_riskscore.csv")
cols <- names(d)[vapply(d, function(x) any(grepl("Non Tumoral", as.character(x), fixed = TRUE)), logical(1))]
if (length(cols)) {
  hit <- Reduce(`|`, lapply(cols, function(cc) grepl("Non Tumoral", as.character(d[[cc]]), fixed = TRUE)))
  d <- d[!hit]
}
d <- d[!is.na(os_time) & !is.na(os_event) & !is.na(riskScore) & cit.molecularsubtype != "N/A"]
d[, cit_f := factor(cit.molecularsubtype)]
age_col <- names(d)[grepl("age", names(d), ignore.case = TRUE)][1]
stage_col <- names(d)[grepl("tnm.stage", names(d), ignore.case = TRUE)][1]
d[, age_n := suppressWarnings(as.numeric(get(age_col)))]
d[, stage_n := factor(suppressWarnings(as.numeric(get(stage_col))))]
d <- d[!is.na(age_n) & !is.na(stage_n)]
log_msg("n=", nrow(d))

models <- list(
  riskScore        = "riskScore",
  CIT              = "cit_f",
  age_stage        = "age_n + stage_n",
  age_stage_CIT    = "age_n + stage_n + cit_f",
  age_stage_CIT_RS = "age_n + stage_n + cit_f + riskScore"
)

fit <- lapply(models, function(rhs) coxph(as.formula(paste("Surv(os_time, os_event) ~", rhs)), data = d, x = TRUE))
cidx <- sapply(fit, function(m) summary(m)$concordance["C"])

# 3 年 time-dependent AUC（连续 marker 用线性预测值）
tt <- 36  # months
aucs <- sapply(fit, function(m) {
  lp <- as.numeric(predict(m, type = "lp"))
  r <- timeROC(T = d$os_time, delta = d$os_event, marker = lp,
               cause = 1, weighting = "aalen", times = tt, ROC = TRUE)
  r$AUC[2]
})

# NRI / IDI（3 年，用两个嵌套模型的线性预测值；连续型 NRI，survival 设置）
base <- fit$age_stage_CIT; full <- fit$age_stage_CIT_RS
# 手工计算：用 coxph 的 baseline survival 把线性预测映射为 3 年发病概率
bh <- basehaz(base, centered = FALSE)
s0 <- exp(-bh$hazard)
s0_t <- approx(c(0, bh$time), c(1, s0), xout = tt, method = "constant")$y
lp0 <- predict(base, type = "lp"); lp1 <- predict(full, type = "lp")
# full 模型的基准生存
bh1 <- basehaz(full, centered = FALSE)
s0_t1 <- approx(c(0, bh1$time), c(1, exp(-bh1$hazard)), xout = tt, method = "constant")$y
risk0 <- 1 - s0_t^exp(lp0)      # 旧模型 3 年风险
risk1 <- 1 - s0_t1^exp(lp1)     # 新模型（加 riskScore）3 年风险
ev <- d$os_event == 1 & d$os_time <= tt   # 3 年内事件
nv <- d$os_time > tt                     # 3 年后仍无事件（含删失在 3 年后）
# 连续 NRI（Pencina 定义，生存数据用 case/control 近似）
nri_event <- mean(risk1[ev] > risk0[ev]) - mean(risk1[ev] <= risk0[ev])
nri_nonev <- mean(risk1[!ev & nv] <= risk0[!ev & nv]) - mean(risk1[!ev & nv] > risk0[!ev & nv])
nri <- nri_event + nri_nonev
idi <- mean(risk1[ev] - risk0[ev]) - mean(risk1[!ev & nv] - risk0[!ev & nv])

out <- data.table(
  model = names(models), cindex = cidx, auc_3y = aucs,
  nri_event = c(rep(NA_real_, 4), nri_event),
  nri_nonevent = c(rep(NA_real_, 4), nri_nonev),
  nri_total = c(rep(NA_real_, 4), nri),
  idi = c(rep(NA_real_, 4), idi))
fwrite(out, "02_results/deepen/cindex_auc_nri_cit.csv")
log_msg("C: ", paste(names(cidx), sprintf("%.3f", cidx), sep = "=", collapse = " "))
log_msg("AUC3y: ", paste(names(aucs), sprintf("%.3f", aucs), sep = "=", collapse = " "))
log_msg("NRI=", round(nri, 4), " IDI=", round(idi, 4))
