#!/usr/bin/env Rscript
# Deepen-05: 多因素 Cox 补协变量 —— 报告只在 TCGA 上放了年龄+分期，
# 这里把外部队列的可用协变量（年龄、分期、性别、MMR、CIT、辅助化疗）逐个加进去，
# 看风险评分的 HR 是否被临床病理因素吃掉。
suppressPackageStartupMessages({
  library(data.table); library(survival)
})
setwd("/root/projects/crc_caf_myeloid")
dir.create("02_results/deepen", recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

num_or_na <- function(x) suppressWarnings(as.numeric(as.character(x)))

mv <- function(tag, d, form_rhs, time_col = "os_time", event_col = "os_event") {
  d <- copy(d)
  d <- d[!is.na(get(time_col)) & !is.na(get(event_col)) & !is.na(riskScore)]
  f_uni <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ riskScore"))
  f_mv  <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ riskScore + ", form_rhs))
  dd <- d[, all.vars(f_mv), with = FALSE]
  # 把因子列转成 factor，数值列保持数值
  for (cc in setdiff(names(dd), c(time_col, event_col, "riskScore"))) {
    if (!is.numeric(dd[[cc]])) dd[, (cc) := factor(as.character(get(cc)))]
  }
  cc_u <- try(coxph(f_uni, data = dd), silent = TRUE)
  cc_m <- try(coxph(f_mv, data = dd), silent = TRUE)
  if (inherits(cc_u, "try-error") || inherits(cc_m, "try-error")) {
    log_msg("skip ", tag, ": model failed"); return(NULL)
  }
  su <- summary(cc_u); sm <- summary(cc_m)
  data.table(cohort = tag, covariates = form_rhs,
             n = cc_m$n, events = sum(dd[[event_col]]),
             hr_uni = su$coefficients[1, 2], lo_uni = su$conf.int[1, 3], hi_uni = su$conf.int[1, 4],
             hr_adj = sm$coefficients["riskScore", 2], lo_adj = sm$conf.int["riskScore", 3],
             hi_adj = sm$conf.int["riskScore", 4], p_adj = sm$coefficients["riskScore", 5],
             delta_pct = 100 * (sm$coefficients["riskScore", 2] - su$coefficients[1, 2]) /
               (su$coefficients[1, 2] - 1))
}

out <- list()

# ---- TCGA ----
dt <- fread("02_results/model/tcga_riskscore.csv")
dt[, os_month := os_time / 30.4375]
dt[, stage_simple := ifelse(grepl("IV", stage), "IV",
                     ifelse(grepl("III", stage), "III",
                     ifelse(grepl("II", stage), "II",
                     ifelse(grepl("I", stage), "I", NA_character_))))]
out[[length(out) + 1]] <- mv("TCGA (train)", dt, "age + stage_simple", "os_month")
out[[length(out) + 1]] <- mv("TCGA (train)", dt, "age")

# ---- GSE39582（肿瘤样本） ----
d3 <- fread("02_results/validation/gse39582_riskscore.csv")
cols <- names(d3)[vapply(d3, function(x) any(grepl("Non Tumoral", as.character(x), fixed = TRUE)), logical(1))]
if (length(cols)) {
  hit <- Reduce(`|`, lapply(cols, function(cc) grepl("Non Tumoral", as.character(d3[[cc]]), fixed = TRUE)))
  d3 <- d3[!hit]
}
nm <- names(d3)
age_col <- nm[grepl("age", nm, ignore.case = TRUE)][1]
stage_col <- nm[grepl("tnm.stage", nm, ignore.case = TRUE)][1]
mmr_col <- nm[grepl("mmr", nm, ignore.case = TRUE)][1]
cit_col <- nm[grepl("cit", nm, ignore.case = TRUE)][1]
chemo_col <- nm[grepl("chemotherapy.adjuvant$", nm, ignore.case = TRUE)][1]
log_msg("GSE39582 cols: age=", age_col, " stage=", stage_col, " mmr=", mmr_col, " cit=", cit_col, " chemo=", chemo_col)
d3[, age_n := num_or_na(get(age_col))]
d3[, stage_n := factor(num_or_na(get(stage_col)))]
d3[, mmr_f := factor(as.character(get(mmr_col)))]
d3[, cit_f := factor(as.character(get(cit_col)))]
d3[, chemo_f := factor(as.character(get(chemo_col)))]
out[[length(out) + 1]] <- mv("GSE39582 (tumor)", d3, "age_n + stage_n")
out[[length(out) + 1]] <- mv("GSE39582 (tumor)", d3, "age_n + stage_n + mmr_f")
out[[length(out) + 1]] <- mv("GSE39582 (tumor)", d3, "age_n + stage_n + mmr_f + cit_f")
out[[length(out) + 1]] <- mv("GSE39582 (tumor)", d3, "age_n + stage_n + mmr_f + chemo_f")

# ---- GSE17536 ----
d1 <- fread("02_results/validation/gse17536_riskscore.csv")
d1[, age_n := num_or_na(age)]
d1[, stage_f := factor(as.character(ajcc_stage))]
d1[, gender_f := factor(as.character(gender))]
out[[length(out) + 1]] <- mv("GSE17536", d1, "age_n + stage_f")
out[[length(out) + 1]] <- mv("GSE17536", d1, "age_n + stage_f + gender_f")

# ---- GSE17537 ----
d7 <- fread("02_results/deepen/gse17537_riskscore.csv")
d7[, age_n := num_or_na(age)]
d7[, stage_f := factor(as.character(ajcc_stage))]
d7[, gender_f := factor(as.character(gender))]
out[[length(out) + 1]] <- mv("GSE17537", d7, "age_n + stage_f")
out[[length(out) + 1]] <- mv("GSE17537", d7, "age_n + stage_f + gender_f")

res <- rbindlist(Filter(Negate(is.null), out))
fwrite(res, "02_results/deepen/multivariate_covariates.csv")
print(res)
log_msg("DONE deepen_05")
