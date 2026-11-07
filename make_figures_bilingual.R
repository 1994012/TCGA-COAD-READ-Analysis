#!/usr/bin/env Rscript
# =============================================================================
# make_figures_bilingual.R —— 生成本报告全部图件的中文 / 英文两套
#
# 设计约束（与报告第 2 章「可复现性约定」一致）：
#   * 只读取 02_分析结果/ 下已落盘的结果表与风险评分表；
#   * 不重新估计任何统计量。唯一例外是 FigE1 / FigE2 的时间依赖 ROC，
#     这两张图在原流程中即由 timeROC 直接读取风险评分表绘制，属绘图环节而非重算；
#   * 随机数种子固定为 20240924（与原 deepen 系列一致）。
#
# 用法：
#   Rscript make_figures_bilingual.R zh 03_图件      # 中文：只写新增的 FigN*.png
#   Rscript make_figures_bilingual.R en 03_图件_en   # 英文：写全部 26 张
#
# 中文模式不覆盖 FigD*/FigE*：那 18 张是服务器上生成、已经过核对的定稿，
# 本机字体与服务器不同，重新渲染只会引入无意义的差异。
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
  library(survival); library(survminer); library(timeROC)
  library(patchwork); library(ggrepel); library(ggridges)
})
set.seed(20240924)

args   <- commandArgs(trailingOnly = TRUE)
LANG   <- if (length(args) >= 1) args[1] else "zh"
FIGDIR <- if (length(args) >= 2) args[2] else "03_图件"
stopifnot(LANG %in% c("zh", "en"))
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

BASE <- file.path(Sys.getenv("HOME"), "Documents", "first-year", "docx", "8.25",
                  "肝癌，肺癌，结直肠癌", "结直肠癌_CAF髓系")
D <- function(...) file.path(BASE, "02_分析结果", ...)
log_msg <- function(...)
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "")

# 中文字体：本机已装 Noto Sans CJK SC（OTF，systemfonts 可直接解析，不走 showtext）
base_family <- if (LANG == "zh") "Noto Sans CJK SC" else "Helvetica Neue"

theme_set(theme_bw(base_size = 11, base_family = base_family) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey92"),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  axis.text = element_text(colour = "black")))

save_png <- function(p, name, w = 7, h = 5) {
  fn <- file.path(FIGDIR, paste0(name, ".png"))
  ggsave(fn, p, width = w, height = h, dpi = 300, type = "cairo", bg = "white")
  log_msg("saved ", basename(fn))
}

pal_caf <- "#2c7bb6"; pal_mye <- "#d7191c"; pal_alt <- "#fdae61"
pal4 <- c("TCGA" = "#2C7BB6", "GSE39582" = "#D7191C",
          "GSE17536" = "#FDAE61", "GSE17537" = "#1A9641")

# =============================================================================
# 文案表
# =============================================================================
T <- if (LANG == "zh") list(
  train = "训练", pooled = "合并（随机效应）", pooled_ext = "合并（仅外部队列）",
  hr_axis = "每单位风险评分的风险比 HR (95% CI)", hr_axis2 = "高 vs 低危组 HR (95% CI)",
  d1 = "图 D1  多队列随机效应 meta 分析",
  d2a = "图 D2  GSE17537 外部验证（第三个 GEO 队列）", d2a_sub = "总生存 OS",
  d2b = "图 D2b  GSE17537 无病生存（探索性）", d2b_sub = "DFS / RFS",
  time_year = "时间（年）", surv = "总生存概率", riskgrp = "风险分组",
  lowN = "低危", highN = "高危", atrisk = "风险人数",
  d3 = "图 D3  切点敏感性分析", d3_sub = "三种切点定义下风险分层的稳定性",
  cut_med = "中位数", cut_ter = "三分位 T1 vs T3", cut_max = "maxstat 最优",
  d4 = "图 D4  空转逐切片：签名与髓系 / CAF 的相关",
  d4_sub = "灰线连接同一切片；两队列共 8 张切片",
  d4_rho = "13 基因签名 vs 细胞评分的 Spearman rho",
  mye = "髓系评分", caf = "CAF 评分",
  d5 = "图 D5  单细胞髓系优势的置换检验",
  d5_sub = "GSE132465 子集；零分布由 2000 次细胞标签重排得到",
  d5_x = "髓系 - T 细胞的平均差值（随机重排零分布）", dens = "密度",
  obs = "观测值 %.3f\nz = %.1f, p = %.4f",
  d6 = "图 D6  CAF -> 髓系 通讯的主要配体-受体对", d6_sub = "按贡献排序的前 15 对",
  d6_y = "CAF 配体平均表达 x 髓系受体平均表达",
  d7 = "图 D7  风险分组通路信号的跨队列一致性", d7_sub = "每个点为一个 Hallmark / KEGG 通路",
  nes_a = "TCGA-COADREAD  NES（高危 vs 低危）", nes_b = "GSE39582  NES（高危 vs 低危）",
  both = "两队列均显著", one = "单队列显著", ns = "不显著",
  e1 = "TCGA-COADREAD 训练集时间依赖 ROC", e1_leg = "随访时点",
  fpr = "假阳性率 (1 - 特异度)", tpr = "真阳性率 (敏感度)", auc_y = "%d 年 (AUC = %.3f)",
  e2 = "时间依赖 AUC 随访轨迹：训练集与三个外部队列",
  e2_x = "随访时间（年）", e2_y = "时间依赖 AUC", cohort_leg = "队列",
  e3 = "四个队列中风险评分的分布", e3_x = "队列内 z-score",
  e4 = "TCGA 高低危组差异表达火山图（805 个显著基因）",
  up_high = "高危组上调", up_low = "低危组上调", notsig = "不显著",
  e5 = "通路富集的跨队列热图", e5_sub = "按平均 |NES| 取前 24 条；* 校正 p<0.05，** <0.01，*** <0.001",
  e6 = "CAF → 髓系前 20 配体-受体对：强度与零模型证据",
  e6_sub = "虚线：z = 1.96（零模型 p ≈ 0.05）",
  e6_x = "配体-受体平均表达乘积（CAF × 髓系）", e6_y = "表达匹配零模型 z 值",
  z_lab = "零模型 z", contrib = "贡献占比 (%)",
  e7 = "全部 380 个 source → target 组合中强度最高的 20 个",
  e7_sub = "红色：CAF（Myofibroblasts / Stromal 1–3）作为信号来源",
  e7_x = "CellChat 通讯聚合强度", caf_src = "CAF 来源", oth_src = "其他来源",
  e8 = "空间转录组逐 spot：CAF 评分与髓系评分的关系",
  e8_sub = "颜色为 13 基因签名评分；虚线为线性趋势", sig13_lab = "13 基因\n签名评分",
  cit = "GSE39582：按 CIT 分子亚型分层的风险评分效应",
  cit_x = "HR（每单位风险评分，对数轴）", overall = "GSE39582 全队列",
  sub_lab = "%s (n=%d, %d 事件)",
  qc = "GSE226997 逐切片质量与构成复核（红：异例切片 P2 / P3）",
  qc_a = "spot 质量", qc_ay = "QC 通过率（%）",
  qc_b = "测序深度", qc_by = "中位 UMI / spot",
  qc_c = "肿瘤区域构成", qc_cy = "肿瘤高表达 spot 占比（%）"
) else list(
  train = "training", pooled = "Pooled (random effects)", pooled_ext = "Pooled (external only)",
  hr_axis = "Hazard ratio per unit risk score, HR (95% CI)",
  hr_axis2 = "Hazard ratio, high vs low risk (95% CI)",
  d1 = "Figure D1  Random-effects meta-analysis across cohorts",
  d2a = "Figure D2  External validation in GSE17537 (third GEO cohort)", d2a_sub = "Overall survival",
  d2b = "Figure D2b  Disease-free survival in GSE17537 (exploratory)", d2b_sub = "DFS / RFS",
  time_year = "Time (years)", surv = "Overall survival probability", riskgrp = "Risk group",
  lowN = "Low", highN = "High", atrisk = "No. at risk",
  d3 = "Figure D3  Cut-point sensitivity analysis",
  d3_sub = "Stability of risk stratification across three cut-point definitions",
  cut_med = "Median", cut_ter = "Tertile T1 vs T3", cut_max = "Optimal (maxstat)",
  d4 = "Figure D4  Per-slice correlation of the signature with myeloid / CAF scores",
  d4_sub = "Grey lines join the two estimates from the same slice; 8 slices from two cohorts",
  d4_rho = "Spearman rho, 13-gene signature vs cell-type score",
  mye = "Myeloid score", caf = "CAF score",
  d5 = "Figure D5  Permutation test for myeloid dominance in single cells",
  d5_sub = "GSE132465 subset; null distribution from 2,000 cell-label permutations",
  d5_x = "Mean difference, myeloids minus T cells (permutation null)", dens = "Density",
  obs = "Observed %.3f\nz = %.1f, p = %.4f",
  d6 = "Figure D6  Top ligand-receptor pairs in CAF-to-myeloid signalling",
  d6_sub = "Top 15 pairs ranked by contribution",
  d6_y = "mean CAF ligand expression x mean myeloid receptor expression",
  d7 = "Figure D7  Cross-cohort concordance of risk-group pathway signals",
  d7_sub = "Each point is one Hallmark / KEGG pathway",
  nes_a = "TCGA-COADREAD NES (high vs low risk)", nes_b = "GSE39582 NES (high vs low risk)",
  both = "Significant in both", one = "Significant in one", ns = "Not significant",
  e1 = "Time-dependent ROC in the TCGA-COADREAD training set", e1_leg = "Time horizon",
  fpr = "False positive rate (1 - specificity)", tpr = "True positive rate (sensitivity)",
  auc_y = "%d years (AUC = %.3f)",
  e2 = "Trajectory of time-dependent AUC: training set and three external cohorts",
  e2_x = "Follow-up time (years)", e2_y = "Time-dependent AUC", cohort_leg = "Cohort",
  e3 = "Distribution of the risk score in four cohorts", e3_x = "Within-cohort z-score",
  e4 = "Volcano plot of differential expression between risk groups in TCGA (805 genes)",
  up_high = "Up in high risk", up_low = "Up in low risk", notsig = "Not significant",
  e5 = "Heatmap of pathway enrichment across cohorts",
  e5_sub = "Top 24 pathways by mean |NES|; * adj. p<0.05, ** <0.01, *** <0.001",
  e6 = "Top 20 CAF-to-myeloid ligand-receptor pairs: strength and null-model evidence",
  e6_sub = "Dashed line: z = 1.96 (null-model p approximately 0.05)",
  e6_x = "Product of mean ligand-receptor expression (CAF x myeloid)",
  e6_y = "z value, expression-matched null model",
  z_lab = "Null-model z", contrib = "Contribution (%)",
  e7 = "Top 20 of all 380 source-to-target combinations by signalling strength",
  e7_sub = "Red: combinations with CAF (myofibroblasts / stromal 1-3) as the source",
  e7_x = "CellChat aggregate signalling strength", caf_src = "CAF source", oth_src = "Other source",
  e8 = "Per-spot relationship between CAF and myeloid scores in spatial transcriptomics",
  e8_sub = "Colour, 13-gene signature score; dashed line, linear trend",
  sig13_lab = "13-gene\nsignature score",
  cit = "GSE39582: effect of the risk score stratified by CIT molecular subtype",
  cit_x = "Hazard ratio per unit risk score (log scale)", overall = "Whole cohort",
  sub_lab = "%s (n=%d, %d events)",
  qc = "Per-slice quality and composition check in GSE226997 (red, outlier slices P2 / P3)",
  qc_a = "Spot quality", qc_ay = "QC pass rate (%)",
  qc_b = "Sequencing depth", qc_by = "Median UMI / spot",
  qc_c = "Tumour composition", qc_cy = "Tumour-high spots (%)"
)

N <- if (LANG == "zh") list(
  n1 = "图 N1  13 基因模型的 Cox 系数", n1_sub = "结果为模型落盘的系数点估计",
  n1_x = "多因素 Cox 系数 beta（每基因 z-score 单位）",
  beta_up = "推高评分", beta_dn = "压低评分",
  n2 = "图 N2  bootstrap 乐观度校正（1000 次重抽样）", n2_y = "取值",
  apparent = "表观值", corrected = "校正后", optband = "乐观度 95% 范围",
  n2_note = "标注为「表观值 / 校正后」；灰色横线为乐观度 95% 范围对应的校正值区间",
  n3 = "图 N3  校正协变量前后的风险评分效应", n3_x = "HR（对数轴，每单位风险评分）",
  unadj = "单因素", adj = "校正后", combo_lab = "队列与协变量组合",
  n4 = "图 N4  meta 分析的留一法敏感性", leave_out = "剔除队列",
  n4_sub = "菱形为显著（p<0.05）；灰色虚线为 HR = 1",
  n5 = "图 N5  相对临床病理指标与 CIT 分型的判别力增量",
  n5_y = "Harrell C 指数 / 3 年 AUC", model = "模型",
  rs_only = "仅风险评分", cit_only = "仅 CIT 分型", as_lab = "年龄 + 分期",
  as_cit = "年龄 + 分期 + CIT", as_cit_rs = "年龄 + 分期 + CIT + 风险评分",
  n6 = "图 N6  两种签名口径的逐切片共定位（同一批 spot）",
  n6_y = "rho（签名, CAF 评分）", sig13n = "13 基因程序", sig17n = "17 基因扰动签名",
  n7 = "图 N7  逐切片：边际相关与偏相关", n7_y = "rho（签名, 髓系评分）",
  marginal = "边际相关", partial = "偏相关（控制 CAF）",
  n7_sub = "同一切片的两种估计以竖线相连；三角为 GSE226997，圆为 GSE267401",
  n8 = "图 N8  单细胞置换检验的全部六个对比", n8_x = "z 值（观测差 / 零分布标准差）",
  prog13 = "13 基因程序", snai1 = "SNAI1 表达", n8_sub = "虚线为 z = ±1.96",
  c_my_t = "髓系 vs T 细胞", c_my_st = "髓系 vs 基质", c_st_t = "基质 vs T 细胞"
) else list(
  n1 = "Figure N1  Cox coefficients of the 13-gene model",
  n1_sub = "Persisted coefficient point estimates from the fitted model",
  n1_x = "Multivariable Cox coefficient beta (per gene z-score unit)",
  beta_up = "Raises the score", beta_dn = "Lowers the score",
  n2 = "Figure N2  Bootstrap optimism correction (1,000 replicates)", n2_y = "Value",
  apparent = "Apparent", corrected = "Corrected", optband = "95% range of optimism",
  n2_note = "Labels read apparent / corrected; grey bars span the corrected values implied by the 95% range of optimism",
  n3 = "Figure N3  Effect of the risk score before and after covariate adjustment",
  n3_x = "Hazard ratio (log scale, per unit risk score)",
  unadj = "Univariable", adj = "Adjusted", combo_lab = "Cohort and covariate set",
  n4 = "Figure N4  Leave-one-out sensitivity of the meta-analysis", leave_out = "Cohort left out",
  n4_sub = "Filled diamonds mark p<0.05; grey dashed line at HR = 1",
  n5 = "Figure N5  Incremental discrimination over clinicopathological variables and CIT subtype",
  n5_y = "Harrell C-index / 3-year AUC", model = "Model",
  rs_only = "Risk score only", cit_only = "CIT subtype only", as_lab = "Age + stage",
  as_cit = "Age + stage + CIT", as_cit_rs = "Age + stage + CIT + risk score",
  n6 = "Figure N6  Per-slice colocalisation under two signature definitions (same spots)",
  n6_y = "rho (signature, CAF score)", sig13n = "13-gene programme", sig17n = "17-gene perturbation signature",
  n7 = "Figure N7  Per-slice marginal and partial correlation", n7_y = "rho (signature, myeloid score)",
  marginal = "Marginal", partial = "Partial (adjusted for CAF)",
  n7_sub = "Vertical lines join the two estimates from the same slice; triangles denote GSE226997, circles GSE267401",
  n8 = "Figure N8  All six contrasts of the single-cell permutation test",
  n8_x = "z value (observed difference / null SD)",
  prog13 = "13-gene programme", snai1 = "SNAI1 expression", n8_sub = "Dashed lines at z = +/-1.96",
  c_my_t = "Myeloid vs T cells", c_my_st = "Myeloid vs stromal", c_st_t = "Stromal vs T cells"
)

cohort_short <- function(x) sub(" \\(tumor\\)", "", x)
tag_train <- function(x) {
  x <- sub(" \\(train\\)", "", cohort_short(x))
  fifelse(grepl("TCGA", x),
          if (LANG == "zh") paste0(x, "（", T$train, "）") else paste0(x, " (", T$train, ")"),
          x)
}

# 后面多处共用的两张表
sp   <- fread(D("deepen", "spatial_permutation_both_signatures.csv"))
sp13 <- sp[grepl("13-gene", signature)]
sp13[, slice_lab := paste0(substr(cohort, 1, 9), " ", gsub("Ajou_Visium_", "", slice))]

WRITE_LEGACY <- (LANG == "en")   # 仅英文模式重绘既有 18 张

# =============================================================================
# 既有 18 张
# =============================================================================
if (WRITE_LEGACY) {

  ## FigD1 meta 森林图
  ct <- fread(D("deepen", "meta_cohort_table.csv"))
  mr <- fread(D("deepen", "meta_random_effects.csv"))
  d1 <- rbind(
    ct[, .(cohort = tag_train(cohort), hr = HR, lo, hi, type = "cohort")],
    data.table(cohort = T$pooled, type = "pooled",
               hr = mr[analysis == "all cohorts (train + external)"]$hr,
               lo = mr[analysis == "all cohorts (train + external)"]$lo,
               hi = mr[analysis == "all cohorts (train + external)"]$hi),
    data.table(cohort = T$pooled_ext, type = "pooled",
               hr = mr[analysis == "external only"]$hr,
               lo = mr[analysis == "external only"]$lo,
               hi = mr[analysis == "external only"]$hi))
  d1[, cohort := factor(cohort, levels = rev(cohort))]
  save_png(ggplot(d1, aes(x = hr, y = cohort, colour = type)) +
             geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
             geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) + geom_point(size = 3) +
             scale_x_log10(breaks = c(0.5, 1, 1.5, 2, 3, 4), labels = number_format(accuracy = 0.01)) +
             scale_colour_manual(values = c(cohort = pal_caf, pooled = pal_mye), guide = "none") +
             labs(x = T$hr_axis, y = NULL, title = T$d1),
           "FigD1_meta_forest", 7.5, 4.2)

  ## FigD2a/b KM
  rs <- fread(D("deepen", "gse17537_riskscore.csv"))
  km_plot <- function(dat, tcol, ecol, title, subtitle) {
    d <- copy(dat)[, `:=`(.t = get(tcol) / 12, .e = get(ecol))]
    d <- d[!is.na(.t) & !is.na(.e) & .t > 0]
    d[, .g := fifelse(riskScore >= median(riskScore, na.rm = TRUE), T$highN, T$lowN)]
    fit <- survfit(Surv(.t, .e) ~ .g, data = d)
    pl <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                     palette = c(pal_mye, pal_caf), legend.title = T$riskgrp,
                     legend.labs = c(T$highN, T$lowN), risk.table.title = T$atrisk,
                     xlab = T$time_year, ylab = T$surv, title = title, subtitle = subtitle,
                     ggtheme = theme_bw(base_size = 11, base_family = base_family))
    pl$plot <- pl$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                               plot.subtitle = element_text(hjust = 0.5))
    pl$plot
  }
  tc <- intersect(names(rs), c("os_time", "os_delay", "os.delay (months)"))[1]
  save_png(km_plot(rs, tc, "os_event", T$d2a, T$d2a_sub), "FigD2a_gse17537_os", 6.5, 5)
  dfc <- intersect(names(rs), c("dfs_time", "dfs_delay", "rfs_delay", "rfs.delay"))[1]
  if (!is.na(dfc)) save_png(km_plot(rs, dfc, "dfs_event", T$d2b, T$d2b_sub),
                            "FigD2b_gse17537_dfs", 6.5, 5)

  ## FigD3 切点敏感性
  cp <- fread(D("deepen", "cutpoint_sensitivity.csv"))
  cp[, cohort2 := tag_train(cohort)]
  cp[, cut_lab := factor(cut, levels = c("median", "tertile T1 vs T3", "maxstat (exploratory)"),
                         labels = c(T$cut_med, T$cut_ter, T$cut_max))]
  save_png(ggplot(cp, aes(x = hr, y = interaction(cut_lab, cohort2), colour = cohort2)) +
             geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
             geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) + geom_point(size = 2.5) +
             scale_x_log10(breaks = c(0.5, 1, 2, 3, 5)) +
             labs(x = T$hr_axis2, y = NULL, title = T$d3, subtitle = T$d3_sub, colour = NULL),
           "FigD3_cutpoint_sensitivity", 7.5, 5)

  ## FigD4 空间逐切片 rho
  longd <- rbind(sp13[, .(slice_lab, cohort, rho = rho_sig_mye, target = T$mye)],
                 sp13[, .(slice_lab, cohort, rho = rho_sig_caf, target = T$caf)])
  save_png(ggplot(longd, aes(x = target, y = rho, fill = target, colour = cohort)) +
             geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
             geom_boxplot(alpha = 0.25, outlier.size = 0.8, width = 0.55) +
             geom_line(aes(group = slice_lab), colour = "grey40", linewidth = 0.35,
                       position = position_dodge(width = 0)) +
             geom_point(size = 2, position = position_jitter(width = 0.08, height = 0)) +
             scale_fill_manual(values = setNames(c(pal_mye, pal_caf), c(T$mye, T$caf)), guide = "none") +
             labs(x = NULL, y = T$d4_rho, title = T$d4, subtitle = T$d4_sub),
           "FigD4_spatial_slice_rho", 6.5, 5)

  ## FigD5 单细胞置换零分布
  pm <- fread(D("deepen", "scrna_permutation.csv"))
  r <- pm[statistic == "13-gene program consensus" & contrast == "Myeloids - T cells"][1]
  nullx <- seq(r$null_2.5 * 1.6, r$null_97.5 * 1.6, length.out = 400)
  dfn <- data.frame(x = nullx, y = dnorm(nullx, mean = r$null_mean, sd = r$null_sd))
  save_png(ggplot(dfn, aes(x = x, y = y)) + geom_area(fill = "grey70", alpha = 0.6) +
             geom_vline(xintercept = r$observed, colour = pal_mye, linewidth = 1.1) +
             geom_vline(xintercept = c(r$null_2.5, r$null_97.5), colour = "grey30", linetype = "dotted") +
             annotate("text", x = r$observed, y = max(dfn$y) * 0.75,
                      label = sprintf(T$obs, r$observed, r$z, r$p_permutation),
                      colour = pal_mye, hjust = -0.05, size = 3.6) +
             labs(x = T$d5_x, y = T$dens, title = T$d5, subtitle = T$d5_sub),
           "FigD5_scrna_permutation_null", 6.5, 4.2)

  ## FigD6 配体-受体
  lr <- fread(D("deepen", "cellchat_lr_top30.csv"))
  lr15 <- head(lr, 15)[, pair := paste0(ligand, " - ", receptor)]
  lr15[, pair := factor(pair, levels = rev(pair))]
  save_png(ggplot(lr15, aes(x = pair, y = product)) + geom_col(fill = pal_caf, width = 0.7) +
             coord_flip() + labs(x = NULL, y = T$d6_y, title = T$d6, subtitle = T$d6_sub),
           "FigD6_lr_pairs", 6.5, 5)

  ## FigD7 通路跨队列散点
  cc <- fread(D("deepen", "riskgroup_pathway_crosscohort.csv"))
  cc[, sig := fifelse(padj_TCGA < 0.05 & padj_GSE39582 < 0.05, T$both,
                      fifelse(padj_TCGA < 0.05 | padj_GSE39582 < 0.05, T$one, T$ns))]
  lab <- head(cc[grepl("HALLMARK", pathway) & padj_TCGA < 0.05][order(-abs(nes_TCGA))], 8)
  save_png(ggplot(cc, aes(x = nes_TCGA, y = nes_GSE39582, colour = sig)) +
             geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
             geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
             geom_point(size = 1.8, alpha = 0.8) +
             geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40") +
             geom_text_repel(aes(label = sub("HALLMARK_", "", pathway)), data = lab,
                             size = 2.8, max.overlaps = 20, segment.size = 0.2) +
             scale_colour_manual(values = setNames(c(pal_mye, pal_alt, "grey70"), c(T$both, T$one, T$ns))) +
             labs(x = T$nes_a, y = T$nes_b, colour = NULL, title = T$d7, subtitle = T$d7_sub) +
             theme(legend.position = "bottom"),
           "FigD7_pathway_crosscohort", 6.5, 5.5)

  ## FigE1-E3
  load_cohort <- function(tag, file, time_div) {
    d <- fread(file)
    tcol <- if ("os_time" %in% names(d)) "os_time" else "os.delay (months)"
    ecol <- if ("os_event" %in% names(d)) "os_event" else "os.event"
    d <- d[!is.na(get(tcol)) & !is.na(get(ecol)) & !is.na(riskScore)]
    data.table(cohort = tag, t = as.numeric(d[[tcol]]) / time_div,
               e = as.numeric(d[[ecol]]), mk = as.numeric(d$riskScore))
  }
  tcga <- load_cohort("TCGA", D("model", "tcga_riskscore.csv"), 365.25)
  coh <- rbind(tcga,
               load_cohort("GSE39582", D("validation", "gse39582_riskscore.csv"), 12),
               load_cohort("GSE17536", D("validation", "gse17536_riskscore.csv"), 12),
               load_cohort("GSE17537", D("deepen", "gse17537_riskscore.csv"), 12))

  mx <- max(tcga$t[tcga$e == 1]); tt <- c(1, 3, 5)[c(1, 3, 5) < mx]
  ro <- timeROC(T = tcga$t, delta = tcga$e, marker = tcga$mk, cause = 1, times = tt, iid = FALSE)
  roc_dt <- rbindlist(lapply(seq_along(tt), function(i) {
    j <- which(abs(as.numeric(ro$times) - tt[i]) < 1e-8)[1]
    if (is.na(j)) return(NULL)
    data.table(year = tt[i], fpr = as.numeric(ro$FP[, j]),
               tpr = as.numeric(ro$TP[, j]), auc = as.numeric(ro$AUC)[j])}))
  save_png(ggplot(roc_dt, aes(fpr, tpr, color = factor(year))) + geom_path(linewidth = 0.9) +
             geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
             scale_color_manual(values = c("1" = "#2C7BB6", "3" = "#D7191C", "5" = "#1A9641"),
                                labels = sapply(split(roc_dt, roc_dt$year),
                                                function(x) sprintf(T$auc_y, unique(x$year), unique(x$auc))),
                                name = T$e1_leg) +
             coord_equal() + lims(x = c(0, 1), y = c(0, 1)) +
             labs(x = T$fpr, y = T$tpr, title = T$e1),
           "FigE1_roc_tcga", 6.2, 5.6)

  auc_at <- function(t, e, mk, times) {
    ok <- !is.na(t) & !is.na(e) & !is.na(mk); t <- t[ok]; e <- e[ok]; mk <- mk[ok]
    ttk <- times[times < max(t[e == 1])]
    if (!length(ttk)) return(rep(NA_real_, length(times)))
    rr <- timeROC(T = t, delta = e, marker = mk, cause = 1, times = ttk, iid = FALSE)
    out <- rep(NA_real_, length(times))
    for (k in seq_along(ttk)) {
      j <- which(abs(as.numeric(rr$times) - ttk[k]) < 1e-8)[1]
      if (!is.na(j)) out[which(abs(times - ttk[k]) < 1e-8)] <- as.numeric(rr$AUC)[j]
    }
    out
  }
  grid <- seq(0.5, 10, by = 0.25)
  auc_dt <- rbindlist(lapply(split(coh, coh$cohort), function(d) {
    tg <- grid[grid < max(d$t[d$e == 1]) * 0.999]
    data.table(cohort = unique(d$cohort), year = tg, auc = auc_at(d$t, d$e, d$mk, tg))}))
  save_png(ggplot(auc_dt[!is.na(auc)], aes(year, auc, color = cohort)) +
             geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60") +
             geom_path(linewidth = 0.85, alpha = 0.9) +
             geom_point(data = auc_dt[year %in% c(1, 3, 5) & !is.na(auc)], size = 1.8) +
             scale_color_manual(values = pal4) +
             scale_y_continuous(limits = c(0.4, 0.85), breaks = seq(0.4, 0.8, 0.1)) +
             labs(x = T$e2_x, y = T$e2_y, color = T$cohort_leg, title = T$e2),
           "FigE2_auc_vs_time", 7.2, 4.6)

  coh[, mk_z := (mk - mean(mk)) / sd(mk), by = cohort]
  save_png(ggplot(coh, aes(x = mk_z, y = reorder(cohort, mk_z, median), fill = cohort)) +
             geom_density_ridges(alpha = 0.75, scale = 1.1, rel_min_height = 0.01) +
             scale_fill_manual(values = pal4) +
             labs(x = T$e3_x, y = NULL, title = T$e3) + theme(legend.position = "none"),
           "FigE3_riskscore_ridge", 6.4, 3.8)

  ## FigE4 火山图
  de <- fread(D("deepen", "riskgroup_de_tcga.csv"))
  de[, sig := fcase(adj_p < 0.05 & logFC >  log2(1.5), T$up_high,
                    adj_p < 0.05 & logFC < -log2(1.5), T$up_low, default = T$notsig)]
  lab <- unique(rbind(de[sig == T$up_high][order(adj_p)][1:6],
                      de[sig == T$up_low][order(adj_p)][1:6],
                      de[gene %in% c("SNAI1","TIMP1","MMP1","IL6","SPP1","VEGFA","FSTL3","COL1A1")][!is.na(gene)]))
  save_png(ggplot(de, aes(logFC, -log10(pmax(adj_p, 1e-300)), color = sig)) +
             geom_point(size = 0.7, alpha = 0.55) +
             geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey55") +
             geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55") +
             geom_point(data = lab, size = 1.1, color = "black") +
             geom_text_repel(data = lab, aes(label = gene), size = 2.6, max.overlaps = 20,
                             segment.size = 0.25, color = "black") +
             scale_color_manual(values = setNames(c("#D7191C", "#2C7BB6", "grey75"),
                                                c(T$up_high, T$up_low, T$notsig)), name = NULL) +
             labs(x = expression(log[2]~"fold change (high / low risk)"),
                  y = expression(-log[10]~"(adjusted "*italic(p)*")"),
                  title = T$e4) +
             theme(legend.position = c(0.08, 0.9), legend.background = element_blank()),
           "FigE4_volcano_de", 7.0, 5.4)

  ## FigE5 通路热图
  cr <- fread(D("deepen", "riskgroup_pathway_crosscohort.csv"))
  cr[, nes_mean := (abs(nes_TCGA) + abs(nes_GSE39582)) / 2]
  top <- cr[order(-nes_mean)][1:24]
  top[, pw := gsub("^(HALLMARK_|KEGG_)", "", pathway)]
  top[, pw := factor(pw, levels = rev(pw))]
  mm <- rbindlist(list(data.table(pw = top$pw, cohort = "TCGA", NES = top$nes_TCGA, padj = top$padj_TCGA),
                       data.table(pw = top$pw, cohort = "GSE39582", NES = top$nes_GSE39582, padj = top$padj_GSE39582)))
  mm[, star := fcase(padj < 0.001, "***", padj < 0.01, "**", padj < 0.05, "*", default = "")]
  save_png(ggplot(mm, aes(cohort, pw, fill = NES)) + geom_tile(color = "white", linewidth = 0.4) +
             geom_text(aes(label = star), size = 3, vjust = 0.75) +
             scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                                  midpoint = 0, limits = c(-3, 3), oob = squish, name = "NES") +
             labs(x = NULL, y = NULL, title = T$e5, subtitle = T$e5_sub) +
             theme(plot.title = element_text(size = 12), axis.text.y = element_text(size = 7.5)),
           "FigE5_pathway_heatmap", 6.2, 7.8)

  ## FigE6 气泡图
  lr20 <- copy(lr)[, pair := paste0(ligand, " \u2192 ", receptor)][order(-product)][1:20]
  save_png(ggplot(lr20, aes(product, z_pair)) +
             geom_hline(yintercept = 1.96, linetype = "dashed", color = "grey55") +
             geom_point(aes(size = contribution_pct, fill = z_pair), shape = 21,
                        color = "grey25", alpha = 0.9) +
             geom_text_repel(aes(label = pair), size = 2.5, max.overlaps = 16, segment.size = 0.25) +
             scale_fill_gradient(low = "#2C7BB6", high = "#D7191C", name = T$z_lab) +
             scale_size_continuous(range = c(2, 9), name = T$contrib) +
             labs(x = T$e6_x, y = T$e6_y, title = T$e6, subtitle = T$e6_sub) +
             theme(legend.box = "vertical"),
           "FigE6_lr_bubble", 7.4, 5.6)

  ## FigE7 棒棒糖
  rk <- fread(D("deepen", "cellchat_rank_specificity.csv"))
  rk[, is_caf := sub(" ->.*", "", pair) %in% c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3")]
  top20 <- rk[order(-agg_sum)][1:20]
  top20[, pair_f := factor(pair, levels = rev(top20$pair))]
  save_png(ggplot(top20, aes(agg_sum, pair_f)) +
             geom_segment(aes(xend = 0, yend = pair_f, color = is_caf), linewidth = 0.8) +
             geom_point(aes(color = is_caf, fill = is_caf), shape = 21, size = 2.6) +
             scale_color_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "grey55"),
                                labels = c(`TRUE` = T$caf_src, `FALSE` = T$oth_src), name = NULL) +
             scale_fill_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "grey55"), guide = "none") +
             labs(x = T$e7_x, y = NULL, title = T$e7, subtitle = T$e7_sub) +
             theme(axis.text.y = element_text(size = 7.5), legend.position = c(0.8, 0.15)),
           "FigE7_cellchat_lollipop", 7.2, 5.2)

  ## FigE8 空间散点
  sps <- fread(D("deepen", "spatial_scores_both_signatures.csv"))
  save_png(ggplot(sps, aes(caf, mye, color = sig13)) + geom_point(size = 0.55, alpha = 0.5) +
             geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5, linetype = "dashed") +
             scale_color_viridis_c(option = "mako", name = T$sig13_lab) +
             facet_wrap(~ cohort, scales = "free") +
             labs(x = T$caf, y = T$mye, title = T$e8, subtitle = T$e8_sub) +
             theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7)),
           "FigE8_spatial_spot", 8.6, 4.4)

  ## FigD5(cit)
  cs <- fread(D("deepen", "cit_stratified.csv"))
  cs[, lab := sprintf(T$sub_lab, subtype, n, events)]
  pf <- rbind(cs[, .(lab, hr, lo, hi, p)],
              data.table(lab = T$overall, hr = 1.20, lo = 1.05, hi = 1.37, p = 0.0065))
  pf[, sig := p < 0.05]
  pf[, lab := factor(lab, levels = rev(pf$lab))]
  save_png(ggplot(pf, aes(x = hr, y = lab)) +
             geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
             geom_pointrange(aes(xmin = lo, xmax = hi, color = sig), size = 0.55, linewidth = 0.7) +
             scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey55"), guide = "none") +
             scale_x_log10() + coord_cartesian(xlim = c(0.4, 3.5)) +
             labs(x = T$cit_x, y = NULL, title = T$cit),
           "FigD5_cit_stratified", 7.2, 4.2)

  ## FigD6 质量三联图
  qc2 <- fread(D("deepen", "spatial_quality_check.csv"))[cohort == "GSE226997"]
  qc2[, slabel := factor(gsub("Ajou_Visium_", "", sample), levels = gsub("Ajou_Visium_", "", sample))]
  fv <- c(P1 = "#4393C3", P2 = "#D6604D", P3 = "#D6604D", P4 = "#4393C3")
  bar <- function(yvar, ylab, ttl, digits) {
    ggplot(qc2, aes(x = slabel, y = get(yvar), fill = slabel)) +
      geom_col(width = 0.62, show.legend = FALSE) +
      geom_text(aes(label = if (digits == 0) sprintf("%d", round(get(yvar)))
                             else sprintf("%.1f%%", get(yvar))), vjust = -0.35, size = 3.2) +
      scale_fill_manual(values = fv) + labs(x = NULL, y = ylab, title = ttl)
  }
  pp <- (bar("qc_pass_pct", T$qc_ay, T$qc_a, 1) |
           bar("median_umi", T$qc_by, T$qc_b, 0) |
           bar("tumor_spot_pct", T$qc_cy, T$qc_c, 1)) +
    plot_annotation(title = T$qc, theme = theme(plot.title = element_text(size = 12, face = "bold")))
  ggsave(file.path(FIGDIR, "FigD6_spatial_quality.png"), pp, width = 9.6, height = 3.6,
         dpi = 300, type = "cairo", bg = "white")
  log_msg("saved FigD6_spatial_quality.png")
}

# =============================================================================
# 新增深化图件 FigN1 - FigN8（两种语言都画）
# =============================================================================

## N1 13 基因模型系数
mb <- fread(D("model", "model_genes_beta.csv"))
mb[, direction := fifelse(beta > 0, N$beta_up, N$beta_dn)]
mb[, gene := factor(gene, levels = gene[order(beta)])]
save_png(ggplot(mb, aes(x = gene, y = beta, fill = direction)) +
           geom_hline(yintercept = 0, colour = "grey40") + geom_col(width = 0.68) +
           coord_flip() +
           geom_text(aes(label = sprintf("%+.3f", beta), y = beta),
                     hjust = fifelse(mb$beta > 0, -0.15, 1.15), size = 2.6) +
           scale_fill_manual(values = setNames(c(pal_mye, pal_caf), c(N$beta_up, N$beta_dn)), name = NULL) +
           scale_y_continuous(limits = c(-0.28, 0.34)) +
           labs(x = NULL, y = N$n1_x, title = N$n1, subtitle = N$n1_sub) +
           theme(legend.position = "bottom"),
         "FigN1_model_coefficients", 6.6, 4.8)

## N2 bootstrap 乐观度校正
op <- fread(D("deepen", "tcga_optimism_bootstrap.csv"))
op[, metric_ord := .I]
op[, `:=`(lo_c = apparent - opt_97.5pct, hi_c = apparent - opt_2.5pct)]
opm <- op[, metric := factor(metric, levels = metric)]
p <- ggplot(op) +
  geom_errorbar(aes(x = metric, ymin = lo_c, ymax = hi_c), colour = "grey60",
                width = 0.2, linewidth = 0.9) +
  geom_point(aes(x = metric, y = apparent), shape = 21, size = 3.6,
             fill = pal_alt, colour = "grey25") +
  geom_point(aes(x = metric, y = optimism_corrected), shape = 21, size = 3.6,
             fill = pal_caf, colour = "grey25") +
  geom_segment(aes(x = metric_ord, xend = metric_ord, y = apparent, yend = optimism_corrected),
               data = op[, metric_ord := as.numeric(metric)], inherit.aes = FALSE,
               arrow = arrow(length = unit(0.12, "cm")), colour = "grey35", linewidth = 0.35)
p <- p + geom_point(data = op, aes(x = metric, y = apparent), shape = 21, size = 3.6,
                    fill = pal_alt, colour = "grey25")
ann <- op[, paste0(sprintf("%.3f", apparent), " / ", sprintf("%.3f", optimism_corrected))]
p <- p + annotate("text", x = op$metric, y = op$hi_c + 0.012, label = ann, size = 2.6) +
  scale_y_continuous(limits = c(0.60, 0.80), breaks = seq(0.60, 0.80, 0.05)) +
  labs(x = NULL, y = N$n2_y, title = N$n2, caption = N$n2_note)
save_png(p, "FigN2_optimism_correction", 7.0, 4.4)

## N3 多因素 Cox 森林图
mv <- fread(D("deepen", "multivariate_covariates.csv"))
mv[, combo := paste(cohort_short(cohort), covariates, sep = " | ")]
mv[, combo := factor(combo, levels = rev(combo))]
mvl <- rbindlist(list(
  mv[, .(combo, hr = hr_uni, lo = lo_uni, hi = hi_uni, what = N$unadj)],
  mv[, .(combo, hr = hr_adj, lo = lo_adj, hi = hi_adj, what = N$adj)]))
mvl[, what := factor(what, levels = c(N$unadj, N$adj))]
save_png(ggplot(mvl, aes(x = hr, y = combo, colour = what)) +
           geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
           geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2,
                          position = position_dodge(width = 0.6)) +
           geom_point(size = 2.1, position = position_dodge(width = 0.6)) +
           scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
           scale_colour_manual(values = setNames(c("grey55", pal_mye), c(N$unadj, N$adj)), name = NULL) +
           labs(x = N$n3_x, y = N$combo_lab, title = N$n3) +
           theme(legend.position = "bottom", axis.text.y = element_text(size = 6.6)),
         "FigN3_multivariate_forest", 7.8, 5.4)

## N4 meta 留一法
loo <- fread(D("deepen", "meta_leave_one_out.csv"))
loo[, y := tag_train(left_out)]
loo[, y := factor(y, levels = rev(y))]
loo[, sig := p < 0.05]
save_png(ggplot(loo, aes(x = hr, y = y)) +
           geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
           geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, colour = "grey45") +
           geom_point(aes(fill = sig), shape = 23, size = 3.8, colour = "grey25") +
           geom_text(aes(label = sprintf("%.2f [%.2f, %.2f]", hr, lo, hi)),
                     vjust = -1.3, size = 2.7, colour = "grey25") +
           scale_x_log10(breaks = c(0.8, 1, 1.2, 1.5, 2, 2.5), limits = c(0.85, 2.6)) +
           scale_fill_manual(values = c(`TRUE` = pal_mye, `FALSE` = "#ffffff"), guide = "none") +
           labs(x = T$hr_axis, y = N$leave_out, title = N$n4, subtitle = N$n4_sub),
         "FigN4_loo_sensitivity", 7.4, 3.9)

## N5 判别力增量
ci <- fread(D("deepen", "cindex_auc_nri_cit.csv"))
cmap <- c(riskScore = N$rs_only, CIT = N$cit_only, age_stage = N$as_lab,
          age_stage_CIT = N$as_cit, age_stage_CIT_RS = N$as_cit_rs)
ci[, lab := factor(cmap[model], levels = cmap)]
ci_l <- rbindlist(list(ci[, .(lab, value = cindex, what = "C-index")],
                       ci[, .(lab, value = auc_3y, what = "AUC-3y")]))
save_png(ggplot(ci_l, aes(x = lab, y = value, fill = what)) +
           geom_col(position = position_dodge(width = 0.72), width = 0.6) +
           geom_text(aes(label = sprintf("%.3f", value)), size = 2.6, vjust = -0.4,
                     position = position_dodge(width = 0.72)) +
           scale_fill_manual(values = c("C-index" = pal_caf, "AUC-3y" = pal_alt), name = NULL) +
           scale_y_continuous(limits = c(0, 0.80), expand = expansion(mult = c(0, 0.03))) +
           labs(x = N$model, y = N$n5_y, title = N$n5) +
           theme(axis.text.x = element_text(size = 7.2, angle = 16, hjust = 1),
                 legend.position = "bottom"),
         "FigN5_discrimination_increment", 7.4, 4.5)

## N6 两种签名口径对比
spp <- copy(sp)
spp[, sig_lab := fifelse(grepl("13-gene", signature), N$sig13n, N$sig17n)]
spp[, sl := paste0(substr(cohort, 1, 9), " ", gsub("Ajou_Visium_", "", slice))]
ord <- spp[sig_lab == N$sig13n][order(rho_sig_caf)]$sl
spp[, sl := factor(sl, levels = ord)]
save_png(ggplot(spp, aes(x = sl, y = rho_sig_caf, fill = sig_lab)) +
           geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
           geom_col(position = position_dodge(width = 0.78), width = 0.62) +
           geom_text(aes(label = sprintf("%.2f", rho_sig_caf)), size = 2.2, vjust = -0.35,
                     position = position_dodge(width = 0.78)) +
           scale_fill_manual(values = setNames(c(pal_caf, pal_alt), c(N$sig13n, N$sig17n)), name = NULL) +
           scale_y_continuous(limits = c(-0.1, 1.0)) +
           labs(x = NULL, y = N$n6_y, title = N$n6) +
           theme(axis.text.x = element_text(size = 6.8, angle = 28, hjust = 1),
                 legend.position = "bottom"),
         "FigN6_spatial_signature_contrast", 7.8, 4.4)

## N7 边际 vs 偏相关
longp <- rbind(
  sp13[, .(slice_lab, cohort, rho = rho_sig_mye, what = N$marginal)],
  sp13[, .(slice_lab, cohort, rho = partial_sig_mye_given_caf, what = N$partial)])
longp[, what := factor(what, levels = c(N$marginal, N$partial))]
ord7 <- sp13[order(rho_sig_mye)]$slice_lab
longp[, slice_lab := factor(slice_lab, levels = ord7)]
save_png(ggplot(longp, aes(x = slice_lab, y = rho, group = slice_lab)) +
           geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
           geom_line(colour = "grey45", linewidth = 0.4, position = position_dodge(width = 0.55)) +
           geom_point(aes(colour = what, shape = cohort), size = 2.6, fill = "white",
                      position = position_dodge(width = 0.55)) +
           scale_colour_manual(values = setNames(c(pal_caf, pal_mye), c(N$marginal, N$partial)), name = NULL) +
           scale_shape_manual(values = c("GSE267401" = 21, "GSE226997" = 24), name = NULL) +
           scale_y_continuous(limits = c(-0.08, 0.6)) +
           labs(x = NULL, y = N$n7_y, title = N$n7, subtitle = N$n7_sub) +
           theme(axis.text.x = element_text(size = 6.8, angle = 28, hjust = 1),
                 legend.position = "bottom", legend.box = "horizontal"),
         "FigN7_spatial_partial_correlation", 7.6, 4.8)

## N8 单细胞置换检验六对比
pm2 <- fread(D("deepen", "scrna_permutation.csv"))
cmap2 <- c("Myeloids - T cells" = N$c_my_t, "Myeloids - Stromal cells" = N$c_my_st,
           "Stromal cells - T cells" = N$c_st_t)
pm2[, contrast_lab := cmap2[contrast]]
pm2[, panel := fifelse(grepl("13-gene", statistic), N$prog13, N$snai1)]
save_png(ggplot(pm2, aes(x = z, y = contrast_lab, fill = z > 0)) +
           geom_vline(xintercept = 0, colour = "grey40") +
           geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey60") +
           geom_col(width = 0.62) +
           geom_text(aes(label = sprintf("%.1f", z)), size = 2.6,
                     hjust = fifelse(pm2$z > 0, -0.15, 1.15)) +
           facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
           scale_fill_manual(values = c(`TRUE` = pal_mye, `FALSE` = pal_caf), guide = "none") +
           scale_x_continuous(limits = c(-95, 95)) +
           labs(x = N$n8_x, y = NULL, title = N$n8, subtitle = N$n8_sub) +
           theme(strip.background = element_rect(fill = "grey92")),
         "FigN8_scrna_permutation_contrasts", 7.4, 4.4)

log_msg("DONE lang=", LANG, " dir=", FIGDIR)
