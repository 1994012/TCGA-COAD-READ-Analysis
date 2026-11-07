#!/usr/bin/env Rscript
# =============================================================================
# deepen_11: 深化分析图件
#
# 产出 03_figures/deepen/ 下的一组图，全部由 02_results/deepen/ 里已落盘的
# 结果表绘制，不重新计算任何统计量（保证图与表严格一致）。
#
#   FigD1  多队列随机效应 meta 森林图（HR / 95%CI，含合并估计）
#   FigD2  GSE17537 外部队列 KM 曲线（OS 与 DFS）
#   FigD3  切点敏感性：三种切点 x 四个队列的 HR 森林图
#   FigD4  空转逐切片：13 基因签名与髓系 / CAF 评分的 Spearman rho 对比
#   FigD5  单细胞置换零分布（13 基因程序：髓系 - T 细胞）
#   FigD6  CellChat CAF -> 髓系 前 15 配体-受体对
#   FigD7  风险分组通路 NES 的跨队列一致性散点
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(survival); library(survminer)
  library(patchwork); library(scales)
})

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
dir.create("03_figures/deepen", showWarnings = FALSE, recursive = TRUE)

theme_set(theme_bw(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  axis.text = element_text(colour = "black")))

ggsave2 <- function(p, name, w = 7, h = 5) {
  fn <- file.path("03_figures/deepen", name)
  ggsave(fn, p, width = w, height = h, dpi = 300)
  log_msg("saved ", fn)
}

# ------------------------------------------------------------------ FigD1
ct <- fread("02_results/deepen/meta_cohort_table.csv")
mr <- fread("02_results/deepen/meta_random_effects.csv")
ct[, cohort := fifelse(grepl("TCGA", cohort), paste0(cohort, " (训练)"), cohort)]

d1 <- rbind(
  ct[, .(cohort, hr = HR, lo, hi, n, events, type = "队列")],
  data.table(cohort = "合并（随机效应）", hr = mr[analysis == "all cohorts (train + external)"]$hr,
             lo = mr[analysis == "all cohorts (train + external)"]$lo,
             hi = mr[analysis == "all cohorts (train + external)"]$hi,
             n = NA_integer_, events = NA_integer_, type = "合并"),
  data.table(cohort = "合并（仅外部队列）", hr = mr[analysis == "external only"]$hr,
             lo = mr[analysis == "external only"]$lo,
             hi = mr[analysis == "external only"]$hi,
             n = NA_integer_, events = NA_integer_, type = "合并")
)
d1[, cohort := factor(cohort, levels = rev(cohort))]

p1 <- ggplot(d1, aes(x = hr, y = cohort, colour = type)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 2, 3, 4), labels = number_format(accuracy = 0.01)) +
  scale_colour_manual(values = c("队列" = "#2c7bb6", "合并" = "#d7191c"), guide = "none") +
  labs(x = "每单位风险评分的风险比 HR (95% CI)", y = NULL,
       title = "图 D1  多队列随机效应 meta 分析") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave2(p1, "FigD1_meta_forest.png", 7.5, 4.2)

# ------------------------------------------------------------------ FigD2
rs <- fread("02_results/deepen/gse17537_riskscore.csv")
log_msg("GSE17537 columns: ", paste(names(rs), collapse = ", "))
km_plot <- function(dat, tcol, ecol, mcol, title, subtitle) {
  d <- copy(dat)
  d[, .t := get(tcol) / 12]                       # GEO 时间为月 -> 年
  d[, .e := get(ecol)]
  d[, .g := fifelse(get(mcol) >= median(get(mcol), na.rm = TRUE), "High", "Low")]
  d <- d[!is.na(.t) & !is.na(.e) & .t > 0]
  fit <- survfit(Surv(.t, .e) ~ .g, data = d)
  p <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE,
                  conf.int = FALSE, palette = c("#d7191c", "#2c7bb6"),
                  legend.title = "风险分组", legend.labs = c("High", "Low"),
                  xlab = "时间（年）", ylab = "总生存概率",
                  title = title, subtitle = subtitle,
                  ggtheme = theme_bw(base_size = 11))
  p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                           plot.subtitle = element_text(hjust = 0.5))
  p
}
os_cols <- intersect(names(rs), c("os_time", "os_delay", "os.delay (months)"))
if (length(os_cols)) {
  k <- km_plot(rs, os_cols[1], "os_event", "riskScore",
               "图 D2  GSE17537 外部验证（第三个 GEO 队列）", "总生存 OS")
  fn <- file.path("03_figures/deepen", "FigD2a_gse17537_os.png")
  ggsave(fn, k$plot, width = 6.5, height = 5, dpi = 300); log_msg("saved ", fn)
}
dfs_cols <- intersect(names(rs), c("dfs_time", "dfs_delay", "rfs_delay", "rfs.delay"))
if (length(dfs_cols)) {
  k <- km_plot(rs, dfs_cols[1], "dfs_event", "riskScore",
               "图 D2b  GSE17537 无病生存（探索性）", "DFS / RFS")
  fn <- file.path("03_figures/deepen", "FigD2b_gse17537_dfs.png")
  ggsave(fn, k$plot, width = 6.5, height = 5, dpi = 300); log_msg("saved ", fn)
}

# ------------------------------------------------------------------ FigD3
cp <- fread("02_results/deepen/cutpoint_sensitivity.csv")
cp[, cohort2 := fifelse(grepl("TCGA", cohort), "TCGA（训练）", cohort)]
cp[, cut := factor(cut, levels = c("median", "tertile T1 vs T3", "maxstat (exploratory)"),
                   labels = c("中位数", "三分位 T1 vs T3", "maxstat 最优"))]
p3 <- ggplot(cp, aes(x = hr, y = interaction(cut, cohort2), colour = cohort2)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = c(0.5, 1, 2, 3, 5)) +
  labs(x = "High vs Low 风险比 HR (95% CI)", y = NULL,
       title = "图 D3  切点敏感性分析",
       subtitle = "三种切点定义下风险分层的稳定性") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave2(p3, "FigD3_cutpoint_sensitivity.png", 7.5, 5)

# ------------------------------------------------------------------ FigD4
sp <- fread("02_results/deepen/spatial_permutation_both_signatures.csv")
sp13 <- sp[grepl("13-gene", signature)]
sp13[, slice_lab := paste0(substr(cohort, 1, 9), " ", slice)]
longd <- rbind(
  sp13[, .(slice_lab, cohort, rho = rho_sig_mye, target = "髓系评分")],
  sp13[, .(slice_lab, cohort, rho = rho_sig_caf, target = "CAF 评分")]
)
p4 <- ggplot(longd, aes(x = target, y = rho, fill = target, colour = cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_boxplot(alpha = 0.25, outlier.size = 0.8, width = 0.55) +
  geom_line(aes(group = slice_lab), colour = "grey40", linewidth = 0.35,
            position = position_dodge(width = 0)) +
  geom_point(size = 2, position = position_jitter(width = 0.08, height = 0)) +
  scale_fill_manual(values = c("髓系评分" = "#d7191c", "CAF 评分" = "#2c7bb6"), guide = "none") +
  labs(x = NULL, y = "13 基因签名 vs 细胞评分的 Spearman rho",
       title = "图 D4  空转逐切片：签名与髓系 / CAF 的相关",
       subtitle = "灰线连接同一切片；两队列共 8 张切片") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave2(p4, "FigD4_spatial_slice_rho.png", 6.5, 5)

# ------------------------------------------------------------------ FigD5
pm <- fread("02_results/deepen/scrna_permutation.csv")
rowsel <- pm[statistic == "13-gene program consensus" & contrast == "Myeloids - T cells"]
if (nrow(rowsel)) {
  r <- rowsel[1]
  nullx <- seq(r$null_2.5 * 1.6, r$null_97.5 * 1.6, length.out = 400)
  dfn <- data.frame(x = nullx,
                    y = dnorm(nullx, mean = r$null_mean, sd = r$null_sd))
  p5 <- ggplot(dfn, aes(x = x, y = y)) +
    geom_area(fill = "grey70", alpha = 0.6) +
    geom_vline(xintercept = r$observed, colour = "#d7191c", linewidth = 1.1) +
    geom_vline(xintercept = c(r$null_2.5, r$null_97.5), colour = "grey30",
               linetype = "dotted") +
    annotate("text", x = r$observed, y = max(dfn$y) * 0.75,
             label = sprintf("观测值 %.3f\nz = %.1f, p = %.4f",
                             r$observed, r$z, r$p_permutation),
             colour = "#d7191c", hjust = -0.05, size = 3.6) +
    labs(x = "髓系 - T 细胞的平均差值（随机重排零分布）", y = "密度",
         title = "图 D5  单细胞髓系优势的置换检验",
         subtitle = "GSE132465 子集；零分布由 2000 次细胞标签重排得到") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave2(p5, "FigD5_scrna_permutation_null.png", 6.5, 4.2)
}

# ------------------------------------------------------------------ FigD6
lr <- fread("02_results/deepen/cellchat_lr_top30.csv")
if (nrow(lr)) {
  lr15 <- head(lr, 15)
  lr15[, pair := paste0(ligand, " - ", receptor)]
  lr15[, pair := factor(pair, levels = rev(pair))]
  p6 <- ggplot(lr15, aes(x = pair, y = product)) +
    geom_col(fill = "#2c7bb6", width = 0.7) +
    coord_flip() +
    labs(x = NULL, y = "CAF 配体平均表达 x 髓系受体平均表达",
         title = "图 D6  CAF -> 髓系 通讯的主要配体-受体对",
         subtitle = "按贡献排序的前 15 对") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  ggsave2(p6, "FigD6_lr_pairs.png", 6.5, 5)
}

# ------------------------------------------------------------------ FigD7
cc <- fread("02_results/deepen/riskgroup_pathway_crosscohort.csv")
if (all(c("nes_TCGA", "nes_GSE39582") %in% names(cc))) {
  cc[, sig := fifelse(padj_TCGA < 0.05 & padj_GSE39582 < 0.05, "两队列均显著",
                      fifelse(padj_TCGA < 0.05 | padj_GSE39582 < 0.05, "单队列显著", "不显著"))]
  lab <- cc[grepl("HALLMARK", pathway) & (padj_TCGA < 0.05)]
  lab <- head(lab[order(-abs(nes_TCGA))], 8)
  p7 <- ggplot(cc, aes(x = nes_TCGA, y = nes_GSE39582, colour = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40") +
    ggrepel::geom_text_repel(aes(label = sub("HALLMARK_", "", pathway)), data = lab,
                             size = 2.8, max.overlaps = 20, segment.size = 0.2) +
    scale_colour_manual(values = c("两队列均显著" = "#d7191c",
                                   "单队列显著" = "#fdae61",
                                   "不显著" = "grey70")) +
    labs(x = "TCGA-COADREAD  NES（高危 vs 低危）",
         y = "GSE39582  NES（高危 vs 低危）",
         colour = NULL,
         title = "图 D7  风险分组通路信号的跨队列一致性",
         subtitle = "每个点为一个 Hallmark / KEGG 通路") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  ggsave2(p7, "FigD7_pathway_crosscohort.png", 6.5, 5.5)
}

log_msg("DONE deepen_11")
