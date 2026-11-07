# =============================================================================
# 补充图件：8 种新图型（ROC / AUC-时间线 / 风险评分分布 / 火山图 / 通路热图 /
#           气泡图 / 棒棒糖图 / 空间二维密度散点）
# 只读取已落盘的结果表与风险评分表，不改动任何已有文件。
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(timeROC); library(survival)
  library(patchwork); library(ggrepel); library(ggridges); library(scales)
})
options(stringsAsFactors = FALSE)
set.seed(20240924)
fig_dir <- "03_figures/deepen"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(sprintf(...), "\n", sep = "")

pal4 <- c("TCGA" = "#2C7BB6", "GSE39582" = "#D7191C",
          "GSE17536" = "#FDAE61", "GSE17537" = "#1A9641")
theme_set(theme_bw(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey92")))

# ---------------------------------------------------------------- 工具函数
# 读取某队列的 (时间[年], 事件, 评分)
load_cohort <- function(tag, file, time_div) {
  d <- fread(file)
  tcol <- if ("os_time" %in% names(d)) "os_time" else "os.delay (months)"
  ecol <- if ("os_event" %in% names(d)) "os_event" else "os.event"
  mcol <- if ("riskScore" %in% names(d)) "riskScore" else "riskScore"
  d <- d[!is.na(get(tcol)) & !is.na(get(ecol)) & !is.na(get(mcol))]
  data.table(cohort = tag,
             t = as.numeric(d[[tcol]]) / time_div,
             e = as.numeric(d[[ecol]]),
             mk = as.numeric(d[[mcol]]))
}
# timeROC 的 AUC 按 times 对齐提取
auc_at <- function(t, e, mk, times) {
  ok <- !is.na(t) & !is.na(e) & !is.na(mk)
  t <- t[ok]; e <- e[ok]; mk <- mk[ok]
  tt <- times[times < max(t[e == 1])]
  if (!length(tt)) return(rep(NA_real_, length(times)))
  r <- timeROC(T = t, delta = e, marker = mk, cause = 1, times = tt, iid = FALSE)
  out <- rep(NA_real_, length(times))
  for (k in seq_along(tt)) {
    j <- which(abs(as.numeric(r$times) - tt[k]) < 1e-8)[1]
    if (!is.na(j)) out[which(abs(times - tt[k]) < 1e-8)] <- as.numeric(r$AUC)[j]
  }
  out
}

# ============================================================ FigE1: ROC 曲线
tcga <- load_cohort("TCGA", "02_results/model/tcga_riskscore.csv", 365.25)
mx <- max(tcga$t[tcga$e == 1])
tt <- c(1, 3, 5)[c(1, 3, 5) < mx]
r <- timeROC(T = tcga$t, delta = tcga$e, marker = tcga$mk, cause = 1,
             times = tt, iid = FALSE)
roc_dt <- rbindlist(lapply(seq_along(tt), function(i) {
  j <- which(abs(as.numeric(r$times) - tt[i]) < 1e-8)[1]
  if (is.na(j)) return(NULL)
  data.table(year = tt[i],
             fpr = as.numeric(r$FP[, j]),
             tpr = as.numeric(r$TP[, j]),
             auc = as.numeric(r$AUC)[j])
}))
p1 <- ggplot(roc_dt, aes(fpr, tpr, color = factor(year))) +
  geom_path(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  scale_color_manual(values = c("1" = "#2C7BB6", "3" = "#D7191C", "5" = "#1A9641"),
                     labels = sapply(split(roc_dt, roc_dt$year),
                                     function(x) sprintf("%d 年 (AUC = %.3f)",
                                                         unique(x$year), unique(x$auc))),
                     name = "随访时点") +
  coord_equal() + lims(x = c(0, 1), y = c(0, 1)) +
  labs(x = "假阳性率 (1 - 特异度)", y = "真阳性率 (敏感度)",
       title = "TCGA-COADREAD 训练集时间依赖 ROC")
ggsave(file.path(fig_dir, "FigE1_roc_tcga.png"), p1,
       width = 6.2, height = 5.6, dpi = 300)

# ==================================================== FigE2: AUC 随时间变化
coh <- rbind(
  tcga,
  load_cohort("GSE39582", "02_results/validation/gse39582_riskscore.csv", 12),
  load_cohort("GSE17536", "02_results/validation/gse17536_riskscore.csv", 12),
  load_cohort("GSE17537", "02_results/deepen/gse17537_riskscore.csv", 12))
grid <- seq(0.5, 10, by = 0.25)
auc_dt <- rbindlist(lapply(split(coh, coh$cohort), function(d) {
  tg <- grid[grid < max(d$t[d$e == 1]) * 0.999]
  a <- auc_at(d$t, d$e, d$mk, tg)
  data.table(cohort = unique(d$cohort), year = tg, auc = a)
}))
p2 <- ggplot(auc_dt[!is.na(auc)], aes(year, auc, color = cohort)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_path(linewidth = 0.85, alpha = 0.9) +
  geom_point(data = auc_dt[year %in% c(1, 3, 5) & !is.na(auc)],
             size = 1.8) +
  scale_color_manual(values = pal4) +
  scale_y_continuous(limits = c(0.4, 0.85), breaks = seq(0.4, 0.8, 0.1)) +
  labs(x = "随访时间（年）", y = "时间依赖 AUC",
       color = "队列", title = "时间依赖 AUC 随访轨迹：训练集与三个外部队列")
ggsave(file.path(fig_dir, "FigE2_auc_vs_time.png"), p2,
       width = 7.2, height = 4.6, dpi = 300)

# ============================================ FigE3: 风险评分分布（山脊图）
coh[, mk_z := (mk - mean(mk)) / sd(mk), by = cohort]
p3 <- ggplot(coh, aes(x = mk_z, y = reorder(cohort, mk_z, median), fill = cohort)) +
  ggridges::geom_density_ridges(alpha = 0.75, scale = 1.1, rel_min_height = 0.01) +
  scale_fill_manual(values = pal4) +
  labs(x = "队列内 z-score", y = NULL, title = "四个队列中风险评分的分布") +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "FigE3_riskscore_ridge.png"), p3,
       width = 6.4, height = 3.8, dpi = 300)

# ================================================== FigE4: 火山图（差异表达）
de <- fread("02_results/deepen/riskgroup_de_tcga.csv")
de[, sig := fcase(adj_p < 0.05 & logFC >  log2(1.5), "高危组上调",
                  adj_p < 0.05 & logFC < -log2(1.5), "低危组上调",
                  default = "不显著")]
lab <- rbind(
  de[sig == "高危组上调"][order(adj_p)][1:6],
  de[sig == "低危组上调"][order(adj_p)][1:6],
  de[gene %in% c("SNAI1","TIMP1","MMP1","IL6","SPP1","VEGFA","FSTL3","COL1A1")][!is.na(gene)])
lab <- unique(lab)
p4 <- ggplot(de, aes(logFC, -log10(pmax(adj_p, 1e-300)), color = sig)) +
  geom_point(size = 0.7, alpha = 0.55) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55") +
  geom_point(data = lab, size = 1.1, color = "black") +
  geom_text_repel(data = lab, aes(label = gene), size = 2.6,
                  max.overlaps = 20, segment.size = 0.25, color = "black") +
  scale_color_manual(values = c("高危组上调" = "#D7191C",
                                "低危组上调" = "#2C7BB6",
                                "不显著" = "grey75"), name = NULL) +
  labs(x = expression(log[2]*" 倍数变化（高危 / 低危）"),
       y = expression(-log[10]*"（校正 "*italic(p)*"）"),
       title = "TCGA 高低危组差异表达火山图（805 个显著基因）") +
  theme(legend.position = c(0.08, 0.9), legend.background = element_blank())
ggsave(file.path(fig_dir, "FigE4_volcano_de.png"), p4,
       width = 7.0, height = 5.4, dpi = 300)

# ==================================================== FigE5: 通路 NES 热图
cr <- fread("02_results/deepen/riskgroup_pathway_crosscohort.csv")
cr[, nes_mean := (abs(nes_TCGA) + abs(nes_GSE39582)) / 2]
top <- cr[order(-nes_mean)][1:24]
lab_clean <- function(x) gsub("^(HALLMARK_|KEGG_)", "", x)
top[, pw := lab_clean(pathway)]
top[, pw := factor(pw, levels = rev(pw))]
melten <- melt(top[, .(pw, TCGA = nes_TCGA, GSE39582 = nes_GSE39582)],
               id.vars = "pw", variable.name = "cohort", value.name = "NES")
sigm <- melt(top[, .(pw, TCGA = padj_TCGA, GSE39582 = padj_GSE39582)],
             id.vars = "pw", variable.name = "cohort", value.name = "padj")
melten <- merge(melten, sigm, by = c("pw", "cohort"))
melten[, star := fcase(padj < 0.001, "***", padj < 0.01, "**",
                       padj < 0.05, "*", default = "")]
p5 <- ggplot(melten, aes(cohort, pw, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 3, vjust = 0.75) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, limits = c(-3, 3), oob = squish,
                       name = "NES") +
  labs(x = NULL, y = NULL,
       title = "通路富集的跨队列热图",
       subtitle = "按平均 |NES| 取前 24 条；* 校正 p<0.05，** <0.01，*** <0.001") +
  theme(plot.title = element_text(size = 12))
  theme(axis.text.y = element_text(size = 7.5))
ggsave(file.path(fig_dir, "FigE5_pathway_heatmap.png"), p5,
       width = 6.2, height = 7.8, dpi = 300)

# ==================================================== FigE6: 配体-受体气泡图
lr <- fread("02_results/deepen/cellchat_lr_top30.csv")
lr[, pair := paste0(ligand, " \u2192 ", receptor)]
lr20 <- lr[order(-product)][1:20]
p6 <- ggplot(lr20, aes(product, z_pair)) +
  geom_hline(yintercept = 1.96, linetype = "dashed", color = "grey55") +
  geom_point(aes(size = contribution_pct, fill = z_pair),
             shape = 21, color = "grey25", alpha = 0.9) +
  geom_text_repel(aes(label = pair), size = 2.5, max.overlaps = 16,
                  segment.size = 0.25) +
  scale_fill_gradient(low = "#2C7BB6", high = "#D7191C", name = "零模型 z") +
  scale_size_continuous(range = c(2, 9), name = "贡献占比 (%)") +
  labs(x = expression("配体-受体平均表达乘积（CAF" %.% "髓系）"),
       y = "表达匹配零模型 z 值",
       title = "CAF \u2192 髓系前 20 配体-受体对：强度与零模型证据",
       subtitle = "虚线：z = 1.96（零模型 p ≈ 0.05）") +
  theme(legend.box = "vertical")
ggsave(file.path(fig_dir, "FigE6_lr_bubble.png"), p6,
       width = 7.4, height = 5.6, dpi = 300)

# ================================================= FigE7: 通讯排位棒棒糖图
rk <- fread("02_results/deepen/cellchat_rank_specificity.csv")
rk[, src := sub(" ->.*", "", pair)]
rk[, is_caf_src := src %in% c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3")]
top20 <- rk[order(-agg_sum)][1:20]
top20[, pair_f := factor(pair, levels = rev(top20$pair))]
p7 <- ggplot(top20, aes(agg_sum, pair_f)) +
  geom_segment(aes(xend = 0, yend = pair_f, color = is_caf_src), linewidth = 0.8) +
  geom_point(aes(color = is_caf_src, fill = is_caf_src), shape = 21, size = 2.6) +
  scale_color_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "grey55"),
                     labels = c(`TRUE` = "CAF 来源", `FALSE` = "其他来源"), name = NULL) +
  scale_fill_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "grey55"), guide = "none") +
  labs(x = "CellChat 通讯聚合强度", y = NULL,
       title = "全部 380 个 source \u2192 target 组合中强度最高的 20 个",
       subtitle = "红色：CAF（Myofibroblasts / Stromal 1\u20133）作为信号来源") +
  theme(axis.text.y = element_text(size = 7.5), legend.position = c(0.8, 0.15))
ggsave(file.path(fig_dir, "FigE7_cellchat_lollipop.png"), p7,
       width = 7.2, height = 5.2, dpi = 300)

# ============================================ FigE8: 空间 spot 二维密度散点
sp <- fread("02_results/deepen/spatial_scores_both_signatures.csv")
sp[, slice_n := paste0(cohort, " ", slice)]
p8 <- ggplot(sp, aes(caf, mye, color = sig13)) +
  geom_point(size = 0.55, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black",
              linewidth = 0.5, linetype = "dashed") +
  scale_color_viridis_c(option = "mako", name = "13 基因\n签名评分") +
  facet_wrap(~ cohort, scales = "free") +
  labs(x = "CAF 评分", y = "髓系评分",
       title = "空间转录组逐 spot：CAF 评分与髓系评分的关系",
       subtitle = "颜色为 13 基因签名评分；虚线为线性趋势") +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))
ggsave(file.path(fig_dir, "FigE8_spatial_spot.png"), p8,
       width = 8.6, height = 4.4, dpi = 300)

log_msg("all 8 extended figures written to %s", fig_dir)
