#!/usr/bin/env Rscript
# =============================================================================
# deepen_09: CAF -> 髓系 配体-受体优先级与零模型检验
#
# 背景：主报告的 CellChat 只给出 CAF->髓系 的通讯强度原始值（如 Myofibroblasts
# -> SPP1+ = 0.0140），没有零模型、没有 p 值、也没有拆分到具体的配体-受体对。
# 只凭原始权重无法回答"这条轴是否真的强于随机期望"。本脚本补三件事：
#
#   1. 把 CAF->髓系 的聚合强度放到全部 source->target 组合里排位，看它是不是
#      全局特异性的（而非所有细胞对都高）。
#   2. 用"表达量匹配的随机基因对"作为零模型，检验 CellChat 收录的 CAF 配体 x
#      髓系受体 共表达统计量是否高于同丰度随机基因对（B = 1000）。
#   3. 拆到具体配体-受体对，给出贡献排序与每对的零模型 z / p。
#
# 输出：
#   02_results/deepen/cellchat_rank_specificity.csv
#   02_results/deepen/cellchat_lr_pairs_caf_myeloid.csv
#   02_results/deepen/cellchat_null_summary.csv
# =============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(data.table)
  library(Matrix)
})

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

dir.create("02_results/deepen", showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- 载入对象
log_msg("loading cellchat object ...")
cc <- readRDS("02_results/scrna/cellchat_obj.rds")
ct_levels <- levels(cc@idents)
log_msg("cell types (", length(ct_levels), "): ", paste(ct_levels, collapse = ", "))

stroma_src <- intersect(c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3"), ct_levels)
myeloid_tgt <- intersect(c("SPP1+", "Pro-inflammatory", "Myeloid_cDC",
                           "Myeloid_Unknown", "Myeloid_Proliferating"), ct_levels)
log_msg("stromal sources: ", paste(stroma_src, collapse = ", "))
log_msg("myeloid targets: ", paste(myeloid_tgt, collapse = ", "))

# ------------------------------------------------- 1) 全矩阵排位特异性检验
prob <- cc@netP$prob          # [source, target, pathway]
pw <- dimnames(prob)[[3]]

# 每个 source->target 的聚合强度（跨通路求和 / 均值）
agg <- function(s, t, fun = sum) fun(prob[s, t, ])
all_pairs <- CJ(source = ct_levels, target = ct_levels, sorted = FALSE)
all_pairs <- all_pairs[source != target]
all_pairs[, agg_sum := mapply(function(s, t) agg(s, t, sum), source, target)]
all_pairs[, agg_mean := mapply(function(s, t) agg(s, t, mean), source, target)]
setorder(all_pairs, -agg_sum)
all_pairs[, rank_sum := seq_len(.N)]

ca2my <- all_pairs[source %in% stroma_src & target %in% myeloid_tgt]
oth <- all_pairs[!(source %in% stroma_src & target %in% myeloid_tgt)]

rank_tab <- data.table(
  pair = paste0(ca2my$source, " -> ", ca2my$target),
  agg_sum = ca2my$agg_sum,
  agg_mean = ca2my$agg_mean,
  rank_sum = ca2my$rank_sum,
  n_total_pairs = nrow(all_pairs),
  percentile = 100 * (1 - (ca2my$rank_sum - 1) / nrow(all_pairs))
)

# CAF->髓系 vs 其余组合的 Mann-Whitney
wt <- wilcox.test(ca2my$agg_sum, oth$agg_sum, alternative = "greater")
spec <- data.table(
  group = c("CAF->myeloid", "all other source->target"),
  n = c(nrow(ca2my), nrow(oth)),
  median_agg_sum = c(median(ca2my$agg_sum), median(oth$agg_sum)),
  mean_agg_sum = c(mean(ca2my$agg_sum), mean(oth$agg_sum)),
  median_agg_mean = c(median(ca2my$agg_mean), median(oth$agg_mean))
)
spec_tab <- rbind(
  data.table(metric = "CAF->myeloid vs others (agg_sum, one-sided MWU)",
             estimate = wt$statistic, p = wt$p.value,
             n_a = nrow(ca2my), n_b = nrow(oth)),
  data.table(metric = "median rank of CAF->myeloid pairs",
             estimate = median(ca2my$rank_sum), p = NA_real_,
             n_a = nrow(ca2my), n_b = nrow(all_pairs))
)

fwrite(rank_tab, "02_results/deepen/cellchat_rank_specificity.csv")
fwrite(spec_tab, "02_results/deepen/cellchat_specificity_test.csv")
fwrite(spec, "02_results/deepen/cellchat_specificity_groups.csv")
log_msg("rank/specificity done; CAF->myeloid n=", nrow(ca2my),
        " median rank=", median(ca2my$rank_sum), "/", nrow(all_pairs),
        " MWU p=", signif(wt$p.value, 3))

# --------------------------------------------- 2) 表达匹配随机基因对零模型
log_msg("building expression-matched null ...")
# 注意：CellChat 是在自己下采样后的 19425 个细胞上构建的（每类上限 2500），
# 与 gse132465_scored_subset.rds 的 15003 个细胞不是同一集合，标签粒度也不同。
# 因此这里统一使用 CellChat 自带的 idents + data.signaling，保证标签与矩阵严格一致。
ct <- as.character(cc@idents)
expr <- as(cc@data.signaling, "dgCMatrix")
stopifnot(length(ct) == ncol(expr))
log_msg("using CellChat's own matrix: ", nrow(expr), " genes x ", ncol(expr), " cells")

# CellChat 内部使用的信号基因（data.signaling 的行名即参与通讯计算的基因）
sig_genes <- rownames(cc@data.signaling)
sig_genes <- intersect(sig_genes, rownames(expr))
log_msg("signaling genes available: ", length(sig_genes))

# 提取 CellChat LR 数据库（兼容新旧版本槽位结构）
lrdb <- NULL
if (!is.null(cc@LR$LRdb) && nrow(cc@LR$LRdb) > 0) {
  lrdb <- as.data.frame(cc@LR$LRdb)
} else {
  lrdb <- tryCatch(CellChatDB.human$complex_input_expanded, error = function(e) NULL)
  if (is.null(lrdb)) {
    lrdb <- tryCatch(CellChat::CellChatDB.human$interaction, error = function(e) NULL)
  }
}
if (is.null(lrdb)) stop("cannot locate CellChat LR database")

col_lig <- names(lrdb)[grep("ligand", names(lrdb), ignore.case = TRUE)][1]
col_rec <- names(lrdb)[grep("receptor", names(lrdb), ignore.case = TRUE)][1]
lrdb$.lig <- as.character(lrdb[[col_lig]])
lrdb$.rec <- as.character(lrdb[[col_rec]])
lrdb$.pair <- paste0(lrdb$.lig, "_", lrdb$.rec)
log_msg("LR db rows: ", nrow(lrdb), " using cols: ", col_lig, " / ", col_rec)

# 只保留单基因配体 / 单基因受体（跳过复合物，避免拆分歧义）
lrdb <- lrdb[lrdb$.lig %in% rownames(expr) & lrdb$.rec %in% rownames(expr), ]
log_msg("single-gene LR pairs measurable: ", nrow(lrdb))

# 预计算一次「基因 x 细胞类型」均值矩阵，后续全部走查表（否则逐对重算会慢到不可行）
cell_idx <- lapply(ct_levels, function(cl) which(ct == cl))
names(cell_idx) <- ct_levels
MEAN_CT <- sapply(cell_idx, function(ix) Matrix::rowMeans(expr[, ix, drop = FALSE]))
rownames(MEAN_CT) <- rownames(expr)
mean_by_ct <- function(genes) MEAN_CT[genes, , drop = FALSE]

# 统计量与零模型 -------------------------------------------------------
# stat = mean_over_pairs( mean_ligand_in_CAF * mean_receptor_in_myeloid )
# 零模型：把配体、受体分别替换成"同表达十分位"的随机信号基因，重复 B 次。
B <- 1000
set.seed(20240924)

calc_stat <- function(lig_genes, rec_genes) {
  Lg <- rowMeans(mean_by_ct(lig_genes)[, stroma_src, drop = FALSE])
  Rg <- rowMeans(mean_by_ct(rec_genes)[, myeloid_tgt, drop = FALSE])
  mean(Lg * Rg)
}

obs <- calc_stat(lrdb$.lig, lrdb$.rec)
log_msg("observed CAF-ligand x myeloid-receptor stat = ", signif(obs, 5))

# 表达十分位分箱：按全细胞平均表达
gene_mean <- Matrix::rowMeans(expr[sig_genes, , drop = FALSE])
names(gene_mean) <- sig_genes
dec <- cut(gene_mean, breaks = quantile(gene_mean, probs = seq(0, 1, 0.1)),
           include.lowest = TRUE, labels = FALSE)
names(dec) <- sig_genes
dec_lig <- dec[lrdb$.lig]; dec_rec <- dec[lrdb$.rec]
pool_by_dec <- split(sig_genes, dec)

nullv <- numeric(B)
for (b in seq_len(B)) {
  lg <- vapply(seq_len(nrow(lrdb)), function(i) {
    pool <- pool_by_dec[[as.character(dec_lig[i])]]
    if (is.null(pool) || !length(pool)) return(NA_character_)
    sample(pool, 1)
  }, character(1))
  rc <- vapply(seq_len(nrow(lrdb)), function(i) {
    pool <- pool_by_dec[[as.character(dec_rec[i])]]
    if (is.null(pool) || !length(pool)) return(NA_character_)
    sample(pool, 1)
  }, character(1))
  keep <- !is.na(lg) & !is.na(rc)
  nullv[b] <- calc_stat(lg[keep], rc[keep])
}

null_mean <- mean(nullv); null_sd <- sd(nullv)
z <- (obs - null_mean) / null_sd
p_perm <- (sum(nullv >= obs) + 1) / (B + 1)
log_msg(sprintf("null mean=%.5f sd=%.5f  z=%.2f  p=%.4g", null_mean, null_sd, z, p_perm))

# ------------------------------------------- 3) 逐配体-受体对的贡献与零模型
Lg_all <- rowMeans(mean_by_ct(lrdb$.lig)[, stroma_src, drop = FALSE])
Rg_all <- rowMeans(mean_by_ct(lrdb$.rec)[, myeloid_tgt, drop = FALSE])
pair_tab <- data.table(
  ligand = lrdb$.lig,
  receptor = lrdb$.rec,
  pair = lrdb$.pair,
  mean_ligand_in_CAF = Lg_all,
  mean_receptor_in_myeloid = Rg_all,
  product = Lg_all * Rg_all
)
setorder(pair_tab, -product)
pair_tab[, contribution_pct := 100 * product / sum(product)]

# 每个配对单独的零模型：把配体/受体分别换成同表达十分位的随机信号基因，
# 重复 1000 次取 product 的零分布（查表实现，成本可忽略）。
top_pairs <- head(pair_tab, 30)
lig_ct <- rowMeans(MEAN_CT[, stroma_src, drop = FALSE])     # 每个基因在 CAF 的平均表达
rec_ct <- rowMeans(MEAN_CT[, myeloid_tgt, drop = FALSE])    # 每个基因在髓系的平均表达
rand_prod <- function(g_lig, g_rec, B = 1000) {
  pooll <- pool_by_dec[[as.character(dec[g_lig])]]
  poolr <- pool_by_dec[[as.character(dec[g_rec])]]
  if (is.null(pooll) || is.null(poolr)) return(c(NA_real_, NA_real_))
  lg <- sample(pooll, B, replace = TRUE); rc <- sample(poolr, B, replace = TRUE)
  v <- lig_ct[lg] * rec_ct[rc]
  c(mean(v), sd(v))
}
nm <- numeric(nrow(top_pairs)); ns <- numeric(nrow(top_pairs))
for (i in seq_len(nrow(top_pairs))) {
  r <- rand_prod(top_pairs$ligand[i], top_pairs$receptor[i])
  nm[i] <- r[1]; ns[i] <- r[2]
}
top_pairs[, null_mean := nm]
top_pairs[, null_sd := ns]
top_pairs[, z_pair := (product - null_mean) / null_sd]
top_pairs[, p_pair := (pnorm(z_pair, lower.tail = FALSE))]
fwrite(pair_tab, "02_results/deepen/cellchat_lr_pairs_caf_myeloid.csv")
fwrite(top_pairs[, .(ligand, receptor, mean_ligand_in_CAF, mean_receptor_in_myeloid,
                     product, contribution_pct, null_mean, null_sd, z_pair, p_pair)],
       "02_results/deepen/cellchat_lr_top30.csv")

fwrite(data.table(
  statistic = c("observed", "null_mean", "null_sd", "z", "p_permutation",
                "B", "n_LR_pairs", "null_2.5pct", "null_97.5pct"),
  value = c(obs, null_mean, null_sd, z, p_perm, B, nrow(lrdb),
            quantile(nullv, 0.025), quantile(nullv, 0.975))
), "02_results/deepen/cellchat_null_summary.csv")

log_msg("DONE deepen_09")
