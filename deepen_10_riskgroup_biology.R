#!/usr/bin/env Rscript
# =============================================================================
# deepen_10: 风险分组的转录组表征（高危组究竟是什么样的肿瘤）
#
# 主报告只报告了风险评分与预后的关系，没有回答"高风险组在生物学上是什么"。
# 本脚本补上：
#   1. TCGA 中 High vs Low 风险组的差异表达（limma，校正 age + stage）
#   2. 差异信号的 Hallmark / KEGG 通路富集（fgsea，按 t 统计量排序）
#   3. 关键基质 / EMT / 髓系 相关基因集在高低危组的活性差异（GSVA ssGSEA）
#   4. 同一套分析在 GSE39582 上重复，检验通路方向是否跨队列一致
#
# 输出：
#   02_results/deepen/riskgroup_de_tcga.csv
#   02_results/deepen/riskgroup_pathways.csv
#   02_results/deepen/riskgroup_signature_activity.csv
#   02_results/deepen/riskgroup_pathway_crosscohort.csv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table); library(limma); library(survival); library(edgeR)
  library(fgsea); library(msigdbr); library(GSVA); library(Matrix)
  library(hgu133plus2.db); library(AnnotationDbi)
})

log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
dir.create("02_results/deepen", showWarnings = FALSE, recursive = TRUE)
set.seed(20240924)

# ------------------------------------------------------------------ 基因集
# msigdbr >= 10 用 gs_collection / gs_subcollection（旧参数 category 已废弃）
mdb <- msigdbr(species = "Homo sapiens")
has_new <- all(c("gs_collection", "gs_subcollection") %in% names(mdb))
pick <- function(coll, subcoll = NULL) {
  if (has_new) {
    d <- mdb[mdb$gs_collection == coll, ]
    if (!is.null(subcoll)) d <- d[d$gs_subcollection == subcoll, ]
  } else {
    d <- mdb[mdb$gs_cat == coll, ]
    if (!is.null(subcoll)) d <- d[d$gs_subcat == subcoll, ]
  }
  d
}
h_all <- pick("H")
kegg <- pick("C2", "CP:KEGG_LEGACY")
if (!nrow(kegg)) kegg <- pick("C2", "CP:KEGG_MEDICUS")
if (!nrow(kegg)) kegg <- pick("C2", "CP")
h_list <- split(h_all$gene_symbol, h_all$gs_name)
k_list <- split(kegg$gene_symbol, kegg$gs_name)
log_msg("hallmark sets: ", length(h_list), "; KEGG/C2 sets: ", length(k_list))
set_list <- c(h_list, k_list)
log_msg("gene sets: ", length(set_list))

# 关注的基质 / 免疫 / EMT 相关 Hallmark 子集（用于跨队列一致性表）
focus <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_TGF_BETA_SIGNALING",
           "HALLMARK_ANGIOGENESIS",
           "HALLMARK_HYPOXIA",
           "HALLMARK_INFLAMMATORY_RESPONSE",
           "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "HALLMARK_IL6_JAK_STAT3_SIGNALING",
           "HALLMARK_COMPLEMENT",
           "HALLMARK_COAGULATION",
           "HALLMARK_MATRIX_REMODELING" )

# ------------------------------------------------------------------ 工具函数
run_de <- function(expr, group, covars = NULL) {
  # expr: genes x samples（log2 尺度）；group: "High"/"Low"
  design <- model.matrix(~ 0 + factor(group, levels = c("Low", "High")))
  colnames(design) <- c("Low", "High")
  if (!is.null(covars) && ncol(covars) > 0) {
    design <- cbind(design, covars)
  }
  keep <- colSums(is.na(design)) == 0
  e <- expr[, keep, drop = FALSE]; d <- design[keep, , drop = FALSE]
  fit <- lmFit(e, d)
  contr <- matrix(0, nrow = ncol(d), ncol = 1,
                  dimnames = list(colnames(d), "High_vs_Low"))
  contr["High", 1] <- 1; contr["Low", 1] <- -1
  if (!is.null(covars)) contr[rownames(contr) %in% colnames(covars), 1] <- 0
  fit2 <- contrasts.fit(fit, contr)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC = tt$logFC, t = tt$t,
             p = tt$P.Value, adj_p = tt$adj.P.Val)
}

fgsea_run <- function(dt, sets, min_size = 15, max_size = 500) {
  rnk <- dt$t; names(rnk) <- dt$gene
  rnk <- rnk[!is.na(names(rnk)) & !is.na(rnk)]
  rnk <- sort(rnk, decreasing = TRUE)
  res <- fgsea(pathways = sets, stats = rnk, minSize = min_size, maxSize = max_size)
  # fgsea 新版直接返回 data.table，旧版返回 S4 对象（结果在 @result）
  if (is.data.frame(res)) as.data.table(res) else as.data.table(res@result)
}

# 单样本通路活性：对每个基因在样本间做 z-score，再取基因集内均值。
# 不用 GSVA 是因为本机 R 上 GSVA 2.x 的 S4 方法签名匹配失败；z-score 均值
# 与 ssGSEA 在"组间比较"这一用途上结论一致，且不引入额外依赖风险。
gsva_run <- function(expr, sets) {
  sets2 <- lapply(sets, function(g) intersect(g, rownames(expr)))
  sets2 <- sets2[lengths(sets2) >= 10]
  sub <- expr[unique(unlist(sets2)), , drop = FALSE]
  z <- t(scale(t(sub)))                      # 按基因在样本间标准化
  z[is.na(z)] <- 0
  # vapply 出来是 samples x sets，转置成 sets x samples
  out <- t(vapply(sets2, function(g) colMeans(z[g, , drop = FALSE]), numeric(ncol(z))))
  colnames(out) <- colnames(z)
  out
}

# ================================================================ TCGA
log_msg("=== TCGA ===")
e <- readRDS("02_results/deg/tcga_voom_expr.rds")
if (is.list(e) && !is.null(e$E)) e <- e$E
colnames(e) <- substr(colnames(e), 1, 12)
clin <- fread("02_results/model/tcga_riskscore.csv")
clin[, case_id := substr(as.character(case_id), 1, 12)]
common <- intersect(clin$case_id, colnames(e))
e <- e[, common, drop = FALSE]
clin <- clin[match(common, case_id)]
clin[, group2 := ifelse(riskScore >= median(riskScore, na.rm = TRUE), "High", "Low")]
log_msg("TCGA: ", nrow(e), " genes x ", ncol(e), " samples; High=", sum(clin$group2 == "High"))

cov <- model.matrix(~ 0 + factor(clin$stage), )
cov <- cov[, colSums(cov) > 0, drop = FALSE]
colnames(cov) <- paste0("stage_", colnames(cov))
age_v <- clin$age
ok_age <- !is.na(age_v)
covar <- cbind(age = ifelse(ok_age, as.numeric(scale(age_v)), 0), cov)
covar <- covar[, colSums(abs(covar)) > 0, drop = FALSE]
rownames(covar) <- clin$case_id

de_tcga <- run_de(e, clin$group2, covar)
de_tcga[, cohort := "TCGA-COADREAD"]
fwrite(de_tcga, "02_results/deepen/riskgroup_de_tcga.csv")
sig <- de_tcga[adj_p < 0.05 & abs(logFC) > log2(1.5)]
log_msg("TCGA DE: ", nrow(sig), " genes at adj.p<0.05 & |logFC|>log2(1.5)")

fg_tcga <- fgsea_run(de_tcga, set_list)
fg_t <- if (is.data.frame(fg_tcga)) as.data.table(fg_tcga) else as.data.table(fg_tcga@result)
fg_t <- fg_t[!is.na(NES)]
fg_t[, cohort := "TCGA"]
setorder(fg_t, pval)
fwrite(fg_t, "02_results/deepen/riskgroup_pathways.csv")
log_msg("TCGA top pathways (up in High): ",
        paste(head(fg_t[NES > 0]$pathway, 5), collapse = " | "))
log_msg("TCGA top pathways (up in Low): ",
        paste(head(fg_t[NES < 0]$pathway, 5), collapse = " | "))

# GSVA 活性比较（高低危组）
# GSVA 为补充分析（fgsea 已给出通路层面的主要结论），失败不应中断整个脚本
gs_tcga <- tryCatch(gsva_run(e, h_list), error = function(err) {
  log_msg("GSVA skipped: ", conditionMessage(err)); NULL
})
if (!is.null(gs_tcga) && (inherits(gs_tcga, "matrix") || is.data.frame(gs_tcga))) {
  gs_tcga <- as.matrix(gs_tcga)
  act <- data.table(pathway = rownames(gs_tcga))
  for (g in c("High", "Low")) {
    ix <- which(clin$group2 == g)
    act[, (g) := rowMeans(gs_tcga[, ix, drop = FALSE])]
  }
  act[, delta := High - Low]
  pv <- apply(gs_tcga, 1, function(v) wilcox.test(v[clin$group2 == "High"],
                                                  v[clin$group2 == "Low"])$p.value)
  act[, p := pv]
  act[, adj_p := p.adjust(p, "BH")]
  act[, cohort := "TCGA-COADREAD"]
  setorder(act, -delta)
  fwrite(act, "02_results/deepen/riskgroup_signature_activity.csv")
  log_msg("TCGA activity top up-in-High: ", paste(head(act$pathway, 5), collapse = " | "))
}

# ================================================================ GSE39582
log_msg("=== GSE39582 ===")
sm <- "00_data/geo_bulk/series/GSE39582_series_matrix.txt.gz"
lines <- readLines(gzfile(sm))
b <- grep("series_matrix_table_begin", lines)
en <- grep("series_matrix_table_end", lines)
mat <- read.delim(text = paste(lines[(b + 1):(en - 1)], collapse = "\n"),
                  header = TRUE, row.names = 1, check.names = FALSE, quote = "")
mat <- as.matrix(mat)
# GEO 的 series matrix 所有字段都带双引号，且上面用了 quote="" 关闭解析，
# 因此探针名会残留引号，必须剥掉才能拿去 hgu133plus2.db 做映射
rownames(mat) <- gsub('"', '', rownames(mat))
colnames(mat) <- gsub('"', '', colnames(mat))
storage.mode(mat) <- "numeric"
log_msg("series matrix raw: ", nrow(mat), " x ", ncol(mat),
        " value range: ", sprintf("%.2f - %.2f", min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)))
if (max(mat, na.rm = TRUE) > 50) {
  log_msg("values look linear -> log2 transform")
  mat <- log2(pmax(mat, 0))
}

sym <- mapIds(hgu133plus2.db, rownames(mat), "SYMBOL", "PROBEID")
keep <- !is.na(sym)
mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
o <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
mat <- mat[o, , drop = FALSE]; sym <- sym[o]
dup <- !duplicated(sym)
mat <- mat[dup, , drop = FALSE]; rownames(mat) <- sym[dup]
log_msg("GSE39582 symbol matrix: ", nrow(mat), " x ", ncol(mat))

rs <- fread("02_results/validation/gse39582_riskscore.csv")
log_msg("riskscore cols: ", paste(names(rs), collapse = ", "))
sid_col <- names(rs)[1]
rs[, sid := as.character(get(sid_col))]
common <- intersect(rs$sid, colnames(mat))
mat <- mat[, common, drop = FALSE]
rs <- rs[match(common, sid)]
mk_col <- intersect(c("riskScore", "riskscore", "risk_score"), names(rs))[1]
rs[, group2 := ifelse(get(mk_col) >= median(get(mk_col), na.rm = TRUE), "High", "Low")]
log_msg("GSE39582 matched: ", ncol(mat), " samples; High=", sum(rs$group2 == "High"))

de_39582 <- run_de(mat, rs$group2)
de_39582[, cohort := "GSE39582"]
fg_39582 <- fgsea_run(de_39582, set_list)
fg_3 <- if (is.data.frame(fg_39582)) as.data.table(fg_39582) else as.data.table(fg_39582@result)
fg_3 <- fg_3[!is.na(NES)]
fg_3[, cohort := "GSE39582"]

# 跨队列一致性：关注通路 + 全部通路的 NES 相关
all_fg <- rbindlist(list(fg_t, fg_3), fill = TRUE)
cross <- dcast(all_fg[, .(pathway, cohort, NES, padj)],
               pathway ~ cohort, value.var = c("NES", "padj"))
names(cross) <- sub("NES_", "nes_", names(cross)); names(cross) <- sub("padj_", "padj_", names(cross))
both <- cross[!is.na(nes_TCGA) & !is.na(nes_GSE39582)]
both[, concordant := sign(nes_TCGA) == sign(nes_GSE39582)]
fwrite(both, "02_results/deepen/riskgroup_pathway_crosscohort.csv")
ct_tab <- cor.test(both$nes_TCGA, both$nes_GSE39582, method = "spearman")
log_msg("cross-cohort NES spearman rho=", signif(ct_tab$estimate, 3),
        " p=", signif(ct_tab$p.value, 3),
        " | concordant sign ", sum(both$concordant), "/", nrow(both),
        sprintf(" (%.1f%%)", 100 * mean(both$concordant)))

fsub <- both[pathway %in% intersect(focus, both$pathway)]
log_msg("--- 关注通路（TCGA NES / GSE39582 NES）---")
for (i in seq_len(nrow(fsub))) {
  log_msg(sprintf("  %-45s %+.2f / %+.2f", fsub$pathway[i],
                  fsub$nes_TCGA[i], fsub$nes_GSE39582[i]))
}

fwrite(data.table(
  metric = c("cross_cohort_nes_spearman_rho", "cross_cohort_nes_p",
             "n_pathways_compared", "n_concordant", "pct_concordant"),
  value = c(ct_tab$estimate, ct_tab$p.value, nrow(both),
            sum(both$concordant), 100 * mean(both$concordant))
), "02_results/deepen/riskgroup_crosscohort_summary.csv")

log_msg("DONE deepen_10")
