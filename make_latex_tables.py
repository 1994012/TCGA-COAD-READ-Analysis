#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 02_分析结果/deepen/ 里已落盘的 CSV 渲染成 LaTeX 表格片段 tables.tex。

设计原则：只做格式转换，不重算任何统计量 —— 保证文档里的数字与结果表严格一致。
"""
import csv
import os

BASE = os.path.expanduser(
    "~/Documents/first-year/docx/8.25/肝癌，肺癌，结直肠癌/结直肠癌_CAF髓系")
DEEP = os.path.join(BASE, "02_分析结果", "deepen")
OUT = os.path.join(BASE, "04_LaTeX", "tables.tex")


def read(name):
    with open(os.path.join(DEEP, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def esc(s):
    """LaTeX 转义；下划线后插入断点，长英文串每 12 字符插一个断点。"""
    import re
    t = (str(s).replace("&", r"\&").replace("%", r"\%")
         .replace("#", r"\#"))
    t = t.replace("_", r"\_\allowbreak{}")

    def brk(m):
        w = m.group(0)
        return r"\allowbreak{}".join(w[i:i + 12] for i in range(0, len(w), 12))

    t = re.sub(r"[A-Za-z0-9]{13,}", brk, t)
    return t


def fmt(x, nd=2):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "--"
    if v != v:
        return "--"
    return f"{v:.{nd}f}"


def pfmt(x):
    """p 值格式化：小于 1e-4 走科学计数法，否则四位小数。"""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "--"
    if v != v:
        return "--"
    if v == 0:
        return "$<1\\times10^{-4}$"
    if v < 1e-4:
        import math
        e = int(math.floor(math.log10(v)))
        m = v / (10 ** e)
        return f"${m:.1f}\\times10^{{{e}}}$"
    return f"{v:.4f}"


def longtable(caption, label, colspec, header, rows, note=None, fontsize="small",
              colsep=None):
    ncols = len(header)
    out = []
    out.append("\\begin{table}[htbp]")
    out.append("\\centering")
    out.append(f"\\{fontsize}")
    if colsep:
        out.append(f"\\setlength{{\\tabcolsep}}{{{colsep}pt}}")
    out.append("\\begin{threeparttable}")
    out.append(f"\\caption{{{caption}}}")
    out.append(f"\\label{{{label}}}")
    out.append(f"\\begin{{tabular}}{{{colspec}}}")
    out.append("\\toprule")
    out.append(" & ".join(header) + " \\\\")
    out.append("\\midrule")
    for r in rows:
        out.append(" & ".join(r) + " \\\\")
    out.append("\\bottomrule")
    out.append("\\end{tabular}")
    if note:
        out.append("\\begin{tablenotes}\\footnotesize")
        out.append("\\item[] " + note)
        out.append("\\end{tablenotes}")
    out.append("\\end{threeparttable}")
    out.append("\\end{table}")
    out.append("")
    return "\n".join(out)


PARTS = []

# ============================================================ 表：GSE17537
rows = read("gse17537_stats.csv")
body = []
lab = {"GSE17537 OS (raw log2)": "OS（原始 log2 口径）",
       "GSE17537 OS (cohort z-score)": "OS（队列内 z-score 口径）",
       "GSE17537 DFS (companion)": "DFS（伴随终点，探索性）"}
for r in rows:
    body.append([
        lab.get(r["config"], esc(r["config"])),
        r["n"], r["events"],
        pfmt(r["km_p"]),
        f'{fmt(r["hr_cont"],2)} ({fmt(r["hr_low"],2)}--{fmt(r["hr_high"],2)})',
        pfmt(r["p_cont"]),
        f'{fmt(r["calib_slope"],2)} ({fmt(r["calib_low"],2)}--{fmt(r["calib_high"],2)})',
        f'{fmt(r["auc_1y"],2)} / {fmt(r["auc_3y"],2)} / {fmt(r["auc_5y"],2)}',
    ])
PARTS.append(longtable(
    "第三个外部队列 GSE17537 的验证结果（Affy HG-U133 Plus 2.0，与 GSE39582 / GSE17536 同平台）",
    "tab:gse17537",
    "p{3.4cm}ccccp{2.0cm}p{2.4cm}p{2.0cm}",
    ["终点", "n", "事件", "KM $p$", "HR (95\\%CI)", "$p$",
     "校准斜率 (95\\%CI)", "AUC 1/3/5 年"],
    body,
    note="HR 为风险评分每升高一个单位的连续型风险比；校准斜率由转移的线性预测子直接进入 Cox 模型估计，"
         "斜率 $<1$ 提示外部队列上存在效应收缩。DFS 为伴随终点，仅作方向性参考。", fontsize="footnotesize", colsep=3))

# ============================================================ 表：四队列汇总
rows = read("meta_cohort_table.csv")
g17537 = {r["config"]: r for r in read("gse17537_stats.csv")}
body = []
for r in rows:
    c = r["cohort"]
    c = c.replace("TCGA-COADREAD (train)", "TCGA-COADREAD（训练）")
    c = c.replace("GSE39582 (tumor)", "GSE39582")
    body.append([
        esc(c), r["n"], r["events"],
        f'{fmt(r["HR"],2)} ({fmt(r["lo"],2)}--{fmt(r["hi"],2)})',
        pfmt(r["p_cont"]), pfmt(r["km_p"]),
        (f'{fmt(r["auc_1y"],3)} / {fmt(r["auc_3y"],3)} / '
         f'{fmt(g17537["GSE17537 OS (raw log2)"]["auc_5y"],3)}'
         if c == "GSE17537" else
         f'{fmt(r["auc_1y"],3)} / {fmt(r["auc_3y"],3)} / {fmt(r["auc_5y"],3)}'),
    ])
PARTS.append(longtable(
    "四个队列的验证指标汇总",
    "tab:cohorts",
    "p{3.6cm}ccp{2.6cm}ccp{2.6cm}",
    ["队列", "n", "事件数", "HR (95\\%CI)", "$p$（连续）", "$p$（KM）",
     "AUC 1/3/5 年"],
    body,
    note="HR 为连续型风险比（每单位风险评分）。AUC 为时间依赖 ROC 曲线下面积。", fontsize="footnotesize", colsep=4))

# ============================================================ 表：meta
rows = read("meta_random_effects.csv")
body = []
lab = {"all cohorts (train + external)": "全部四个队列",
       "external only": "仅外部队列"}
for r in rows:
    body.append([
        lab.get(r["analysis"], esc(r["analysis"])), r["k"],
        f'{fmt(r["hr"],2)} ({fmt(r["lo"],2)}--{fmt(r["hi"],2)})',
        pfmt(r["p"]), f'{fmt(r["I2"],1)}\\%', fmt(r["tau2"], 4),
        pfmt(r["Q_p"]), f'{fmt(r["pred_lo"],2)}--{fmt(r["pred_hi"],2)}',
    ])
PARTS.append(longtable(
    "随机效应 meta 分析结果",
    "tab:meta",
    "p{3.0cm}cp{2.6cm}ccccp{1.8cm}",
    ["分析范围", "$k$", "合并 HR (95\\%CI)", "$p$", "$I^2$", "$\\tau^2$",
     "$Q$ 检验 $p$", "预测区间"],
    body,
    note="预测区间给出的是未来一个新的相似队列中真实效应量的可能范围，通常明显宽于合并估计的置信区间。"))

rows = read("meta_leave_one_out.csv")
body = []
for r in rows:
    c = r["left_out"].replace("TCGA-COADREAD (train)", "TCGA-COADREAD（训练）")
    c = c.replace("GSE39582 (tumor)", "GSE39582")
    body.append([f"剔除 {esc(c)}", f'{fmt(r["hr"],2)} ({fmt(r["lo"],2)}--{fmt(r["hi"],2)})',
                 pfmt(r["p"]), f'{fmt(r["I2"],1)}\\%', fmt(r["tau2"], 4)])
PARTS.append(longtable(
    "meta 分析的留一法敏感性",
    "tab:loo",
    "p{4.2cm}p{2.8cm}ccc",
    ["操作", "合并 HR (95\\%CI)", "$p$", "$I^2$", "$\\tau^2$"],
    body,
    note="逐个剔除队列后重新合并。剔除训练集一行等价于「仅外部队列」的合并结果。"))

# ============================================================ 表：乐观度
rows = read("tcga_optimism_bootstrap.csv")
body = []
lab = {"Harrell C": "Harrell $C$ 指数", "AUC 1y": "AUC（1 年）",
       "AUC 3y": "AUC（3 年）", "AUC 5y": "AUC（5 年）"}
for r in rows:
    body.append([
        lab.get(r["metric"], esc(r["metric"])),
        fmt(r["apparent"], 3), fmt(r["optimism"], 3), fmt(r["optimism_corrected"], 3),
        f'{fmt(r["opt_2.5pct"],3)} -- {fmt(r["opt_97.5pct"],3)}',
        r["n_replicates_used"],
    ])
PARTS.append(longtable(
    "TCGA 训练集上 1000 次 bootstrap 的乐观度校正",
    "tab:optimism",
    "p{2.8cm}ccccc",
    ["指标", "表观值", "乐观度", "校正后", "乐观度 95\\% 范围", "有效重复数"],
    body,
    note="每次 bootstrap 均重新拟合 13 基因 Cox 模型后再评价，因此表观值与原报告直接落盘的 "
         "\\texttt{riskScore} 评价结果略有差异（本表基于表达矩阵与临床表可匹配的 591 例，"
         "原报告为 604 例）。本表的使用价值在于乐观度的量级，而非复现表观数值。"))

# ============================================================ 表：切点
rows = read("cutpoint_sensitivity.csv")
body = []
lab = {"median": "中位数", "tertile T1 vs T3": "三分位 T1 vs T3",
       "maxstat (exploratory)": "maxstat 最优"}
for r in rows:
    c = r["cohort"].replace("TCGA (train)", "TCGA（训练）")
    c = c.replace("GSE39582 (tumor)", "GSE39582")
    body.append([
        esc(c), lab.get(r["cut"], esc(r["cut"])),
        f'{r["n_high"]}/{r["n_low"]}',
        f'{fmt(r["hr"],2)} ({fmt(r["lo"],2)}--{fmt(r["hi"],2)})',
        pfmt(r["p"]), pfmt(r["km_p"]),
        pfmt(r["p_maxstat_adjusted"]) if r["p_maxstat_adjusted"] else "--",
    ])
PARTS.append(longtable(
    "切点敏感性分析：三种切点定义 $\\times$ 四个队列",
    "tab:cutpoint",
    "p{2.3cm}p{2.4cm}cp{2.4cm}ccc",
    ["队列", "切点定义", "高危/低危", "HR (95\\%CI)", "$p$", "KM $p$",
     "maxstat 校正 $p$"],
    body,
    note="maxstat 校正 $p$ 为条件 Monte Carlo 重抽样（$B=9999$）后对最优切点选择过程的校正，"
         "仅在探索性切点一行给出。12 个组合的 HR 全部 $>1$", fontsize="footnotesize", colsep=4))

# ============================================================ 表：多因素
rows = read("multivariate_covariates.csv")
body = []
for r in rows:
    c = r["cohort"].replace("TCGA (train)", "TCGA（训练）")
    c = c.replace("GSE39582 (tumor)", "GSE39582")
    body.append([
        esc(c), esc(r["covariates"]), r["n"], r["events"],
        fmt(r["hr_uni"], 2),
        f'{fmt(r["hr_adj"],2)} ({fmt(r["lo_adj"],2)}--{fmt(r["hi_adj"],2)})',
        pfmt(r["p_adj"]), f'{fmt(r["delta_pct"],1)}\\%',
    ])
PARTS.append(longtable(
    "扩展协变量后的多因素 Cox 结果",
    "tab:multivariate",
    "p{2.2cm}p{2.6cm}cccp{2.4cm}cc",
    ["队列", "校正协变量", "n", "事件", "单因素 HR",
     "校正后 HR (95\\%CI)", "校正后 $p$", "效应衰减"],
    body,
    note="效应衰减按 $\\log$HR 的相对变化百分比计算。GSE39582 的 \\texttt{cit\\_f} 为 CIT 分子分型，"
         "加入后效应衰减 56\\% 且 $p$ 值升至 0.097，提示风险评分与 CIT 分型存在明显共线性", fontsize="footnotesize", colsep=4))

# ============================================================ 表：单细胞置换
rows = read("scrna_permutation.csv")
body = []
for r in rows:
    body.append([
        esc(r["statistic"]), esc(r["contrast"]),
        fmt(r["observed"], 3), fmt(r["null_mean"], 5), fmt(r["null_sd"], 4),
        fmt(r["z"], 1), pfmt(r["p_permutation"]),
        f'{r["n_a"]}/{r["n_b"]}',
    ])
PARTS.append(longtable(
    "单细胞数据的置换检验（GSE132465 子集，2000 次细胞标签重排）",
    "tab:scperm",
    "p{3.2cm}p{3.0cm}cccccc",
    ["统计量", "对比", "观测值", "零分布均值", "零分布 SD", "$z$", "$p$", "$n_a/n_b$"],
    body,
    note="零分布由 2000 次随机重排细胞标签得到，$p=(k+1)/(B+1)$，故最小可报告值为 "
         "$1/2001=5.0\\times10^{-4}$", fontsize="footnotesize", colsep=4))

rows = read("scrna_subtype_contrast.csv")
body = []
for r in rows:
    body.append([esc(r["contrast"]), f'{r["n_a"]}/{r["n_b"]}',
                 fmt(r["median_a"], 3), fmt(r["median_b"], 3),
                 pfmt(r["wilcox_p"]), fmt(r["cliffs_delta"], 3)])
PARTS.append(longtable(
    "髓系亚型内部的程序活性对比",
    "tab:subtype",
    "p{4.0cm}ccccc",
    ["对比", "n\\_a/n\\_b", "中位数 a", "中位数 b", "Wilcoxon $p$", "Cliff's $\\delta$"],
    body,
    note="Cliff's $\\delta$ 为效应量，$|\\delta|<0.147$ 通常被视为可忽略量级的差异。"))

# ============================================================ 表：CellChat
rows = read("cellchat_rank_specificity.csv")
body = []
for r in rows[:10]:
    body.append([esc(r["pair"]), fmt(r["agg_sum"], 5), fmt(r["agg_mean"], 6),
                 r["rank_sum"], r["n_total_pairs"], f'{fmt(r["percentile"],1)}\\%'])
PARTS.append(longtable(
    "CellChat 通讯强度的全局排位（CAF $\\to$ 髓系）",
    "tab:cellchat_rank",
    "p{5.0cm}ccccc",
    ["source $\\to$ target", "聚合强度（求和）", "聚合强度（均值）",
     "排位", "总组合数", "百分位"],
    body,
    note="在 $20\\times19=380$ 个 source $\\to$ target 组合中按聚合强度降序排位，"
         "此处列出 CAF $\\to$ 髓系的前 10 个组合", colsep=4))

rows = read("cellchat_null_summary.csv")
d = {r["statistic"]: r["value"] for r in rows}
body = [
    ["观测统计量", fmt(d["observed"], 4)],
    ["零分布均值", fmt(d["null_mean"], 4)],
    ["零分布标准差", fmt(d["null_sd"], 4)],
    ["$z$ 值", fmt(d["z"], 2)],
    ["置换 $p$ 值", pfmt(d["p_permutation"])],
    ["重复次数 $B$", d["B"]],
    ["纳入的配体-受体对数", d["n_LR_pairs"]],
    ["零分布 95\\% 区间", f'{fmt(d["null_2.5pct"],4)} -- {fmt(d["null_97.5pct"],4)}'],
]
PARTS.append(longtable(
    "表达匹配随机基因对零模型的检验结果",
    "tab:cellchat_null",
    "p{5.0cm}c",
    ["项目", "取值"],
    body,
    note="统计量定义为 $\\frac{1}{P}\\sum_{p}\\left(\\bar{x}^{\\mathrm{CAF}}_{\\mathrm{lig},p}"
         "\\times\\bar{x}^{\\mathrm{mye}}_{\\mathrm{rec},p}\\right)$；零模型把配体与受体分别替换为"
         "同表达十分位内的随机信号基因，重复 1000 次。"))

rows = read("cellchat_lr_top30.csv")
body = []
for r in rows[:15]:
    body.append([
        f'{esc(r["ligand"])} -- {esc(r["receptor"])}',
        fmt(r["mean_ligand_in_CAF"], 2), fmt(r["mean_receptor_in_myeloid"], 2),
        fmt(r["product"], 2), f'{fmt(r["contribution_pct"],1)}\\%',
        fmt(r["z_pair"], 2), pfmt(r["p_pair"]),
    ])
PARTS.append(longtable(
    "CAF $\\to$ 髓系通讯的主要配体-受体对（前 15）",
    "tab:lr",
    "p{3.0cm}cccccc",
    ["配体 -- 受体", "CAF 配体均值", "髓系受体均值", "乘积", "贡献占比",
     "$z$", "$p$"],
    body,
    note="$z$ 与 $p$ 来自该配对各自的表达匹配零模型（$B=1000$）。"
         "贡献占比按乘积占全部 1115 个可测单基因配对总和的比例计算。"))

# ============================================================ 表：空转
rows = read("spatial_pooled_both_signatures.csv")
body = []
for r in rows:
    body.append([
        esc(r["cohort"]),
        "13 基因程序（报告口径）" if r["signature"] == "sig13" else "17 基因扰动签名（落盘口径）",
        r["n"], fmt(r["pooled_rho_sig_mye"], 3), fmt(r["pooled_rho_sig_caf"], 3),
        fmt(r["pooled_rho_caf_mye"], 3),
    ])
PARTS.append(longtable(
    "空间转录组：两种签名口径下的合并 Spearman 相关",
    "tab:spatial_pooled",
    "p{2.4cm}p{4.0cm}cccc",
    ["队列", "签名口径", "spot 数", "$\\rho$(签名, 髓系)", "$\\rho$(签名, CAF)",
     "$\\rho$(CAF, 髓系)"],
    body,
    note="13 基因口径下 GSE267401 的 $\\rho$(签名, CAF)$=0.4347$，与原报告正文的 0.434 吻合，"
         "确认报告正文使用的是 13 基因程序；而落盘文件中的 \\texttt{sig} 列实为 17 基因扰动签名", colsep=4))

rows = read("spatial_permutation_both_signatures.csv")
rows13 = [r for r in rows if r["signature"].startswith("13-gene")]
body = []
for r in rows13:
    body.append([
        esc(r["cohort"]), esc(r["slice"]), r["n_spots"],
        fmt(r["rho_sig_mye"], 3), fmt(r["rho_sig_caf"], 3),
        f'{fmt(r["diff_mye_minus_caf"],3)}\\newline ({fmt(r["diff_lo"],3)}, {fmt(r["diff_hi"],3)})',
        fmt(r["partial_sig_mye_given_caf"], 3), fmt(r["partial_sig_caf_given_mye"], 3),
    ])
PARTS.append(longtable(
    "空间转录组逐切片结果（13 基因程序口径）",
    "tab:spatial_slice",
    "p{1.8cm}p{1.9cm}ccccp{2.7cm}cc",
    ["队列", "切片", "n", "$\\rho$(签名,髓系)", "$\\rho$(签名,CAF)",
     "差值 (95\\%CI)", "偏相关$_{\\text{髓}|\\text{CAF}}$",
     "偏相关$_{\\text{CAF}|\\text{髓}}$"],
    body,
    note="差值为 $\\rho$(签名,髓系)$-\\rho$(签名,CAF)，括号内为其 95\\% 置信区间，"
         "配对置换检验 $p$ 除 GSE226997 的 P2（$p=0.983$）外均为 $5.0\\times10^{-4}$。"
         "8 张切片中 7 张差值为正；P2（两者均无关联）与 P3（反转为 CAF 占优）为异例。",
    fontsize="scriptsize", colsep=2))

# ============================================================ 表：通路
rows = read("riskgroup_pathways.csv")
tc = [r for r in rows if r["cohort"] == "TCGA"]
up = sorted([r for r in tc if float(r["NES"]) > 0], key=lambda r: float(r["padj"]))[:12]
dn = sorted([r for r in tc if float(r["NES"]) < 0], key=lambda r: float(r["padj"]))[:12]
body = []
for r in up:
    body.append([esc(r["pathway"]), r["size"], fmt(r["NES"], 2), pfmt(r["padj"]), "高危组上调"])
for r in dn:
    body.append([esc(r["pathway"]), r["size"], fmt(r["NES"], 2), pfmt(r["padj"]), "低危组上调"])
PARTS.append(longtable(
    "TCGA-COADREAD 高危组 vs 低危组的通路富集（fgsea，按校正 $p$ 排序）",
    "tab:pathway",
    "p{7.6cm}ccc p{2.0cm}",
    ["通路", "基因集大小", "NES", "校正 $p$", "方向"],
    body,
    note="排序统计量为 limma 的 $t$ 值（模型校正 age + stage）。"
         "上 12 行为高危组上调，下 12 行为低危组上调。", colsep=4))

rows = read("riskgroup_pathway_crosscohort.csv")
focus = ["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING",
         "HALLMARK_ANGIOGENESIS", "HALLMARK_HYPOXIA", "HALLMARK_INFLAMMATORY_RESPONSE",
         "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
         "HALLMARK_COMPLEMENT", "HALLMARK_COAGULATION"]
body = []
for name in focus:
    for r in rows:
        if r["pathway"] == name:
            body.append([esc(name.replace("HALLMARK_", "")),
                         fmt(r["nes_TCGA"], 2), pfmt(r["padj_TCGA"]),
                         fmt(r["nes_GSE39582"], 2), pfmt(r["padj_GSE39582"]),
                         "同向" if float(r["nes_TCGA"]) * float(r["nes_GSE39582"]) > 0 else "反向"])
PARTS.append(longtable(
    "基质 / 炎症相关通路在两个队列中的 NES 对比",
    "tab:focus",
    "p{5.0cm}ccccc",
    ["通路", "TCGA NES", "校正 $p$", "GSE39582 NES", "校正 $p$", "方向"],
    body,
    note="九条通路在两个队列中全部同向，且均为正值（高危组上调）。"))

rows = read("riskgroup_signature_activity.csv")
body = []
for r in rows[:12]:
    body.append([esc(r["pathway"].replace("HALLMARK_", "")),
                 fmt(r["High"], 3), fmt(r["Low"], 3), fmt(r["delta"], 3),
                 pfmt(r["p"]), pfmt(r["adj_p"])])
PARTS.append(longtable(
    "TCGA 高危组与低危组的通路活性比较（基因 z-score 均值法）",
    "tab:activity",
    "p{6.4cm}ccccc",
    ["通路", "高危组", "低危组", "差值", "$p$", "校正 $p$"],
    body,
    note="对每个基因在样本间做 $z$-score 标准化后，取基因集内均值作为单样本活性，"
         "再用 Wilcoxon 秩和检验比较两组。此处列出高危组活性最高的 12 条通路。"))

# ============================================================ 表：跨队列汇总
rows = read("riskgroup_crosscohort_summary.csv")
d = {r["metric"]: r["value"] for r in rows}
body = [
    ["跨队列 NES Spearman $\\rho$", fmt(d["cross_cohort_nes_spearman_rho"], 3)],
    ["相关检验 $p$", pfmt(d["cross_cohort_nes_p"])],
    ["纳入比较的通路数", d["n_pathways_compared"]],
    ["方向一致的通路数", d["n_concordant"]],
    ["方向一致占比", f'{fmt(d["pct_concordant"],1)}\\%'],
]
PARTS.append(longtable(
    "风险分组通路信号的跨队列一致性汇总",
    "tab:cross",
    "p{6.0cm}c",
    ["指标", "取值"],
    body))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
import re
n = 0
for content in PARTS:
    m = re.search(r"\\label\{tab:([^}]*)\}", content)
    key = m.group(1) if m else f"t{n:02d}"
    fn = os.path.join(os.path.dirname(OUT), f"tab_{key}.tex")
    with open(fn, "w", encoding="utf-8") as f:
        f.write("% 自动生成，请勿手工编辑。数值来自 02_分析结果/deepen/ 下的 CSV。\n")
        f.write(content + "\n")
    n += 1
print("written", n, "table files into", os.path.dirname(OUT))
