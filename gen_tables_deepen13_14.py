#!/usr/bin/env python3
# 生成 deepen_13 / deepen_14 结果的两张 LaTeX 表（风格与 make_latex_tables.py 一致）。
# 输入：02_分析结果/deepen/cit_stratified.csv, cit_stratified_meta.csv, spatial_quality_check.csv
# 输出：04_LaTeX/tab_cit_strat.tex, tab_spatial_quality.tex
import csv, os

ROOT = os.path.expanduser("~/Documents/first-year/docx/8.25/肝癌，肺癌，结直肠癌/结直肠癌_CAF髓系")
DEEP = os.path.join(ROOT, "02_分析结果", "deepen")
OUT = os.path.join(ROOT, "04_LaTeX")

def esc(s):
    return str(s).replace("_", r"\_").replace("%", r"\%")

def fmt_p(p):
    if p is None: return "--"
    if p < 1e-4: return rf"$ {p:.1e} $".replace("e-0", "\\times10^{-").replace("e-", "\\times10^{-") + "}$" if False else ""
    return f"{p:.4f}"

# 数值 p：科学计数法
def p_latex(p):
    if p is None or p != p:
        return "--"
    if p < 1e-4:
        mant, exp = f"{p:.1e}".split("e")
        return f"${mant}\\times10^{{{int(exp)}}}$"
    return f"{p:.4f}"

# ---------------- 表 1：CIT 分层 ----------------
with open(os.path.join(DEEP, "cit_stratified.csv")) as f:
    rows = list(csv.DictReader(f))
with open(os.path.join(DEEP, "cit_stratified_meta.csv")) as f:
    meta = list(csv.DictReader(f))[0]

lines = [
    "% 自动生成，请勿手工编辑。数值来自 02_分析结果/deepen/cit_stratified*.csv。",
    r"\begin{table}[htbp]", r"\centering", r"\small",
    r"\begin{threeparttable}",
    r"\caption{GSE39582 按 CIT 分子亚型分层的风险评分内部验证}",
    r"\label{tab:cit_strat}",
    r"\begin{tabular}{lcccccc}", r"\toprule",
    r"CIT 亚型 & $n$ & 事件数 & HR (95\%CI) & $p$ & 高低危 HR & KM $p$ \\",
    r"\midrule",
]
for r in rows:
    st = esc(r["subtype"])
    n = r["n"]; ev = r["events"]
    hr = f"{float(r['hr']):.2f}"; lo = f"{float(r['lo']):.2f}"; hi = f"{float(r['hi']):.2f}"
    p = p_latex(float(r["p"]))
    hrh = f"{float(r['hr_high']):.2f}"; ph = p_latex(float(r["p_logrank"]))
    lines.append(f"{st} & {n} & {ev} & {hr} ({lo}--{hi}) & {p} & {hrh} & {ph} \\\\")
ilp = float(meta["interaction_lrt_p"])
lines += [r"\bottomrule", r"\end{tabular}",
    r"\begin{tablenotes}\footnotesize",
    rf"\item[] HR 为亚型内部连续风险评分（每单位）的 Cox 回归；高低危 HR 为中位数二分后的 Cox 回归；KM $p$ 为 log-rank 检验。",
    rf"\item[] 风险评分 $\times$ CIT 亚型交互检验（似然比）：$p = {ilp:.3f}$。",
    r"\end{tablenotes}", r"\end{threeparttable}", r"\end{table}"]
with open(os.path.join(OUT, "tab_cit_strat.tex"), "w") as f:
    f.write("\n".join(lines) + "\n")

# ---------------- 表 2：空间质量复核 ----------------
name_map = {"Ajou_Visium_P1": "P1", "Ajou_Visium_P2": "P2（异例）",
            "Ajou_Visium_P3": "P3（异例）", "Ajou_Visium_P4": "P4",
            "CTC17M": "CTC17M", "CTC17P": "CTC17P", "CTC21M": "CTC21M", "CTC21P": "CTC21P"}
with open(os.path.join(DEEP, "spatial_quality_check.csv")) as f:
    rows = list(csv.DictReader(f))
lines = [
    "% 自动生成，请勿手工编辑。数值来自 02_分析结果/deepen/spatial_quality_check.csv。",
    r"\begin{table}[htbp]", r"\centering", r"\footnotesize",
    r"\setlength{\tabcolsep}{4pt}",
    r"\begin{threeparttable}",
    r"\caption{空间转录组逐切片质量与构成复核}",
    r"\label{tab:spatial_quality}",
    r"\begin{tabular}{llccccc}", r"\toprule",
    r"队列 & 切片 & spot 数 & QC（\%） & 中位 UMI & 中位基因 & 肿瘤 spot（\%）\\",
    r"\midrule",
]
for r in rows:
    coh = esc(r["cohort"]); sl = esc(name_map.get(r["sample"], r["sample"]))
    sq = f"{int(r['spots_qc']):,}".replace(",", " ")
    qc = f"{float(r['qc_pass_pct']):.1f}"
    umi = f"{int(float(r['median_umi'])):,}".replace(",", " ")
    gn = f"{float(r['median_genes']):.0f}"
    tp = f"{float(r['tumor_spot_pct']):.1f}"
    lines.append(f"{coh} & {sl} & {sq} & {qc} & {umi} & {gn} & {tp} \\\\")
lines += [r"\bottomrule", r"\end{tabular}",
    r"\begin{tablenotes}\footnotesize",
    r"\item[] QC 阈值与主分析一致：spot 基因数 $>200$ 且 UMI 数 $>500$。",
    r"\item[] 肿瘤高表达 spot 为上皮标志物 EPCAM、KRT8/18/19/20、CEACAM5、VIL1 的",
    r"log 归一化平均表达高于全切片均值加一倍标准差的 spot，为构成占比的代理指标，",
    r"非病理医师 annotation。",
    r"\end{tablenotes}", r"\end{threeparttable}", r"\end{table}"]
with open(os.path.join(OUT, "tab_spatial_quality.tex"), "w") as f:
    f.write("\n".join(lines) + "\n")

print("written: tab_cit_strat.tex, tab_spatial_quality.tex")

# ---------------- 表 3：C 指数 / AUC / NRI 对比 ----------------
with open(os.path.join(DEEP, "cindex_auc_nri_cit.csv")) as f:
    rows = list(csv.DictReader(f))
name_map = {
    "riskScore": "风险评分",
    "CIT": "CIT 分型",
    "age_stage": "年龄 + 分期",
    "age_stage_CIT": "年龄 + 分期 + CIT",
    "age_stage_CIT_RS": "年龄 + 分期 + CIT + 风险评分",
}
lines = [
    "% 自动生成，请勿手工编辑。数值来自 02_分析结果/deepen/cindex_auc_nri_cit.csv。",
    r"\begin{table}[htbp]", r"\centering", r"\small",
    r"\begin{threeparttable}",
    r"\caption{GSE39582 上风险评分与 CIT 分型的预测能力对比}",
    r"\label{tab:cindex_cit}",
    r"\begin{tabular}{lcc}", r"\toprule",
    r"模型 & Harrell $C$ & 3 年 time-dependent AUC \\",
    r"\midrule",
]
for r in rows:
    lines.append(f"{name_map.get(r['model'], r['model'])} & {float(r['cindex']):.3f} & {float(r['auc_3y']):.3f} \\\\")
nri = float(rows[-1]["nri_total"]); idi = float(rows[-1]["idi"])
lines += [r"\bottomrule", r"\end{tabular}",
    r"\begin{tablenotes}\footnotesize",
    r"\item[] 全部模型基于同一批 $n=555$ 的完整病例。3 年 AUC 为 time-dependent AUC（Aalen 加权）。",
    rf"\item[] 在年龄 + 分期 + CIT 基础上加入风险评分：3 年 NRI $= {nri:.3f}$，IDI $= {idi:.4f}$。",
    r"\end{tablenotes}", r"\end{threeparttable}", r"\end{table}"]
with open(os.path.join(OUT, "tab_cindex_cit.tex"), "w") as f:
    f.write("\n".join(lines) + "\n")
print("written: tab_cindex_cit.tex")
