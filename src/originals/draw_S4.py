# -*- coding: utf-8 -*-
"""Merged robustness plate (5 panels) = old S2 + old S3.

a  nine prespecified perturbations of the primary correlation
b  five adjustment sets (rank partial correlation)
c  on-column loading versus abundance
d  per-microlitre versus per-microgram normalisation
e  RNA integrity and sequencing batch versus the programme score

Data and statistics are read from the same data.json / stats.json as the rest of
the synovial-fluid set; nothing is recomputed by hand.
"""
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np

SCRIPTS = Path("project_root/"
               "01_正文图/20261107_投稿图版可编辑包/最小复现总包/"
               "24_图版01-07_最小复现总包_20261106/07_SynovialFluid/scripts")
sys.path.insert(0, str(SCRIPTS))

import make_figures_portable as m   # noqa: E402
import journal_style as J           # noqa: E402

S, PR, NP = m.S, m.PR, m.NP
H_MM = 158.0
OUT = Path(__file__).resolve().parent

# $\rho$ / $R^2$ 若走默认 dejavusans，会把 DejaVu 塞进 PDF（混排字重）。
# 全图只用 Arial，包括数学斜体。
mpl.rcParams.update({
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial",
    "mathtext.tt": "Arial",
    "mathtext.cal": "Arial",
    "mathtext.default": "regular",
})


def yf(mm):
    return mm / H_MM


def hf(mm):
    return mm / H_MM


fig = J.plate(J.W2, H_MM)
PW = J.plate_mm(fig)[0]

# a 与 b 并排后每块只剩约 51 / 62 mm 画幅，左侧标签栏必须压到 ~24 mm；
# 这里给中文标签一条固定短译，不用 m.en() 的长句，否则标签会吃掉半个面板。
SHORT = {
    "全部可配对病例（主分析）": "All pairable (primary)",
    "排除 < LLOQ": "Excl. < LLOQ",
    "仅无血液污染": "No blood only",
    "仅单次冻融": "Single freeze\u2013thaw",
    "排除近期关节腔注射": "No recent injection",
    "仅测序批次 B1": "Batch B1 only",
    "仅 KL IV": "KL IV only",
    "仅肽 P1 定量": "P1 only",
    "仅肽 P2 定量": "P2 only",
    "仅肽 P3（确认肽）": "P3 (confirmatory)",
    "年龄 + 性别 + BMI": "Age, sex, BMI",
    "滑膜炎症（Krenn）": "Krenn inflammation",
    "免疫细胞评分": "Immune-cell score",
    "泛 ECM 评分": "Global ECM score",
    "滑液血液污染分级": "Fluid blood grade",
}

# ---------------------------------------------------------------- a  sensitivity
axA = fig.add_axes([0.168, yf(86), 0.287, hf(48)])
rows = [r for r in S["sens"] if not r.get("na")]
labels = [SHORT.get(r["label"], m.en(r["label"])) for r in rows]
vals = [r["rho"] for r in rows]
los = [r["ci"][0] for r in rows]
his = [r["ci"][1] for r in rows]
y = np.arange(len(rows))[::-1]
for i in range(len(rows)):
    col = J.OA if i == 0 else J.GRAY
    axA.plot([los[i], his[i]], [y[i], y[i]], color=col, lw=0.7, zorder=2)
    axA.scatter([vals[i]], [y[i]], s=17 if i == 0 else 12,
                facecolor=col, edgecolor="white", linewidth=0.35, zorder=3)
axA.axvline(PR["rho"], color=J.OA, lw=0.5, ls=(0, (4, 2)), zorder=1)
axA.axvline(0, color=J.HAIRLINE, lw=0.55, zorder=1)
axA.set_yticks(y)
axA.set_yticklabels(labels, fontsize=6.2)
axA.set_xlim(-0.42, 1.02)
axA.set_xticks([-0.4, 0, 0.4, 0.8])
axA.set_ylim(-0.7, len(rows) - 0.3)
axA.set_xlabel("Spearman $\\rho$ (95% CI)")
J.style_axis(axA)
J.tag(fig, axA, "a", "Every perturbation keeps the sign",
      f"{len(rows)} analyses: primary + {len(rows) - 1} perturbations",
      max_mm=46.0)

# ---------------------------------------------------------------- b  partial
axB = fig.add_axes([0.632, yf(86), 0.343, hf(48)])
parts = S["partial"]
labels = [SHORT.get(p["label"], m.en(p["label"])) for p in parts]
y = np.arange(len(parts))[::-1]
for i, p in enumerate(parts):
    axB.plot([p["ci"][0], p["ci"][1]], [y[i], y[i]], color=J.GRAY, lw=0.7, zorder=2)
    axB.scatter([p["rho"]], [y[i]], s=12, facecolor=J.GRAY,
                edgecolor="white", linewidth=0.35, zorder=3)
axB.axvline(PR["rho"], color=J.OA, lw=0.5, ls=(0, (4, 2)), zorder=1)
axB.axvline(0, color=J.HAIRLINE, lw=0.55, zorder=1)
axB.set_yticks(y)
axB.set_yticklabels(labels, fontsize=6.2)
axB.set_xlim(-0.14, 0.92)
axB.set_xticks([-0.1, 0.1, 0.3, 0.5, 0.7, 0.9])
axB.set_ylim(-0.7, len(parts) - 0.3)
axB.set_xlabel("Rank partial $\\rho$ (95% CI)")
J.style_axis(axB)
J.tag(fig, axB, "b", "Adjustment attenuates but never removes it",
      "dashed line is the unadjusted value", max_mm=58.0)

# ---------------------------------------------------------------- c  loading
NM = S["normalise"]
axC = fig.add_axes([0.072, yf(28), 0.232, hf(38)])
axC.scatter(m.TP, m.SF, s=15, facecolor=J.OA, edgecolor="white", linewidth=0.35,
            zorder=3)
m.logy(axC, (0.2, 0.5, 1, 2, 4))
axC.set_ylim(0.18, 5.2)
axC.set_xlabel("Total protein loaded (µg)")
axC.set_ylabel("Fluid PCDH9\n(fmol µL$^{-1}$)")
J.style_axis(axC)
axC.text(0.04, 0.96,
         f"$\\rho$ = {NM['sf_tp']['rho']:+.2f}\n{m.pstr(NM['sf_tp']['p'])}",
         transform=axC.transAxes, va="top", ha="left", fontsize=6.2,
         linespacing=1.45,
         bbox=dict(facecolor="white", alpha=0.86, edgecolor="none", pad=1.4))
J.tag(fig, axC, "c", "Abundance is not a loading artefact",
      "on-column mass varies with starting volume", max_mm=46.0)

# ---------------------------------------------------------------- d  normalisation
axD = fig.add_axes([0.372, yf(28), 0.180, hf(38)])
names = ["Per µL\nof fluid", "Per µg\nof protein"]
v = [PR["rho"], NM["sf_ug_rho"]["rho"]]
lo = [PR["ci"][0], NM["sf_ug_ci"][0]]
hi = [PR["ci"][1], NM["sf_ug_ci"][1]]
x = np.arange(2)
for i in range(2):
    axD.plot([x[i], x[i]], [lo[i], hi[i]], color=J.GRAY, lw=0.8, zorder=2)
    axD.scatter([x[i]], [v[i]], s=22, facecolor=J.OA, edgecolor="white",
                linewidth=0.4, zorder=3)
    axD.text(x[i] + 0.11, v[i], f"{v[i]:.2f}", fontsize=6.3, va="center", ha="left")
axD.axhline(0, color=J.HAIRLINE, lw=0.55, zorder=1)
axD.set_xticks(x)
axD.set_xticklabels(names, fontsize=6.2)
axD.set_xlim(-0.5, 1.62)
axD.set_ylim(-0.08, 0.95)
axD.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
axD.set_ylabel("$\\rho$ with CORE20 (95% CI)")
J.style_axis(axD)
J.tag(fig, axD, "d", "Normalisation leaves the endpoint unchanged",
      "same patients, two reporting units", max_mm=40.0)

# ---------------------------------------------------------------- e  RNA QC
RQ = S["rnaqc"]
axE = fig.add_axes([0.678, yf(28), 0.296, hf(38)])
for b, col, mk in (("B1", J.HEALTHY, "o"), ("B2", J.SAND, "s")):
    k = np.array([bb == b for bb in m.BATCH])
    axE.scatter(m.RIN[k], m.CORE[k], s=15, marker=mk, facecolor=col,
                edgecolor="white", linewidth=0.35, zorder=3, label=b)
axE.axhline(0, color=J.HAIRLINE, lw=0.55, zorder=1)
axE.set_xticks([6, 7, 8, 9, 10])
axE.set_xlim(5.9, 10.3)
axE.set_xlabel("RNA integrity number (RIN)")
axE.set_ylabel("CORE20 score (z)")
J.style_axis(axE)
axE.legend(loc="upper right", fontsize=6.2, handletextpad=0.4, borderpad=0.2,
           labelspacing=0.35, frameon=False)
axE.text(0.03, 0.04,
         f"$\\rho$ = {RQ['core_rin']['rho']:+.2f}\nbatch {m.pstr(RQ['batch_p'])}",
         transform=axE.transAxes, va="bottom", ha="left", fontsize=6.2,
         linespacing=1.45,
         bbox=dict(facecolor="white", alpha=0.9, edgecolor="none", pad=1.4))
J.tag(fig, axE, "e", "Neither RNA quality nor batch drives the score",
      f"B1 n = {RQ['b1_n']}, B2 n = {RQ['b2_n']}", max_mm=50.0)

# ---------------------------------------------------------------- frame
J.header(fig, "Figure S4",
         "The primary association survives every prespecified perturbation and "
         "every technical check",
         f"One cohort, {NP} paired patients; nine perturbations, five adjustment "
         "sets, three technical controls")
J.footer(fig,
         f"Analytic unit is the patient (n = {NP} paired). Intervals are bootstrap "
         "percentile intervals from 2,000 resamples; partial correlations "
         "rank-transform all variables before residualising; total protein is the "
         "on-column mass measured by BCA. The figure does not show that the "
         f"association is causal or specific to PCDH9: every analysis re-uses the "
         f"same {NP} patients, so the analyses are correlated with one another and "
         "none of them is independent evidence. Nor does it prove that joint "
         "dilution is irrelevant in general \u2014 total protein happens to be "
         "uncorrelated with CORE20 here, and a cohort in which it were correlated "
         "would need the normalised value as primary.")
J.save(fig, str(OUT), "FigureS4_robustness_and_technical_controls", strict=True)
print("S4 done")
