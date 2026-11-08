# -*- coding: utf-8 -*-
"""Figure S6 — genetic context, four panels (a b c d).

a  FinnGen R12 look-up at rs1413573 across 12 pain and joint endpoints
b  regional association of the PCDH9 locus for pain pleiotropy
c  tissue-specific IVW estimates across 15 GTEx tissues
d  per-instrument Wald ratios for the fibroblast instrument set

Everything is read from the released reproduction tables:
  gwas_reproduced/tables/10_selected_forest_source.csv
  gwas_reproduced/tables/05_Mocci_region_plot_source.tsv
  source/53_图版_顶刊交付_20261024/mr_independent_summary.json
No statistic is recomputed by hand.
"""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys

SCRIPTS = Path("project_root/"
               "01_正文图/20261107_投稿图版可编辑包/最小复现总包/"
               "24_图版01-07_最小复现总包_20261106/07_SynovialFluid/scripts")
sys.path.insert(0, str(SCRIPTS))
import journal_style as J   # noqa: E402

REPRO = Path("project_root/"
             "01_正文图/20261107_投稿图版可编辑包/最小复现总包/"
             "24_图版01-07_最小复现总包_20261106/04_Tissue/"
             "Figure4_data_package_20261106/data/_external/"
             "PCDH9_integrated_reproduction_20261024")
FOREST = REPRO / "gwas_reproduced/tables/10_selected_forest_source.csv"
REGION = REPRO / "gwas_reproduced/tables/05_Mocci_region_plot_source.tsv"
MRJSON = REPRO / "source/53_图版_顶刊交付_20261024/mr_independent_summary.json"

OUT = Path(__file__).resolve().parent
H_MM = 192.0

# $−log_{10}$, $P$, $^{-8}$ 若走默认 dejavusans，会把 DejaVu 塞进 PDF。
# 全图只用 Arial，包括数学斜体。
matplotlib.rcParams.update({
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


# ---------------------------------------------------------------- data
rows = list(csv.DictReader(open(FOREST, encoding="utf-8-sig")))
traits = [(r["label"], float(r["OR_A"]), float(r["OR_low"]), float(r["OR_high"]),
           float(r["p_wald"]), float(r["p_Holm_85"])) for r in rows]

reg = list(csv.DictReader(open(REGION, encoding="utf-8-sig"), delimiter="\t"))
pos = np.array([float(r["base_pair_location"]) for r in reg]) / 1e6
mlp = -np.log10(np.array([float(r["p_value"]) for r in reg]))
lead = int(np.argmin([float(r["p_value"]) for r in reg]))
lead_rs = reg[lead]["rsid"]
lead_p = float(reg[lead]["p_value"])
lead_mb = float(reg[lead]["base_pair_location"]) / 1e6

mr = json.load(open(MRJSON, encoding="utf-8"))
tis = [(v["label"], v["k"], v["ivw_fe"]["beta"], v["ivw_fe"]["se"])
       for k, v in mr.items() if k.startswith("oa_") and k != "oa_pooled"]
tis.sort(key=lambda t: -t[2])
FIB = mr["oa_Cells_Cultured_fibroblasts"]
INST = FIB["instruments"]
IVW = FIB["ivw_fe"]

fig = J.plate(J.W2, H_MM)

# ---------------------------------------------------------------- a  FinnGen
axA = fig.add_axes([0.232, yf(124), 0.345, hf(46)])
y = np.arange(len(traits))[::-1]
for i, (lab, orr, lo, hi, p, holm) in enumerate(traits):
    if holm < 0.05:
        col = J.OA
    elif p < 0.05:
        col = J.HEALTHY
    else:
        col = J.GRAY
    axA.plot([lo, hi], [y[i], y[i]], color=col, lw=0.7, zorder=2)
    axA.scatter([orr], [y[i]], s=15, facecolor=col, edgecolor="white",
                linewidth=0.35, zorder=3)
    if holm < 0.05:
        axA.text(hi + 0.0035, y[i], "*", fontsize=7.4, va="center", ha="left",
                 color=J.OA)
axA.axvline(1.0, color=J.HAIRLINE, lw=0.6, ls=(0, (4, 2)), zorder=1)
axA.set_yticks(y)
axA.set_yticklabels([t[0] for t in traits], fontsize=6.2)
axA.set_xlim(0.945, 1.135)
axA.set_xticks([0.95, 1.00, 1.05, 1.10])
axA.set_ylim(-0.7, len(traits) - 0.3)
axA.set_xlabel("Odds ratio per A allele (95% CI)")
J.style_axis(axA)
J.tag(fig, axA, "a", "The pain signal is not an osteoarthritis signal",
      "FinnGen R12, rs1413573. Red, survives Holm across 85 endpoints (*); "
      "blue, nominal P < 0.05 only; grey, not significant",
      max_mm=62.0)

# ---------------------------------------------------------------- b  region
axB = fig.add_axes([0.665, yf(124), 0.310, hf(46)])
axB.scatter(pos, mlp, s=7, facecolor="#7FA3C4", edgecolor="none", alpha=0.9,
            zorder=2)
axB.axhline(-np.log10(5e-8), color=J.GRAY, lw=0.5, ls=(0, (1.6, 1.8)), zorder=1)
axB.text(67.33, -np.log10(5e-8) + 0.30, "5 \u00d7 10$^{-8}$", fontsize=5.8,
         color=J.MUTED, ha="right", va="bottom")
axB.scatter([lead_mb], [mlp[lead]], s=26, marker="D", facecolor=J.ACCENT,
            edgecolor="white", linewidth=0.5, zorder=5)
axB.annotate(f"{lead_rs}\n$P$ = {lead_p:.2e}".replace("e-09", "\u00d7 10$^{-9}$")
             .replace("e-08", "\u00d7 10$^{-8}$"),
             xy=(lead_mb, mlp[lead]), xytext=(lead_mb + 0.085, mlp[lead] + 0.02),
             fontsize=6.0, color=J.ACCENT, ha="left", va="center", linespacing=1.35,
             arrowprops=dict(arrowstyle="-", lw=0.5, color=J.ACCENT,
                             shrinkA=1, shrinkB=2))
oa_p = 0.9577579311445938
axB.scatter([lead_mb], [-np.log10(oa_p)], s=22, marker="o",
            facecolor="white", edgecolor=J.GRAY, linewidth=0.8, zorder=5)
axB.annotate(f"Knee OA, same variant\n$P$ = {oa_p:.2f}",
             xy=(lead_mb, -np.log10(oa_p)), xytext=(lead_mb + 0.10, 3.4),
             fontsize=6.0, color=J.INK, ha="left", va="center", linespacing=1.35,
             bbox=dict(facecolor="white", alpha=0.86, edgecolor="none", pad=1.4),
             arrowprops=dict(arrowstyle="-", lw=0.5, color=J.GRAY,
                             shrinkA=2, shrinkB=3))
axB.set_xlim(66.15, 67.35)
axB.set_ylim(-0.4, 10.1)
axB.set_xticks([66.2, 66.4, 66.6, 66.8, 67.0, 67.2])
axB.set_yticks([0, 2, 4, 6, 8, 10])
axB.set_xlabel("Chromosome 13 (GRCh38, Mb)")
axB.set_ylabel("\u2212log$_{10}$(P), pain pleiotropy")
J.style_axis(axB)
J.tag(fig, axB, "b", "The locus peaks for pain and not for knee OA",
      f"{len(reg)} variants, 66.20\u201367.30 Mb", max_mm=55.0)

# ---------------------------------------------------------------- c  tissue IVW
axC = fig.add_axes([0.232, yf(40), 0.345, hf(64)])
y = np.arange(len(tis))[::-1]
for i, (lab, k, b, se) in enumerate(tis):
    lo, hi = b - 1.96 * se, b + 1.96 * se
    sig = (lo > 0) or (hi < 0)
    col = (J.OA if b > 0 else J.HEALTHY) if sig else J.GRAY
    axC.plot([lo, hi], [y[i], y[i]], color=col, lw=0.7, zorder=2)
    axC.scatter([b], [y[i]], s=13, facecolor=col if sig else "white",
                edgecolor=col, linewidth=0.7, zorder=3)
axC.axvline(0, color=J.HAIRLINE, lw=0.6, zorder=1)
NICE = {"Cells_Cultured_fibroblasts": "Fibroblasts",
        "Skin_Sun_Exposed_Lower_leg": "Skin lower leg",
        "Skin_Not_Sun_Exposed_Suprapubic": "Skin not sun-exp.",
        "Muscle_Skeletal": "Muscle skeletal",
        "Nerve_Tibial": "Nerve tibial",
        "Artery_Tibial": "Artery tibial",
        "Artery_Aorta": "Artery aorta",
        "Heart_Left_Ventricle": "Heart LV",
        "Adipose_Subcutaneous": "Adipose subcut.",
        "Colon_Sigmoid": "Colon sigmoid",
        "Thyroid": "Thyroid", "Lung": "Lung", "Prostate": "Prostate",
        "Stomach": "Stomach", "Testis": "Testis"}
axC.set_yticks(y)
axC.set_yticklabels([f"{NICE.get(t[0], t[0])} (k={t[1]})" for t in tis],
                    fontsize=6.2)
axC.set_xlim(-0.095, 0.115)
axC.set_xticks([-0.075, -0.050, -0.025, 0, 0.025, 0.050, 0.075, 0.100])
axC.set_ylim(-0.7, len(tis) - 0.3)
axC.set_xlabel("IVW log OR per SD of expression (95% CI)")
J.style_axis(axC)
J.tag(fig, axC, "c", "Seven tissues carry a positive, two a negative estimate",
      "15 GTEx tissues, OA as outcome; filled, interval excludes zero",
      max_mm=62.0)

# ---------------------------------------------------------------- d  instruments
axD = fig.add_axes([0.700, yf(40), 0.275, hf(64)])
xg = np.array([i["beta_gx"] for i in INST])
yg = np.array([i["beta_gy"] for i in INST])
sx = np.array([i["se_gx"] for i in INST])
sy = np.array([i["se_gy"] for i in INST])
axD.axhline(0, color=J.HAIRLINE, lw=0.55, ls=(0, (4, 2)), zorder=1)
axD.errorbar(xg, yg, xerr=1.96 * sx, yerr=1.96 * sy, fmt="none",
             ecolor="#CFCFCF", elinewidth=0.6, capsize=0, zorder=2)
xs = np.linspace(xg.min() - 0.06, xg.max() + 0.04, 20)
axD.plot(xs, IVW["beta"] * xs, color=J.INK, lw=0.8, zorder=3)
axD.scatter(xg, yg, s=24, facecolor=J.OA, edgecolor="white", linewidth=0.5,
            zorder=4)
axD.axvline(0, color=J.HAIRLINE, lw=0.55, zorder=1)
axD.set_xlim(-0.42, 0.62)
axD.set_ylim(-0.035, 0.062)
axD.set_xticks([-0.25, 0, 0.25, 0.50])
axD.set_yticks([-0.02, 0, 0.02, 0.04, 0.06])
axD.set_xlabel("SNP\u2013expression association (SD)")
axD.set_ylabel("SNP\u2013OA association (log OR)")
J.style_axis(axD)
axD.text(0.04, 0.96,
         f"k = {IVW['k']} | IVW = {IVW['beta']:+.3f}\n"
         f"IVW P = {IVW['p']:.1e}\n"
         f"Egger intercept P = {FIB['egger']['intercept_p']:.2f}",
         transform=axD.transAxes, va="top", ha="left", fontsize=6.2,
         linespacing=1.45,
         bbox=dict(facecolor="white", alpha=0.88, edgecolor="none", pad=1.6))
J.tag(fig, axD, "d", "The fibroblast estimate rests on four variants",
      "bars, 95% CI; line, fixed-effect IVW", max_mm=43.0)

# ---------------------------------------------------------------- frame
J.header(fig, "Figure S6",
         "The genetic signal at PCDH9 is pleiotropic for pain and weak for "
         "osteoarthritis itself",
         "FinnGen R12 look-up, regional association, and two-sample Mendelian "
         "randomisation across 15 GTEx tissues")
J.footer(fig,
         "Analytic unit is the variant for panels a and b, the tissue for panel c "
         "and the instrument for panel d. FinnGen R12 effect sizes are log OR per "
         "A allele; Holm correction is across 85 pain and joint endpoints. "
         "Mendelian randomisation uses GTEx v8 cis-eQTL instruments and FinnGen "
         "knee-OA summary statistics, fixed-effect IVW with the Egger intercept as "
         "a directional-pleiotropy check. The figure does not establish that PCDH9 "
         "causes osteoarthritis: the pooled estimate is small, several tissue "
         "estimates rest on two variants, and panel b shows the same locus carrying "
         "a much stronger pain than joint signal. Panel b plots the pain series "
         "only; the knee-OA regional series is not reproduced here, and the knee-OA "
         "result at the lead variant is the single point marked in that panel.")
J.save(fig, str(OUT), "FigureS6_genetic_context", strict=True)
print("S6 done")
