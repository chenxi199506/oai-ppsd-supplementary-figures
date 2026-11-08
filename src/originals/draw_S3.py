# -*- coding: utf-8 -*-
"""Figure S3 — targeted PRM validation across two matrices.

Merges what used to be two separate plates (fluid assay validation + tissue
protein quantification) into one methods figure, because both rest on the same
four-peptide inclusion list and the same instrument method.  The two matrices
are kept visually and textually separate: each block has its own header band,
its own reporting unit and its own n, and the legend says so.

  a   the shared peptide panel (one inclusion list, roles fixed before acquisition)
  b   synovial fluid  — matrix-matched calibration
  c   synovial fluid  — back-calculated accuracy
  d   synovial fluid  — agreement between the two quantifier peptides
  e   synovial fluid  — pooled-QC drift across the injection series
  f   synovial fluid  — pre-analytical recovery and carry-over
  g   both matrices   — cytoplasmic-peptide domain control
  h   synovial tissue — agreement between the same two quantifiers

Layout is computed top-down in millimetres: every panel reserves the exact
height its own tag block needs, so no tag can collide with the row above.
"""
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
from scipy import stats

PKG = Path("project_root/"
           "01_正文图/20261107_投稿图版可编辑包/最小复现总包/"
           "24_图版01-07_最小复现总包_20261106")
FLUID = PKG / "07_SynovialFluid" / "scripts"
sys.path.insert(0, str(FLUID))
import make_figures_portable as m    # noqa: E402
import journal_style as J            # noqa: E402

S, D, P, NP, NPRM = m.S, m.D, m.P, m.NP, m.NPRM
OUT = Path(__file__).resolve().parent

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

# ---------------------------------------------------------------- tissue data
TD = json.load(open(PKG / "06_Experiments" / "data" / "data.json", encoding="utf-8"))
PEP_T = {p["id"]: p for p in TD["pep"]}
IDS_T = [p["id"] for p in TD["patients"]]
NT = len(IDS_T)
PROT_T = TD["prot"]
QC_T = TD["qc"]
RT_T = TD["rt"]

P1T = np.array([PROT_T[i]["P1"] for i in IDS_T])
P2T = np.array([PROT_T[i]["P2"] for i in IDS_T])
MEAN_T = np.array([PROT_T[i]["mean"] for i in IDS_T])
R_T = stats.pearsonr(P1T, P2T)
K_T, B_T = np.polyfit(P1T, P2T, 1)
DEV_T = (P2T - P1T) / P1T * 100.0
BIAS_T = float(np.mean(P2T - P1T))
SD_T = float(np.std(P2T - P1T, ddof=1))
LOA_T = (BIAS_T - 1.96 * SD_T, BIAS_T + 1.96 * SD_T)
BIAS_T_PCT = float(np.mean(DEV_T))
SPAN_T = MEAN_T.max() / MEAN_T.min()

FLUID_C = J.HEALTHY
TISSUE_C = J.SAGE

W = 180.0

# ---------------------------------------------------------------- tag geometry
_DUMMY = J.plate(W, 120.0)
FS_TITLE, FS_NOTE = 7.3, 6.3


def tspace(title, note, max_mm, gap=1.8):
    """Exact height the tag block occupies above an axes, from journal_style.tag."""
    nt = (J.wrap_mm(_DUMMY, title, max_mm, FS_TITLE).count("\n") + 1) if title else 0
    nn = (J.wrap_mm(_DUMMY, note, max_mm, FS_NOTE).count("\n") + 1) if note else 0
    return (nt * FS_TITLE * 1.35 * 0.3528 + 0.9
            + nn * FS_NOTE * 1.35 * 0.3528 + gap)


# ---------------------------------------------------------------- panel copy
T = {}   # letter -> (title, note, max_mm)

T["a"] = ("One inclusion list serves both matrices",
          "P1/P2 quantifiers, P3 confirmation, P4 domain control; L/H against a 0.50 fmol heavy "
          "standard added after digestion", 170.0)
T["b"] = ("Calibration is linear over four orders", "1/x² weighted, seven levels in triplicate", 41.0)
T["c"] = ("Accuracy stays inside ±15% at every level", "dashed line: ±15% acceptance", 41.0)
T["d"] = ("The two quantifier peptides agree", "Bland–Altman, n = 31", 39.0)
T["e"] = ("Pooled QC holds across the run", "12 injections, one every third sample", 41.0)
T["f"] = ("Recovery clears every limit", "spike-in recovery, single batch", 41.0)
T["g"] = ("The cytoplasmic peptide stays intracellular", "0/33 fluid vs 5/5 tissue lysate", 39.0)
T["h"] = ("The same quantifiers agree in tissue",
          "level offset, not noise; per-donor −23% to −4%", 165.0)

H_TABLE, H_R1, H_R2, H_H = 22.0, 26.0, 24.0, 28.0
BAND_H = 5.6
HEADER = 11.6
XL = 7.2        # tick labels + axis title below a plotted row
GAP = 3.4

FOOT = (
    f"Analytic unit is the injection for b, c, e and f, the sample for d (n = {NP} paired fluid samples) "
    f"and the donor for h (n = {NT}). Calibration uses light/heavy area ratios with 1/x² weighting and "
    "accuracy is back-calculated from the fitted line. Fluid is reported per microlitre of neat fluid and "
    "tissue per microgram of digested protein, so the two y axes are not on the same scale and the two "
    "matrices must not be pooled. The figure does not establish clinical performance: precision and "
    "recovery are single-batch, and the absence of the cytoplasmic peptide in fluid is consistent with "
    "shedding but does not prove it. Tissue levels are apparent, not absolute, because the heavy standard "
    "is added after digestion and recovery is not corrected."
)
_foot_w = W - 2 * 0.008 * W - 2.0
_foot_lines = J.wrap_mm(_DUMMY, FOOT, _foot_w, 6.0).count("\n") + 1
FOOT_H = _foot_lines * 6.0 * 1.55 * 0.3528

H = (HEADER
     + tspace(*T["a"][:2], T["a"][2]) + H_TABLE + 3.0
     + BAND_H + 3.0
     + tspace(*T["b"][:2], T["b"][2]) + H_R1 + XL + GAP
     + tspace(*T["e"][:2], T["e"][2]) + H_R2 + XL + GAP
     + BAND_H + 2.0
     + tspace(*T["h"][:2], T["h"][2]) + H_H + XL + GAP
     + FOOT_H + 3.4)

print(f"H = {H:.1f} mm   footer {_foot_lines} lines = {FOOT_H:.1f} mm")

fig = J.plate(W, H)


def yf(mm):
    return mm / H


def hf(mm):
    return mm / H


def panel(x_mm, w_mm, h_mm, letter, y_top, xlab=0.0):
    title, note, max_mm = T[letter]
    ts = tspace(title, note, max_mm)
    ax = fig.add_axes([x_mm / W, yf(y_top - ts - h_mm), w_mm / W, hf(h_mm)])
    J.tag(fig, ax, letter, title, note, max_mm=max_mm)
    return ax, y_top - ts - h_mm - xlab


def band(y_top, title, note, colour):
    a = fig.add_axes([0.008, yf(y_top - BAND_H), 0.984, hf(BAND_H)])
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    a.axis("off")
    a.add_patch(mpl.patches.Rectangle((0, 0), 1, 1, facecolor=colour,
                                      edgecolor="none", alpha=0.16, zorder=1))
    a.plot([0, 0], [0, 1], color=colour, lw=2.4, solid_capstyle="butt", zorder=2)
    a.text(0.007, 0.5, title, fontsize=7.3, fontweight="bold", color=J.INK,
           va="center", ha="left", zorder=3)
    a.text(0.995, 0.5, note, fontsize=6.4, color=J.MUTED, va="center",
           ha="right", zorder=3)
    return y_top - BAND_H


# ================================================================ a  peptide panel
y = H - HEADER
axA, y = panel(0.012 * W, 0.976 * W, H_TABLE, "a", y)
J.no_frame(axA)

ROLE = {"P1": "primary", "P2": "primary", "P3": "confirm", "P4": "auxiliary"}
rows = []
for pk in ["P1", "P2", "P3", "P4"]:
    vals = [PROT_T[i].get(pk) for i in IDS_T]
    vals = [v for v in vals if v is not None]
    cvs = [q["cv"] * 100 for q in QC_T if q["pep"] == pk]
    lh = f"{np.mean(vals):.2f} ({min(vals):.2f}\u2013{max(vals):.2f})" if vals else "\u2014"
    cv_s = f"{np.mean(cvs):.1f}" if cvs else "\u2014"
    rows.append([pk, ROLE[pk], PEP_T[pk]["seq"],
                 f'{PEP_T[pk]["start"]}\u2013{PEP_T[pk]["end"]}',
                 f'{PEP_T[pk]["prec2"]:.4f}', f'{RT_T[pk]:.2f}', lh, cv_s])

col = ["Peptide", "Role", "Sequence", "Residues", "m/z (2+)", "RT (min)",
       "L/H, mean (range)", "Tech. CV (%)"]
tb = axA.table(cellText=rows, colLabels=col, loc="upper center", cellLoc="center",
               bbox=[0.0, 0.02, 1.0, 0.90],
               colWidths=[0.075, 0.105, 0.250, 0.090, 0.135, 0.080, 0.175, 0.090])
tb.auto_set_font_size(False)
tb.set_fontsize(6.7)
for (r_, c_), cell in tb.get_celld().items():
    cell.set_linewidth(0.5)
    cell.set_edgecolor("#CCCCCC")
    if r_ == 0:
        cell.set_facecolor(J.INK)
        cell.get_text().set_color("white")
        cell.get_text().set_fontweight("bold")
        cell.get_text().set_fontsize(6.3)
    else:
        cell.set_facecolor("white" if r_ % 2 else "#F5F6F8")
        if c_ == 0:
            cell.get_text().set_fontweight("bold")
        if c_ == 1 and cell.get_text().get_text() == "primary":
            cell.get_text().set_fontweight("bold")

# ================================================================ fluid block
y = band(y - 3.0, "Synovial fluid",
         f"matrix-matched calibration \u00b7 {NP} paired samples \u00b7 "
         f"fmol \u00b5L$^{{-1}}$ of neat fluid", FLUID_C) - 3.0

W3 = 0.262 * W
XS = (0.078 * W, 0.408 * W, 0.738 * W)
WS = (0.262 * W, 0.262 * W, 0.246 * W)

# ---------------------------------------------------------------- b  calibration
axB, _ = panel(XS[0], WS[0], H_R1, "b", y, xlab=XL)
calrows = list(csv.DictReader(open(os.path.join(m.DATA, "PRM_校准曲线.csv"),
                                   encoding="utf-8-sig")))
cal = {}
for r in calrows:
    cal.setdefault(r["Peptide"], []).append(
        (float(r["Nominal_fmol_per_uL"]), float(r["Area_ratio"])))
for pep, colr, mk in (("P1", J.OA, "o"), ("P2", J.HEALTHY, "s")):
    x = np.array([v[0] for v in cal[pep]])
    yy = np.array([v[1] for v in cal[pep]])
    axB.scatter(x, yy, s=9, marker=mk, facecolor=colr, edgecolor="white",
                linewidth=0.3, zorder=3, label=pep)
    xs = np.logspace(np.log10(0.04), np.log10(12), 40)
    k = S["cal"][pep]
    axB.plot(xs, k["slope"] * xs + k["intercept"], color=colr, lw=0.7, zorder=2)
axB.set_xscale("log")
axB.set_yscale("log")
axB.set_xlim(0.03, 14)
axB.set_ylim(0.02, 14)
axB.set_xticks([0.05, 0.2, 1, 5, 10])
axB.set_xticklabels(["0.05", "0.2", "1", "5", "10"])
axB.set_yticks([0.05, 0.2, 1, 5, 10])
axB.set_yticklabels(["0.05", "0.2", "1", "5", "10"])
axB.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
axB.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
axB.set_xlabel("Nominal (fmol µL$^{-1}$)")
axB.set_ylabel("Light/heavy area ratio")
J.style_axis(axB)
axB.text(0.04, 0.96,
         f"P1  slope {S['cal']['P1']['slope']:.3f}, $R^2$ {S['cal']['P1']['r2']:.3f}\n"
         f"P2  slope {S['cal']['P2']['slope']:.3f}, $R^2$ {S['cal']['P2']['r2']:.3f}",
         transform=axB.transAxes, va="top", ha="left", fontsize=6.2, linespacing=1.45)

# ---------------------------------------------------------------- c  accuracy
axC, _ = panel(XS[1], WS[1], H_R1, "c", y, xlab=XL)
lv = S["validation"]
conc = sorted({v["conc"] for v in lv})
for pep, colr, mk in (("P1", J.OA, "o"), ("P2", J.HEALTHY, "s")):
    xs = [v["conc"] for v in lv if v["pep"] == pep]
    ys = [v["acc"] for v in lv if v["pep"] == pep]
    axC.scatter(xs, ys, s=12, marker=mk, facecolor=colr, edgecolor="white",
                linewidth=0.3, zorder=3)
axC.axhline(100, color=J.HAIRLINE, lw=0.55, ls=(0, (4, 2)), zorder=1)
for yv in (85, 115):
    axC.axhline(yv, color=J.GRAY, lw=0.5, ls=(0, (1.5, 2)), zorder=1)
axC.set_xscale("log")
axC.set_xticks(conc)
axC.set_xticklabels([f"{c:g}" for c in conc])
axC.set_xlim(conc[0] * 0.7, conc[-1] * 1.4)
axC.set_ylim(78, 122)
axC.set_xlabel("Nominal (fmol µL$^{-1}$)")
axC.set_ylabel("Back-calculated accuracy (%)")
J.style_axis(axC)

# ---------------------------------------------------------------- d  agreement (fluid)
axD, _ = panel(XS[2], WS[2], H_R1, "d", y, xlab=XL)
P1 = np.array([D["prm"][f"{p}|P1"] for p in P])
P2 = np.array([D["prm"][f"{p}|P2"] for p in P])
diff = (P2 - P1) / P1 * 100
mean = (P1 + P2) / 2
pa = S["pep_agree"]
axD.scatter(mean, diff, s=13, facecolor=J.OA, edgecolor="white",
            linewidth=0.35, zorder=3)
axD.axhline(pa["p12_bias"], color=J.OA, lw=0.6, zorder=2)
for yv in pa["p12_loa"]:
    axD.axhline(yv, color=J.GRAY, lw=0.5, ls=(0, (4, 2)), zorder=1)
axD.set_xscale("log")
axD.set_xlim(0.4, 6)
axD.set_xticks([0.5, 1, 2, 4])
axD.set_xticklabels(["0.5", "1", "2", "4"])
axD.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
axD.set_ylim(-40, 30)
axD.set_xlabel("Mean of P1 and P2 (fmol µL$^{-1}$)")
axD.set_ylabel("P2 − P1 (% of mean)")
J.style_axis(axD, zero=False)
axD.text(0.04, 0.96,
         f"bias {pa['p12_bias']:+.1f}%\n"
         f"limits {pa['p12_loa'][0]:+.0f}% to {pa['p12_loa'][1]:+.0f}%\n"
         f"$\\rho$ = {pa['p12_rho']:.3f}",
         transform=axD.transAxes, va="top", ha="left", fontsize=6.2,
         linespacing=1.45,
         bbox=dict(facecolor="white", alpha=1.0, edgecolor="none", pad=1.8))

# ---------------------------------------------------------------- e  pool QC
y2 = y - tspace(*T["b"][:2], T["b"][2]) - H_R1 - XL - GAP
axE, _ = panel(XS[0], WS[0], H_R2, "e", y2, xlab=XL)
q = S["poolqc"]
xq = np.arange(1, len(q["series1"]) + 1)
for ser, colr, mk in ((q["series1"], J.OA, "o"), (q["series2"], J.HEALTHY, "s")):
    base = np.mean(ser)
    axE.plot(xq, (np.array(ser) / base - 1) * 100, color=colr, lw=0.7,
             marker=mk, ms=2.6, zorder=3)
axE.axhline(0, color=J.HAIRLINE, lw=0.55, ls=(0, (4, 2)), zorder=1)
axE.axhspan(-15, 15, color="#F6F6F6", zorder=0)
axE.set_xlim(0.5, len(xq) + 0.5)
axE.set_ylim(-18, 18)
axE.set_xticks([1, 4, 8, 12])
axE.set_xlabel("Pooled-QC injection")
axE.set_ylabel("Deviation from mean (%)")
J.style_axis(axE)
axE.text(0.04, 0.95, f"CV  P1 {q['cv1']:.1f}%\n      P2 {q['cv2']:.1f}%",
         transform=axE.transAxes, va="top", ha="left", fontsize=6.2,
         linespacing=1.45,
         bbox=dict(facecolor="white", alpha=0.86, edgecolor="none", pad=1.4))

# ---------------------------------------------------------------- f  recovery
axF, _ = panel(XS[1], WS[1], H_R2, "f", y2, xlab=XL)
a = S["assay"]
names = ["Digestion", "Bench 4 h", "3 freeze–thaw", "Carry-over"]
vals = [a["digest"], a["bench"], a["ft3"]]
ypos = np.arange(4)[::-1]
for i, v_ in enumerate(vals):
    yv = ypos[i]
    axF.plot([0, v_], [yv, yv], color=J.LIGHT, lw=2.0, solid_capstyle="butt", zorder=2)
    axF.scatter([v_], [yv], s=14, facecolor=J.OA, edgecolor="white",
                linewidth=0.35, zorder=4)
    axF.text(v_ + 2.5, yv, f"{v_:.1f}%", fontsize=6.2, va="center", ha="left",
             color=J.INK)
axF.scatter([a["carryover"]], [ypos[3]], s=14, facecolor=J.ACCENT,
            edgecolor="white", linewidth=0.35, zorder=4)
axF.text(a["carryover"] + 3.0, ypos[3], f"{a['carryover']:.2f}%  (limit 1.0%)",
         fontsize=6.2, va="center", ha="left", color=J.ACCENT)
axF.set_yticks(ypos)
axF.set_yticklabels(names)
axF.set_xlim(0, 108)
axF.set_ylim(-0.7, 3.6)
axF.set_xticks([0, 25, 50, 75, 100])
axF.set_xlabel("Recovery (%)")
J.style_axis(axF)

# ---------------------------------------------------------------- g  domain control
axG, _ = panel(XS[2], WS[2], H_R2, "g", y2, xlab=XL)
axG.bar([0, 1], [0, 100], width=0.46, color=[J.LIGHT, TISSUE_C], linewidth=0)
axG.set_xticks([0, 1])
axG.set_xticklabels([f"Synovial fluid\n(0/{NPRM})", "Tissue lysate\n(5/5)"], fontsize=6.2)
axG.set_xlim(-0.6, 1.6)
axG.set_ylim(0, 118)
axG.set_ylabel("Samples with P4 detected (%)")
axG.text(0, 3, "not detected", ha="center", va="bottom", fontsize=6.0, color=J.ACCENT)
axG.text(1, 103, "detected", ha="center", va="bottom", fontsize=6.0, color=J.INK)
J.style_axis(axG)

# ================================================================ tissue block
y3 = band(y2 - tspace(*T["e"][:2], T["e"][2]) - H_R2 - XL - GAP, "Synovial tissue",
          f"same inclusion list, different matrix \u00b7 {NT} OA donors \u00d7 3 digestions \u00b7 "
          f"fmol \u00b5g$^{{-1}}$ of protein", TISSUE_C) - 2.0

axH1, _ = panel(0.105 * W, 0.335 * W, H_H, "h", y3, xlab=XL)
axH1.plot([0, 1.55], [0, 1.55], ls=(0, (4, 2)), lw=0.55, color=J.GRAY, zorder=1)
xx = np.linspace(0.0, 1.55, 50)
axH1.plot(xx, K_T * xx + B_T, lw=0.7, color=J.ACCENT, zorder=2)
axH1.scatter(P1T, P2T, s=22, color=J.OA, edgecolor="white", linewidth=0.6, zorder=3)
axH1.set_xlim(0, 1.55)
axH1.set_ylim(0, 1.55)
axH1.set_xticks([0, 0.5, 1.0, 1.5])
axH1.set_yticks([0, 0.5, 1.0, 1.5])
J.style_axis(axH1, "P1 (fmol/µg)", "P2 (fmol/µg)")
axH1.text(0.06, 0.94, f"r = {R_T[0]:.3f}\nslope = {K_T:.2f}\nb = {B_T:+.2f}",
          transform=axH1.transAxes, fontsize=6.0, va="top", color=J.INK,
          linespacing=1.45)

axH2 = fig.add_axes([0.565, axH1.get_position().y0, 0.335, hf(H_H)])
axH2.axhline(BIAS_T, lw=0.7, color=J.ACCENT, zorder=2)
for _y in LOA_T:
    axH2.axhline(_y, lw=0.55, ls=(0, (4, 2)), color=J.OA, zorder=2)
axH2.scatter((P1T + P2T) / 2, P2T - P1T, s=22, color=J.OA, edgecolor="white",
             linewidth=0.6, zorder=3)
J.style_axis(axH2, "Mean of P1 and P2 (fmol/µg)", "P2 \u2212 P1 (fmol/µg)")
axH2.set_xlim(0.10, 1.45)
axH2.set_ylim(-0.195, 0.038)
axH2.set_xticks([0.2, 0.6, 1.0, 1.4])
axH2.set_yticks([-0.15, -0.10, -0.05])
axH2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
axH2.tick_params(axis="y", pad=1.5)

# ================================================================ header / footer
J.header(fig, "Figure S3",
         "Targeted PRM validation in two matrices: synovial fluid and synovial tissue",
         f"One inclusion list and one instrument method throughout; fluid n = {NP} paired samples, "
         f"tissue n = {NT} donors; the two matrices are not interchangeable")
J.footer(fig, FOOT)

mpl.rcParams["savefig.bbox"] = "standard"
mpl.rcParams["savefig.pad_inches"] = 0.0
fig.savefig(OUT / "FigureS3_PRM_two_matrices.png", dpi=600, facecolor="white")
fig.savefig(OUT / "FigureS3_PRM_two_matrices.pdf", facecolor="white")
print(f"S3 done  ({W} x {H:.1f} mm)")
