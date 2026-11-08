# -*- coding: utf-8 -*-
"""Figure S2 — covariate balance behind the weighted comparisons of Figures 1 and 2.

Main Figures 1 and 2 each used to carry a fourth panel holding a balance
ledger.  Both panels were squeezed into a 13-16 mm strip, their column labels
ran into the panel letter of the row below, and each carried its own colour
ramp and its own column meanings.  They are merged here into a single matrix:

    one row per covariate, one column per weighting stage,
    one 0-1 colour scale shared by the two cohorts.

That shared scale is the reason to merge rather than stack: a covariate
measured in both cohorts now sits on one line, and the OAI imbalance that
overlap weighting had to remove (baseline expected pain, |SMD| = 0.98) can be
read against the much smaller TKA imbalance it faced (anxiety, 0.24).

No value is re-estimated.  OAI values are those of 01_OAI/data/balance_values.csv;
TKA values are panel_heat() of 02_TKA/scripts/redraw_master.py recomputed from
the source workbook - the recomputation reproduces every printed value of the
original panel to two decimals.  The TKA signed differences and the z-scored
stratum means, which the merged matrix does not carry, are written to
S2_balance_values.csv next to the figure so that nothing is lost.
"""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

SCRIPTS = Path("project_root/"
               "01_正文图/20261107_投稿图版可编辑包/最小复现总包/"
               "24_图版01-07_最小复现总包_20261106/07_SynovialFluid/scripts")
sys.path.insert(0, str(SCRIPTS))
import journal_style as J   # noqa: E402

OUT = Path(__file__).resolve().parent

mpl.rcParams.update({
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial",
    "mathtext.tt": "Arial",
    "mathtext.cal": "Arial",
})

# ---------------------------------------------------------------- data
# (block, label, OAI |SMD| before, OAI |SMD| after,
#  TKA z in high stratum, TKA z in the remaining stratum,
#  TKA SMD before, TKA SMD after)
#
# OAI: 104 knees (52 lower vs 52 upper early PPSD decile), overlap weighting
# TKA: 453 knees (115 persistently high PPSD vs 338), ATT weighting
ROWS = [
    ("both", "Age",
     0.422900, 6.36e-09, -0.010624, 0.003615, -0.014166, 0.005964),
    ("both", "Female proportion",
     0.165632, 5.57e-09, 0.0, 0.0, 0.000000, 0.000000),
    ("both", "BMI",
     0.771044, 2.53e-09, 0.123165, -0.041905, 0.168739, -0.033066),
    ("both", "Medial JSN grade",
     0.205696, 1.72e-08, -0.134375, 0.045719, -0.186280, -0.003935),
    ("both", "Lateral JSN grade",
     0.113975, 2.50e-09, 0.085687, -0.029154, 0.115587, 0.005800),
    ("both", "Medial osteophyte grade",
     0.430175, 9.54e-08, 0.074310, -0.025283, 0.104744, -0.009388),
    ("oai", "KL grade 3 proportion",
     0.076750, 8.13e-09, None, None, None, None),
    ("oai", "Lateral osteophyte grade",
     0.672802, 9.93e-09, None, None, None, None),
    ("oai", "Baseline medial JSW",
     0.110711, 3.87e-09, None, None, None, None),
    ("oai", "Baseline expected pain",
     0.980525, 0.038148, None, None, None, None),
    ("tka", "KL grade",
     None, None, 0.043654, -0.014853, 0.059386, -0.016867),
    ("tka", "Comorbidity",
     None, None, 0.059171, -0.020132, 0.078390, 0.037313),
    ("tka", "ASA class",
     None, None, 0.008970, -0.003061, 0.012404, 0.019472),
    ("tka", "Anxiety",
     None, None, 0.183183, -0.062326, 0.241661, 0.010971),
    ("tka", "Depression",
     None, None, 0.114901, -0.039094, 0.153873, -0.001906),
]

BLOCKS = [("both", "Both cohorts"),
          ("oai", "OAI only"),
          ("tka", "TKA only")]

GROUPS = [
    ("OAI \u00b7 extreme PPSD deciles", "n = 104; overlap weighting", 0),
    ("TKA \u00b7 high PPSD vs rest", "n = 453; inverse-probability weighting", 2),
]

# ---------------------------------------------------------------- geometry (mm)
W = 180.0
H = 160.0
LEFT = 6.0
LAB_W = 32.0
COL_W = 27.0
COL_GAP = 1.2
GROUP_GAP = 6.0
CB_W = 5.0
CB_GAP = 4.0
RH = 6.3          # data row height
BH = 5.4          # block header height
MATRIX_TOP = H - 26.0

X = [LEFT + LAB_W + 4.0]
X.append(X[0] + COL_W + COL_GAP)
X.append(X[1] + COL_W + GROUP_GAP)
X.append(X[2] + COL_W + COL_GAP)
X_END = X[3] + COL_W
CB_X = X_END + CB_GAP

CMAP = LinearSegmentedColormap.from_list(
    "s2_balance", ["#EAF1F7", "#D8E7F1", "#BCD3E4", "#8FB2CC", "#5A87A8", "#2F5D7C"])
NORM = Normalize(0.0, 1.0)
THRESH = 0.10


def fmt(v):
    """Two decimals, with everything under 0.005 shown as <0.01."""
    a = abs(v)
    if a == 0.0:
        return "0.00"
    if a < 0.005:
        return "<0.01"
    return f"{a:.2f}"


# ---------------------------------------------------------------- canvas
fig = J.plate(W, H)
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.set_xlim(0.0, W)
ax.set_ylim(0.0, H)
J.no_frame(ax)

# ---------------------------------------------------------------- column head
for title, sub, ci in GROUPS:
    xc = X[ci] + (COL_W * 2 + COL_GAP) / 2.0
    ax.text(xc, MATRIX_TOP + 8.4, title, fontsize=6.3, color=J.INK,
            va="bottom", ha="center", fontweight="bold")
    ax.text(xc, MATRIX_TOP + 4.6, sub, fontsize=5.9, color=J.MUTED,
            va="bottom", ha="center")
    for k, lab in ((ci, "before"), (ci + 1, "after")):
        ax.text(X[k] + COL_W / 2, MATRIX_TOP + 1.2, lab, fontsize=6.1,
                color=J.INK, va="bottom", ha="center")

ax.plot([LEFT, X_END], [MATRIX_TOP, MATRIX_TOP], color=J.RULE, lw=0.7,
        solid_capstyle="butt", zorder=4)

# ---------------------------------------------------------------- matrix
y = MATRIX_TOP
for key, label in BLOCKS:
    # block header: hairline + label
    if key != "both":
        ax.plot([LEFT, X_END], [y, y], color=J.HAIRLINE, lw=0.5,
                solid_capstyle="butt", zorder=4)
    ax.text(LEFT, y - BH / 2.0, label, fontsize=5.9, color=J.MUTED,
            va="center", ha="left", fontweight="bold")
    y -= BH

    for row in ROWS:
        if row[0] != key:
            continue
        _, cov, oai_b, oai_a, _zh, _zr, tka_b, tka_a = row
        vals = [oai_b, oai_a, tka_b, tka_a]
        ax.text(LEFT, y - RH / 2.0, cov, fontsize=6.4, color=J.INK,
                va="center", ha="left")

        for k, v in enumerate(vals):
            xy = (X[k], y - RH)
            if v is None:
                ax.add_patch(Rectangle(xy, COL_W, RH, facecolor="#F0F0F0",
                                       edgecolor="white", lw=1.1, zorder=1))
                continue
            a = min(abs(v), 1.0)
            fc = CMAP(NORM(a))
            ax.add_patch(Rectangle(xy, COL_W, RH, facecolor=fc,
                                   edgecolor="white", lw=1.1, zorder=1))
            if a >= THRESH:
                ax.add_patch(Rectangle(xy, COL_W, RH, facecolor="none",
                                       edgecolor=J.OA, lw=0.7, zorder=3))
            ink = "white" if a > 0.52 else J.INK
            ax.text(X[k] + COL_W / 2.0, y - RH / 2.0, fmt(v), fontsize=6.4,
                    color=ink, va="center", ha="center", zorder=4)
        y -= RH

ax.plot([LEFT, X_END], [y, y], color=J.RULE, lw=0.7, solid_capstyle="butt",
        zorder=4)

# ---------------------------------------------------------------- colour bar
CB_H = 52.0
CB_TOP = MATRIX_TOP - 6.0
cax = fig.add_axes([CB_X / W, (CB_TOP - CB_H) / H, CB_W / W, CB_H / H])
sm = mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP)
sm.set_array([])
cb = fig.colorbar(sm, cax=cax)
cb.outline.set_visible(True)
cb.outline.set_linewidth(0.4)
cb.outline.set_edgecolor(J.HAIRLINE)
cb.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
cb.ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1.0"])
cb.ax.tick_params(length=2.0, width=0.5, labelsize=5.9, pad=1.8, colors=J.RULE)
for t in cb.ax.get_yticklabels():
    t.set_color(J.INK)
ax.text(CB_X, CB_TOP + 1.8, "|SMD|", fontsize=5.9, color=J.INK,
        va="bottom", ha="left")

# legend swatch for the outline
sw_y = CB_TOP - CB_H - 9.0
ax.add_patch(Rectangle((CB_X, sw_y), 5.0, 3.2, facecolor="#DCE8F2",
                       edgecolor=J.OA, lw=0.7, zorder=2))
ax.text(CB_X + 6.6, sw_y + 1.6, "\u2265 0.1", fontsize=5.9, color=J.INK,
        va="center", ha="left")

# ---------------------------------------------------------------- text
J.header(fig, "Figure S2",
         "Covariate balance behind the weighted comparisons in Figures 1 and 2",
         "Absolute standardised differences, one row per covariate, one 0\u20131 "
         "colour scale shared by both cohorts")
J.footer(fig,
         "Cells are absolute standardised mean differences (group mean difference over "
         "the pooled SD), shaded on one 0\u20131 scale shared by both cohorts; values below "
         "0.005 print as <0.01, cells at or above the prespecified 0.1 threshold are "
         "outlined, and blank cells mark covariates not recorded in that cohort. OAI, 104 "
         "knees in the extreme PPSD deciles (52 lower vs 52 upper), before and after "
         "overlap weighting; TKA, 453 knees (115 persistently high PPSD vs 338), before "
         "and after inverse-probability-of-treatment weighting. Seven TKA covariates "
         "removed by the |SMD| \u2264 0.25 screen are not shown, and the two cohorts "
         "parameterise Kellgren\u2013Lawrence differently, so it occupies one row in each "
         "block. Balance is marginal and covariate by covariate; it does not establish "
         "exchangeability on unmeasured variables, and the two cohorts are not pooled.")

# ---------------------------------------------------------------- side table
with open(OUT.parent.parent / "04_S2平衡账数值.csv", "w", encoding="utf-8-sig",
          newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["block", "covariate", "oai_smd_before", "oai_smd_after",
                "tka_z_high", "tka_z_rest", "tka_smd_before", "tka_smd_after"])
    for row in ROWS:
        w.writerow(["" if v is None else (repr(v) if isinstance(v, float) else v)
                    for v in row])

mpl.rcParams["savefig.bbox"] = "standard"
mpl.rcParams["savefig.pad_inches"] = 0.0
J.save(fig, str(OUT), "FigureS2_covariate_balance", strict=True)
print(f"S2 done  ({W} x {H} mm)")
