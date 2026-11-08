# -*- coding: utf-8 -*-
"""Figure S1 — study design, population selection and analytical framework.

The supplied artwork was a 1536 x 1024 raster (217 dpi at 180 mm, below the
300 dpi floor), so it is reconstructed here as vector from its own text and box
geometry.  Content, numbers and topology are unchanged; typography and the
colour key are normalised to the rest of the supplement.
"""
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
from matplotlib.patches import FancyBboxPatch, Rectangle

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

# ---------------------------------------------------------------- palette
EDGE = "#2F5D7C"      # borders, connectors, cohort bars
TINT = "#D9E5EF"      # card title strips
GREY_E = "#8C8C8C"    # exclusion branch
GREY_T = "#EDEDED"

# ---------------------------------------------------------------- metrics (mm)
W = 180.0
PT, PB = 2.8, 2.8     # padding inside the cohort frame
GAP_R = 4.4           # vertical gap, right column
GAP_EX = 5.0          # database -> excluded card
GAP_CC = 4.6          # excluded card -> complete-case card
HBAR = 8.6
HEADER_ZONE = 11.5
GAP_BAND = 4.2
MARGIN_BOT = 2.6
BAND_PAD = 1.95
BAND_LEAD = 3.05
PAD_X = 1.2

LH_T, LH_B = 2.55, 2.80
PAD_T, PAD_TB, PAD_B1, PAD_B2 = 1.30, 1.00, 1.10, 1.20
FS_T, FS_B = 6.6, 6.2
FS_H, FS_H2 = 7.4, 6.6
FS_BAND = 6.05


def card_h(n_title, n_body):
    return PAD_T + n_title * LH_T + PAD_TB + PAD_B1 + n_body * LH_B + PAD_B2


# ---------------------------------------------------------------- content
L_DB = (["OAI public database (PSAP 002)"],
        ["4,796 participants (2004\u20132006, 4 centers)",
         "Baseline and 0\u201396 m follow-up"])
L_INC = (["Inclusion criteria"],
         [("\u2022  KL \u2265 2", 3.0),
          ("\u2022  Complete baseline radiographic", 3.0),
          ("details, WOMAC pain,", 5.7),
          ("age, BMI, and sex", 5.7),
          ("\u2022  One knee per participant", 3.0)])
L_SET = (["OAI analysis set"], ["512 knees"])
L_KIDS = [
    (6.0, 31.0, ["Prediction model"], ["Pain: 512 knees", "Function: 505 knees"]),
    (33.0, 58.0, ["Baseline balance", "(extreme PPSD)"],
     ["104 knees", "(Lower 10%: $n$ = 52;", "Upper 10%: $n$ = 52)"]),
    (60.0, 85.0, ["Longitudinal", "trajectory (36\u201348 m)"],
     ["505 knees", "(High burden: $n$ = 49;", "Others: $n$ = 456)"]),
]

R_DB = (["Existing institutional TKA database"],
        ["453 knees from 453 unique patients"])
R_EX = (["Excluded from complete-case", "analysis"],
        ["$n$ = 27",
         "Missing one or more required",
         "preoperative clinical, radiographic,",
         "or postoperative outcome variables"])
R_CC = (["Complete-case analytical cohort"],
        ["426 knees",
         "Preoperative WOMAC, preoperative knee radiographs, PPSD,",
         "and postoperative WOMAC at approximately 6 and 12 months"])
R_K1 = (["Continuous PPSD analysis"], ["426 knees", "Primary prognostic analyses"])
R_K2 = (["Extreme-phenotype selection", "within the 426-knee cohort"], ["182 knees"])
R_LH = [(["Low PPSD"], ["$n$ = 91", "Lower end of the", "PPSD distribution"]),
        (["High PPSD"], ["$n$ = 91", "Upper end of the", "PPSD distribution"])]

# ---------------------------------------------------------------- heights
def _h(c):
    return card_h(len(c[0]), len(c[1]))


hA, hInc, hC = _h(L_DB), _h(L_INC), _h(L_SET)
hK = max(card_h(len(k[2]), len(k[3])) for k in L_KIDS)
hDb, hEx, hCc = _h(R_DB), _h(R_EX), _h(R_CC)
hKrow = max(card_h(len(R_K1[0]), len(R_K1[1])), card_h(len(R_K2[0]), len(R_K2[1])))
hLh = max(card_h(len(t), len(b)) for t, b in R_LH)

left_fixed = PT + hA + hInc + hC + hK + PB
right_fixed = (PT + hDb + GAP_EX + hEx + GAP_CC + hCc + GAP_R
               + hKrow + GAP_R + hLh + PB)
PANEL_H = right_fixed
GAP_L = (PANEL_H - left_fixed) / 3.0

# (starts a bullet, bold lead-in, remainder)
RUNS = [
    (True, "Upstream (OAI):", " characterize structure\u2013pain mismatch, derive PPSD, and evaluate"),
    (False, None, "its longitudinal association with pain/function."),
    (True, "Downstream (TKA):", " validate the clinical relevance of PPSD in surgical outcomes and"),
    (False, None, "select extreme-phenotype knees for tissue analyses."),
    (True, "Two independent cohorts:", " no participant overlap; not combined; not interchangeable."),
]
BAND_H = 2 * BAND_PAD + len(RUNS) * BAND_LEAD

H = HEADER_ZONE + HBAR + PANEL_H + GAP_BAND + BAND_H + MARGIN_BOT
print(f"panel {PANEL_H:.2f}  gapL {GAP_L:.2f}  band {BAND_H:.2f}  H {H:.2f}")

fig = J.plate(W, H)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")

BAD = []


def fits(text, w_mm, size, weight="normal"):
    got = J.measure_mm(None, text, size, weight)
    if got > w_mm:
        BAD.append((got, w_mm, text))
        print(f"  !! {got:6.1f} > {w_mm:6.1f} mm  {text[:56]!r}")


def card(x, y_top, w, title, body, h=None, dash=False, tint=TINT):
    n_t, n_b = len(title), len(body)
    hn = card_h(n_t, n_b)
    h = hn if h is None else max(h, hn)
    extra = (h - hn) / 2.0
    y0 = y_top - h
    ec = GREY_E if dash else EDGE
    ls = (0, (2.0, 1.3)) if dash else "-"
    box = FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0,rounding_size=1.1",
                         fc="white", ec=ec, lw=0.85, ls=ls, zorder=2,
                         joinstyle="round")
    ax.add_patch(box)

    th = PAD_T + n_t * LH_T + PAD_TB
    strip = Rectangle((x, y0 + h - th), w, th, fc=tint, ec="none", zorder=3)
    ax.add_patch(strip)
    strip.set_clip_path(box)
    ax.plot([x, x + w], [y0 + h - th] * 2, lw=0.85, color=ec, zorder=4,
            solid_capstyle="butt")

    inner = w - 2 * PAD_X
    yt = y0 + h - PAD_T
    for i, line_t in enumerate(title):
        fits(line_t, inner, FS_T, "bold")
        ax.text(x + w / 2, yt - LH_T * (i + 0.5), line_t, ha="center", va="center",
                fontsize=FS_T, fontweight="bold", color=J.INK, zorder=5)
    yb = y0 + h - th - PAD_B1 - extra
    for i, line_b in enumerate(body):
        dx = None
        if isinstance(line_b, tuple):
            line_b, dx = line_b
        fits(line_b, (w - dx - PAD_X) if dx is not None else inner + 0.4, FS_B)
        ax.text(x + dx if dx is not None else x + w / 2, yb - LH_B * (i + 0.5),
                line_b, ha="left" if dx is not None else "center", va="center",
                fontsize=FS_B, color=J.INK, zorder=5)
    return y0, h


def arrow(x0, y0, x1, y1, color=EDGE, lw=0.85, head=5.2):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", lw=lw, color=color,
                                shrinkA=0, shrinkB=0, mutation_scale=head),
                zorder=6)


def line(xs, ys, color=EDGE, lw=0.85):
    ax.plot(xs, ys, lw=lw, color=color, solid_capstyle="butt",
            solid_joinstyle="miter", zorder=6)


def bus(x_parent, y_from, xs, y_to, color=EDGE):
    y_bus = (y_from + y_to) / 2.0
    line([x_parent, x_parent], [y_from, y_bus], color=color)
    line([min(xs), max(xs)], [y_bus, y_bus], color=color)
    for x in xs:
        arrow(x, y_bus, x, y_to, color=color)


# ---------------------------------------------------------------- frame
J.header(fig, "Figure S1",
         "Study design, population selection, and analytical framework",
         "Two independent cohorts with no participant overlap: not combined, "
         "not interchangeable")

PB_TOP = H - HEADER_ZONE - HBAR
FRAME_BOT = PB_TOP - PANEL_H

L0, L1 = 3.0, 88.0
R0, R1 = 92.0, 177.0
LX, RX = 45.0, 134.0

for x0, x1 in ((L0, L1), (R0, R1)):
    ax.add_patch(FancyBboxPatch((x0, FRAME_BOT), x1 - x0, PB_TOP - FRAME_BOT,
                                boxstyle="round,pad=0,rounding_size=1.6",
                                fc="none", ec=EDGE, lw=0.9, zorder=1))

for x0, x1, t1, t2 in [
    (L0, L1, "Osteoarthritis Initiative (OAI) cohort", "(upstream clinical characterization)"),
    (R0, R1, "Institutional TKA cohort", "(tissue and mechanistic study)"),
]:
    ax.add_patch(FancyBboxPatch((x0, PB_TOP), x1 - x0, HBAR,
                                boxstyle="round,pad=0,rounding_size=1.6",
                                fc=EDGE, ec=EDGE, lw=0.9, zorder=3))
    ax.text((x0 + x1) / 2, PB_TOP + HBAR - 2.85, t1, ha="center", va="center",
            fontsize=FS_H, fontweight="bold", color="white", zorder=4)
    ax.text((x0 + x1) / 2, PB_TOP + 2.45, t2, ha="center", va="center",
            fontsize=FS_H2, color="white", zorder=4)

# ================================================================ LEFT
AX0, AX1 = 20.0, 70.0
IX0, IX1 = 5.0, 42.5

y = PB_TOP - PT
y, _ = card(AX0, y, AX1 - AX0, *L_DB)
y_db = y

y_inc = y - GAP_L
y, h_inc = card(IX0, y_inc, IX1 - IX0, *L_INC)
y_tee = y_inc - h_inc / 2.0

y_c = y - GAP_L
line([LX, LX], [y_db, y_c + 2.6])
arrow(LX, y_c + 2.6, LX, y_c)
arrow(LX, y_tee, IX1, y_tee)
y, _ = card(AX0, y_c, AX1 - AX0, *L_SET)
y_set = y

y_k = y_set - GAP_L
bus(LX, y_set, [18.5, 45.5, 72.5], y_k)
for x0, x1, t, b in L_KIDS:
    card(x0, y_k, x1 - x0, t, b, h=hK)
y_left = y_k - hK

# ================================================================ RIGHT
BX0, BX1 = 105.0, 163.0
CX0, CX1 = 100.0, 168.0
EX0, EX1 = 93.0, 130.0

y = PB_TOP - PT
y, _ = card(BX0, y, BX1 - BX0, *R_DB)
y_rdb = y

y_ex = y - GAP_EX
y, h_ex = card(EX0, y_ex, EX1 - EX0, *R_EX, dash=True, tint=GREY_T)
y_tee = y_ex - h_ex / 2.0

y_cc = y - GAP_CC
line([RX, RX], [y_rdb, y_cc + 2.6])
arrow(RX, y_cc + 2.6, RX, y_cc)
arrow(RX, y_tee, EX1, y_tee, color=GREY_E, head=4.6)
y, _ = card(CX0, y_cc, CX1 - CX0, *R_CC)
y_ccb = y

y_k2 = y_ccb - GAP_R
bus(RX, y_ccb, [112.5, 153.5], y_k2)
card(95.0, y_k2, 35.0, *R_K1, h=hKrow)
card(133.0, y_k2, 41.0, *R_K2, h=hKrow)
y_ext = y_k2 - hKrow

y_lh = y_ext - GAP_R
bus(153.5, y_ext, [141.75, 164.25], y_lh)
for (x0, x1), (t, b) in zip([(132.0, 151.5), (154.5, 174.0)], R_LH):
    card(x0, y_lh, x1 - x0, t, b)
y_right = y_lh - hLh

# ================================================================ bottom band
BAND_TOP = FRAME_BOT - GAP_BAND
band = FancyBboxPatch((L0, BAND_TOP - BAND_H), R1 - L0, BAND_H,
                      boxstyle="round,pad=0,rounding_size=1.6",
                      fc="none", ec=EDGE, lw=0.9, zorder=1)
ax.add_patch(band)
CELL = 30.0
cell = Rectangle((L0, BAND_TOP - BAND_H), CELL, BAND_H, fc=TINT, ec="none", zorder=2)
ax.add_patch(cell)
cell.set_clip_path(band)
ax.plot([L0 + CELL] * 2, [BAND_TOP - BAND_H, BAND_TOP], lw=0.9, color=EDGE, zorder=3)
ax.text(L0 + CELL / 2, BAND_TOP - BAND_H / 2, "Overall\nanalytical\nframework",
        ha="center", va="center", fontsize=6.8, fontweight="bold",
        color=J.INK, linespacing=1.26, zorder=4)

TX = L0 + CELL + 3.2
for i, (bullet, bold, rest) in enumerate(RUNS):
    yy = BAND_TOP - BAND_PAD - BAND_LEAD * (i + 0.5)
    if bullet:
        ax.text(TX, yy, "\u2022", ha="center", va="center", fontsize=FS_BAND,
                color=J.INK, zorder=4)
    xx = TX + 2.3
    if bold:
        ax.text(xx, yy, bold, ha="left", va="center", fontsize=FS_BAND,
                fontweight="bold", color=J.INK, zorder=4)
        xx += J.measure_mm(None, bold, FS_BAND, "bold")
    ax.text(xx, yy, rest, ha="left", va="center", fontsize=FS_BAND,
            color=J.INK, zorder=4)
    fits(rest, R1 - xx - 1.0, FS_BAND)

# ================================================================ export
print(f"frame {FRAME_BOT:.2f}..{PB_TOP:.2f}   left end {y_left:.2f}   "
      f"right end {y_right:.2f}   band {BAND_TOP - BAND_H:.2f}..{BAND_TOP:.2f}")
if BAD:
    print(f"WARNING: {len(BAD)} strings overflow")
mpl.rcParams["savefig.bbox"] = "standard"
mpl.rcParams["savefig.pad_inches"] = 0.0
fig.savefig(OUT / "FigureS1_study_design.png", dpi=600, facecolor="white")
fig.savefig(OUT / "FigureS1_study_design.pdf", facecolor="white")
print("S1 done")
