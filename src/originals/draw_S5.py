# -*- coding: utf-8 -*-
"""Figure S5 — donor-level consistency and rescue of neuronal responses.

Only the three panels that are not already in the main text are kept:

  a  the design (five donors, four fibroblast states, one read-out per well)
  b  per-donor trajectories, each donor normalised to its own NT-CM
  c  per-donor effect size and rescue recovery

The mean firing rate, the inter-spike interval and the network-burst count are
shown in main Figure 6 F–H and are deliberately not repeated here.  What the
main figure cannot show is whether the effect depends on one or two donors, and
that is the whole point of b and c.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl

AUDIT = Path("audit_root/editable")
SCRIPTS = Path("project_root/"
               "01_正文图/20261107_投稿图版可编辑包/最小复现总包/"
               "24_图版01-07_最小复现总包_20261106/07_SynovialFluid/scripts")
sys.path.insert(0, str(AUDIT))
sys.path.insert(0, str(SCRIPTS))
import fig5                      # noqa: E402
import style as S                # noqa: E402
import journal_style as J        # noqa: E402

OUT = Path(__file__).resolve().parent
W, H = 180.0, 107.0

# 面板字母重排：源图的 A/E/F -> a/b/c
_REMAP = {"A": "a", "E": "b", "F": "c"}
_head = S.panel_head


def _panel_head(ax, letter, *a, **k):
    return _head(ax, _REMAP.get(letter, letter), *a, **k)


S.panel_head = _panel_head

fig = J.plate(W, H)


def yf(mm):
    return mm / H


def hf(mm):
    return mm / H


TOP = H - 14.0
ax_a = fig.add_axes([0.018, yf(TOP - 26.0), 0.964, hf(26.0)])
ax_b = fig.add_axes([0.072, yf(22.0), 0.352, hf(31.0)])
ax_c = fig.add_axes([0.588, yf(22.0), 0.366, hf(31.0)])

fig5.panel_a(ax_a)
fig5.panel_e(ax_b)
fig5.panel_f(ax_c)

# S1–S4 全部是左+下轴线（J.style_axis 约定），这里把 b / c 的闭合框去掉，
# 四张补充图放在一起时轴线风格一致。
for _ax in (ax_b, ax_c):
    S.strip_spines(_ax, keep=("left", "bottom"))

J.header(fig, "Figure S5",
         "Donor-level consistency and rescue of neuronal responses",
         "PCDH9 knockdown conditioned medium; five independent OA donors; the donor, not the well, "
         "is the analytic unit")
J.footer(fig,
         "Each donor's wells are averaged before plotting, so every point in b and c is one donor and the "
         "lines connect that donor's own NT-CM, PCDH9-knockdown and rescue conditions; n = 5 throughout. "
         "The mean firing rate, inter-spike interval and network-burst count are shown in main Figure 6 F–H "
         "and are not repeated here. The figure does not show that the knockdown effect is equivalent to "
         "vehicle \u2014 the design was not powered for equivalence \u2014 and the P01 recovery above 100% is a "
         "denominator artefact of a small knockdown effect in that donor, not supra-normal rescue.")

mpl.rcParams["savefig.bbox"] = "standard"
mpl.rcParams["savefig.pad_inches"] = 0.0
fig.savefig(OUT / "FigureS5_MEA_donor_level.png", dpi=600, facecolor="white")
fig.savefig(OUT / "FigureS5_MEA_donor_level.pdf", facecolor="white")
print(f"S5 done  ({W} x {H} mm)")
