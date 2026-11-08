"""
fig5.py — Figure 5（PCDH9 缺失降低感觉神经元兴奋性）A–F 六个面板。

运行：python fig5.py
输出：out/Figure5.pdf + out/Figure5.png

数据来源：
    B/C/D  ← plot_data/fig5_well_level_84.csv（84 孔）
    E/F    ← plot_data/fig5_donor_trajectory.csv（5 供体）
    A      ← 无数据，纯示意图，坐标全在本文件
"""

from __future__ import annotations

import os

import numpy as np

import data as D
import style as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

COL = {g: S.C[D.M5_KEY[g]] for g in D.M5}

NB_CAP = 1.75          # 网络爆发轴截断值（原件设定）


# ------------------------------------------------------------------ 工具
def _sem(v):
    v = np.asarray(v, float)
    return v.std(ddof=1) / np.sqrt(len(v))


def bar_wells(ax, groups, table, metric, width=0.58, seed=1, cap=None):
    """柱=孔均值，误差线=SEM，点=单孔（横向抖动）。cap 给定时超出部分画空心三角。"""
    rng = np.random.default_rng(seed)
    for i, g in enumerate(groups):
        vals = np.array(table[g][metric], float)
        m, se = vals.mean(), _sem(vals)
        col = COL[g]
        ax.bar(i, m, width, facecolor=col, alpha=S.BAR_ALPHA,
               edgecolor=col, linewidth=0.8, zorder=1)
        ax.errorbar(i, m, yerr=se, fmt="none", ecolor=col,
                    elinewidth=0.9, capsize=2.2, capthick=0.9, zorder=3)

        off = rng.uniform(-0.16, 0.16, len(vals))
        if cap is None:
            ax.scatter(np.full(len(vals), i) + off, vals, s=10,
                       facecolor=col, edgecolor="none", alpha=S.DOT_ALPHA, zorder=4)
        else:
            lo = vals <= cap
            ax.scatter((np.full(lo.sum(), i) + off[lo]), vals[lo], s=10,
                       facecolor=col, edgecolor="none", alpha=S.DOT_ALPHA, zorder=4)
            hi = ~lo
            k = int(hi.sum())
            if k == 0:
                continue
            spread = np.linspace(-0.13, 0.13, k) if k > 1 else np.zeros(1)
            ax.scatter(i + spread, np.full(k, cap + 0.02),
                       s=17, facecolor="none", edgecolor=col, linewidth=0.8,
                       marker="^", zorder=5, clip_on=True)

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=45, ha="right", rotation_mode="anchor",
                       fontsize=S.FS["tick"])
    ax.set_xlim(-0.62, len(groups) - 0.38)


# ------------------------------------------------------------------ 面板
def panel_a(ax):
    S.blank(ax)
    S.panel_head(ax, "A", "Five donors, four medium states, one read-out per well",
                 pad=2)

    boxes = [
        ("OA synovium", "5 donors\nindependent"),
        ("Synovial fibroblasts", "P2–P3\npaired prep."),
        ("Four fibroblast states", "NT siRNA\nPCDH9 siRNA\nsiRNA + Rescue / Mock"),
        ("Conditioned medium", "24 h\n1:1, per donor"),
        ("Sensory neurons", "iPSC-derived\n2 lots / 2 plates"),
        ("MEA", "48 wells\n16 elec."),
    ]
    w, gap, y, h = 0.146, 0.023, 0.16, 0.56
    for i, (title, sub) in enumerate(boxes):
        x = 0.010 + i * (w + gap)
        S.rounded_box(ax, x, y, w, h, f"{title}\n\n{sub}",
                      fs=S.FS["small"] - 0.8, fill="#F2F5F8", edge="#9AA5B1")
        if i < len(boxes) - 1:
            S.arrow(ax, x + w + 0.002, y + h / 2, x + w + gap - 0.002, y + h / 2,
                    lw=0.7)


def panel_b(ax):
    t = D.fig5_wells()
    bar_wells(ax, D.M5, t, "firing")

    ax.set_ylim(0, 5.3)
    ax.set_yticks([0, 1, 2, 3, 4, 5])

    S.bracket(ax, 0, 1, 3.60, "***", h=0.22)
    S.bracket(ax, 1, 2, 4.02, "***", h=0.22)
    S.bracket(ax, 0, 3, 4.46, "n.s.", h=0.22, fs=S.FS["annot"])
    S.bracket(ax, 1, 3, 4.92, "**", h=0.22)

    S.style_axis(ax, "Mean firing rate (Hz)")
    S.panel_head(ax, "B", "Firing falls below vehicle",
                 "84 wells; LMM, donor random", pad=3)


def panel_c(ax):
    t = D.fig5_wells()
    bar_wells(ax, D.M5, t, "isi", seed=2)

    ax.set_ylim(0, 1120)
    ax.set_yticks([0, 200, 400, 600, 800, 1000])

    S.bracket(ax, 0, 1, 855, "***", h=38)
    S.bracket(ax, 0, 3, 935, "***", h=38)
    S.bracket(ax, 1, 2, 1015, "***", h=38)

    S.style_axis(ax, "Mean inter-spike interval (ms)")
    S.panel_head(ax, "C", "Intervals lengthen past vehicle",
                 "inverse of B; same model", pad=3)


def panel_d(ax):
    t = D.fig5_wells()
    bar_wells(ax, D.M5, t, "nburst", seed=3, cap=NB_CAP)

    ax.set_ylim(0, NB_CAP + 0.07)
    ax.set_yticks([0.0, 0.5, 1.0, 1.5])
    ax.set_yticklabels(["0.0", "0.5", "1.0", "1.5"])

    S.bracket(ax, 0, 1, 1.28, "**", h=0.06)
    S.bracket(ax, 0, 2, 1.47, "*", h=0.06)
    S.bracket(ax, 1, 2, 1.09, "n.s.", h=0.06, fs=S.FS["annot"])

    S.style_axis(ax, "Network bursts (per min)")
    S.panel_head(ax, "D", "Network bursts reach vehicle",
                 "25% of electrodes in 100 ms", pad=3)


def panel_e(ax):
    traj = D.fig5_traj()
    palette = {"P01": "#4C78A8", "P02": "#B4482F", "P03": "#6B8E6F",
               "P04": "#D9A441", "P05": "#8E6C8A"}
    xs = [0, 1, 2]
    lab = []

    for donor, nt, kd, rs in traj:
        ys = [100.0, kd, rs]
        c = palette.get(donor, S.C["faint"])
        ax.plot(xs, ys, color=c, lw=1.1, zorder=3, solid_capstyle="round")
        ax.scatter(xs, ys, s=14, facecolor=c, edgecolor="none", zorder=4)
        lab.append((rs, donor, c))

    # 右侧供体标签按 y 排开，避免重叠
    lab.sort()
    gap = 6.2
    ys_lab = [lab[0][0]]
    for k in range(1, len(lab)):
        ys_lab.append(max(lab[k][0], ys_lab[-1] + gap))
    for (y0, donor, c), yl in zip(lab, ys_lab):
        if abs(yl - y0) > 0.8:
            ax.plot([2.0, 2.12], [y0, yl], color=c, lw=0.45, zorder=2)
        ax.text(2.15, yl, donor, color=c, fontsize=S.FS["annot"],
                va="center", ha="left")

    ax.axhline(100, color=S.C["faint"], lw=0.7, ls=(0, (4, 3)), zorder=0)
    ax.set_ylim(30, 112)
    ax.set_yticks([40, 60, 80, 100])
    ax.set_xticks(xs)
    ax.set_xticklabels(["NT-CM", "KD-CM", "Rescue"], fontsize=S.FS["tick"])
    ax.set_xlim(-0.35, 2.75)

    S.style_axis(ax, "% of the donor's own NT-CM")
    S.panel_head(ax, "E", "Every donor moves the same way",
                 "5 donors; 3–4 wells each", pad=3)


def panel_f(ax):
    traj = D.fig5_traj()
    donors = [t[0] for t in traj]
    effect = [t[2] - 100.0 for t in traj]                       # KD 效应
    recovery = [(t[3] - t[2]) / (100.0 - t[2]) * 100.0 for t in traj]  # (G3-G2)/(G1-G2)

    x = np.arange(len(donors))
    w = 0.36
    ax.bar(x - w / 2 - 0.01, effect, w, color=S.C["kd"], alpha=0.82,
           edgecolor="none", label="KD effect", zorder=2)
    ax.bar(x + w / 2 + 0.01, recovery, w, color=S.C["rescue"], alpha=0.82,
           edgecolor="none", label="Rescue recovery", zorder=2)

    ax.axhline(0, color=S.C["axis"], lw=0.7, zorder=3)
    ax.axhline(100, color=S.C["faint"], lw=0.7, ls=(0, (4, 3)), zorder=1)

    ax.set_ylim(-62, 172)
    ax.set_yticks([-50, 0, 50, 100, 150])
    ax.set_xticks(x)
    ax.set_xticklabels(donors, fontsize=S.FS["tick"])
    ax.set_xlim(-0.6, len(donors) - 0.4)

    ax.legend(loc="upper right", handlelength=0.9, handletextpad=0.5,
              borderaxespad=0.1, labelspacing=0.3)

    S.style_axis(ax, "Change from NT-CM (%)")
    S.panel_head(ax, "F", "Magnitude, not direction, varies",
                 "recovery = (G3-G2)/(G1-G2); P01 >100%", pad=3)


# ------------------------------------------------------------------ 组装
def build():
    S.apply()
    fig = S.figure(180, 152)

    ax_a = fig.add_axes([0.018, 0.700, 0.964, 0.170])
    ax_b = fig.add_axes([0.060, 0.400, 0.230, 0.215])
    ax_c = fig.add_axes([0.390, 0.400, 0.230, 0.215])
    ax_d = fig.add_axes([0.715, 0.400, 0.230, 0.215])
    ax_e = fig.add_axes([0.070, 0.060, 0.360, 0.205])
    ax_f = fig.add_axes([0.585, 0.060, 0.370, 0.205])

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)
    panel_d(ax_d)
    panel_e(ax_e)
    panel_f(ax_f)

    S.stamp(fig)
    return fig


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    f = build()
    f.savefig(os.path.join(OUT, "Figure5.pdf"))
    f.savefig(os.path.join(OUT, "Figure5.png"))
    print("saved:", os.path.join(OUT, "Figure5.pdf"))
