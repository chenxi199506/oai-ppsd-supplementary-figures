# -*- coding: utf-8 -*-
"""生成补充图版 S1–S6（统一编号、180 mm、600 dpi、PNG + 矢量 PDF）。

每张图版一个独立驱动脚本，按顺序在子进程里跑，互不污染模块状态：

  draw_S1.py  研究设计与人群筛选流程（流程图，原 Figure 1 位图重绘为矢量）
  draw_S2.py  协变量平衡账（原正文 Fig 1d + Fig 2d 两个热条合并为一张热图）
  draw_S3.py  两种基质的靶向 PRM 验证（a–h，原 S2 滑液 + 原 S6 滑膜合并）
  draw_S4.py  稳健性 + 三项技术对照（a–e）
  draw_S5.py  MEA 供体级一致性与回补（a–c，原 S7 只留 A/E/F）
  draw_S6.py  PCDH9 遗传背景（a–d，按作者要求保留）

本文件就住在 src/originals/ 里，脚本与成品同源：直接 `python3 build_supp.py`
即从当前目录出图、把 FigureS* 拷进 ../figures/，再生成三份文档。

用法：cd src/originals && python3 build_supp.py
"""
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent.parent
FIG = OUT / "figures"
ARCH = OUT / "归档"

PKG = Path("project_root/"
           "01_正文图/20261107_投稿图版可编辑包/最小复现总包/"
           "24_图版01-07_最小复现总包_20261106")
AUDIT = Path("audit_root/editable")
PREV = Path("project_root/"
            "99_归档/20261108_补充图版_S1-S7_旧版第3版")
REVIEW = Path("project_root/"
              "01_正文图/20261108_NC正文图审阅与图注")

SCRIPTS = ["draw_S1.py", "draw_S2.py", "draw_S3.py", "draw_S4.py",
           "draw_S5.py", "draw_S6.py"]

SOURCES = [
    (PKG / "07_SynovialFluid" / "scripts" / "make_figures_portable.py",
     "make_figures_portable.py"),
    (PKG / "07_SynovialFluid" / "scripts" / "journal_style.py", "journal_style.py"),
    (PKG / "06_Experiments" / "scripts" / "make_figure_S1_portable.py",
     "make_figure_S1_portable.py"),
    (AUDIT / "fig5.py", "fig5.py"),
    (AUDIT / "style.py", "style.py"),
    (AUDIT / "data.py", "data.py"),
]

# 移出本套的图版：整份归档，不删除
ARCHIVE = [
    (PREV / "figures" / "FigureS4_cutpoint_and_outcome_sensitivity.png",
     "旧 S4_切点与增量AUC.png"),
    (PREV / "figures" / "FigureS4_cutpoint_and_outcome_sensitivity.pdf",
     "旧 S4_切点与增量AUC.pdf"),
    (REVIEW / "整改产物/FigS6/FigureS6_covariate_balance.png",
     "旧 FigureS6_合并前两联热图.png"),
    (REVIEW / "整改产物/FigS6/FigureS6_covariate_balance.pdf",
     "旧 FigureS6_合并前两联热图.pdf"),
]


def build():
    for d in (FIG, ARCH):
        d.mkdir(parents=True, exist_ok=True)
    for s in SCRIPTS:
        print(f"\n=== {s} ===", flush=True)
        subprocess.run([sys.executable, str(HERE / s)], check=True, cwd=str(HERE))

    for p in sorted(HERE.glob("FigureS*")):
        if p.is_file():
            shutil.copy2(p, FIG / p.name)

    # 投稿系统对 alpha 通道挑剔（正文 Fig 6 就因为 RGBA 被点过一次）：
    # 六张 PNG 统一拍平到 RGB，画布本来就是纯白，转换无信息损失。
    from PIL import Image

    for p in sorted(FIG.glob("*.png")):
        im = Image.open(p)
        if im.mode != "RGB":
            im.convert("RGB").save(p, dpi=(600, 600))

    # 依赖副本：只刷新，不动正在用的脚本（脚本本身就在原地）
    for src, name in SOURCES:
        if src.exists() and src.resolve() != (HERE / name).resolve():
            shutil.copy2(src, HERE / name)

    for old, new in ARCHIVE:
        if old.exists():
            shutil.copy2(old, ARCH / new)

    src = PREV / "src" / "originals" / "draw_S4.py"
    if src.exists():
        shutil.copy2(src, ARCH / "旧 S4_切点与增量AUC.py")

    # make_docs.py is an internal legend/preview generator and is not part of
    # this public code archive; skip it when the file is absent.
    docs = HERE / "make_docs.py"
    if docs.exists():
        subprocess.run([sys.executable, str(docs)], check=True, cwd=str(HERE))
    else:
        print("skip docs: make_docs.py not distributed", flush=True)
    print("\nALL DONE ->", OUT, flush=True)


if __name__ == "__main__":
    build()
