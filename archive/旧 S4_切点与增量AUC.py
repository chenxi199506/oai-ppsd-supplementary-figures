# -*- coding: utf-8 -*-
"""Figure S4 — cut-point and outcome sensitivity (source Figure S4, same content)."""
import sys
from pathlib import Path

SCRIPTS = Path("project_root/"
               "20261107_投稿图版可编辑包/最小复现总包/"
               "24_图版01-07_最小复现总包_20261106/07_SynovialFluid/scripts")
sys.path.insert(0, str(SCRIPTS))
import make_figures_portable as m   # noqa: E402

m.FIG = str(Path(__file__).resolve().parent)
m.figureS4()
print("S4 done")
