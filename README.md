# Supplementary Figures S1–S6 — Plotting Code

Plotting source code for the supplementary figure plates S1–S6 of the
OAI–PPSD / PCDH9 manuscript (180 mm wide, 600 dpi PNG + vector PDF).

**This repository contains plotting code only.** Source datasets, rendered
figures, figure legends and manuscript text are not included.

> The associated manuscript has not yet been published. No license is granted;
> all rights are reserved. Please contact the authors before any reuse.

## Contents

| Path | Role |
|---|---|
| `src/originals/draw_S1.py` … `draw_S6.py` | One driver script per supplementary plate |
| `src/originals/build_supp.py` | Runs the six drivers in isolated subprocesses, flattens PNGs to RGB, collects outputs |
| `src/originals/journal_style.py` | Shared typography, geometry and layout helpers (Arial, 180 mm, PDF Type 42) |
| `src/originals/make_figures_portable.py`, `make_figure_S1_portable.py` | Panel-level drawing functions reused by S1 / S3 / S4 |
| `src/originals/fig5.py`, `data.py`, `style.py` | Donor-level MEA panels, their CSV readers, and the palette |
| `archive/旧 S4_切点与增量AUC.py` | Retired cut-point / ΔAUC plate driver, kept for reference |

`make_docs.py` (an internal legend and preview-page generator) is deliberately
**not** distributed. `build_supp.py` skips that step when the file is absent;
figure rendering is unaffected.

## Requirements

- Python 3.9+
- `matplotlib`, `numpy`, `scipy`, `Pillow`
- Arial (falls back to DejaVu Sans when unavailable)

## Configuration before running

The scripts were written on the authors' workstation and originally referenced
absolute local paths. In this public archive those prefixes are replaced by the
placeholders `project_root/` and `audit_root/`, e.g.

```python
SCRIPTS = Path("project_root/"
               "01_正文图/...")
```

Point them at your own copy of the figure-source tree. The input data are
**not** bundled here and must be supplied separately: `data.json`, `stats.json`,
`plot_data/*.csv` (e.g. `fig5_well_level_84.csv`, `fig5_donor_trajectory.csv`,
`figS1_tidy.csv`), the `PRM_校准曲线.csv` calibration table, the covariate-balance
tables, and the GWAS / MR source tables referenced by `draw_S6.py`. Some panels
(e.g. S2) carry their summary values inline in the script.

Because the data are not bundled, treat this repository as a **source-code
record**, not a turnkey reproduction package.

## Reproducing

```bash
cd src/originals
python3 build_supp.py
```

Or render a single plate: `python3 draw_S2.py`.
