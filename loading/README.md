# Loading Assay Quantification Scripts

Scripts for quantifying MCM helicase loading from gel images processed in ImageJ.

## Experiment Context

Loading assays measure how much MCM helicase is loaded onto chromatin under
different conditions. ORC4R is an ORC4 subunit with an R-to-A mutation that
severely impairs loading. "sofa" (suppressor of four R/A) mutations are
second-site mutations in other ORC subunits (orc1, orc3, orc4, orc5, orc6)
that partially rescue ORC4R loading activity.

## Scripts

### 260331_loading-quantification_wt-4r-supps_350mM-KGlut.R

Plots MCM loading at the 350 mM KGlut salt condition for WT, ORC4R, and all
sofa suppressors. Input data (Excel) is already normalized to WT = 100% per
replicate. Produces a bar chart with error bars and individual replicate
overlay, plus summary and full data CSVs.

- **Input**: `260331_aggregate-analysis_load_wt-4r-supps_350mM-KGlut.xlsx`
  (in the experiment's `analysis/` directory)
- **Output**: PDF bar chart, summary CSV, full data CSV (same directory)
- **Conditions**: WT, ORC4R, +1sofa, +3sofa, +4sofa, +5sofa, +6sofa
- **Replicates**: 3 biological replicates per condition
- **Metrics computed**: Percent Wildtype, percent difference (from WT and
  ORC4R), percent change (from WT and ORC4R), fold change (from ORC4R).
  See `METRICS_REFERENCE.md` for definitions and usage guidance.

### quantification_kgluttitr_wt-4r-ps.R

Plots MCM loading across three KGlut salt concentrations (250, 300, 350 mM)
for a subset of conditions (WT, ORC4R, +4sofa). Input data requires
background subtraction, input correction, and WT normalization from raw
ImageJ intensity measurements. Produces multiple plot variants; the primary
output is the faceted-by-KGlut bar chart (`faceted_by_kglut_plot`). Other
plots in the script are exploratory and kept for reference.

- **Input**: `Analysis.xlsx` (in the experiment directory, sheets 2-4)
- **Output**: PDFs and CSVs in `~/data/loading_analysis/`
- **Conditions plotted**: WT, ORC4R, +4sofa (others excluded in filtering)
- **Salt concentrations**: 250, 300, 350 mM KGlut

### inspect_excel_file.R

Utility script for inspecting Excel files before analysis. Writes sheet
names, dimensions, column names/types, and row previews to a text file.
Not experiment-specific.

## Data Location

All experiment data lives under:
```
MIT Dropbox/Luis Martinez/Lab/Experiments/Loading/
  2022_12_18 Loading Assays Repeats for publication/
```

Scripts reference this path via either the `MC_DROPBOX_PATH` environment
variable (multi-salt script) or a hardcoded WSL mount path (350 mM script).

## Dependencies

- `readxl` (350 mM script)
- `xlsx` (multi-salt script)
- `tidyverse` (both)
