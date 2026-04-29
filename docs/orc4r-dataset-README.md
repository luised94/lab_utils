# Dataset: A lethal ORC ATPase mutation is suppressed by alterations in ORC subunits and RNA Pol II transcription components

Martinez-Rodriguez L. and Bell S.P., 2026

DOI: [To be assigned upon Zenodo upload]

## Overview

This dataset contains raw images, processed data, analysis scripts, and
sequencing alignments supporting all figures and supplementary figures in
the publication. Files are named with a figure-number prefix indicating
which panel they support.

## License

This dataset is released under the Creative Commons Attribution 4.0
International License (CC BY 4.0). You are free to share and adapt this
material for any purpose, provided appropriate credit is given.

## Directory Structure

```
dataset/
|-- README.md
|-- bam_files/                     Pooled segregant analysis alignments
|   |-- Exp_N_ctrl_sorted.bam      Control pool (no phenotype)
|   |-- Exp_N_ctrl_sorted.bam.bai  BAM index
|   |-- Exp_N_hit_sorted.bam       Selected pool (screened phenotype)
|   +-- Exp_N_hit_sorted.bam.bai   BAM index
|-- fig1_*                         ORC4R suppressor screen plates and analysis
|-- fig2_*                         Protein structure visualization
|-- fig3_*                         ATPase activity, EMSA, and MCM loading
|-- fig4_*                         Genetic interaction plates
|-- fig5_*                         GAL-CDC6 system plates
|-- fig6_*                         General transcription factor allele plates
|-- sfig1_*                        Drug sensitivity plates
|-- sfig2_*                        Multiple sequence alignment pipeline
|-- sfig3_*                        EMSA and MCM loading supplementary
|-- sfig4_*                        MCM loading supplementary
+-- sfig5_*                        General transcription factor allele plates (additional)
```

## BAM Files

Nine pooled segregant analysis experiments. Each experiment contains two
samples: a "hit" pool (cells displaying the screened phenotype) and a
"ctrl" pool (cells without the phenotype, used for filtering variant
calls). All files are coordinate-sorted with BAI indices.

Reads were aligned to the S. cerevisiae reference genome R64
(GCA_000146045.2) using BWA.

| File pattern | Description |
|---|---|
| `Exp_N_hit_sorted.bam` | Sequencing of pooled cells with selected phenotype |
| `Exp_N_ctrl_sorted.bam` | Sequencing of pooled cells without phenotype (filter control) |

## Figure 1 - ORC4R Suppressor Screen

| File | Description |
|---|---|
| `fig1_2022_12_08 yCL4 WT-4R_Supps YPD 30C_d2_005.tif` | Plate image: WT, ORC4R, and suppressor strains on YPD at 30°C, day 2 |
| `fig1_2022_12_09 yCL4 WT-4R_Supps FOA 30C_d3_010.tif` | Plate image: WT, ORC4R, and suppressor strains on FOA at 30°C, day 3 |
| `fig1_orc4r-screen_01_prepare-reference.py` | Python script: download and process yeast gene coordinates from SGD (SGD_features.tab). Filters for gene-like features on nuclear chromosomes, maps chromosome identifiers to chrI-chrXVI format, and saves as TSV with provenance metadata. Output is used by the BAM analysis script for gene annotation. |
| `fig1_orc4r-screen_02_analyze-bam.py` | Python script: analyze BAM files from pooled segregant analysis for variant calling and gene annotation |

## Figure 2 - Protein Structure Visualization

| File | Description |
|---|---|
| `fig2_7mca_visualization.py` | Python script: visualize ORC-Cdc6 structure (PDB 7MCA) |
| `fig2b_session_prepare_1nh2_tfiia-tbp-dna.py` | Python script: visualize TFIIA-TBP-DNA complex (PDB 1NH2) |

## Figure 3 - ATPase Activity, EMSA, and MCM Loading

| File | Description |
|---|---|
| `fig3a_atpase-timecourse_processed-data.csv` | Processed ATPase timecourse data |
| `fig3a_plot-sofa-atpase.R` | R script: plot ATPase timecourse data |
| `fig3b1_orc5,6_20220307-162606-[Phosphor].tif` | Phosphor image: EMSA with ORC5/ORC6 suppressor mutants |
| `fig3b2_20210908-orc1,3,4-suppressor effect and ATP dependence-[Phosphor].tif` | Phosphor image: EMSA with ORC1/ORC3/ORC4 suppressors, ATP dependence |
| `fig3c_lemr load wt,4r,supps 300mM kgl 2023.04.06_20.01.22_Fl-Green.tif` | Fluorescence gel image: MCM loading on DNA (Krypton stain), WT/ORC4R/suppressors at 300 mM potassium glutamate |

## Figure 4 - Genetic Interactions

| File | Description |
|---|---|
| `fig4a1_2022_12_22_AIAy19_WT-KT-R6_Supps_YPD_30C_d2_009.tif` | Plate image: ORC1 plasmid shuffle strains (WT ORC1, K485T, R616A) with suppressors on YPD at 30�C, day 2 |
| `fig4a2_2022_12_22_AIAy19_WT-KT-R6_Supps_FOA_30C_d3_018.tif` | Plate image: same strains on FOA at 30°C, day 3 |
| `fig4b1_260410_exp-014_plates-{01-2}_day-04_time-1058am_rawscan-black.tiff` | Plate scan: orc2-1 temperature-sensitive strains, plates 1-2 |
| `fig4b2_260410_exp-014_plates-{03}_day-04_time-1058am_rawscan-black.tiff` | Plate scan: orc2-1 temperature-sensitive strains, plate 3 |
| `fig4b3_2022_12_19_O21-O1161_YPD_25C_d3_007.tif` | Plate image: orc1-161 temperature-sensitive strains on YPD at 25°C, day 3 |
| `fig4b4_2022_12_17_O21-O1161_YPD_37C_d2_006.tif` | Plate image: orc1-161 temperature-sensitive strains on YPD at 37°C, day 2 |

## Figure 5 - GAL-CDC6 System

| File | Description |
|---|---|
| `fig5a_galON_ypd_glu-gal_30_HU000_HUcontrols.tiff` | Plate scan: GAL-CDC6 strains on YPD, glucose vs galactose at 30°C, HU controls |
| `fig5b_galON_foa_glu-gal_30_HU000_HUcontrols.tiff` | Plate scan: GAL-CDC6 strains on FOA, glucose vs galactose at 30°C, HU controls |

## Figure 6 - General Transcription Factor Alleles

| File | Description |
|---|---|
| `fig6_20260310_exp-006_plates-01,02_day-05_time-211pm_rawscan-black.tiff` | Plate scan: general transcription factor allele strains, plates 1-2, day 5 |

## Supplementary Figure 1 - Drug Sensitivity

| File | Description |
|---|---|
| `sfig1a_2022_12_17_yCL4_WT-4R_Supps_FOA_02MMS_30C_d3_018.tif` | Plate image: WT/ORC4R/suppressors on FOA + 0.2 mM MMS at 30°C, day 3 |
| `sfig1b_2022_12_17_yCL4_WT-4R_Supps_FOA_50mMHU_30C_d3_016.tif` | Plate image: same strains on FOA + 50 mM HU at 30°C, day 3 |
| `sfig1c_2022_12_17_yCL4_WT-4R_Supps_YPD_02MMS_30C_d3_017.tif` | Plate image: same strains on YPD + 0.2 mM MMS at 30°C, day 3 |
| `sfig1d_2022_12_17_yCL4_WT-4R_Supps_YPD_50mMHU_30C_d3_015.tif` | Plate image: same strains on YPD + 50 mM HU at 30°C, day 3 |

## Supplementary Figure 2 - Multiple Sequence Alignment Pipeline

See `sfig2_orc4r-screen_00_README.md` for detailed pipeline documentation.

| File | Description |
|---|---|
| `sfig2_orc4r-screen_00_README.md` | Pipeline documentation and execution instructions |
| `sfig2_orc4r-screen_01_protein-accessions.tsv` | Input accession table: 94 accessions across 14 organisms, 7 genes |
| `sfig2_orc4r-screen_02_download-blosum62.R` | R script: download and validate BLOSUM62 substitution matrix |
| `sfig2_orc4r-screen_03_fetch-and-align.R` | R script: fetch sequences from UniProt/NCBI, align with DECIPHER, compute conservation |
| `sfig2_orc4r-screen_04_msa-visualization.R` | R script: generate zoomed MSA plots at positions of interest |
| `sfig2_orc4r-screen_05_renv.lock` | R environment lockfile for reproducibility |

## Supplementary Figure 3 - EMSA and MCM Loading Supplementary

| File | Description |
|---|---|
| `sfig3a_20210927-Seq Spec and ATP titr repeat-[Phosphor].tif` | Phosphor image: EMSA, sequence specificity and ATP titration |
| `sfig3b_20211208-ATPgS and Cdc6 repeat-[Phosphor].tif` | Phosphor image: EMSA, ATPgammaS and Cdc6 effects |
| `sfig3c_2021.02.23 orc titration of 4r,kt .2mM.tif` | Fluorescence gel image: MCM loading assay (Krypton stain), ORC titration of ORC4R and K485T mutants |

## Supplementary Figure 4 - MCM Loading Supplementary

| File | Description |
|---|---|
| `sfig4a_lemr load,input 350KGlut wt-4r supps.jpg` | Gel image: MCM loading, WT/ORC4R/suppressors at 350 mM potassium glutamate |
| `sfig4b_260331_loading-summary_wt-4r-supps_350mM-KGlut.csv` | Quantification: MCM loading summary at 350 mM potassium glutamate |
| `sfig4c_lemr load wt,4r,4rps arsseq kgluttit.jpg` | Gel image: MCM loading with ARS sequence variants across potassium glutamate titration |
| `sfig4d_loading-wt,4r,4ps-kglut-titr_summary-statistics.csv` | Quantification: MCM loading summary statistics across potassium glutamate titration |

## Supplementary Figure 5 - General Transcription Factor Alleles (Additional)

| File | Description |
|---|---|
| `sfig5a_20260310_exp-006_plates-03,04_day-05_time-211pm_rawscan-black.tiff` | Plate scan: general transcription factor allele strains, plates 3-4, day 5 |
| `sfig5b_20260310_exp-006_plates-05,06_day-05_time-211pm_rawscan-black.tiff` | Plate scan: general transcription factor allele strains, plates 5-6, day 5 |

## Software Requirements

### Python scripts (Figures 1, 2)

Managed with [uv](https://docs.astral.sh/uv/). Dependencies are declared
via inline script metadata (`# /// script` blocks) in each file and are
installed automatically by `uv run`.

**`fig1_orc4r-screen_01_prepare-reference.py`:**
- Python >= 3.10
- pandas

**`fig1_orc4r-screen_02_analyze-bam.py`:**
- Python >= 3.10
- numpy, pandas, pysam

**`fig2_*.py`:** Requires [PyMOL](https://pymol.org/) (>= 3.1.3.1). See script headers for additional dependencies.

### R scripts

**Figure 3a** (`fig3a_plot-sofa-atpase.R`):
- R >= 4.4.3
- tidyverse, ggplot2

**Supplementary Figure 2** (MSA pipeline):
- R >= 4.4.3
- Bioconductor: Biostrings (>= 2.74.1), DECIPHER (>= 3.2.0)
- GitHub: YuLab-SMU/ggtree (>= 4.1.2), YuLab-SMU/ggmsa (>= 1.17.0)
- CRAN: ggplot2 (>= 4.0.2)
- vctrs >= 0.7.1 (required by ggmsa)
- See `sfig2_orc4r-screen_00_README.md` for full instructions
- A pinned renv.lock is provided as `sfig2_orc4r-screen_05_renv.lock`

## Image File Formats

| Extension | Description |
|---|---|
| `.tif` / `.tiff` | Uncompressed plate photographs, gel scans, and phosphor images |
| `.jpg` | Gel images |
| `.csv` | Comma-separated processed data and quantifications |
| `.tsv` | Tab-separated data tables |
| `.bam` | Binary alignment/map files (coordinate-sorted, BWA to R64) |
| `.bam.bai` | BAM index files |

## Reproducing Analyses

**Pooled segregant analysis (Figure 1):**
1. Run `fig1_orc4r-screen_01_prepare-reference.py` to download and process SGD gene coordinates
2. Run `fig1_orc4r-screen_02_analyze-bam.py` on the BAM files in `bam_files/`

```bash
uv run fig1_orc4r-screen_01_prepare-reference.py
uv run fig1_orc4r-screen_02_analyze-bam.py
```

**Protein structure visualization (Figure 2):**
1. Download PDB structures 7MCA and 1NH2
2. Run the corresponding Python visualization scripts

**ATPase timecourse plot (Figure 3a):**
1. Run `fig3a_plot-sofa-atpase.R` with `fig3a_atpase-timecourse_processed-data.csv` as input

**Multiple sequence alignments (Supplementary Figure 2):**
1. See `sfig2_orc4r-screen_00_README.md` for step-by-step instructions

## Contact

Luis Martinez-Rodriguez - luised94@mit.edu
Stephen P. Bell (corresponding author) - spbell@mit.edu
