# CosMx FOV QC Script

Automated FOV quality control for CosMx SMI flat file exports. Identifies FOVs with low overall signal or biased reporter cycle expression compared to spatially similar regions, and generates a PDF report.

## What it does

The script applies the FOV QC approach from [Danaher et al. (NanoString Biostats)](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/fov-qc/index.html). For each FOV, a 7x7 grid is placed across the tissue. Each grid square is compared to the 10 most similar squares in other FOVs (one neighbor per FOV). FOVs are then flagged for:

1. **Signal loss**: >75% of grid squares have total counts below the threshold compared to comparators (default threshold: 60% loss).
2. **Reporter cycle bias**: Systematic underexpression of genes sharing a reporter cycle, indicating potential imaging or chemistry issues.

## Exporting data from AtoMx

Before running this script, you need to export two flat files from AtoMx SIP:

1. Log in to [AtoMx SIP](https://app.atomx.io).
2. Navigate to your study and select the experiment/slide you want to QC.
3. Go to **Export** (or **Download Data**).
4. Select **Flat file** export format.
5. Export the following two files:
   - **Expression matrix** (`*_exprMat_file.csv` or `.csv.gz`): Cell-by-gene counts matrix containing `fov`, `cell_ID`, and gene expression columns.
   - **Cell metadata** (`*_metadata_file.csv` or `.csv.gz`): Per-cell metadata containing `fov`, `cell_ID`, spatial coordinates (`CenterX_global_px`/`CenterY_global_px` or `x_slide_mm`/`y_slide_mm`), and cell morphology metrics.
6. Save both files to a local folder. Note the full file paths for use in the script.

**Important**: Export one slide/tissue at a time. This script should be run per-slide, not across multiple slides in a single run.

## Quick start

1. Open `CosMx_FOV_QC.R` in RStudio.
2. Edit the **4 lines** in Section 1:

```r
exprmat_path  <- "C:/path/to/ExperimentName_exprMat_file.csv.gz"   # full path to exprMat
metadata_path <- "C:/path/to/ExperimentName_metadata_file.csv.gz"  # full path to metadata
panel         <- "Hs_6k"                                           # see panel options below
output_dir    <- "C:/Users/YourName/Desktop"                       # where the PDF goes
```

3. Select all (Ctrl+A) and run (Ctrl+Enter).

The script extracts the experiment name from the exprMat filename and writes a PDF report named `{ExperimentName}_FOV_QC_Report.pdf` to `output_dir`.

## Requirements

R 4.0+ with the following packages (auto-installed on first run):

- data.table
- Matrix
- FNN
- pheatmap
- scales

An internet connection is required on first run to download the FOV QC utility functions and barcode maps from GitHub.

## Panel options

| Panel ID   | Description                                  |
|------------|----------------------------------------------|
| `Hs_IO`    | Human Immuno-Oncology (100-plex)             |
| `Hs_UCC`   | Human Universal Cell Characterization (1K)   |
| `Hs_6k`    | Human 6K Discovery                           |
| `Mm_Neuro` | Mouse Neuroscience (1K)                      |
| `Mm_UCC`   | Mouse Universal Cell Characterization (1K)   |

**Note**: The Human Whole Transcriptome (18K) panel and custom/RBS panels are not currently supported by the upstream barcode map file. Contact Bioinformatics (Patrick Danaher) for the appropriate barcode map if using those panels.

## Input files

The script expects two files from the standard CosMx/AtoMx flat file export:

| File | Required columns |
|------|-----------------|
| `*_exprMat_file.csv(.gz)` | `fov`, `cell_ID`, plus one column per gene |
| `*_metadata_file.csv(.gz)` | `fov`, `cell_ID`, and either `x_slide_mm`/`y_slide_mm` or `CenterX_global_px`/`CenterY_global_px` |

Both compressed (.csv.gz) and uncompressed (.csv) files are supported. The script auto-detects which coordinate system is available in the metadata.

## Output

A multi-page PDF report:

| Page | Contents |
|------|----------|
| 1 | Summary: experiment name, panel, total FOVs, flagged count/percentage, flagged FOV IDs by failure type |
| 2 | Flagged FOV map: spatial layout with flagged FOVs in red, passing FOVs in blue |
| 3 | Signal loss plot: log2 fold-change in total counts per grid square vs. comparable regions |
| 4 | Flagged FOV heatmap: reporter cycle bias for flagged FOVs only |
| 5 | Full heatmap: reporter cycle bias across all FOVs |

If any FOVs are flagged for reporter cycle bias, a CSV of flagged FOV/gene pairs (`{ExperimentName}_flagged_fov_gene_pairs.csv`) is also saved to `output_dir`.

## QC thresholds

Default thresholds are set in Section 1b of the script and are recommended for most experiments:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `max_prop_loss` | 0.7 | Flags FOVs with reporter bias exceeding 60% signal loss |
| `max_totalcounts_loss` | 0.7 | Flags FOVs with total counts below 40% of comparable regions |

Lower values are stricter (flag more FOVs). Only adjust if you have a specific reason to do so.

## Interpreting results

- **All-white bias heatmap for flagged FOVs**: Flags are driven entirely by total signal loss, not reporter cycle issues. Likely cause is tissue or section quality (thin sections, damage, detachment, high autofluorescence). Recommendation: exclude those FOVs from downstream analysis.
- **Colored bands in the heatmap**: Specific reporter cycles are underperforming in those FOVs. May indicate imaging or chemistry issues. Check if affected cycles share a common fluorescent channel.
- **Low failure rate (<5%)**: Typical for a good experiment. Exclude flagged FOVs and proceed.
- **High failure rate (>10%)**: May indicate a systemic issue. Review slide preparation and imaging conditions before proceeding with analysis.

## Running multiple experiments

To process multiple experiments, update the 4 paths in Section 1 for each experiment and re-source the script. The output PDF is auto-named per experiment, so files will not overwrite each other.

## References

- [FOV QC blog post (concept and approach)](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/fov-qc/index.html)
- [FOV QC vignette (worked example)](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/_code/FOV%20QC/FOV-QC-vignette.html)
- [FOV QC source code (utility functions)](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/FOV%20QC)
