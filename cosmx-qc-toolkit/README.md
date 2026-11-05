# CosMx QC Toolkit

GUI-based quality control toolkit for CosMx SMI flat file exports. Wraps published bioinformatics best practices ([Danaher et al., NanoString Biostats](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/fov-qc/index.html)) into a Shiny app that any FAS or analyst can run without writing R code.

## What it does

Upload two flat files from AtoMx, select your panel, click Run. The toolkit produces an interactive report plus a downloadable PDF and CSV.

### Run Summary
A Bruker-style descriptive summary of the run at a glance: FOV count, tissue area, cell counts, transcript stats, control probe stats, and SNR. Toggle between Pre-QC (raw run) and Post-QC (final analysis) views. Exportable to CSV for downstream reporting.

### FOV-Level QC
Identifies FOVs with low overall signal or biased reporter cycle expression compared to spatially similar regions. For each FOV, a 7x7 grid is placed across the tissue. Each grid square is compared to the 10 most similar squares in other FOVs. FOVs are flagged for:

1. **Signal loss**: >75% of grid squares have total counts below the threshold compared to comparators.
2. **Reporter cycle bias**: Systematic underexpression of genes sharing a reporter cycle, indicating potential imaging or chemistry issues.

### Cell-Level QC
Flags individual cells that should be excluded from downstream analysis:

1. **Low count cells**: Cells below a minimum total count threshold (panel-specific defaults, auto-capped at the 10th percentile to prevent over-filtering).
2. **Oversized cells**: Cells exceeding a maximum area threshold, likely resulting from segmentation errors.

Cell QC thresholds update reactively. Adjust them after the initial run and see results instantly without re-running the FOV QC computation.

### Combined QC Summary
Merges FOV-level and cell-level flags into a single report with color-coded severity badges, spatial flag maps, and a total counts-over-space visualization.

## Quick start

1. Open `app.R` in RStudio.
2. Click **Run App** (or run `shiny::runApp("path/to/this/folder")`).
3. Upload your expression matrix and metadata files.
4. Select your panel and adjust thresholds if needed.
5. Click **Run QC**.
6. Review results across the four tabs: Run Summary, FOV QC, Cell QC, QC Summary.
7. Click **Save PDF to Desktop** for a shareable report, or use **Download CSV** on the Run Summary tab for the metrics table.

## Exporting data from AtoMx

Before running the toolkit, export two flat files from AtoMx SIP:

1. Log in to [AtoMx SIP](https://app.atomx.io).
2. Navigate to your study and select the experiment/slide to QC.
3. Go to **Export** and select **Flat file** format.
4. Export:
   - **Expression matrix** (`*_exprMat_file.csv` or `.csv.gz`)
   - **Cell metadata** (`*_metadata_file.csv` or `.csv.gz`)
5. Save both files locally.

**Important**: Export one slide at a time. This toolkit should be run per-slide.

## Requirements

- R 4.0+ with RStudio
- Internet connection on first run (downloads FOV QC functions and barcode maps from GitHub)
- RAM: 16GB minimum for 1K panels, 32GB recommended for 6K panels, 64GB recommended for WTX

The following packages are auto-installed on first run:

- shiny
- data.table
- Matrix
- FNN
- pheatmap
- scales
- grDevices

## Panel options

| Panel ID   | Description                                | Default count threshold |
|------------|--------------------------------------------|------------------------|
| `Hs_IO`    | Human Immuno-Oncology (100-plex)           | 20                     |
| `Hs_UCC`   | Human Universal Cell Characterization (1K) | 20                     |
| `Hs_6k`    | Human 6K Discovery                         | 50                     |
| `Hs_WTX`   | Human Whole Transcriptome (~18K)           | 100                    |
| `Mm_Neuro` | Mouse Neuroscience (1K)                    | 20                     |
| `Mm_UCC`   | Mouse Universal Cell Characterization (1K) | 20                     |
| `Mm_WTX`   | Mouse Whole Transcriptome (~21K)           | 100                    |

**Note**: Custom and RBS panels are not currently supported by the upstream barcode map file. Contact Bioinformatics for the appropriate barcode map if using those panels.

**Note on WTX runtime**: WTX panels (Hs_WTX and Mm_WTX) have ~3x the gene count of 6K panels and can take significantly longer to load and compute. Loading phase alone can be 20-40 minutes on a FAS laptop; full FOV QC can push past an hour. Progress messages now show elapsed time so you can distinguish "still working" from "hung."

## Input files

| File | Required columns |
|------|-----------------|
| `*_exprMat_file.csv(.gz)` | `fov`, `cell_ID`, plus one column per gene. Negative probes and SystemControl columns are used for Run Summary metrics. |
| `*_metadata_file.csv(.gz)` | `fov`, `cell_ID`, spatial coordinates (`x_slide_mm`/`y_slide_mm` or `CenterX_global_px`/`CenterY_global_px`), `nCount_RNA` (optional), `nFeature_RNA` (optional), `Area` (optional) |

Both compressed (.csv.gz) and uncompressed (.csv) files are supported. The toolkit auto-detects coordinate systems and handles missing `nCount_RNA` (computed from expression matrix), missing `nFeature_RNA` (computed from expression matrix), and missing `Area` (area filtering skipped) gracefully.

## Output

### Interactive (in-app)

Four tabs:

**Run Summary tab:**
- Bruker-style descriptive metrics table
- Pre-QC / Post-QC toggle
- Download CSV button

**FOV QC tab:**
- Flagged FOV spatial map
- Signal loss plot
- Bias heatmap (flagged FOVs)
- Bias heatmap (all FOVs)

**Cell QC tab:**
- Total counts histogram with threshold
- Cell area histogram with threshold
- Flagged cells spatial map (color-coded by flag reason)

**QC Summary tab:**
- Color-coded severity badges (green/yellow/red)
- Combined text summary
- All-flags spatial map
- Counts-over-space heatmap (viridis)

### PDF report (9 pages)

| Page | Contents |
|------|----------|
| 1 | Run Summary metrics table (Bruker-style) |
| 2 | Combined summary: experiment, panel, FOV flags, cell flags, overall pass rate |
| 3 | Flagged FOV map |
| 4 | Signal loss plot |
| 5 | Flagged FOV bias heatmap |
| 6 | All-FOV bias heatmap |
| 7 | Cell QC: count distribution histogram |
| 8 | Cell QC: area distribution histogram |
| 9 | Counts over space (viridis heatmap) |

## Run Summary metrics

| Metric | Definition |
|--------|-----------|
| Number of FOVs | Unique FOVs in the (scoped) cell set |
| Total Tissue Area | Convex hull area of cell centroids (mm^2 or px^2) |
| Mean Single Cell Size | Mean of the `Area` metadata column |
| Number of Cells in Final Analysis | Cell count under the current scope |
| Total Transcripts Assigned to Cells | Sum of `nCount_RNA` across cells |
| Mean / Max Transcripts per Cell | Central tendency of `nCount_RNA` |
| Mean Unique Genes per Cell | Nonzero gene columns per cell |
| Mean Transcripts per Area | Total transcripts / tissue area |
| Mean Negative per Plex per Cell | Mean negative-probe count per cell / number of negative probes |
| Mean False Codes per Plex per Cell | Same, for SystemControl probes |
| SNR | Mean gene signal per probe / mean negative signal per probe |
| QC Pass Rate | Fraction of cells passing all Cell and FOV QC flags |

**Note on tissue area**: Convex hull is a reasonable approximation but will differ from methods based on bounding boxes or FOV-area sums. Numbers may not exactly match AtoMx-reported area figures.

**Note on SNR**: Uses a standard signal-to-noise ratio. Bruker's own SNR formula is not published, so exact matches to Bruker reports are not guaranteed.

## QC thresholds

### FOV QC

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `max_prop_loss` | 0.7 | Flags FOVs with reporter bias exceeding signal loss threshold |
| `max_totalcounts_loss` | 0.7 | Flags FOVs with total counts below threshold vs comparable regions |

Lower values are stricter (flag more FOVs).

### Cell QC

| Parameter | Default (6K) | Meaning |
|-----------|-------------|---------|
| Min counts per cell | 50 | Cells below this are flagged. Auto-capped at 10th percentile. |
| Max cell area | 30,000 px | Cells above this are flagged as likely segmentation errors. |

Count threshold defaults are panel-specific and update automatically when the panel is changed. The 10th percentile cap prevents over-filtering in datasets with genuinely lower signal.

**Note on area threshold**: Skip or increase for tissues with biologically large cells (cardiomyocytes, Purkinje neurons).

## Interpreting results

- **FOV flag rate <5%**: Typical for a good experiment. Exclude flagged FOVs and proceed.
- **FOV flag rate >10%**: May indicate a systemic issue. Review slide preparation and imaging conditions.
- **Cell flag rate <10%**: Normal range. Proceed with filtered dataset.
- **Cell flag rate 10-20%**: Worth investigating. Check if flags are spatially clustered (technical) or dispersed (biological).
- **Cell flag rate >20%**: Likely a tissue quality or segmentation issue. Review in AtoMx tissue viewer.
- **Counts-over-space plot**: Look for sharp FOV-boundary artifacts (technical) vs smooth spatial gradients (biological). Smooth gradients are expected.
- **10th percentile cap triggered**: Indicates the default threshold would flag >10% of cells. The dataset may have generally lower signal. Check tissue quality.

## Changelog

### v0.3.2
- Added an explicit WTX runtime-expectations warning. Selecting Hs_WTX or Mm_WTX and clicking Run QC now opens a modal describing expected runtimes (60-120 min load, 30-60 min FOV QC) and RAM requirements (32 GB minimum, 64 GB preferred). The warning fires once per session per WTX panel selection. Non-WTX panels are unaffected.

### v0.3.1
- Added Mm_WTX (Mouse Whole Transcriptome, ~21K-plex) as a selectable panel with count threshold default of 100.
- Progress messages during loading and FOV QC now include elapsed time, so long-running WTX jobs can be distinguished from a hung session at a glance.

### v0.3.0
- Added Run Summary tab as the first tab. Provides a Bruker-style descriptive metrics table with Pre-QC / Post-QC toggle and CSV export.
- Chunked loader now preserves per-cell sums across Negative probe and SystemControl columns, enabling probe-level Run Summary metrics without noticeable memory impact.
- Added Run Summary as page 1 of the PDF export.

### v0.2.x
- Chunked expression matrix loading. Files load in 5000-cell chunks, reducing peak memory usage on FAS laptops. Progress bar reports chunk number and running cell count during load.
- Replaced viridis dependency with scales::viridis_pal() for the same colormap. Resolves S7 install failures on some machines.
- Restored Hs_WTX (Human Whole Transcriptome, ~18K) as a selectable panel with count threshold default of 100.
- Increased max upload size from 500MB to 8GB to handle WTX-sized exprMat files.
- FOV QC threshold defaults remain at 0.7 per Bioinformatics guidance.

## Roadmap

- [x] Phase 1: FOV QC (signal loss + reporter bias)
- [x] Phase 2: Cell-level QC (count + area thresholds)
- [x] Phase 3: Run Summary (Bruker-style descriptive metrics)
- [ ] Phase 4: Normalization preview

## References

- [FOV QC method (Danaher et al.)](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/fov-qc/index.html)
- [Basic analysis vignette (updated Dec 2025)](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/vignette-basic-analysis-updated/01_preprocessing.html)
- [CosMx Analysis Scratch Space](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/about.html)
- [FOV QC source code](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/FOV%20QC)
