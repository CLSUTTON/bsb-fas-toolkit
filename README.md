# BSB FAS Toolkit

A collection of field application tools for Bruker Spatial Biology platforms. Each tool wraps published bioinformatics methods into accessible interfaces that FAS teams and customers can run without writing code.

## Tools

| Tool | Platform | Description | Status |
|------|----------|-------------|--------|
| [CosMx QC Toolkit](cosmx-qc-toolkit/) | CosMx SMI | Shiny GUI for FOV-level and cell-level QC on flat file exports. Flags bad FOVs (signal loss, reporter bias) and bad cells (low counts, segmentation errors). Generates interactive reports and downloadable PDFs. | Active |
| [CosMx FOV QC Script](cosmx-fov-qc/) | CosMx SMI | Standalone R script for FOV QC. Edit 4 lines, source, get a PDF report. For users who prefer scripts over GUIs. | Active |

## Who this is for

- **FAS teams** running QC during customer engagements
- **Customers** performing initial QC on their own data
- **RAMs/FSEs** who need a quick data quality check before escalating to bioinformatics

## Requirements

- R 4.0+ (RStudio recommended)
- Internet connection on first run (downloads functions and barcode maps from GitHub)
- See individual tool READMEs for RAM requirements

## References

- [CosMx Analysis Scratch Space](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/about.html)
- [FOV QC method — Danaher et al.](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/fov-qc/index.html)
- [Basic analysis vignette (updated Dec 2025)](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/vignette-basic-analysis-updated/01_preprocessing.html)
