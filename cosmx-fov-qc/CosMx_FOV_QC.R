## ============================================================================
## CosMx FOV QC Script
## Outputs: PDF report with flagged FOV map, signal loss plot, bias heatmap,
##          and a text summary of results.
##
## USAGE: Edit the 4 variables in Section 1 below, then source the entire script.
## ============================================================================

## -- 0. Packages -------------------------------------------------------------
required_pkgs <- c("data.table", "Matrix", "FNN", "pheatmap", "scales")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

## -- 1. User settings (EDIT THESE 4 LINES) -----------------------------------

# Full path to the exprMat flat file (.csv or .csv.gz):
exprmat_path  <- "C:/path/to/ExperimentName_exprMat_file.csv.gz"

# Full path to the metadata flat file (.csv or .csv.gz):
metadata_path <- "C:/path/to/ExperimentName_metadata_file.csv.gz"

# Which panel barcode map to use:
#   Hs_IO     = Human Immuno-Oncology (100-plex)
#   Hs_UCC    = Human Universal Cell Characterization (1K)
#   Hs_6k     = Human 6K Discovery
#   Mm_Neuro  = Mouse Neuroscience (1K)
#   Mm_UCC    = Mouse Universal Cell Characterization (1K)
panel <- "Hs_6k"

# Output directory for the PDF report and any CSVs:
output_dir <- "C:/Users/YourName/Desktop"

## -- 1b. Derived settings (no edits needed) ----------------------------------

# Experiment name from filename:
expt_name <- sub("_exprMat_file\\.csv(\\.gz)?$", "", basename(exprmat_path))

# Output PDF path:
output_pdf <- file.path(output_dir, paste0(expt_name, "_FOV_QC_Report.pdf"))

# QC thresholds (defaults recommended; lower = stricter):
max_prop_loss        <- 0.6
max_totalcounts_loss <- 0.6

cat("Experiment:  ", expt_name, "\n")
cat("exprMat:     ", exprmat_path, "\n")
cat("metadata:    ", metadata_path, "\n")
cat("Panel:       ", panel, "\n")
cat("Output PDF:  ", output_pdf, "\n\n")

## -- 2. Source FOV QC functions and load barcode map --------------------------

source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/Main/_code/FOV%20QC/FOV%20QC%20utils.R")

allbarcodes <- readRDS(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))
barcodemap  <- allbarcodes[[panel]]

if (is.null(barcodemap)) stop("Panel '", panel, "' not found. Available: ", paste(names(allbarcodes), collapse = ", "))
cat("Functions sourced. Barcode map loaded for panel:", panel, "\n")

## -- 3. Load and prepare counts matrix ---------------------------------------

cat("Loading expression matrix...\n")
counts       <- as.data.frame(fread(exprmat_path))
counts$cell  <- paste0("c_1_", counts$fov, "_", counts$cell_ID)
rownames(counts) <- counts$cell

# Remove negative probes and system controls:
counts <- counts[, !grepl("Negative", colnames(counts)) &
                   !grepl("SystemControl", colnames(counts))]

# Remove non-gene columns and convert to sparse matrix:
counts2 <- counts[, -which(names(counts) %in% c("fov", "cell_ID", "cell"))]
counts2 <- as.matrix(counts2[order(row.names(counts2)), ])
counts.final <- Matrix::Matrix(counts2, sparse = TRUE)

cat("Counts matrix ready:", nrow(counts.final), "cells x", ncol(counts.final), "genes\n")
rm(counts, counts2)

## -- 4. Load and prepare metadata --------------------------------------------

cat("Loading metadata...\n")
obs <- as.data.frame(fread(metadata_path))

# Use existing cell column if present, otherwise construct it:
if (!"cell" %in% colnames(obs)) {
  obs$cell <- paste0("c_1_", obs$fov, "_", obs$cell_ID)
}
rownames(obs) <- obs$cell
obs <- obs[order(row.names(obs)), ]

# Determine XY column names (slide_mm preferred, global_px as fallback):
if (all(c("x_slide_mm", "y_slide_mm") %in% colnames(obs))) {
  xy <- obs[c("x_slide_mm", "y_slide_mm")]
} else if (all(c("CenterX_global_px", "CenterY_global_px") %in% colnames(obs))) {
  xy <- obs[c("CenterX_global_px", "CenterY_global_px")]
  colnames(xy) <- c("x_slide_mm", "y_slide_mm")
  cat("Note: Using global pixel coordinates (no slide_mm columns found).\n")
} else {
  stop("Could not find suitable XY coordinate columns in metadata.")
}

fov <- obs$fov

cat("Metadata ready:", nrow(obs), "cells,", length(unique(fov)), "FOVs\n")

## -- 5. Run FOV QC -----------------------------------------------------------

cat("Running FOV QC (this may take a few minutes)...\n")
res <- runFOVQC(counts     = counts.final,
                xy          = xy,
                fov         = fov,
                barcodemap  = barcodemap,
                max_prop_loss        = max_prop_loss,
                max_totalcounts_loss = max_totalcounts_loss)

## -- 6. Summary --------------------------------------------------------------

total_fovs   <- length(unique(fov))
n_flagged    <- length(res$flaggedfovs)
pct_flagged  <- round(100 * n_flagged / total_fovs, 1)

cat("\n")
cat("====== FOV QC SUMMARY ======\n")
cat("Experiment:              ", expt_name, "\n")
cat("Panel:                   ", panel, "\n")
cat("Total FOVs:              ", total_fovs, "\n")
cat("Flagged FOVs:            ", n_flagged, " (", pct_flagged, "%)\n")
cat("Flagged for signal loss: ", ifelse(length(res$flaggedfovs_fortotalcounts) == 0, "None",
                                         paste(res$flaggedfovs_fortotalcounts, collapse = ", ")), "\n")
cat("Flagged for bias:        ", ifelse(length(res$flaggedfovs_forbias) == 0, "None",
                                         paste(res$flaggedfovs_forbias, collapse = ", ")), "\n")
cat("All flagged FOVs:        ", ifelse(length(res$flaggedfovs) == 0, "None",
                                         paste(res$flaggedfovs, collapse = ", ")), "\n")
cat("============================\n\n")

## -- 7. Generate PDF report --------------------------------------------------

cat("Generating PDF report:", output_pdf, "\n")
pdf(output_pdf, width = 14, height = 10)

# Page 1: Summary text
plot.new()
text(0.5, 0.92, "CosMx FOV QC Report", cex = 2, font = 2)
text(0.5, 0.82, paste("Experiment:", expt_name), cex = 1.4, font = 2)
text(0.5, 0.73, paste("Panel:", panel), cex = 1.2)
text(0.5, 0.65, paste("Total FOVs:", total_fovs), cex = 1.2)
text(0.5, 0.57, paste("Flagged FOVs:", n_flagged, "(", pct_flagged, "%)"), cex = 1.2)
text(0.5, 0.49, paste("Thresholds: max_prop_loss =", max_prop_loss,
                       " max_totalcounts_loss =", max_totalcounts_loss), cex = 1)
text(0.5, 0.37, "FOVs flagged for signal loss:", cex = 1, font = 2)
text(0.5, 0.30, ifelse(length(res$flaggedfovs_fortotalcounts) == 0, "None",
                        paste(res$flaggedfovs_fortotalcounts, collapse = ", ")), cex = 1)
text(0.5, 0.22, "FOVs flagged for reporter bias:", cex = 1, font = 2)
text(0.5, 0.15, ifelse(length(res$flaggedfovs_forbias) == 0, "None",
                        paste(res$flaggedfovs_forbias, collapse = ", ")), cex = 1)
text(0.5, 0.05, "FOV QC functions: Danaher et al., NanoString Biostats", cex = 0.8, col = "grey50")

# Page 2: Flagged FOV map
mapFlaggedFOVs(res)

# Page 3: Signal loss spatial plot
FOVSignalLossSpatialPlot(res)

# Page 4: Heatmap (flagged FOVs only)
if (n_flagged > 0) {
  mat <- res$fovstats$bias[res$flaggedfovs, , drop = FALSE] *
         res$fovstats$flag[res$flaggedfovs, , drop = FALSE]
  pheatmap::pheatmap(mat,
                     col = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(100),
                     breaks = seq(-2, 2, length.out = 101),
                     cluster_rows = (nrow(mat) > 1),
                     main = "Flagged FOVs: bias by reporter cycle")
} else {
  plot.new()
  text(0.5, 0.5, "No FOVs flagged - heatmap not applicable", cex = 1.5)
}

# Page 5: Full heatmap (all FOVs)
FOVEffectsHeatmap(res)

dev.off()
cat("PDF report saved to:", output_pdf, "\n")

## -- 8. FOVs to exclude (for downstream use) ---------------------------------

cat("\nFOVs to exclude from downstream analysis:\n")
if (length(res$flaggedfovs) == 0) {
  cat("None - all FOVs passed QC\n")
} else {
  cat(paste(res$flaggedfovs, collapse = ", "), "\n")
}

# Save flagged FOV/gene pairs if any exist:
if (!is.null(res$flagged_fov_x_gene) && nrow(res$flagged_fov_x_gene) > 0) {
  flagged_csv <- file.path(output_dir, paste0(expt_name, "_flagged_fov_gene_pairs.csv"))
  write.csv(res$flagged_fov_x_gene, flagged_csv, row.names = FALSE)
  cat("Flagged FOV x gene pairs saved to:", flagged_csv, "\n")
}

cat("\nDone.\n")
