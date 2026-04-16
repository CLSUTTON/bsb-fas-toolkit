## ============================================================================
## CosMx QC Toolkit - Shiny App
## Phase 1: FOV QC (signal loss + reporter bias)
## Phase 2: Cell-level QC (count threshold + area threshold)
## Phase 3: Smoothed background visualization + filtered data export
##
## Launch with: shiny::runApp("path/to/this/folder")
## Or just open this file in RStudio and click "Run App"
## ============================================================================

# -- Max file upload size (500MB to handle large exprMat files) ---------------
options(shiny.maxRequestSize = 500 * 1024^2)

# -- Install/load packages ----------------------------------------------------
required_pkgs <- c("shiny", "data.table", "Matrix", "FNN", "pheatmap", 
                   "scales", "grDevices", "viridis")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(shiny)
library(data.table)
library(Matrix)
library(pheatmap)
library(scales)
library(viridis)
library(FNN)

# -- Helper: get a writable PDF path, falling back to timestamp if locked -----
# Solves the Windows file-lock problem when the destination PDF is already
# open in RStudio's viewer or another PDF reader. Trying to open a file for
# writing will fail fast if the lock is held, so we detect that BEFORE
# calling pdf() and fall back to a timestamped filename instead.
get_writable_pdf_path <- function(preferred_path) {
  con <- try(file(preferred_path, open = "wb"), silent = TRUE)
  if (!inherits(con, "try-error")) {
    close(con)
    return(list(path = preferred_path, fallback_used = FALSE))
  }
  dir_part  <- dirname(preferred_path)
  base_part <- tools::file_path_sans_ext(basename(preferred_path))
  ext_part  <- tools::file_ext(preferred_path)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  fallback  <- file.path(dir_part, paste0(base_part, "_", timestamp, ".", ext_part))
  return(list(path = fallback, fallback_used = TRUE))
}

# -- Source FOV QC functions from GitHub --------------------------------------
source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/Main/_code/FOV%20QC/FOV%20QC%20utils.R")

# -- Load barcode maps from GitHub --------------------------------------------
allbarcodes <- readRDS(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))

# -- Panel-specific default count thresholds ----------------------------------
# From updated vignette (Danaher, Dec 2025):
#   1000plex: 20 (permissive) - 50 (conservative)
#   6000plex: 50 (permissive) - 100 (conservative)
#   WTX:     100 (permissive) - 200 (conservative)
panel_count_defaults <- list(
  Hs_IO    = 20,   # 100-plex
  Hs_UCC   = 20,   # 1K
  Hs_6k    = 50,   # 6K
  Mm_Neuro = 20,   # 1K
  Mm_UCC   = 20    # 1K
)

# ============================================================================
# UI
# ============================================================================
ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      body {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        background-color: #f7f8fa;
      }
      .title-bar {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        color: white;
        padding: 20px 30px;
        margin: -15px -15px 20px -15px;
        border-radius: 0 0 8px 8px;
      }
      .title-bar h2 {
        margin: 0 0 4px 0;
        font-weight: 600;
        font-size: 24px;
      }
      .title-bar p {
        margin: 0;
        opacity: 0.7;
        font-size: 13px;
      }
      .card {
        background: white;
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 16px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.08);
        border: 1px solid #e8eaed;
      }
      .card h4 {
        margin-top: 0;
        color: #1a1a2e;
        font-weight: 600;
        border-bottom: 2px solid #0f3460;
        padding-bottom: 8px;
      }
      .btn-run {
        background: #0f3460;
        color: white;
        border: none;
        font-size: 16px;
        padding: 12px 40px;
        border-radius: 6px;
        width: 100%;
        margin-top: 10px;
      }
      .btn-run:hover {
        background: #1a1a2e;
        color: white;
      }
      .summary-box {
        background: #eef2f7;
        border-left: 4px solid #0f3460;
        padding: 15px 20px;
        border-radius: 0 6px 6px 0;
        margin-bottom: 12px;
        font-family: 'Consolas', 'Courier New', monospace;
        font-size: 13px;
        white-space: pre-wrap;
      }
      .plot-tab .nav-tabs > li > a {
        color: #1a1a2e;
        font-weight: 500;
      }
      .progress .progress-bar {
        background-color: #0f3460;
      }
      .status-msg {
        color: #555;
        font-style: italic;
        padding: 10px 0;
      }
      .flag-stat {
        display: inline-block;
        padding: 8px 16px;
        margin: 4px;
        border-radius: 6px;
        font-weight: 600;
        font-size: 14px;
      }
      .flag-pass { background: #d4edda; color: #155724; }
      .flag-warn { background: #fff3cd; color: #856404; }
      .flag-fail { background: #f8d7da; color: #721c24; }
    "))
  ),
  
  div(class = "title-bar",
      h2("CosMx QC Toolkit"),
      p("Automated quality control for CosMx SMI flat file exports \u2014 FOV-level and cell-level QC")
  ),
  
  fluidRow(
    
    # ================================================================
    # LEFT PANEL: inputs
    # ================================================================
    column(4,
      
      # -- File inputs --
      div(class = "card",
        h4("1. Select Input Files"),
        fileInput("exprmat_file", "Expression Matrix (.csv or .csv.gz)",
                  accept = c(".csv", ".csv.gz", ".gz")),
        fileInput("metadata_file", "Cell Metadata (.csv or .csv.gz)",
                  accept = c(".csv", ".csv.gz", ".gz")),
        
        h4("2. Select Panel"),
        selectInput("panel", "Barcode Panel:",
                    choices = c("Human 6K Discovery" = "Hs_6k",
                                "Human Universal Cell Char (1K)" = "Hs_UCC",
                                "Human Immuno-Oncology (100-plex)" = "Hs_IO",
                                "Mouse Neuroscience (1K)" = "Mm_Neuro",
                                "Mouse Universal Cell Char (1K)" = "Mm_UCC"),
                    selected = "Hs_6k")
      ),
      
      # -- FOV QC thresholds --
      div(class = "card",
        h4("3. FOV QC Thresholds"),
        sliderInput("max_totalcounts_loss", "Max total counts loss:",
                    min = 0.1, max = 0.9, value = 0.7, step = 0.1),
        sliderInput("max_prop_loss", "Max reporter bias loss:",
                    min = 0.1, max = 0.9, value = 0.7, step = 0.1)
      ),
      
      # -- Cell QC thresholds (Phase 2) --
      div(class = "card",
        h4("4. Cell QC Thresholds"),
        numericInput("count_threshold", 
                     "Min counts per cell:",
                     value = 50, min = 1, max = 1000, step = 10),
        helpText("Recommended: 20 for 1K panels, 50 for 6K, 100-200 for WTX. ",
                 "Auto-capped at 10th percentile to prevent over-filtering."),
        numericInput("area_threshold",
                     "Max cell area (pixels):",
                     value = 30000, min = 1000, max = 100000, step = 1000),
        helpText("Flags oversized cells likely from segmentation errors. ",
                 "Skip for tissues with large cells (cardiomyocytes, Purkinje).")
      ),
      
      # -- Run button --
      div(class = "card",
        actionButton("run_qc", "Run QC", class = "btn-run"),
        br(), br(),
        actionButton("save_pdf", "Save PDF to Desktop", 
                     style = "width:100%; margin-top:5px; background:#2d6a4f; color:white; border:none; padding:12px; border-radius:6px; font-size:14px;"),
        uiOutput("pdf_status"),
        br(),
        actionButton("export_filtered", "Export Filtered Data",
                     style = "width:100%; margin-top:5px; background:#52796f; color:white; border:none; padding:12px; border-radius:6px; font-size:14px;"),
        helpText("Writes QC-passing cells to Desktop as counts_filtered.csv.gz and metadata_filtered.csv.gz.",
                 style = "font-size:11px; margin-top:4px;"),
        uiOutput("export_status")
      )
    ),
    
    # ================================================================
    # RIGHT PANEL: results
    # ================================================================
    column(8,
      
      uiOutput("status"),
      
      # -- Top-level result tabs --
      tabsetPanel(id = "main_tabs", type = "pills",
        
        # ==============================================================
        # TAB: FOV QC
        # ==============================================================
        tabPanel("FOV QC",
          div(class = "card", style = "margin-top: 12px;",
            h4("FOV QC Results"),
            verbatimTextOutput("fov_summary_text"),
            tabsetPanel(id = "fov_plot_tabs", type = "tabs",
              tabPanel("Flagged FOV Map",   plotOutput("plot_flagged",          height = "500px")),
              tabPanel("Signal Loss",       plotOutput("plot_signal",           height = "500px")),
              tabPanel("Bias Heatmap (Flagged)", plotOutput("plot_heatmap_flagged", height = "500px")),
              tabPanel("Bias Heatmap (All)",     plotOutput("plot_heatmap_all",     height = "500px"))
            )
          )
        ),
        
        # ==============================================================
        # TAB: Cell QC (Phase 2)
        # ==============================================================
        tabPanel("Cell QC",
          div(class = "card", style = "margin-top: 12px;",
            h4("Cell-Level QC Results"),
            verbatimTextOutput("cell_summary_text"),
            tabsetPanel(id = "cell_plot_tabs", type = "tabs",
              tabPanel("Count Distribution",  plotOutput("plot_count_hist",  height = "500px")),
              tabPanel("Area Distribution",   plotOutput("plot_area_hist",   height = "500px")),
              tabPanel("Flagged Cells (Spatial)", plotOutput("plot_cell_flags_spatial", height = "500px"))
            )
          )
        ),
        
        # ==============================================================
        # TAB: QC Summary (combined)
        # ==============================================================
        tabPanel("QC Summary",
          div(class = "card", style = "margin-top: 12px;",
            h4("Combined QC Summary"),
            uiOutput("combined_flag_badges"),
            verbatimTextOutput("combined_summary_text"),
            tabsetPanel(id = "summary_plot_tabs", type = "tabs",
              tabPanel("All Flags (Spatial)", plotOutput("plot_combined_spatial", height = "500px")),
              tabPanel("Counts Over Space",   plotOutput("plot_counts_spatial",  height = "500px")),
              tabPanel("Smoothed Background", plotOutput("plot_smoothed_bg",     height = "500px"))
            )
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================
server <- function(input, output, session) {
  
  # -- Reactive values ---------------------------------------------------------
  rv <- reactiveValues(
    # FOV QC results (Phase 1)
    res         = NULL,
    expt_name   = NULL,
    ready       = FALSE,
    
    # Loaded data retained for Cell QC (Phase 2)
    obs         = NULL,     # full metadata data.frame
    xy          = NULL,     # xy coordinates
    fov         = NULL,     # fov vector
    nCount_RNA  = NULL,     # total counts per cell
    cell_area   = NULL,     # cell areas
    total_cells = NULL,     # total cell count
    
    # Phase 3: Smoothed background + filtered export
    nCount_negprobes = NULL,  # total negprobe counts per cell (for smoothed background)
    smoothed_bg      = NULL,  # cached smoothed background vector
    counts_final     = NULL   # retained sparse counts matrix for filtered export
  )
  
  # -- Update default count threshold when panel changes -----------------------
  observeEvent(input$panel, {
    default_val <- panel_count_defaults[[input$panel]]
    if (!is.null(default_val)) {
      updateNumericInput(session, "count_threshold", value = default_val)
    }
  })
  
  # -- Run QC (both FOV and cell data loading) ---------------------------------
  observeEvent(input$run_qc, {
    
    req(input$exprmat_file, input$metadata_file)
    
    rv$ready <- FALSE
    rv$res   <- NULL
    rv$obs   <- NULL
    
    withProgress(message = "Running QC Toolkit...", value = 0, {
      
      tryCatch({
        
        # --- Load barcode map ---
        incProgress(0.05, detail = "Loading barcode map...")
        barcodemap <- allbarcodes[[input$panel]]
        if (is.null(barcodemap)) {
          stop(paste0("Panel '", input$panel, "' not found."))
        }
        
        # --- Load expression matrix ---
        incProgress(0.1, detail = "Loading expression matrix (this may take a minute)...")
        counts <- as.data.frame(fread(input$exprmat_file$datapath))
        
        # Extract experiment name
        rv$expt_name <- sub("_exprMat_file\\.csv(\\.gz)?$", "", input$exprmat_file$name)
        
        counts$cell <- paste0("c_1_", counts$fov, "_", counts$cell_ID)
        rownames(counts) <- counts$cell
        counts <- counts[, !grepl("Negative", colnames(counts)) &
                           !grepl("SystemControl", colnames(counts))]
        counts2 <- counts[, -which(names(counts) %in% c("fov", "cell_ID", "cell"))]
        counts2 <- as.matrix(counts2[order(row.names(counts2)), ])
        counts.final <- Matrix::Matrix(counts2, sparse = TRUE)
        rm(counts, counts2)
        
        # --- Load metadata ---
        incProgress(0.15, detail = "Loading metadata...")
        obs <- as.data.frame(fread(input$metadata_file$datapath))
        
        if (!"cell" %in% colnames(obs)) {
          obs$cell <- paste0("c_1_", obs$fov, "_", obs$cell_ID)
        }
        rownames(obs) <- obs$cell
        obs <- obs[order(row.names(obs)), ]
        
        # Determine XY columns
        if (all(c("x_slide_mm", "y_slide_mm") %in% colnames(obs))) {
          xy <- obs[c("x_slide_mm", "y_slide_mm")]
        } else if (all(c("CenterX_global_px", "CenterY_global_px") %in% colnames(obs))) {
          xy <- obs[c("CenterX_global_px", "CenterY_global_px")]
          colnames(xy) <- c("x_slide_mm", "y_slide_mm")
        } else {
          stop("No suitable XY coordinate columns found in metadata.")
        }
        
        fov <- obs$fov
        
        # ============================================================
        # PHASE 2: Store metadata for cell-level QC
        # ============================================================
        rv$obs         <- obs
        rv$xy          <- xy
        rv$fov         <- fov
        rv$total_cells <- nrow(obs)
        
        # Extract nCount_RNA: use metadata column if present, else compute from counts
        if ("nCount_RNA" %in% colnames(obs)) {
          rv$nCount_RNA <- obs$nCount_RNA
        } else {
          rv$nCount_RNA <- Matrix::rowSums(counts.final)
        }
        
        # Extract cell area if present
        if ("Area" %in% colnames(obs)) {
          rv$cell_area <- obs$Area
        } else {
          rv$cell_area <- NULL
        }
        
        # Phase 3: Extract negprobe counts if present (for smoothed background plot)
        if ("nCount_negprobes" %in% colnames(obs)) {
          rv$nCount_negprobes <- obs$nCount_negprobes
        } else {
          rv$nCount_negprobes <- NULL
        }
        
        # Phase 3: Retain counts matrix for filtered export
        rv$counts_final <- counts.final
        
        # Reset smoothed background cache (recomputed on demand)
        rv$smoothed_bg <- NULL
        
        # --- Run FOV QC ---
        incProgress(0.2, detail = "Running FOV QC algorithm (this takes a few minutes)...")
        res <- runFOVQC(counts     = counts.final,
                        xy          = xy,
                        fov         = fov,
                        barcodemap  = barcodemap,
                        max_prop_loss        = input$max_prop_loss,
                        max_totalcounts_loss = input$max_totalcounts_loss)
        
        incProgress(0.45, detail = "Generating results...")
        
        rv$res   <- res
        rv$ready <- TRUE
        
        incProgress(0.05, detail = "Done!")
        
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
      })
    })
  })
  
  # ============================================================================
  # REACTIVE: Cell QC flags (recompute when thresholds change - instant)
  # ============================================================================
  cell_qc <- reactive({
    req(rv$ready, rv$nCount_RNA)
    
    n_cells <- rv$total_cells
    
    # -- Count threshold with 10th percentile cap (per updated vignette) --
    count_thresh <- min(input$count_threshold, 
                        quantile(rv$nCount_RNA, 0.10, na.rm = TRUE))
    flag_counts <- rv$nCount_RNA < count_thresh
    
    # -- Area threshold --
    if (!is.null(rv$cell_area)) {
      flag_area <- rv$cell_area > input$area_threshold
    } else {
      flag_area <- rep(FALSE, n_cells)
    }
    
    # -- FOV flags (from Phase 1) --
    flag_fov <- rv$fov %in% rv$res$flaggedfovs
    
    # -- Combined --
    flag_any <- flag_counts | flag_area | flag_fov
    
    list(
      count_thresh_used = count_thresh,
      flag_counts       = flag_counts,
      flag_area         = flag_area,
      flag_fov          = flag_fov,
      flag_any          = flag_any,
      n_flag_counts     = sum(flag_counts),
      n_flag_area       = sum(flag_area),
      n_flag_fov        = sum(flag_fov),
      n_flag_any        = sum(flag_any),
      n_cells           = n_cells
    )
  })
  
  # ============================================================================
  # STATUS MESSAGE
  # ============================================================================
  output$status <- renderUI({
    if (is.null(input$exprmat_file) || is.null(input$metadata_file)) {
      div(class = "status-msg", "Upload both files and click 'Run QC' to begin.")
    } else if (!rv$ready) {
      div(class = "status-msg", "Files selected. Click 'Run QC' to begin.")
    } else {
      NULL
    }
  })
  
  # ============================================================================
  # FOV QC TAB OUTPUTS (Phase 1 - unchanged logic)
  # ============================================================================
  
  output$fov_summary_text <- renderPrint({
    req(rv$ready, rv$res)
    res <- rv$res
    total_fovs  <- length(unique(rv$fov))
    n_flagged   <- length(res$flaggedfovs)
    pct_flagged <- round(100 * n_flagged / total_fovs, 1)
    
    cat(paste0(
      "====== FOV QC SUMMARY ======\n",
      "Experiment:              ", rv$expt_name, "\n",
      "Panel:                   ", input$panel, "\n",
      "Total FOVs:              ", total_fovs, "\n",
      "Flagged FOVs:            ", n_flagged, " (", pct_flagged, "%)\n",
      "Flagged for signal loss: ", ifelse(length(res$flaggedfovs_fortotalcounts) == 0, "None",
                                          paste(res$flaggedfovs_fortotalcounts, collapse = ", ")), "\n",
      "Flagged for bias:        ", ifelse(length(res$flaggedfovs_forbias) == 0, "None",
                                          paste(res$flaggedfovs_forbias, collapse = ", ")), "\n",
      "All flagged FOVs:        ", ifelse(length(res$flaggedfovs) == 0, "None",
                                          paste(res$flaggedfovs, collapse = ", ")), "\n",
      "============================\n"
    ))
  })
  
  output$plot_flagged <- renderPlot({
    req(rv$ready, rv$res)
    mapFlaggedFOVs(rv$res)
  })
  
  output$plot_signal <- renderPlot({
    req(rv$ready, rv$res)
    FOVSignalLossSpatialPlot(rv$res)
  })
  
  output$plot_heatmap_flagged <- renderPlot({
    req(rv$ready, rv$res)
    n_flagged <- length(rv$res$flaggedfovs)
    if (n_flagged > 0) {
      mat <- rv$res$fovstats$bias[rv$res$flaggedfovs, , drop = FALSE] *
             rv$res$fovstats$flag[rv$res$flaggedfovs, , drop = FALSE]
      pheatmap::pheatmap(mat,
                         col = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(100),
                         breaks = seq(-2, 2, length.out = 101),
                         cluster_rows = (nrow(mat) > 1),
                         main = "Flagged FOVs: bias by reporter cycle")
    } else {
      plot.new()
      text(0.5, 0.5, "No FOVs flagged - heatmap not applicable", cex = 1.5)
    }
  })
  
  output$plot_heatmap_all <- renderPlot({
    req(rv$ready, rv$res)
    FOVEffectsHeatmap(rv$res)
  })
  
  # ============================================================================
  # CELL QC TAB OUTPUTS (Phase 2)
  # ============================================================================
  
  # -- Cell QC summary text --
  output$cell_summary_text <- renderPrint({
    req(rv$ready)
    qc <- cell_qc()
    
    cat(paste0(
      "====== CELL QC SUMMARY ======\n",
      "Total cells:               ", qc$n_cells, "\n",
      "Count threshold (used):    ", qc$count_thresh_used,
      ifelse(qc$count_thresh_used < input$count_threshold,
             paste0("  (capped from ", input$count_threshold, " at 10th pctile)"), ""), "\n",
      "Flagged (low counts):      ", qc$n_flag_counts, 
      " (", round(100 * qc$n_flag_counts / qc$n_cells, 1), "%)\n",
      "Area threshold:            ", input$area_threshold, " pixels\n",
      "Flagged (oversized):       ", qc$n_flag_area, 
      " (", round(100 * qc$n_flag_area / qc$n_cells, 1), "%)\n",
      "=============================\n"
    ))
  })
  
  # -- Plot: Count distribution histogram --
  output$plot_count_hist <- renderPlot({
    req(rv$ready, rv$nCount_RNA)
    qc <- cell_qc()
    
    log_counts <- log2(pmax(rv$nCount_RNA, 1))
    thresh_log <- log2(max(qc$count_thresh_used, 1))
    
    par(mar = c(5, 4, 4, 2))
    hist(log_counts, breaks = 100, 
         main = "Total Counts per Cell",
         xlab = "Log2(counts per cell)", 
         ylab = "Number of cells",
         col = "grey80", border = "grey60")
    abline(v = thresh_log, col = "red", lwd = 2, lty = 2)
    legend("topright", 
           legend = c(
             paste0("Threshold: ", qc$count_thresh_used, " counts"),
             paste0("Flagged: ", qc$n_flag_counts, " cells (",
                    round(100 * qc$n_flag_counts / qc$n_cells, 1), "%)")
           ),
           col = c("red", NA), lwd = c(2, NA), lty = c(2, NA),
           bty = "n", cex = 1.1)
  })
  
  # -- Plot: Area distribution histogram --
  output$plot_area_hist <- renderPlot({
    req(rv$ready)
    qc <- cell_qc()
    
    if (is.null(rv$cell_area)) {
      plot.new()
      text(0.5, 0.5, "No 'Area' column found in metadata", cex = 1.5)
      return()
    }
    
    par(mar = c(5, 4, 4, 2))
    hist(rv$cell_area, breaks = 100, 
         main = "Cell Area Distribution",
         xlab = "Cell area (pixels)", 
         ylab = "Number of cells",
         col = "grey80", border = "grey60")
    abline(v = input$area_threshold, col = "red", lwd = 2, lty = 2)
    legend("topright",
           legend = c(
             paste0("Threshold: ", input$area_threshold, " px"),
             paste0("Flagged: ", qc$n_flag_area, " cells (",
                    round(100 * qc$n_flag_area / qc$n_cells, 2), "%)")
           ),
           col = c("red", NA), lwd = c(2, NA), lty = c(2, NA),
           bty = "n", cex = 1.1)
  })
  
  # -- Plot: Flagged cells spatial map --
  output$plot_cell_flags_spatial <- renderPlot({
    req(rv$ready, rv$xy)
    qc <- cell_qc()
    
    xy <- rv$xy
    
    # Color: grey = pass, orange = low counts, red = oversized, purple = both
    cell_col <- rep("grey80", qc$n_cells)
    cell_col[qc$flag_counts] <- "#e6550d"  # orange for low counts
    cell_col[qc$flag_area]   <- "#d62728"  # red for oversized
    cell_col[qc$flag_counts & qc$flag_area] <- "#7b2d8e"  # purple for both
    
    par(mar = c(1, 1, 3, 1))
    # Plot passing cells first (background), then flagged on top
    pass_idx <- !qc$flag_counts & !qc$flag_area
    plot(xy[pass_idx, 1], xy[pass_idx, 2], 
         pch = 16, cex = 0.1, asp = 1,
         col = cell_col[pass_idx],
         xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         main = "Cell-Level QC Flags (Spatial)")
    # Overlay flagged cells
    flag_idx <- qc$flag_counts | qc$flag_area
    if (any(flag_idx)) {
      points(xy[flag_idx, 1], xy[flag_idx, 2],
             pch = 16, cex = 0.3, col = cell_col[flag_idx])
    }
    legend("bottomright", pch = 16, cex = 1,
           col = c("grey80", "#e6550d", "#d62728", "#7b2d8e"),
           legend = c("Pass", 
                      paste0("Low counts (", qc$n_flag_counts, ")"),
                      paste0("Oversized (", qc$n_flag_area, ")"),
                      "Both"),
           bg = "white")
  })
  
  # ============================================================================
  # QC SUMMARY TAB OUTPUTS (Combined)
  # ============================================================================
  
  # -- Combined flag badges --
  output$combined_flag_badges <- renderUI({
    req(rv$ready)
    qc <- cell_qc()
    res <- rv$res
    
    total_fovs  <- length(unique(rv$fov))
    n_fov_flag  <- length(res$flaggedfovs)
    pct_fov     <- round(100 * n_fov_flag / total_fovs, 1)
    pct_cells   <- round(100 * qc$n_flag_any / qc$n_cells, 1)
    
    # Determine severity classes
    fov_class  <- if (pct_fov < 5) "flag-pass" else if (pct_fov < 10) "flag-warn" else "flag-fail"
    cell_class <- if (pct_cells < 5) "flag-pass" else if (pct_cells < 15) "flag-warn" else "flag-fail"
    
    div(style = "margin-bottom: 16px;",
      span(class = paste("flag-stat", fov_class),
           paste0("FOVs flagged: ", n_fov_flag, "/", total_fovs, " (", pct_fov, "%)")),
      span(class = paste("flag-stat", cell_class),
           paste0("Cells flagged: ", qc$n_flag_any, "/", qc$n_cells, " (", pct_cells, "%)"))
    )
  })
  
  # -- Combined summary text --
  output$combined_summary_text <- renderPrint({
    req(rv$ready)
    qc  <- cell_qc()
    res <- rv$res
    
    total_fovs   <- length(unique(rv$fov))
    n_fov_flag   <- length(res$flaggedfovs)
    pct_fov      <- round(100 * n_fov_flag / total_fovs, 1)
    n_pass       <- sum(!qc$flag_any)
    pct_pass     <- round(100 * n_pass / qc$n_cells, 1)
    
    cat(paste0(
      "========== COMBINED QC SUMMARY ==========\n",
      "Experiment:                ", rv$expt_name, "\n",
      "Panel:                     ", input$panel, "\n",
      "\n",
      "--- FOV-Level ---\n",
      "Total FOVs:                ", total_fovs, "\n",
      "Flagged FOVs:              ", n_fov_flag, " (", pct_fov, "%)\n",
      "  Signal loss:             ", ifelse(length(res$flaggedfovs_fortotalcounts) == 0, "None",
                                            paste(res$flaggedfovs_fortotalcounts, collapse = ", ")), "\n",
      "  Reporter bias:           ", ifelse(length(res$flaggedfovs_forbias) == 0, "None",
                                            paste(res$flaggedfovs_forbias, collapse = ", ")), "\n",
      "\n",
      "--- Cell-Level ---\n",
      "Total cells:               ", qc$n_cells, "\n",
      "Count threshold:           ", qc$count_thresh_used, "\n",
      "  Flagged (low counts):    ", qc$n_flag_counts, 
      " (", round(100 * qc$n_flag_counts / qc$n_cells, 1), "%)\n",
      "Area threshold:            ", input$area_threshold, " px\n",
      "  Flagged (oversized):     ", qc$n_flag_area, 
      " (", round(100 * qc$n_flag_area / qc$n_cells, 1), "%)\n",
      "  Flagged (in bad FOVs):   ", qc$n_flag_fov, 
      " (", round(100 * qc$n_flag_fov / qc$n_cells, 1), "%)\n",
      "\n",
      "--- Overall ---\n",
      "Total cells flagged:       ", qc$n_flag_any, 
      " (", round(100 * qc$n_flag_any / qc$n_cells, 1), "%)\n",
      "Cells passing all QC:      ", n_pass, " (", pct_pass, "%)\n",
      "=========================================\n",
      "\nFOVs to exclude:\n",
      ifelse(n_fov_flag == 0, "None - all FOVs passed QC",
             paste(res$flaggedfovs, collapse = ", ")), "\n"
    ))
  })
  
  # -- Plot: All flags spatial --
  output$plot_combined_spatial <- renderPlot({
    req(rv$ready, rv$xy)
    qc <- cell_qc()
    xy <- rv$xy
    
    # Color coding: grey=pass, orange=cell flag only, blue=FOV flag, red=both
    cell_col <- rep("grey85", qc$n_cells)
    cell_flag_only <- (qc$flag_counts | qc$flag_area) & !qc$flag_fov
    fov_flag_only  <- qc$flag_fov & !qc$flag_counts & !qc$flag_area
    both_flags     <- qc$flag_fov & (qc$flag_counts | qc$flag_area)
    
    cell_col[cell_flag_only] <- "#e6550d"  # orange
    cell_col[fov_flag_only]  <- "#3182bd"  # blue
    cell_col[both_flags]     <- "#d62728"  # red
    
    par(mar = c(1, 1, 3, 1))
    pass_idx <- !qc$flag_any
    plot(xy[pass_idx, 1], xy[pass_idx, 2], 
         pch = 16, cex = 0.1, asp = 1,
         col = cell_col[pass_idx],
         xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         main = "Combined QC Flags (Spatial)")
    flag_idx <- qc$flag_any
    if (any(flag_idx)) {
      points(xy[flag_idx, 1], xy[flag_idx, 2],
             pch = 16, cex = 0.3, col = cell_col[flag_idx])
    }
    legend("bottomright", pch = 16, cex = 1,
           col = c("grey85", "#e6550d", "#3182bd", "#d62728"),
           legend = c("Pass", "Cell QC flag", "FOV QC flag", "Both"),
           bg = "white")
  })
  
  # -- Plot: Counts over space (descriptive QC from updated vignette) --
  output$plot_counts_spatial <- renderPlot({
    req(rv$ready, rv$xy, rv$nCount_RNA)
    
    xy <- rv$xy
    counts_vec <- rv$nCount_RNA
    
    # Log2 color scale clamped to [10, 5000] per vignette approach
    log_min <- log2(10)
    log_max <- log2(5000)
    log_counts <- log2(pmax(pmin(counts_vec, 5000), 10))
    scaled <- round(1 + 100 * (log_counts - log_min) / (log_max - log_min))
    scaled <- pmax(1, pmin(101, scaled))
    
    pal <- viridis::viridis(101, option = "B")
    
    par(mar = c(1, 1, 3, 1))
    plot(xy[, 1], xy[, 2], pch = 16, cex = 0.1, asp = 1,
         col = pal[scaled],
         xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         main = "Total Counts per Cell (Spatial)")
    
    legend_vals <- c(10, 100, 500, 1000, 5000)
    legend_scaled <- round(1 + 100 * (log2(legend_vals) - log_min) / (log_max - log_min))
    legend("bottomright", pch = 16, cex = 1,
           col = c("white", pal[legend_scaled]),
           legend = c("Counts:", legend_vals),
           bg = "white")
  })
  
  # ============================================================================
  # PHASE 3: Smoothed background over space (descriptive QC from vignette §4.3)
  # Uses FNN::get.knn to compute 50 nearest neighbors per cell, then averages
  # the totalcount-normalized negprobe rate across each cell's neighborhood.
  # Results are cached in rv$smoothed_bg so the expensive KNN lookup only runs
  # once per dataset rather than on every tab click.
  # ============================================================================
  smoothed_background <- reactive({
    req(rv$ready, rv$xy, rv$nCount_RNA)
    
    # Return cached result if available
    if (!is.null(rv$smoothed_bg)) return(rv$smoothed_bg)
    
    # Can't compute without negprobe counts
    if (is.null(rv$nCount_negprobes)) return(NULL)
    
    # Compute per-cell background rate, normalized by total RNA counts
    # Guard against division by zero
    total_counts_safe <- pmax(rv$nCount_RNA, 1)
    bg_rate <- rv$nCount_negprobes / total_counts_safe
    
    # Compute 50 nearest neighbors in xy space
    xy_mat <- as.matrix(rv$xy)
    knn <- FNN::get.knn(xy_mat, k = 50)
    
    # For each cell, average bg_rate over its 50 neighbors (plus itself)
    # knn$nn.index is N x 50 matrix of neighbor indices
    neighbor_bg <- matrix(bg_rate[knn$nn.index], 
                          nrow = length(bg_rate), 
                          ncol = 50)
    smoothed <- (rowMeans(neighbor_bg) * 50 + bg_rate) / 51
    
    # Cache and return
    rv$smoothed_bg <- smoothed
    smoothed
  })
  
  output$plot_smoothed_bg <- renderPlot({
    req(rv$ready, rv$xy)
    
    smoothed <- smoothed_background()
    
    if (is.null(smoothed)) {
      plot.new()
      text(0.5, 0.55, "No 'nCount_negprobes' column found in metadata",
           cex = 1.4, font = 2)
      text(0.5, 0.45, "Smoothed background plot requires negprobe counts.",
           cex = 1.0, col = "grey40")
      text(0.5, 0.38, "This column is included in standard AtoMx exports.",
           cex = 1.0, col = "grey40")
      return()
    }
    
    xy <- rv$xy
    
    # Scale to 0.9995 quantile per vignette approach to clip outliers
    q_high <- quantile(smoothed, 0.9995, na.rm = TRUE)
    scaled <- round(1 + 100 * pmin(smoothed / q_high, 1))
    scaled <- pmax(1, pmin(101, scaled))
    
    pal <- viridis::viridis(101, option = "B")
    
    par(mar = c(1, 1, 3, 1))
    plot(xy[, 1], xy[, 2], pch = 16, cex = 0.1, asp = 1,
         col = pal[scaled],
         xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         main = "Spatially Smoothed Background (negprobe / total counts)")
    
    # Legend shows multiples of median for interpretability
    med_bg <- median(smoothed, na.rm = TRUE)
    legend_mults <- c(0.5, 1, 2, 5, 10)
    legend_vals <- signif(med_bg * legend_mults, 2)
    legend_scaled <- pmax(1, pmin(101, round(1 + 100 * pmin(legend_vals / q_high, 1))))
    legend("bottomright", pch = 16, cex = 0.9,
           col = c("white", pal[legend_scaled]),
           legend = c("x median background:", paste0(legend_mults, "x")),
           bg = "white")
  })
  
  # ============================================================================
  # SAVE PDF REPORT (updated for Phase 2)
  # ============================================================================
  observeEvent(input$save_pdf, {
    req(rv$ready, rv$res)
    
    res <- rv$res
    qc  <- cell_qc()
    total_fovs  <- length(unique(rv$fov))
    n_fov_flag  <- length(res$flaggedfovs)
    pct_fov     <- round(100 * n_fov_flag / total_fovs, 1)
    expt_name   <- rv$expt_name
    panel       <- input$panel
    
    # Write to user's Desktop
    desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
    if (!dir.exists(desktop_path)) desktop_path <- file.path(path.expand("~"), "Desktop")
    if (!dir.exists(desktop_path)) desktop_path <- path.expand("~")
    preferred_pdf <- file.path(desktop_path, paste0(expt_name, "_QC_Report.pdf"))
    
    # Phase 3: check for file lock (e.g. PDF open in RStudio viewer) BEFORE
    # calling pdf(). If locked, fall back to a timestamped filename so we
    # never fail silently or produce a corrupt PDF.
    path_info <- get_writable_pdf_path(preferred_pdf)
    output_pdf    <- path_info$path
    fallback_used <- path_info$fallback_used
    
    tryCatch({
      pdf(output_pdf, width = 14, height = 10)
      
      # ---- Page 1: Combined Summary ----
      plot.new()
      text(0.5, 0.95, "CosMx QC Toolkit Report", cex = 2.2, font = 2)
      text(0.5, 0.87, paste("Experiment:", expt_name), cex = 1.4, font = 2)
      text(0.5, 0.80, paste("Panel:", panel), cex = 1.2)
      
      # FOV summary
      text(0.5, 0.70, "--- FOV QC ---", cex = 1.2, font = 2)
      text(0.5, 0.64, paste("Total FOVs:", total_fovs, "  |  Flagged:", 
                             n_fov_flag, "(", pct_fov, "%)"), cex = 1.1)
      text(0.5, 0.58, paste("Thresholds: max_prop_loss =", input$max_prop_loss,
                             " max_totalcounts_loss =", input$max_totalcounts_loss), cex = 0.9)
      text(0.5, 0.52, paste("Signal loss:", ifelse(length(res$flaggedfovs_fortotalcounts) == 0, "None",
                              paste(res$flaggedfovs_fortotalcounts, collapse = ", "))), cex = 0.9)
      text(0.5, 0.46, paste("Reporter bias:", ifelse(length(res$flaggedfovs_forbias) == 0, "None",
                              paste(res$flaggedfovs_forbias, collapse = ", "))), cex = 0.9)
      
      # Cell summary
      text(0.5, 0.36, "--- Cell QC ---", cex = 1.2, font = 2)
      text(0.5, 0.30, paste("Total cells:", qc$n_cells, "  |  Total flagged:", 
                             qc$n_flag_any, "(", round(100 * qc$n_flag_any / qc$n_cells, 1), "%)"), cex = 1.1)
      text(0.5, 0.24, paste("Low counts (<", qc$count_thresh_used, "):", 
                             qc$n_flag_counts, "  |  Oversized (>", input$area_threshold, "px):", 
                             qc$n_flag_area), cex = 0.9)
      text(0.5, 0.18, paste("Cells passing all QC:", sum(!qc$flag_any)), cex = 1.1, font = 2)
      
      text(0.5, 0.05, "CosMx QC Toolkit | FOV QC: Danaher et al., NanoString Biostats", 
           cex = 0.8, col = "grey50")
      
      # ---- Page 2: FOV Flagged map ----
      mapFlaggedFOVs(res)
      
      # ---- Page 3: FOV Signal loss ----
      FOVSignalLossSpatialPlot(res)
      
      # ---- Page 4: FOV Flagged heatmap ----
      if (n_fov_flag > 0) {
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
      
      # ---- Page 5: FOV Full heatmap ----
      FOVEffectsHeatmap(res)
      
      # ---- Page 6: Cell QC - Count histogram ----
      log_counts <- log2(pmax(rv$nCount_RNA, 1))
      thresh_log <- log2(max(qc$count_thresh_used, 1))
      hist(log_counts, breaks = 100, 
           main = "Cell QC: Total Counts per Cell",
           xlab = "Log2(counts per cell)", ylab = "Number of cells",
           col = "grey80", border = "grey60")
      abline(v = thresh_log, col = "red", lwd = 2, lty = 2)
      legend("topright", 
             legend = c(paste0("Threshold: ", qc$count_thresh_used),
                        paste0("Flagged: ", qc$n_flag_counts, " (",
                               round(100 * qc$n_flag_counts / qc$n_cells, 1), "%)")),
             col = c("red", NA), lwd = c(2, NA), lty = c(2, NA),
             bty = "n", cex = 1.1)
      
      # ---- Page 7: Cell QC - Area histogram ----
      if (!is.null(rv$cell_area)) {
        hist(rv$cell_area, breaks = 100, 
             main = "Cell QC: Cell Area Distribution",
             xlab = "Cell area (pixels)", ylab = "Number of cells",
             col = "grey80", border = "grey60")
        abline(v = input$area_threshold, col = "red", lwd = 2, lty = 2)
        legend("topright",
               legend = c(paste0("Threshold: ", input$area_threshold, " px"),
                          paste0("Flagged: ", qc$n_flag_area, " (",
                                 round(100 * qc$n_flag_area / qc$n_cells, 2), "%)")),
               col = c("red", NA), lwd = c(2, NA), lty = c(2, NA),
               bty = "n", cex = 1.1)
      }
      
      # ---- Page 8: Counts over space ----
      xy <- rv$xy
      log_min <- log2(10); log_max <- log2(5000)
      log_c <- log2(pmax(pmin(rv$nCount_RNA, 5000), 10))
      sc <- pmax(1, pmin(101, round(1 + 100 * (log_c - log_min) / (log_max - log_min))))
      pal <- viridis::viridis(101, option = "B")
      par(mar = c(1, 1, 3, 1))
      plot(xy[, 1], xy[, 2], pch = 16, cex = 0.1, asp = 1,
           col = pal[sc], xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           main = "Total Counts per Cell (Spatial)")
      lv <- c(10, 100, 500, 1000, 5000)
      ls <- round(1 + 100 * (log2(lv) - log_min) / (log_max - log_min))
      legend("bottomright", pch = 16, col = c("white", pal[ls]),
             legend = c("Counts:", lv), bg = "white")
      
      # ---- Page 9: Smoothed background over space (Phase 3) ----
      smoothed <- smoothed_background()
      if (!is.null(smoothed)) {
        q_high <- quantile(smoothed, 0.9995, na.rm = TRUE)
        sc_bg  <- pmax(1, pmin(101, round(1 + 100 * pmin(smoothed / q_high, 1))))
        par(mar = c(1, 1, 3, 1))
        plot(xy[, 1], xy[, 2], pch = 16, cex = 0.1, asp = 1,
             col = pal[sc_bg], xlab = "", ylab = "", xaxt = "n", yaxt = "n",
             main = "Spatially Smoothed Background (negprobe / total counts)")
        med_bg <- median(smoothed, na.rm = TRUE)
        lm_mults <- c(0.5, 1, 2, 5, 10)
        lm_vals <- signif(med_bg * lm_mults, 2)
        lm_sc <- pmax(1, pmin(101, round(1 + 100 * pmin(lm_vals / q_high, 1))))
        legend("bottomright", pch = 16, cex = 0.9,
               col = c("white", pal[lm_sc]),
               legend = c("x median background:", paste0(lm_mults, "x")),
               bg = "white")
      } else {
        plot.new()
        text(0.5, 0.5, "Smoothed background unavailable (no nCount_negprobes column)",
             cex = 1.3)
      }
      
      dev.off()
      
      output$pdf_status <- renderUI({
        if (fallback_used) {
          div(style = "color:#d68910; font-weight:bold; padding:8px 0; font-size:13px;",
              HTML(paste0("Original PDF was open in another program.<br>",
                          "Saved with timestamp instead:<br>",
                          "<span style='font-weight:normal;'>", output_pdf, "</span>")))
        } else {
          div(style = "color:#2d6a4f; font-weight:bold; padding:8px 0; font-size:13px;",
              paste("PDF saved to:", output_pdf))
        }
      })
      
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      output$pdf_status <- renderUI({
        div(style = "color:red; padding:8px 0; font-size:13px;",
            paste("Error saving PDF:", e$message))
      })
    })
  })
  
  # ============================================================================
  # PHASE 3: Export filtered data (counts + metadata) as gzipped CSVs
  # Writes QC-passing cells to the user's Desktop, matching the same flat-file
  # format as AtoMx exports. This is the hand-off point from QC to downstream
  # analysis (Seurat, Giotto, AtoMx secondary analysis, etc).
  # ============================================================================
  observeEvent(input$export_filtered, {
    req(rv$ready, rv$counts_final, rv$obs)
    
    qc <- cell_qc()
    keep_idx <- !qc$flag_any
    n_keep <- sum(keep_idx)
    
    if (n_keep == 0) {
      output$export_status <- renderUI({
        div(style = "color:red; padding:8px 0; font-size:13px;",
            "Cannot export: all cells flagged. Loosen thresholds.")
      })
      return()
    }
    
    # Resolve Desktop path (same logic as PDF save)
    desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
    if (!dir.exists(desktop_path)) desktop_path <- file.path(path.expand("~"), "Desktop")
    if (!dir.exists(desktop_path)) desktop_path <- path.expand("~")
    
    expt_name <- rv$expt_name
    counts_out_preferred <- file.path(
      desktop_path,
      paste0(expt_name, "_exprMat_filtered.csv.gz"))
    meta_out_preferred <- file.path(
      desktop_path,
      paste0(expt_name, "_metadata_filtered.csv.gz"))
    
    # Apply the same file-lock safety check we use for PDFs
    counts_info <- get_writable_pdf_path(counts_out_preferred)
    meta_info   <- get_writable_pdf_path(meta_out_preferred)
    counts_out  <- counts_info$path
    meta_out    <- meta_info$path
    fallback_used <- counts_info$fallback_used || meta_info$fallback_used
    
    withProgress(message = "Exporting filtered data...", value = 0, {
      tryCatch({
        incProgress(0.2, detail = "Subsetting counts matrix...")
        
        # Subset counts - retain original fov/cell_ID columns by rebuilding
        # from metadata rather than the sparse matrix (which dropped them)
        counts_sparse_kept <- rv$counts_final[keep_idx, , drop = FALSE]
        
        incProgress(0.3, detail = "Converting to dense for export...")
        # Densify only the kept cells - memory hit is bounded by filtered size
        counts_dense <- as.matrix(counts_sparse_kept)
        
        # Prepend fov + cell_ID columns to match AtoMx flat-file format
        meta_kept <- rv$obs[keep_idx, , drop = FALSE]
        counts_df <- data.frame(
          fov     = meta_kept$fov,
          cell_ID = meta_kept$cell_ID,
          counts_dense,
          check.names = FALSE
        )
        
        incProgress(0.3, detail = "Writing gzipped CSVs to Desktop...")
        data.table::fwrite(counts_df, counts_out, compress = "gzip")
        data.table::fwrite(meta_kept, meta_out, compress = "gzip")
        
        incProgress(0.2, detail = "Done!")
        
        n_total  <- length(keep_idx)
        n_flag   <- sum(!keep_idx)
        pct_kept <- round(100 * n_keep / n_total, 1)
        
        output$export_status <- renderUI({
          prefix <- if (fallback_used) {
            "Files were open elsewhere - saved with timestamp instead.<br>"
          } else ""
          div(style = "color:#2d6a4f; font-weight:bold; padding:8px 0; font-size:13px;",
              HTML(paste0(
                prefix,
                "Exported ", format(n_keep, big.mark = ","), " of ",
                format(n_total, big.mark = ","), " cells (", pct_kept, "% kept).<br>",
                "<span style='font-weight:normal;'>",
                basename(counts_out), "<br>", basename(meta_out),
                "</span>"
              )))
        })
        
      }, error = function(e) {
        output$export_status <- renderUI({
          div(style = "color:red; padding:8px 0; font-size:13px;",
              paste("Export error:", e$message))
        })
      })
    })
  })
}

# ============================================================================
# LAUNCH
# ============================================================================
shinyApp(ui = ui, server = server)
