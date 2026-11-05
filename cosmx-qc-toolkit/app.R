## ============================================================================
## CosMx QC Toolkit - Shiny App
## Phase 1: FOV QC (signal loss + reporter bias)
## Phase 2: Cell-level QC (count threshold + area threshold)
##
## Launch with: shiny::runApp("path/to/this/folder")
## Or just open this file in RStudio and click "Run App"
## ============================================================================

# -- Max file upload size (8GB to handle WTX-sized exprMat files) ------------
options(shiny.maxRequestSize = 8 * 1024^3)

# -- Install/load packages ----------------------------------------------------
required_pkgs <- c("shiny", "data.table", "Matrix", "FNN", "pheatmap", 
                   "scales", "grDevices")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(shiny)
library(data.table)
library(Matrix)
library(pheatmap)
library(scales)

# -- Source FOV QC functions from GitHub --------------------------------------
source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/Main/_code/FOV%20QC/FOV%20QC%20utils.R")

# -- Load barcode maps from GitHub --------------------------------------------
allbarcodes <- readRDS(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))

# -- Chunked exprMat loader ---------------------------------------------------
# Reads large exprMat files in row chunks to avoid RAM spikes on FAS laptops.
# Always used (single code path); small files load as a single chunk.
# Peak memory = chunk_size * n_genes dense, instead of n_cells * n_genes dense.
load_exprmat_chunked <- function(path, chunk_size = 5000, progress_fn = NULL) {
  
  # Peek at header to discover columns without loading data
  header <- data.table::fread(path, nrows = 0)
  col_names <- names(header)
  
  # Validate required columns
  required_cols <- c("fov", "cell_ID")
  missing_cols <- setdiff(required_cols, col_names)
  if (length(missing_cols) > 0) {
    stop("exprMat file is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Identify column categories:
  #   - gene_cols: real biological probes (used for FOV/Cell QC)
  #   - neg_cols:  negative probes (used for Run Summary metrics)
  #   - fc_cols:   SystemControl / false-code columns (used for Run Summary metrics)
  neg_mask  <- grepl("Negative",      col_names)
  fc_mask   <- grepl("SystemControl", col_names)
  meta_or_special <- c("fov", "cell_ID", "cell")
  gene_cols <- col_names[!neg_mask & !fc_mask & !(col_names %in% meta_or_special)]
  neg_cols  <- col_names[neg_mask]
  fc_cols   <- col_names[fc_mask]
  
  if (length(gene_cols) == 0) {
    stop("No gene columns detected in exprMat file.")
  }
  
  # Read and sparsify in chunks
  sparse_chunks <- list()
  neg_sums_list <- list()   # per-cell sum across all negative probes
  fc_sums_list  <- list()   # per-cell sum across all SystemControl probes
  chunk_idx   <- 1
  skip_rows   <- 1  # account for header row
  total_cells <- 0
  
  repeat {
    if (!is.null(progress_fn)) {
      progress_fn(chunk_idx, total_cells)
    }
    
    chunk <- tryCatch(
      data.table::fread(path,
                        skip      = skip_rows,
                        nrows     = chunk_size,
                        header    = FALSE,
                        col.names = col_names),
      error = function(e) NULL
    )
    
    # Stop when the file is exhausted
    if (is.null(chunk) || nrow(chunk) == 0) break
    
    n_read <- nrow(chunk)
    
    # Build cell IDs and extract gene matrix
    cell_ids <- paste0("c_1_", chunk$fov, "_", chunk$cell_ID)
    gene_matrix <- as.matrix(chunk[, gene_cols, with = FALSE])
    rownames(gene_matrix) <- cell_ids
    
    # Convert to sparse immediately, then release the dense chunk
    sparse_chunks[[chunk_idx]] <- Matrix::Matrix(gene_matrix, sparse = TRUE)
    
    # Extract per-cell control sums (kept as named numeric vectors, tiny memory)
    if (length(neg_cols) > 0) {
      neg_mat <- as.matrix(chunk[, neg_cols, with = FALSE])
      neg_sums <- rowSums(neg_mat)
      names(neg_sums) <- cell_ids
      neg_sums_list[[chunk_idx]] <- neg_sums
      rm(neg_mat, neg_sums)
    }
    if (length(fc_cols) > 0) {
      fc_mat <- as.matrix(chunk[, fc_cols, with = FALSE])
      fc_sums <- rowSums(fc_mat)
      names(fc_sums) <- cell_ids
      fc_sums_list[[chunk_idx]] <- fc_sums
      rm(fc_mat, fc_sums)
    }
    
    rm(chunk, gene_matrix, cell_ids)
    gc(verbose = FALSE)
    
    total_cells <- total_cells + n_read
    
    # If we got fewer rows than requested, we've hit the end of the file
    if (n_read < chunk_size) break
    
    skip_rows <- skip_rows + chunk_size
    chunk_idx <- chunk_idx + 1
  }
  
  if (length(sparse_chunks) == 0) {
    stop("No cells were read from the exprMat file.")
  }
  
  # Combine sparse chunks and sort rows to match metadata alignment convention
  counts_final <- do.call(rbind, sparse_chunks)
  rm(sparse_chunks)
  gc(verbose = FALSE)
  
  ord <- order(rownames(counts_final))
  counts_final <- counts_final[ord, , drop = FALSE]
  
  # Combine control sums, aligned to the same cell order
  neg_per_cell <- if (length(neg_sums_list) > 0) {
    v <- unlist(neg_sums_list); v[rownames(counts_final)]
  } else NULL
  fc_per_cell  <- if (length(fc_sums_list) > 0) {
    v <- unlist(fc_sums_list); v[rownames(counts_final)]
  } else NULL
  
  list(
    counts       = counts_final,
    neg_per_cell = neg_per_cell,
    fc_per_cell  = fc_per_cell,
    n_neg_probes = length(neg_cols),
    n_fc_probes  = length(fc_cols),
    n_gene_cols  = length(gene_cols)
  )
}

# -- Panel-specific default count thresholds ----------------------------------
# From updated vignette (Danaher, Dec 2025):
#   1000plex: 20 (permissive) - 50 (conservative)
#   6000plex: 50 (permissive) - 100 (conservative)
#   WTX:     100 (permissive) - 200 (conservative)
panel_count_defaults <- list(
  Hs_IO    = 20,   # 100-plex
  Hs_UCC   = 20,   # 1K
  Hs_6k    = 50,   # 6K
  Hs_WTX   = 100,  # WTX (Human Whole Transcriptome, 18K-plex)
  Mm_Neuro = 20,   # 1K
  Mm_UCC   = 20,   # 1K
  Mm_WTX   = 100   # WTX (Mouse Whole Transcriptome, ~21K-plex)
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
                                "Human Whole Transcriptome (18K/WTX)" = "Hs_WTX",
                                "Human Universal Cell Char (1K)" = "Hs_UCC",
                                "Human Immuno-Oncology (100-plex)" = "Hs_IO",
                                "Mouse Neuroscience (1K)" = "Mm_Neuro",
                                "Mouse Universal Cell Char (1K)" = "Mm_UCC",
                                "Mouse Whole Transcriptome (WTX)" = "Mm_WTX"),
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
        uiOutput("pdf_status")
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
        # TAB: Run Summary (Phase 3 - Bruker-style descriptive metrics)
        # ==============================================================
        tabPanel("Run Summary",
          div(class = "card", style = "margin-top: 12px;",
            h4("Run Summary"),
            fluidRow(
              column(6,
                radioButtons("run_summary_scope", "Scope:",
                             choices  = c("Post-QC (final analysis)" = "post",
                                          "Pre-QC (raw run)"          = "pre"),
                             selected = "post", inline = TRUE)
              ),
              column(6, style = "text-align: right; padding-top: 8px;",
                downloadButton("download_summary_csv", "Download CSV",
                               style = "background:#0f3460; color:white; border:none;")
              )
            ),
            br(),
            tableOutput("run_summary_table"),
            helpText("Metrics computed from the loaded expression matrix, metadata, ",
                     "and (if present) negative probe / SystemControl columns. Tissue area ",
                     "is estimated as the convex hull of cell centroids.")
          )
        ),
        
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
              tabPanel("Counts Over Space",   plotOutput("plot_counts_spatial",  height = "500px"))
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
    
    # Control-column counts for Run Summary (Phase 3)
    neg_per_cell = NULL,    # per-cell sum across all Negative probes
    fc_per_cell  = NULL,    # per-cell sum across all SystemControl probes
    n_neg_probes = NULL,    # number of negative probes in exprMat
    n_fc_probes  = NULL,    # number of SystemControl probes in exprMat
    n_gene_cols  = NULL,    # number of gene columns loaded
    
    # Gene-nonzero counts per cell (for "unique genes per cell" metric)
    n_features   = NULL
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
        
        # --- Load expression matrix (chunked) ---
        incProgress(0.1, detail = "Loading expression matrix in chunks...")
        
        # Extract experiment name from uploaded filename
        rv$expt_name <- sub("_exprMat_file\\.csv(\\.gz)?$", "", input$exprmat_file$name)
        
        # Helper: format elapsed seconds as "Xm Ys" (used by long-running WTX runs
        # so users can distinguish "still working" from "hung")
        fmt_elapsed <- function(start_time) {
          secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
          mins <- floor(secs / 60)
          rem  <- round(secs - mins * 60)
          if (mins > 0) paste0(mins, "m ", rem, "s") else paste0(rem, "s")
        }
        
        load_start <- Sys.time()
        
        # Chunked load: reads in 5000-cell chunks, converts each to sparse
        # immediately so peak memory is chunk_size * n_genes instead of
        # n_cells * n_genes dense. Same result as a single fread + sparsify.
        loader_out <- load_exprmat_chunked(
          path       = input$exprmat_file$datapath,
          chunk_size = 5000,
          progress_fn = function(chunk_idx, cells_so_far) {
            setProgress(value  = NULL,
                        detail = paste0("Loading expression matrix: chunk ",
                                        chunk_idx, " (",
                                        format(cells_so_far, big.mark = ","),
                                        " cells, ", fmt_elapsed(load_start),
                                        " elapsed)..."))
          }
        )
        counts.final     <- loader_out$counts
        rv$neg_per_cell  <- loader_out$neg_per_cell
        rv$fc_per_cell   <- loader_out$fc_per_cell
        rv$n_neg_probes  <- loader_out$n_neg_probes
        rv$n_fc_probes   <- loader_out$n_fc_probes
        rv$n_gene_cols   <- loader_out$n_gene_cols
        rm(loader_out); gc(verbose = FALSE)
        
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
        
        # Unique genes per cell: nonzero column count per row of the sparse matrix.
        # Uses nFeature_RNA from metadata if present, else computes from counts.
        if ("nFeature_RNA" %in% colnames(obs)) {
          rv$n_features <- obs$nFeature_RNA
        } else {
          rv$n_features <- Matrix::rowSums(counts.final > 0)
        }
        
        # Extract cell area if present
        if ("Area" %in% colnames(obs)) {
          rv$cell_area <- obs$Area
        } else {
          rv$cell_area <- NULL
        }
        
        # --- Run FOV QC ---
        fovqc_start <- Sys.time()
        incProgress(0.2, detail = paste0(
          "Running FOV QC algorithm (WTX runs can take 30+ min; ",
          "load phase took ", fmt_elapsed(load_start), ")..."
        ))
        res <- runFOVQC(counts     = counts.final,
                        xy          = xy,
                        fov         = fov,
                        barcodemap  = barcodemap,
                        max_prop_loss        = input$max_prop_loss,
                        max_totalcounts_loss = input$max_totalcounts_loss)
        
        incProgress(0.45, detail = paste0("Generating results (FOV QC took ",
                                           fmt_elapsed(fovqc_start), ")..."))
        
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
  # REACTIVE: Run Summary metrics (Bruker-style descriptive table)
  # Recomputes when cell QC thresholds or the pre/post-QC toggle change.
  # ============================================================================
  run_summary <- reactive({
    req(rv$ready, rv$nCount_RNA, rv$xy)
    
    qc <- cell_qc()
    
    # Cell inclusion mask depends on toggle
    use_post_qc <- (isTRUE(input$run_summary_scope == "post"))
    keep_mask <- if (use_post_qc) !qc$flag_any else rep(TRUE, qc$n_cells)
    
    counts_vec   <- rv$nCount_RNA[keep_mask]
    features_vec <- if (!is.null(rv$n_features)) rv$n_features[keep_mask] else NA_real_
    area_vec     <- if (!is.null(rv$cell_area))  rv$cell_area[keep_mask]  else NA_real_
    neg_vec      <- if (!is.null(rv$neg_per_cell)) rv$neg_per_cell[keep_mask] else NA_real_
    fc_vec       <- if (!is.null(rv$fc_per_cell))  rv$fc_per_cell[keep_mask]  else NA_real_
    xy_kept      <- rv$xy[keep_mask, , drop = FALSE]
    fov_kept     <- rv$fov[keep_mask]
    
    n_fovs <- length(unique(fov_kept))
    
    # Tissue area via convex hull of cell centroids.
    # Coordinates are stored under x_slide_mm / y_slide_mm even when the source
    # was pixel coordinates (aliased at load time). We report the raw hull area
    # and label the units accordingly.
    if (all(c("x_slide_mm", "y_slide_mm") %in% colnames(rv$obs))) {
      area_unit <- "mm^2"
      unit_label <- "mm2"
    } else {
      # Coordinates were global pixel; hull area will be in px^2, not mm^2
      area_unit <- "px^2"
      unit_label <- "px2"
    }
    
    hull_idx <- tryCatch(
      grDevices::chull(xy_kept[, 1], xy_kept[, 2]),
      error = function(e) integer(0)
    )
    
    if (length(hull_idx) >= 3) {
      hx <- xy_kept[hull_idx, 1]
      hy <- xy_kept[hull_idx, 2]
      # Shoelace formula for polygon area
      tissue_area <- 0.5 * abs(sum(hx * c(hy[-1], hy[1]) - c(hx[-1], hx[1]) * hy))
    } else {
      tissue_area <- NA_real_
    }
    
    # Aggregate metrics (with defensive NA handling)
    n_cells_final       <- length(counts_vec)
    total_transcripts   <- sum(counts_vec, na.rm = TRUE)
    mean_tx_per_cell    <- if (n_cells_final > 0) mean(counts_vec, na.rm = TRUE) else NA_real_
    max_tx_per_cell     <- if (n_cells_final > 0) max(counts_vec, na.rm = TRUE)  else NA_real_
    mean_genes_per_cell <- if (n_cells_final > 0 && !all(is.na(features_vec))) {
                             mean(features_vec, na.rm = TRUE)
                           } else NA_real_
    mean_cell_size      <- if (n_cells_final > 0 && !all(is.na(area_vec))) {
                             mean(area_vec, na.rm = TRUE)
                           } else NA_real_
    mean_tx_per_area    <- if (!is.na(tissue_area) && tissue_area > 0) {
                             total_transcripts / tissue_area
                           } else NA_real_
    
    # Mean Negative per Plex per Cell:
    #   average negative-probe count per cell, divided by number of negative probes
    mean_neg_per_plex <- if (!is.null(rv$neg_per_cell) && rv$n_neg_probes > 0 && n_cells_final > 0) {
      mean(neg_vec, na.rm = TRUE) / rv$n_neg_probes
    } else NA_real_
    
    # Mean False Codes per Plex per Cell:
    #   same idea with SystemControl probes
    mean_fc_per_plex <- if (!is.null(rv$fc_per_cell) && rv$n_fc_probes > 0 && n_cells_final > 0) {
      mean(fc_vec, na.rm = TRUE) / rv$n_fc_probes
    } else NA_real_
    
    # SNR: mean transcripts per cell over mean negative-probe-per-plex per cell.
    # Common Bruker-style SNR definition.
    snr <- if (!is.na(mean_neg_per_plex) && mean_neg_per_plex > 0) {
      mean_tx_per_cell / (mean_neg_per_plex * rv$n_gene_cols)
    } else if (!is.null(rv$neg_per_cell) && n_cells_final > 0) {
      # Fallback: mean gene signal per probe over mean negative signal per probe
      mean_gene_per_plex <- mean_tx_per_cell / max(rv$n_gene_cols, 1)
      if (!is.na(mean_gene_per_plex) && mean_neg_per_plex > 0) mean_gene_per_plex / mean_neg_per_plex else NA_real_
    } else NA_real_
    
    # QC pass/fail bonus row
    pct_pass <- if (qc$n_cells > 0) 100 * sum(!qc$flag_any) / qc$n_cells else NA_real_
    
    list(
      scope             = if (use_post_qc) "Post-QC" else "Pre-QC",
      expt_name         = rv$expt_name,
      panel             = input$panel,
      n_fovs            = n_fovs,
      tissue_area       = tissue_area,
      tissue_area_unit  = area_unit,
      tissue_area_label = unit_label,
      mean_cell_size    = mean_cell_size,
      n_cells_final     = n_cells_final,
      total_transcripts = total_transcripts,
      mean_tx_per_cell  = mean_tx_per_cell,
      max_tx_per_cell   = max_tx_per_cell,
      mean_genes_per_cell = mean_genes_per_cell,
      mean_tx_per_area  = mean_tx_per_area,
      mean_neg_per_plex = mean_neg_per_plex,
      mean_fc_per_plex  = mean_fc_per_plex,
      snr               = snr,
      pct_pass          = pct_pass
    )
  })
  
  # Helper to build a data.frame view of the summary for rendering and export
  run_summary_df <- reactive({
    rs <- run_summary()
    fmt_int <- function(x) if (is.na(x)) "NA" else formatC(round(x), format = "d", big.mark = ",")
    fmt_num <- function(x, digits = 2) if (is.na(x)) "NA" else formatC(x, format = "f", digits = digits, big.mark = ",")
    fmt_pct <- function(x) if (is.na(x)) "NA" else paste0(formatC(x, format = "f", digits = 1), "%")
    
    data.frame(
      Metric = c(
        paste0("Experiment (", rs$scope, ")"),
        "Panel",
        "Number of FOVs",
        paste0("Total Tissue Area (", rs$tissue_area_unit, ")"),
        "Mean Single Cell Size (px^2)",
        "Number of Cells in Final Analysis",
        "Total Transcripts Assigned to Cells",
        "Mean Transcripts per Cell",
        "Maximum Transcripts per Cell",
        "Mean Unique Genes per Cell",
        paste0("Mean Transcripts per ", rs$tissue_area_unit),
        "Mean Negative per Plex per Cell",
        "Mean False Codes per Plex per Cell",
        "SNR",
        "QC Pass Rate"
      ),
      Value = c(
        as.character(rs$expt_name),
        as.character(rs$panel),
        fmt_int(rs$n_fovs),
        fmt_num(rs$tissue_area, digits = 2),
        fmt_num(rs$mean_cell_size, digits = 2),
        fmt_int(rs$n_cells_final),
        fmt_int(rs$total_transcripts),
        fmt_int(rs$mean_tx_per_cell),
        fmt_int(rs$max_tx_per_cell),
        fmt_int(rs$mean_genes_per_cell),
        fmt_num(rs$mean_tx_per_area, digits = 2),
        fmt_num(rs$mean_neg_per_plex, digits = 4),
        fmt_num(rs$mean_fc_per_plex,  digits = 4),
        fmt_num(rs$snr, digits = 2),
        fmt_pct(rs$pct_pass)
      ),
      stringsAsFactors = FALSE
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
  # RUN SUMMARY TAB OUTPUTS (Phase 3 - Bruker-style descriptive metrics)
  # ============================================================================
  
  output$run_summary_table <- renderTable({
    req(rv$ready)
    df <- run_summary_df()
    df
  }, striped = TRUE, hover = TRUE, width = "100%", align = "lr")
  
  output$download_summary_csv <- downloadHandler(
    filename = function() {
      base <- if (!is.null(rv$expt_name)) rv$expt_name else "run"
      paste0(base, "_RunSummary_",
             if (isTRUE(input$run_summary_scope == "post")) "postQC" else "preQC",
             ".csv")
    },
    content = function(file) {
      write.csv(run_summary_df(), file, row.names = FALSE)
    }
  )
  
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
    
    pal <- scales::viridis_pal(option = "B")(101)
    
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
    output_pdf <- file.path(desktop_path, paste0(expt_name, "_QC_Report.pdf"))
    
    tryCatch({
      pdf(output_pdf, width = 14, height = 10)
      
      # ---- Page 1: Run Summary (Bruker-style descriptive table) ----
      rs_df <- run_summary_df()
      rs    <- run_summary()
      
      plot.new()
      text(0.5, 0.95, "Run Summary", cex = 2.2, font = 2)
      text(0.5, 0.89, paste("Experiment:", expt_name), cex = 1.3, font = 2)
      text(0.5, 0.85, paste("Panel:", panel, "   |   Scope:", rs$scope), cex = 1.0)
      
      # Render the table as text lines (avoids gridExtra dependency)
      n_rows <- nrow(rs_df)
      y_top  <- 0.78
      y_bot  <- 0.08
      y_step <- (y_top - y_bot) / max(n_rows - 1, 1)
      for (i in seq_len(n_rows)) {
        y_i <- y_top - (i - 1) * y_step
        text(0.28, y_i, rs_df$Metric[i], cex = 0.95, adj = c(1, 0.5))
        text(0.32, y_i, rs_df$Value[i],  cex = 0.95, adj = c(0, 0.5), font = 2)
      }
      text(0.5, 0.03, "Run Summary metrics: derived from expression matrix + metadata + control probes",
           cex = 0.75, col = "grey50")
      
      # ---- Page 2: Combined Summary ----
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
      pal <- scales::viridis_pal(option = "B")(101)
      par(mar = c(1, 1, 3, 1))
      plot(xy[, 1], xy[, 2], pch = 16, cex = 0.1, asp = 1,
           col = pal[sc], xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           main = "Total Counts per Cell (Spatial)")
      lv <- c(10, 100, 500, 1000, 5000)
      ls <- round(1 + 100 * (log2(lv) - log_min) / (log_max - log_min))
      legend("bottomright", pch = 16, col = c("white", pal[ls]),
             legend = c("Counts:", lv), bg = "white")
      
      dev.off()
      
      output$pdf_status <- renderUI({
        div(style = "color:#2d6a4f; font-weight:bold; padding:8px 0; font-size:13px;",
            paste("PDF saved to:", output_pdf))
      })
      
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      output$pdf_status <- renderUI({
        div(style = "color:red; padding:8px 0; font-size:13px;",
            paste("Error saving PDF:", e$message))
      })
    })
  })
}

# ============================================================================
# LAUNCH
# ============================================================================
shinyApp(ui = ui, server = server)
