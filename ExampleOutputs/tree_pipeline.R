#!/usr/bin/env Rscript
# =============================================================
# tree_pipeline.R
# Reads methods.txt (key=value format) and runs NJ tree analysis
# with noise-based support estimation.
#
# methods.txt fields:
#   Filename        – Input CSV (rows=samples, cols=features)
#   Input_Type      – log10 | Count | Frequency
#   Distance_Metric – Non-phylo: Bray | Euclidean
#                     Phylo:     UniFrac | UniFracW | MPD | MPDw |
#                                MNTD | MNTDw
#   PhyloTree       – path to .nwk or None (required for phylo metrics)
#   Normalization   – MinMax | Zscore | TSS | Hellinger |
#                     CLR | Rarefaction | Frequency
#   SD              – Non-phylo: 0 = Control only | >0 = also run Test
#                     Phylo: NA (not applicable -- see Support below)
#   Support         – (optional, default 100) number of replicates.
#                     Non-phylo: SD-noise replicates (unchanged).
#                     Phylo: raw-count multinomial bootstrap replicates
#                     -- the established support method for
#                     UniFrac/MPD/MNTD, since perturbing feature
#                     values with synthetic noise doesn't apply to a
#                     distance defined against a fixed reference tree.
#
# Derived automatically from Input_Type (non-phylo methods only):
#   log10     -> Gaussian noise,   limma
#   Count     -> Poisson noise,    MaAsLin2
#   Frequency -> Dirichlet noise,  MaAsLin2
#
# Output structure (folder named after Filename):
#   {FileBase}/{Normalization}/{Distance}/SD0/            (Control, or
#                                                            the single
#                                                            bootstrap
#                                                            run for
#                                                            phylo methods)
#   {FileBase}/{Normalization}/{Distance}/SD{x}/           (Test, non-phylo only)
# =============================================================

n_total <- 101   # 1 main tree + 100 noise replicates

suppressPackageStartupMessages({
  libs <- c("ape","vegan","picante","phyloseq","dplyr","stringr",
            "tools","progress","ggplot2","ggtree","ggnewscale","gtools",
            "reshape2","limma","ggrepel","uwot")
  cran_libs <- setdiff(libs, c("limma","Maaslin2"))
  miss_cran  <- cran_libs[!(cran_libs %in% installed.packages()[,"Package"])]
  if (length(miss_cran)) {
    message("Installing CRAN packages: ", paste(miss_cran, collapse=", "))
    install.packages(miss_cran, repos="https://cran.rstudio.com/")
  }
  bioc_libs <- c("limma")
  miss_bioc  <- bioc_libs[!(bioc_libs %in% installed.packages()[,"Package"])]
  if (length(miss_bioc)) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager", repos="https://cran.rstudio.com/")
    BiocManager::install(miss_bioc, ask=FALSE)
  }
  lapply(libs, require, character.only=TRUE)
})
message("Packages loaded.")

# =============================================================
# 1. Helpers
# =============================================================
get_script_dir <- function() {
  ca       <- commandArgs(trailingOnly=FALSE)
  file_arg <- grep("^--file=", ca, value=TRUE)
  if (length(file_arg))
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  if (requireNamespace("rstudioapi", quietly=TRUE) &&
      rstudioapi::isAvailable() &&
      nchar(rstudioapi::getSourceEditorContext()$path) > 0)
    return(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)))
  getwd()
}
script_dir <- get_script_dir()

ensure_dir <- function(d)
  if (!dir.exists(d)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Parse key=value methods.txt, stripping ## comments
read_methods <- function(path) {
  if (!file.exists(path)) stop("methods.txt not found: ", path)
  lines <- readLines(path, warn=FALSE, encoding="UTF-8")
  lines <- sub("\\s*##.*$", "", lines)
  lines <- gsub("\\r", "", lines)
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0]
  cfg   <- list()
  for (ln in lines) {
    kv  <- strsplit(ln, "=", fixed=TRUE)[[1]]
    if (length(kv) < 2) next
    key <- trimws(kv[1])
    val <- trimws(paste(kv[-1], collapse="="))
    cfg[[key]] <- val
  }
  cfg
}

# =============================================================
# 2. Normalization
# =============================================================
apply_normalization <- function(mat, method) {
  mat <- as.matrix(mat)
  switch(method,
    None = mat,

    MinMax = {
      t(apply(mat, 1, function(x) {
        r <- range(x, na.rm = TRUE)
        if (diff(r) == 0) return(rep(0, length(x)))
        (x - r[1]) / diff(r)
      }))
    },

    Zscore = {
      apply(mat, 2, function(x) {
        s <- sd(x, na.rm=TRUE)
        if (is.na(s) || s == 0) return(rep(0, length(x)))
        (x - mean(x, na.rm=TRUE)) / s
      })
    },

    TSS = {
      rs <- rowSums(mat, na.rm=TRUE)
      mat / ifelse(rs == 0, 1, rs)
    },

    Hellinger = {
      rs  <- rowSums(mat, na.rm=TRUE)
      tss <- mat / ifelse(rs == 0, 1, rs)
      sqrt(tss)
    },

    CLR = {
      mat[mat <= 0] <- 0.5
      t(apply(mat, 1, function(x) log(x) - mean(log(x))))
    },

    Rarefaction = {
      mat <- round(mat); mat[mat < 0] <- 0
      min_depth <- min(rowSums(mat))
      if (min_depth < 1) {
        warning("Rarefaction: some rows sum to 0; skipping rarefaction.")
        return(mat)
      }
      feat_names <- colnames(mat)
      samp_names <- rownames(mat)
      result <- t(apply(mat, 1, function(x) {
        tot <- sum(x)
        if (tot <= min_depth) return(x)
        as.numeric(rmultinom(1, min_depth, x / tot))
      }))
      colnames(result) <- feat_names
      rownames(result) <- samp_names
      result
    },

    Frequency = mat,   # already proportions

    stop("Unknown normalization: ", method)
  )
}

# =============================================================
# 3. Noise injection
# =============================================================
inject_gaussian <- function(mat, sd) {
  if (sd == 0) return(mat)
  noise <- matrix(rnorm(length(mat), mean=0, sd=sd),
                  nrow=nrow(mat), ncol=ncol(mat))
  out <- mat + noise
  out[out < 0] <- 0
  out
}

inject_dirichlet <- function(mat, sd) {
  if (sd == 0) return(mat)
  out <- mat
  for (i in seq_len(nrow(mat))) {
    row <- mat[i, ]
    tot <- sum(row, na.rm=TRUE)
    if (tot <= 0) next
    # Concentration: smaller sd = more concentrated (less noise)
    alpha        <- pmax(row / tot / (sd^2), 1e-6)
    samp         <- as.numeric(gtools::rdirichlet(1, alpha))
    out[i, ]     <- samp * tot
  }
  out
}

inject_poisson <- function(mat, sd) {
  if (sd == 0) return(mat)
  # Rescale counts by 1/sd^2 before Poisson sampling, then rescale back.
  # Smaller sd -> larger scale -> higher effective counts -> less relative
  # sampling noise; larger sd -> more noise. Keeps SD semantics consistent
  # with the Gaussian and Dirichlet noise models.
  scale   <- 1 / (sd^2)
  lambda  <- as.numeric(mat) * scale
  lambda[lambda < 0] <- 0
  samp    <- rpois(length(lambda), lambda) / scale
  out     <- matrix(samp, nrow=nrow(mat), ncol=ncol(mat))
  dimnames(out) <- dimnames(mat)
  out
}

# Dispatcher: Input_Type drives noise model
inject_noise <- function(mat, sd, input_type) {
  it <- tolower(input_type)
  if (it == "log10")      inject_gaussian(mat, sd)
  else if (it == "count") inject_poisson(mat, sd)
  else                    inject_dirichlet(mat, sd)
}

# =============================================================
# 3b. Raw-count bootstrap (phylo methods only)
# =============================================================
# Phylo distances (UniFrac/UniFracW/MPD/MPDw/MNTD/MNTDw) are computed
# against an already-established reference tree, so perturbing the
# feature VALUES with synthetic Gaussian/Poisson/Dirichlet noise (the
# SD mechanism used for Bray/Euclidean) isn't the standard way to get
# support for them. Instead: resample each sample's raw read counts
# via multinomial bootstrap (same total depth, same relative
# proportions) and rebuild the tree -- the standard microbiome
# bootstrap-support approach, matching the strategy already used
# elsewhere for these exact metrics. SD is not applicable here (it's
# reported as NA for phylo methods -- see the main loop below).
bootstrap_reads <- function(mat) {
  boot <- mat
  for (i in seq_len(nrow(mat))) {
    reads <- as.numeric(mat[i, ])
    tot   <- sum(reads, na.rm = TRUE)
    if (tot > 0) {
      prob      <- reads / tot
      boot[i, ] <- as.numeric(rmultinom(1, round(tot), prob))
    }
  }
  boot
}

# =============================================================
# 4. Distance computation
# =============================================================
phylo_metrics <- c("UniFrac","UniFracW","MPD","MPDw","MNTD","MNTDw")

compute_distance <- function(mat, dist_method, ref_tree=NULL, coph=NULL) {
  switch(dist_method,
    Bray      = as.matrix(vegan::vegdist(mat, method="bray", na.rm=TRUE)),
    Euclidean = as.matrix(dist(mat)),
    UniFrac   = {
                  otu <- phyloseq::otu_table(t(mat), taxa_are_rows=TRUE)
                  ps  <- phyloseq::phyloseq(otu, phyloseq::phy_tree(ref_tree))
                  as.matrix(phyloseq::distance(ps, method="unifrac"))
                },
    UniFracW  = {
                  otu <- phyloseq::otu_table(t(mat), taxa_are_rows=TRUE)
                  ps  <- phyloseq::phyloseq(otu, phyloseq::phy_tree(ref_tree))
                  as.matrix(phyloseq::distance(ps, method="wunifrac"))
                },
    MPD  = as.matrix(picante::comdist(mat,  coph, abundance.weighted=FALSE)),
    MPDw = as.matrix(picante::comdist(mat,  coph, abundance.weighted=TRUE)),
    MNTD = as.matrix(picante::comdistnt(mat, coph, abundance.weighted=FALSE)),
    MNTDw= as.matrix(picante::comdistnt(mat, coph, abundance.weighted=TRUE)),
    stop("Unknown distance metric: ", dist_method)
  )
}

# =============================================================
# 5. Tree builder
# =============================================================
make_tree <- function(norm_mat, dist_method, ref_tree=NULL, coph=NULL) {
  d  <- compute_distance(norm_mat, dist_method, ref_tree, coph)
  tr <- ape::nj(as.dist(d))
  tr$edge.length[tr$edge.length < 0] <- 0
  tr
}

# =============================================================
# 6. Align matrix columns to phylo tree tips
# =============================================================
align_to_tree <- function(mat, ref_tree, label="matrix") {
  common <- intersect(colnames(mat), ref_tree$tip.label)
  if (!length(common))
    stop("No overlapping taxa between ", label, " columns and tree tips.")
  dropped_mat  <- setdiff(colnames(mat),     common)
  dropped_tree <- setdiff(ref_tree$tip.label, common)
  if (length(dropped_mat))
    message("  Dropping ", length(dropped_mat), " columns not in tree.")
  if (length(dropped_tree))
    message("  Pruning ", length(dropped_tree), " tree tips not in matrix.")
  mat2  <- mat[, common, drop=FALSE]
  tree2 <- ape::keep.tip(ref_tree, common)
  mat2  <- mat2[, tree2$tip.label, drop=FALSE]
  list(mat=mat2, tree=tree2)
}

# =============================================================
# Patient color palette (IDs 1–11)
# =============================================================
# Fixed colors for tip label text (status)
STATUS_COLORS <- c(
  "progressing" = "#E91E8C",   # pink
  "stable"      = "#26C6DA"    # cyan
)

PATIENT_COLORS <- c(
  "1"  = "#76D1C3",
  "2"  = "#FFF6A0",
  "3"  = "#C7B7E8",
  "4"  = "#FF8C7A",
  "5"  = "#90CAF9",
  "6"  = "#FFB74D",
  "7"  = "#AED581",
  "8"  = "#F8BBD0",
  "9"  = "#BDBDBD",
  "10" = "#BA68C8",
  "11" = "#C5E1A5"
)

# =============================================================
# 7. ggtree visualization — tips coloured by Patient ID
# =============================================================
plot_ggtree <- function(tr, outdir, title, target_yes=character(0),
                        sample_to_patient=NULL,
                        patient_colors=PATIENT_COLORS) {
  ensure_dir(outdir)

  tips    <- tr$tip.label
  n_tips  <- length(tips)
  node_df <- data.frame(
    node      = seq(n_tips + 1, n_tips + tr$Nnode),
    bootstrap = suppressWarnings(as.numeric(tr$node.label)),
    stringsAsFactors=FALSE)

  # --- Patient colouring (primary) ---
  if (!is.null(sample_to_patient)) {
    tip_patient <- as.character(sample_to_patient[tips])
    tip_patient[is.na(tip_patient)] <- "unknown"
    tip_color   <- ifelse(tip_patient %in% names(patient_colors),
                          patient_colors[tip_patient], "#AAAAAA")

    # Extract "stable" or "progressing" from sample name (pattern: patient_status_timepoint)
    tip_status       <- sapply(strsplit(tips, "_"), function(x) if (length(x) >= 2) x[2] else "unknown")
    tip_status_color <- ifelse(tip_status %in% names(STATUS_COLORS),
                               STATUS_COLORS[tip_status], "#555555")

    tip_df <- data.frame(label=tips, status=tip_status, patient=tip_patient,
                         color=tip_color, status_color=tip_status_color,
                         stringsAsFactors=FALSE)

    present_raw <- unique(tip_patient[tip_patient != "unknown"])
    present     <- present_raw[order(suppressWarnings(as.integer(present_raw)))]
    leg_colors  <- setNames(
      ifelse(present %in% names(patient_colors),
             patient_colors[present], "#AAAAAA"),
      paste0("Patient ", present))

    p <- ggtree(tr, layout="rectangular", branch.length="branch.length",
                color="gray30", linewidth=0.7) %<+% tip_df %<+% node_df +
      geom_treescale(x=NULL, y=-2, offset=0.5, fontsize=4.5,
                     linesize=0.6, color="grey40") +
      geom_nodepoint(aes(fill=bootstrap),
                     shape=22, size=5, color="gray25",
                     stroke=0.2, na.rm=TRUE) +
      scale_fill_gradient(
        name   = "Bootstrap",
        low    = "#B3C6FF", high="#1A237E",
        limits = c(0, 100), na.value="white",
        guide  = guide_colorbar(order=2, barwidth=1, barheight=8)) +
      ggnewscale::new_scale_fill() +
      geom_tippoint(aes(fill=color), shape=21, size=4,
                    color="gray25", stroke=0.3) +
      scale_fill_identity(
        name   = "Patient",
        labels = names(leg_colors),
        breaks = leg_colors,
        guide  = guide_legend(order=1,
                              override.aes=list(size=5, shape=21,
                                                color="gray25", stroke=0.3))) +
      geom_tiplab(aes(label=status, color=status_color), size=4,
                  hjust=-0.1, align=FALSE) +
      scale_color_identity()

  } else {
    # Fallback: Y/N target colouring
    status  <- ifelse(tips %in% target_yes, "Y", "N")
    tip_df  <- data.frame(label=tips, status=status, stringsAsFactors=FALSE)
    target_cols <- c("Y"="#E41A1C", "N"="#4878CF")

    p <- ggtree(tr, layout="rectangular", branch.length="branch.length",
                color="gray30", linewidth=0.7) %<+% tip_df %<+% node_df +
      geom_treescale(x=NULL, y=-2, offset=0.5, fontsize=4.5,
                     linesize=0.6, color="grey40") +
      geom_nodepoint(aes(fill=bootstrap),
                     shape=22, size=5, color="gray25",
                     stroke=0.2, na.rm=TRUE) +
      scale_fill_gradient(
        name   = "Bootstrap",
        low    = "#B3C6FF", high="#1A237E",
        limits = c(0, 100), na.value="white",
        guide  = guide_colorbar(order=2, barwidth=1, barheight=8)) +
      ggnewscale::new_scale_fill() +
      geom_tippoint(aes(fill=status), shape=21, size=4,
                    color="gray25", stroke=0.3) +
      scale_fill_manual(
        name   = "Target",
        values = target_cols,
        labels = c("Y"="Y (target)", "N"="N (other)"),
        guide  = guide_legend(order=1, override.aes=list(size=5))) +
      geom_tiplab(aes(label=label), size=4, hjust=-0.1,
                  color="gray15", align=FALSE)
  }

  p <- p +
    hexpand(0.45, direction=1) +
    ggtitle(title) +
    theme_void() +
    theme(
      plot.background   = element_rect(fill="white", color=NA),
      panel.background  = element_rect(fill="white", color=NA),
      legend.background = element_rect(fill="white", color=NA),
      legend.key        = element_rect(fill="white", color=NA),
      legend.position   = "right",
      legend.direction  = "vertical",
      legend.box        = "vertical",
      legend.spacing.y  = unit(8, "pt"),
      legend.margin     = margin(0, 0, 0, 10),
      plot.title        = element_text(size=15, hjust=0.5, face="bold",
                                       color="gray15", margin=margin(b=12)),
      legend.text       = element_text(size=10),
      legend.title      = element_text(size=11, face="bold"),
      plot.margin       = margin(20, 30, 20, 20))

  n_tips   <- length(tr$tip.label)
  fig_h    <- max(20, n_tips * 0.20)   # ~0.20 in per tip, minimum 20 in
  fig_w    <- 14

  ggsave(file.path(outdir, "support_tree.pdf"),
         plot=p, width=fig_w, height=fig_h, units="in",
         bg="white", limitsize=FALSE)
  ggsave(file.path(outdir, "support_tree.png"),
         plot=p, width=fig_w, height=fig_h, units="in",
         dpi=150, bg="white", limitsize=FALSE)

  # tip_render_order.csv + patient_samples.pdf/png -- REMOVED. Both were
  # patient-identity-oriented (patient_samples colors by Patient; the tip
  # order CSV only existed to keep a separate per-patient renderer in
  # sync), which no longer belongs in a pipeline that treats patient
  # identity as unknown for target samples. support_tree.pdf/png above
  # (the plain tree, no patient coloring) is unaffected.

  invisible(NULL)
}

# =============================================================
# 7c. Write run log
# =============================================================
write_run_log <- function(outdir, params, res, clade_info=NULL,
                          rf_vs_ctrl=NA, coph_cor=NA) {
  ensure_dir(outdir)
  log_lines <- c(
    paste0("SampleTree Pipeline Log"),
    paste0("Generated : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0(strrep("-", 50)),
    paste0("Filename       : ", params$filename),
    paste0("Input_Type     : ", params$input_type),
    paste0("Distance       : ", params$dist_method),
    paste0("Normalization  : ", params$norm),
    paste0("SD             : ", params$sd),
    paste0("Noise_Type     : ", params$noise_type),
    paste0("Library        : ", params$library_name),
    paste0("Replicates     : ", params$n_total,
           "  (1 main + ", params$n_total - 1, " noise)"),
    paste0(strrep("-", 50)),
    paste0("Samples        : ", params$n_samples),
    paste0("Features       : ", params$n_features),
    paste0("Target (Y)     : ", params$n_target),
    paste0(strrep("-", 50)),
    paste0("Mean Bootstrap : ", round(res$mean_sup, 2), "%"),
    paste0("Reps Kept      : ", length(res$rep_trees_kept))
  )
  if (!is.na(rf_vs_ctrl))
    log_lines <- c(log_lines, paste0("Mean RF vs Ctrl: ", round(rf_vs_ctrl, 4)))
  if (!is.na(coph_cor))
    log_lines <- c(log_lines, paste0("Mean Coph Corr : ", round(coph_cor, 4)))
  if (!is.null(clade_info)) {
    n_clades <- if (is.null(clade_info$clade_df)) 0 else nrow(clade_info$clade_df)
    log_lines <- c(log_lines,
                   paste0(strrep("-", 50)),
                   paste0("Clades Found   : ", n_clades))
    if (n_clades > 0) {
      log_lines <- c(log_lines,
        paste0("  ", paste(clade_info$clade_df$clade_name, collapse=", ")))
    }
  }
  writeLines(log_lines, file.path(outdir, "log.txt"))
  invisible(NULL)
}

# =============================================================
# 8. Run condition (Control SD=0 or Test SD>0)
# =============================================================
run_condition <- function(mat, sd, outdir, norm, dist_method,
                          use_phy, ref_tree=NULL, coph=NULL,
                          input_type="log10", target_yes=character(0),
                          sample_to_patient=NULL) {
  ensure_dir(outdir)
  n_reps <- n_total - 1   # 100

  # Main tree: always built without noise
  prep_main <- apply_normalization(mat, norm)
  main_tr   <- make_tree(prep_main, dist_method,
                         if (use_phy) ref_tree else NULL,
                         if (use_phy) coph     else NULL)

  # Noise replicates
  rep_trees <- vector("list", n_reps)
  pb <- progress::progress_bar$new(
    total=n_reps,
    format="  Rep :current/:total [:bar] :percent | elapsed :elapsed | remaining :eta",
    clear=FALSE, width=80)

  for (i in seq_len(n_reps)) {
    pb$tick()
    # Phylo methods: raw-count multinomial bootstrap (established
    # method for UniFrac/MPD/MNTD support), not SD noise injection.
    # Non-phylo methods: unchanged SD noise injection.
    noisy <- if (use_phy) bootstrap_reads(mat) else inject_noise(mat, sd, input_type)
    np    <- apply_normalization(noisy, norm)
    rep_trees[[i]] <- tryCatch(
      make_tree(np, dist_method,
                if (use_phy) ref_tree else NULL,
                if (use_phy) coph     else NULL),
      error=function(e) NULL)
  }
  good <- rep_trees[!vapply(rep_trees, is.null, logical(1))]

  # Node support = % of replicates agreeing with main topology
  sup <- NULL
  if (length(good) > 0) {
    sup <- ape::prop.clades(main_tr, good) / length(good) * 100
    main_tr$node.label <- round(sup)
  } else {
    main_tr$node.label <- rep(NA_real_, main_tr$Nnode)
  }

  # RF distances (main vs each replicate)
  rf_main <- numeric(0)
  if (length(good) > 0) {
    all_trees <- structure(c(list(main_tr), good), class="multiPhylo")
    rf_main   <- suppressWarnings(as.matrix(ape::dist.topo(all_trees, method="PH85")))[1, -1]
  }

  # Phylo methods: drop Normalization/SD from the figure titles too --
  # same reasoning as the folder-nesting change (no SD sweep, single
  # fixed normalization per job, so both are just clutter on the plot).
  ggtree_title <- if (use_phy) {
    paste0("NJ | ", dist_method, "  n=", n_total, " bootstrap")
  } else {
    paste0("NJ | ", dist_method, " | ", norm, "  SD=", sd, "  n=", n_total)
  }

  ape::write.tree(main_tr, file.path(outdir, "support_tree.nwk"))
  plot_ggtree(
    main_tr, outdir,
    title=ggtree_title,
    target_yes=target_yes,
    sample_to_patient=sample_to_patient,
    patient_colors=PATIENT_COLORS)

  # Ordination plots (PCoA/NMDS/UMAP) -- REMOVED. These were generic
  # exploratory views unrelated to the Y-target clade question this
  # pipeline actually answers; only the tree (support_tree.*), the
  # Y-target highlight (support_tree_highlighted.jpeg), the detected
  # clades (clades/), and the TREND differential-analysis output
  # (trend_outputs/) matter here.

  list(
    sd=sd, norm=norm, dist_method=dist_method,
    main_tr=main_tr, sup=sup, rep_trees_kept=good,
    mean_sup   = if (!is.null(sup)) mean(sup,             na.rm=TRUE) else NA,
    median_sup = if (!is.null(sup)) median(sup,           na.rm=TRUE) else NA,
    pct_ge70   = if (!is.null(sup)) mean(sup >= 70, na.rm=TRUE) * 100 else NA,
    pct_ge90   = if (!is.null(sup)) mean(sup >= 90, na.rm=TRUE) * 100 else NA,
    pct_same   = if (length(rf_main)) mean(rf_main == 0) * 100        else NA,
    mean_rf    = if (length(rf_main)) mean(rf_main)                   else NA,
    n_valid    = length(good))
}


# =============================================================
# 11. Highlighted tree (base R, red=Y targets, blue=contamination)
# =============================================================
plot_highlighted_tree <- function(tree, targeted_set, contamination_set,
                                  outfile, tag, method,
                                  clade_member_targets=character(0)) {
  # Okabe-Ito colorblind-safe palette -- the standard categorical palette
  # for scientific figures, in place of raw base-R primaries (red/blue/
  # purple/black), which aren't colorblind-accessible and read as
  # informal in a publication figure.
  COL_TARGET <- "#0072B2"  # blue -- Target (Y), any sample with Target=Y
  COL_CONTAM <- "#D55E00"  # vermillion -- Contamination (non-target, admitted into clade)
  COL_OTHER  <- "grey30"   # neutral, matches this file's other grey tip/label text

  tips <- tree$tip.label
  cols <- setNames(rep(COL_OTHER, length(tips)), tips)
  cols[intersect(names(cols), contamination_set)] <- COL_CONTAM
  cols[intersect(names(cols), targeted_set)]      <- COL_TARGET

  n_tips  <- length(tips)

  # Bootstrap values stored in node.label (may be character or numeric)
  boot_vals <- suppressWarnings(as.numeric(tree$node.label))

  ok <- FALSE
  try({
    tip_df2  <- data.frame(label=names(cols), color=cols, stringsAsFactors=FALSE)
    node_df2 <- data.frame(
      node      = seq(n_tips + 1, n_tips + tree$Nnode),
      bootstrap = boot_vals,
      stringsAsFactors=FALSE)

    p <- ggtree(tree, branch.length="branch.length") %<+% tip_df2 %<+% node_df2 +

      # Branch-length scale bar
      geom_treescale(x=NULL, y=-2, offset=0.5, fontsize=4.5, linesize=0.6, color="grey40") +

      # Bootstrap as text on all internal nodes
      geom_nodelab(aes(label=ifelse(!is.na(bootstrap),
                                    as.character(round(bootstrap)), "")),
                   size=2.2, hjust=1.2, vjust=-0.4, color="grey30") +

      # Tip colours (Okabe-Ito colorblind-safe palette) — no legend, identity scale
      geom_tippoint(aes(color=color), size=2.5) +
      geom_tiplab(aes(color=color, label=label), size=4, offset=0.001) +
      scale_color_identity() +

      ggtitle(paste0("SampleTree  ", tag, "  (", method, ")\n",
                     "Blue = Target (Y)  |  Vermillion = Contamination  |  Grey = Other\n",
                     "Node labels = Support/bootstrap values")) +
      theme_tree2() +
      theme(
        plot.title      = element_text(hjust=0.5, size=10),
        legend.position = "none"
      )

    ggsave(outfile, p, width=14, height=max(12, n_tips * 0.12),
           units="in", dpi=300, limitsize=FALSE)
    ok <- TRUE
  }, silent=FALSE)

  if (!ok) {
    try({
      grDevices::jpeg(outfile, width=12, height=max(12, n_tips * 0.12),
                      units="in", res=200)
      par(mar=c(4,1,3,1))
      plot(tree, cex=0.5, tip.color=cols[tree$tip.label],
           show.node.label=TRUE, use.edge.length=TRUE,
           main=paste("SampleTree", tag, "(", method, ")"))
      ape::add.scale.bar(cex=0.7, col="grey40")
      legend("topleft", legend=c("Target (Y)","Contamination","Other"),
             col=c(COL_TARGET, COL_CONTAM, COL_OTHER), pch=16, bty="n", cex=0.8)
      grDevices::dev.off()
    }, silent=TRUE)
  }
  invisible(NULL)
}

# =============================================================
# 12. Clade detection helpers
# =============================================================
get_tip_descendants <- function(tree, node) {
  if (node <= length(tree$tip.label)) return(tree$tip.label[node])
  desc     <- c(node)
  to_visit <- node
  while (length(to_visit) > 0) {
    cur      <- to_visit[1]
    kids     <- tree$edge[tree$edge[,1] == cur, 2]
    desc     <- c(desc, kids)
    to_visit <- c(to_visit[-1], kids[kids > length(tree$tip.label)])
  }
  tree$tip.label[sort(desc[desc <= length(tree$tip.label)])]
}

find_clades_global <- function(tree, tag, method, outdir,
                               min_targeted, max_clade, max_others,
                               targeted_set, assign_policy="best") {
  clade_dir <- file.path(outdir, "clades")
  ensure_dir(clade_dir)
  clade_csv <- file.path(clade_dir, paste0("clades_", tag, "_", method, ".csv"))

  tips_all     <- tree$tip.label
  targeted_set <- intersect(tips_all, targeted_set)
  others_all   <- setdiff(tips_all, targeted_set)
  nTips        <- length(tips_all)
  internal     <- (nTips + 1):(nTips + tree$Nnode)

  cand <- list()
  for (nd in internal) {
    desc  <- unique(get_tip_descendants(tree, nd))
    n_tgt <- length(intersect(desc, targeted_set))
    n_oth <- length(intersect(desc, others_all))
    n_tot <- length(desc)
    if (n_tot  < min_targeted || n_tot > max_clade) next
    if (n_tgt  < min_targeted) next
    if (n_oth  > max_others)   next
    cand[[length(cand)+1]] <- list(node=nd, samples=desc,
      tgt=intersect(desc, targeted_set), oth=intersect(desc, others_all),
      n_tips=n_tot, n_tgt=n_tgt, n_oth=n_oth)
  }

  if (!length(cand)) {
    empty <- data.frame(Tag=character(0), Method=character(0),
      Clade_Node=integer(0), Clade_Name=character(0),
      N_Tips=integer(0), Targeted_Count=integer(0), Others_Count=integer(0),
      Targeted_Samples=character(0), Others=character(0),
      Contamination_Flag=character(0), stringsAsFactors=FALSE)
    utils::write.csv(empty, clade_csv, row.names=FALSE)
    message("  No clades meeting thresholds.")
    return(list(clades=empty, contam_tips=character(0), kept=list()))
  }

  ord <- switch(tolower(assign_policy),
    "first"    = seq_along(cand),
    "largest"  = order(vapply(cand, function(x) -x$n_tips, numeric(1))),
    "smallest" = order(vapply(cand, function(x)  x$n_tips, numeric(1))),
    # "best" is default: most targeted -> fewest others -> smallest total
    order(vapply(cand, function(x) -x$n_tgt,  numeric(1)),
          vapply(cand, function(x)  x$n_oth,  numeric(1)),
          vapply(cand, function(x)  x$n_tips, numeric(1)))
  )
  assigned <- character(0)
  kept     <- list()
  for (i in ord) {
    ci <- cand[[i]]
    if (!length(intersect(ci$samples, assigned))) {
      kept     <- c(kept, list(ci))
      assigned <- union(assigned, ci$samples)
    }
  }

  contam_union <- character(0)
  clade_df <- dplyr::bind_rows(lapply(kept, function(ci) {
    if (length(ci$oth)) contam_union <<- union(contam_union, ci$oth)
    data.frame(Tag=tag, Method=method, Clade_Node=ci$node,
      Clade_Name=paste0("node", ci$node),
      N_Tips=ci$n_tips, Targeted_Count=ci$n_tgt, Others_Count=ci$n_oth,
      Targeted_Samples=paste(sort(ci$tgt), collapse=","),
      Others=paste(sort(ci$oth), collapse=","),
      Contamination_Flag=ifelse(length(ci$oth)>0,"Yes","No"),
      stringsAsFactors=FALSE)
  }))

  utils::write.csv(clade_df, clade_csv, row.names=FALSE)
  message("  Clades found: ", nrow(clade_df), " -> ", clade_csv)
  list(clades=clade_df, contam_tips=contam_union, kept=kept)
}

# =============================================================
# 13. TREND helpers
# =============================================================
plot_group_box <- function(Mat, feats, in_group_vec, group_label, out_png) {
  if (!length(feats)) return(invisible(NULL))
  long_df <- reshape2::melt(as.matrix(Mat[, feats, drop=FALSE]),
                            varnames=c("SampleID","Feature"), value.name="Value")
  long_df$Group <- ifelse(long_df$SampleID %in% in_group_vec, group_label, "other")
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x=Group, y=Value, fill=Group)) +
    ggplot2::geom_boxplot(outlier.shape=NA) +
    ggplot2::geom_jitter(width=0.2, size=1, alpha=0.6) +
    ggplot2::facet_wrap(~Feature, scales="free_y") +
    ggplot2::theme_minimal(base_size=12) +
    ggplot2::labs(title=paste0("Boxplots: ", group_label, " vs other"),
                  x="Group", y="Value") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
  ggplot2::ggsave(out_png, plot=p, width=12, height=7, dpi=300)
  invisible(TRUE)
}

run_trend_for_job <- function(otu, clade_info, outdir, tag, method,
                               matrix_label, library_name, target_yes,
                               fdr_cutoff=0.05, min_in_group=2, min_in_other=2) {
  if (!length(clade_info$kept)) {
    message("  (TREND) No clades; skipping.")
    return(invisible(NULL))
  }
  Mat        <- otu
  trend_root <- file.path(outdir, "trend_outputs")
  ensure_dir(trend_root)
  all_samples <- rownames(Mat)
  lib_lower   <- tolower(library_name)

  for (ci in clade_info$kept) {
    clade_name       <- paste0("node", ci$node)
    in_clade_ids     <- intersect(all_samples, ci$samples)
    # In-group = every sample IN the highlighted clade -- targeted (Y) AND
    # the tolerated contamination (N, admitted by max_others) alike -- not
    # just the Y-targeted subset. "Other" = everyone NOT in this clade at
    # all. (Was: only in_target_ids counted as in-group, so a clade's own
    # contamination fell into "other" instead of the clade it's visually
    # highlighted as part of.)
    grp              <- ifelse(rownames(Mat) %in% in_clade_ids, clade_name, "other")
    names(grp)       <- rownames(Mat)
    n_g <- sum(grp == clade_name); n_o <- sum(grp == "other")
    if (n_g < min_in_group || n_o < min_in_other) {
      message("  (TREND) Skipping ", clade_name, ": too few samples.")
      next
    }

    # tag and matrix_label are the same value (the file_base) for every
    # job in this pipeline, so including both duplicated it in every
    # folder/file name (e.g. "Cell_Bray_node201_Cell"). Only append
    # matrix_label when it actually differs from tag.
    g_name_parts <- c(tag, method, clade_name)
    if (!identical(matrix_label, tag)) g_name_parts <- c(g_name_parts, matrix_label)
    g_safe <- make.names(paste(g_name_parts, collapse="_"))
    g_dir  <- file.path(trend_root, paste0("group_", g_safe))
    ensure_dir(g_dir)
    message("  (TREND) Clade ", clade_name, " | ", library_name)

    if (lib_lower == "limma") {
      fgrp   <- factor(grp, levels=c(clade_name,"other"))
      design <- model.matrix(~0 + fgrp)
      colnames(design) <- c(make.names(clade_name), "other")
      fit    <- limma::lmFit(t(Mat), design)
      contr  <- limma::makeContrasts(
        contrasts=paste0(make.names(clade_name), " - other"), levels=design)
      fit2   <- limma::eBayes(limma::contrasts.fit(fit, contr))
      res    <- limma::topTable(fit2, number=Inf, adjust.method="fdr")
      res$Feature <- rownames(res)
      sig    <- dplyr::filter(res, adj.P.Val < fdr_cutoff)
      if (!nrow(sig)) { message("    No sig features (limma)."); next }
      utils::write.csv(sig,
        file.path(g_dir, paste0("sig_features_", matrix_label, "_limma.csv")),
        row.names=FALSE)
      trend_df <- data.frame(Element=sig$Feature,
        `Observed Shift`=dplyr::case_when(sig$logFC>0~1L, sig$logFC<0~-1L, TRUE~0L),
        check.names=FALSE)
      utils::write.csv(trend_df,
        file.path(g_dir, paste0("Input_", matrix_label, "_limma.csv")),
        row.names=FALSE)
      feats <- intersect(trend_df$Element, colnames(Mat))
      if (length(feats))
        plot_group_box(Mat, feats, in_clade_ids, paste0("Clade_",clade_name),
                       file.path(g_dir, paste0(g_safe, "_limma_box.png")))

    } else if (lib_lower %in% c("maaslin2","maaslin")) {
      if (!requireNamespace("Maaslin2", quietly=TRUE)) {
        message("    Maaslin2 not installed; skipping."); next }
      metadata <- data.frame(SampleID=rownames(Mat), Group=grp,
                              stringsAsFactors=FALSE)
      rownames(metadata) <- metadata$SampleID
      ma_dir <- file.path(g_dir, "Maaslin2_out")
      if (dir.exists(ma_dir)) unlink(ma_dir, recursive=TRUE, force=TRUE)
      fit <- tryCatch(
        Maaslin2::Maaslin2(input_data=t(Mat), input_metadata=metadata,
          output=ma_dir, fixed_effects="Group", random_effects=NULL,
          normalization="NONE", transform="NONE", standardize=FALSE),
        error=function(e) { message("    Maaslin2 error: ", e$message); NULL })
      if (is.null(fit)) next
      res_tab <- if (!is.null(fit$results)) fit$results else {
        rf <- file.path(ma_dir,"all_results.tsv")
        if (file.exists(rf)) tryCatch(read.delim(rf,stringsAsFactors=FALSE),
                                      error=function(e) NULL) else NULL }
      if (is.null(res_tab) || !"qval" %in% names(res_tab)) {
        message("    No Maaslin2 results; skipping."); next }
      res_tab <- res_tab[res_tab$metadata=="Group",,drop=FALSE]
      sig     <- res_tab[!is.na(res_tab$qval) & res_tab$qval < fdr_cutoff,,drop=FALSE]
      if (!nrow(sig)) { message("    No sig features (MaAsLin2)."); next }
      utils::write.csv(sig,
        file.path(g_dir, paste0("sig_features_", matrix_label, "_Maaslin2.csv")),
        row.names=FALSE)
      trend_df <- data.frame(Element=sig$feature,
        `Observed Shift`=ifelse(sig$value>0,1L,ifelse(sig$value<0,-1L,0L)),
        check.names=FALSE)
      utils::write.csv(trend_df,
        file.path(g_dir, paste0("Input_", matrix_label, "_Maaslin2.csv")),
        row.names=FALSE)
      feats <- intersect(trend_df$Element, colnames(Mat))
      if (length(feats))
        plot_group_box(Mat, feats, in_clade_ids, paste0("Clade_",clade_name),
                       file.path(g_dir, paste0(g_safe, "_MaAsLin2_box.png")))
    }
  }
  invisible(TRUE)
}

# =============================================================
# 15. MAIN
# =============================================================
methods_path <- file.path(script_dir, "methods.txt")
cfg          <- read_methods(methods_path)

# --- CLI overrides: SampleBioShift.py passes per-block parameters ---
# Supported: --csv-file= --input-type= --distance= --norm= --sd= --phylo-tree=
cli_args <- commandArgs(trailingOnly=TRUE)
cli      <- list()
for (.a in cli_args) {
  if (grepl("^--", .a) && grepl("=", .a)) {
    .parts      <- strsplit(.a, "=", fixed=TRUE)[[1]]
    .key        <- sub("^--", "", .parts[1])
    .val        <- paste(.parts[-1], collapse="=")
    cli[[.key]] <- trimws(.val)
  }
}

filename    <- if (nchar(cli[["csv-file"]]    %||% "") > 0) cli[["csv-file"]]    else cfg[["Filename"]]
input_type  <- if (nchar(cli[["input-type"]]  %||% "") > 0) cli[["input-type"]]  else cfg[["Input_Type"]]
dist_method <- if (nchar(cli[["distance"]]    %||% "") > 0) cli[["distance"]]    else cfg[["Distance_Metric"]]
phylo_field <- if (nchar(cli[["phylo-tree"]]  %||% "") > 0) cli[["phylo-tree"]]  else cfg[["PhyloTree"]]
norm        <- if (nchar(cli[["norm"]]        %||% "") > 0) cli[["norm"]]        else cfg[["Normalization"]]
sd_val      <- suppressWarnings(as.numeric(
                 if (nchar(cli[["sd"]] %||% "") > 0) cli[["sd"]] else cfg[["SD"]]))

# --- Support/bootstrap replicate count (optional; defaults to the
#     original 100) -- shared knob for both mechanisms: drives the
#     number of SD-noise replicates for non-phylo methods, and the
#     number of raw-count bootstrap replicates for phylo methods. ---
support_field <- cli[["support"]] %||% cfg[["Support"]] %||% cfg[["Bootstrap"]] %||% cfg[["Support/bootstrap"]] %||% NA
support_reps  <- suppressWarnings(as.integer(support_field))
if (is.na(support_reps) || support_reps < 1) support_reps <- n_total - 1
n_total <- support_reps + 1L

# --- Differential analysis on/off gate (methods.txt "Differential
# analysis" field) -- a real gate again: even when target.txt's
# Y-targeted samples DO form a valid clade, "Differential analysis = No"
# skips run_trend_for_job() (the actual TREND/df step) for this
# block/job. Clade detection + the highlighted tree still run either way
# -- only the differential-analysis computation itself is gated. Defaults
# to TRUE (runs) if the field is absent, so blocks that don't set it
# behave as before.
.diff_flag_raw <- cfg[["Differential analysis"]]
# This line uses a single "#" for its trailing comment (e.g.
# "Yes # Chose Yes| No"), unlike the "##" the rest of methods.txt uses,
# so strip that off before comparing the value.
if (!is.null(.diff_flag_raw)) .diff_flag_raw <- trimws(sub("#.*$", "", .diff_flag_raw))
run_diff_analysis <- if (is.null(.diff_flag_raw) || !nchar(.diff_flag_raw)) {
  TRUE
} else {
  tolower(.diff_flag_raw) %in% c("yes","y","true","1")
}

# --- Required field checks ---
if (is.null(filename))    stop(
  "methods.txt missing: Filename\n",
  "  Fix: add  Filename = yourfile.csv")
if (is.null(input_type))  stop(
  "methods.txt missing: Input_Type\n",
  "  Fix: add  Input_Type = log10 | Count | Frequency")
if (is.null(dist_method)) stop(
  "methods.txt missing: Distance_Metric\n",
  "  Fix: add  Distance_Metric = Bray | Euclidean | UniFrac | UniFracW | MPD | MPDw | MNTD | MNTDw")
if (is.null(norm))        stop(
  "methods.txt missing: Normalization\n",
  "  Fix: add  Normalization = MinMax | Zscore | TSS | Hellinger | CLR | Rarefaction | None")
.is_phylo_dist <- !is.null(dist_method) && tolower(trimws(dist_method)) %in% tolower(phylo_metrics)
if (is.na(sd_val) && !.is_phylo_dist)        stop(
  "methods.txt: SD must be numeric (e.g. SD = 0 or SD = 0.01)\n",
  "  SD = 0   -> Control only (no noise)\n",
  "  SD = 0.01 -> Control + noise test\n",
  "  (SD is not applicable for phylogenetic Distance_Metric values --\n",
  "   set SD = NA there; support comes from raw-count bootstrap instead)")
if (is.na(sd_val) && .is_phylo_dist) sd_val <- 0   # placeholder only; unused for phylo (see below)

# --- Valid value checks ---
.valid_input_types  <- c("log10","count","frequency")
.valid_norms        <- c("minmax","zscore","tss","hellinger","clr","rarefaction","none","frequency")
.valid_nonphy_dists <- c("bray","euclidean")
.valid_phy_dists    <- c("unifrac","unifracw","mpd","mpdw","mntd","mntdw")
.valid_dists        <- c(.valid_nonphy_dists, .valid_phy_dists)

.input_lc <- tolower(trimws(input_type))
.norm_lc  <- tolower(trimws(norm))
.dist_lc  <- tolower(trimws(dist_method))

if (!.input_lc %in% .valid_input_types)
  stop("Unknown Input_Type: '", input_type, "'\n",
       "  Valid values: log10 | Count | Frequency")

if (!.norm_lc %in% .valid_norms)
  stop("Unknown Normalization: '", norm, "'\n",
       "  Valid values: MinMax | Zscore | TSS | Hellinger | CLR | Rarefaction | None")

if (!.dist_lc %in% .valid_dists)
  stop("Unknown Distance_Metric: '", dist_method, "'\n",
       "  Non-phylo : Bray | Euclidean\n",
       "  Phylo     : UniFrac | UniFracW | MPD | MPDw | MNTD | MNTDw")

# --- Normalization × Input_Type compatibility ---
if (.input_lc == "log10" && !.norm_lc %in% c("minmax","zscore","none"))
  stop("Normalization='", norm, "' is not valid for Input_Type=log10\n",
       "  Valid for log10: MinMax | Zscore | None\n",
       "  Note: Zscore requires Distance_Metric = Euclidean")

if (.input_lc == "count" && !.norm_lc %in% c("tss","hellinger","clr","rarefaction","none"))
  stop("Normalization='", norm, "' is not valid for Input_Type=Count\n",
       "  Valid for Count: TSS | Hellinger | CLR | Rarefaction | None\n",
       "  Note: CLR requires Distance_Metric = Euclidean")

if (.input_lc == "frequency" && !.norm_lc %in% c("none","hellinger","zscore","frequency"))
  stop("Normalization='", norm, "' is not valid for Input_Type=Frequency\n",
       "  Valid for Frequency: None (Bray) | Hellinger (Bray) | Zscore (Euclidean only)")

# --- Normalization × Distance_Metric compatibility ---
if (.dist_lc == "bray" && .norm_lc %in% c("clr","zscore"))
  stop("Invalid combination: Normalization=", norm, " + Distance_Metric=Bray\n",
       "  ", norm, " produces negative values — Bray-Curtis requires values ≥ 0\n",
       "  Fix A: change Distance_Metric = Euclidean  (keep ", norm, ")\n",
       "  Fix B: change Normalization = TSS | Hellinger | None  (keep Bray)")

if (.dist_lc == "bray" && .input_lc == "log10" && .norm_lc == "minmax") {
  # Valid — per-sample MinMax on log10 + Bray is the intended combo; no warning
}

if (.dist_lc == "euclidean" && .norm_lc %in% c("tss","hellinger","rarefaction"))
  warning("Normalization=", norm, " with Euclidean distance may give poor results.\n",
          "  Recommended for Euclidean: CLR (Count) | Zscore (log10/Frequency)")

if (.dist_lc == "bray" && .norm_lc == "zscore" && .input_lc == "frequency")
  stop("Invalid combination: Zscore + Bray for Frequency data\n",
       "  Zscore produces negatives — Bray requires values ≥ 0\n",
       "  Fix: use Zscore + Euclidean  OR  None/Hellinger + Bray")

# --- Phylo distance requires a tree ---
if (.dist_lc %in% .valid_phy_dists) {
  .phylo_val_check <- tolower(trimws(phylo_field %||% "none"))
  if (.phylo_val_check %in% c("none","0","na",""))
    stop("Distance_Metric=", dist_method, " is a phylogenetic method and requires a tree.\n",
         "  Fix: set PhyloTree = yourfile.nwk in methods.txt")
}

# Derive noise type and differential library from Input_Type
input_lc     <- tolower(trimws(input_type))
noise_type   <- if (input_lc == "log10") "Gaussian" else if (input_lc == "count") "Poisson" else "Dirichlet"
library_name <- if (input_lc == "log10") "limma"     else "MaAsLin2"

message("\n=== SampleTree Pipeline ===")
message("Filename     : ", filename)
message("Input_Type   : ", input_type)
message("Distance     : ", dist_method)
message("Normalization: ", norm)
message("PhyloTree    : ", phylo_field)
if (.is_phylo_dist) {
  message("SD           : NA (phylogenetic method -- support from raw-count bootstrap, not SD noise)")
} else {
  message("SD           : ", sd_val)
}
message("Noise        : ", if (.is_phylo_dist) "None (raw-count bootstrap)" else noise_type,
        "  [derived from Input_Type]")
message("Library      : ", library_name, "  [derived from Input_Type]")
message("Replicates   : ", n_total, "  (1 main + ", n_total - 1,
        if (.is_phylo_dist) " bootstrap)" else " noise)")

# --- Differential analysis (TREND) is the clade-based path below:
# whenever target.txt's Y-targeted samples form a real clade
# (min_targeted/max_others from sampletree_control.txt) AND methods.txt's
# "Differential analysis" field isn't explicitly "No" (run_diff_analysis,
# computed above), run_trend_for_job() runs automatically on the
# original/unnormalized input. "Differential analysis = No" stops the
# pipeline after the plain sample tree -- no clade detection, no
# highlighted tree, no TREND. "Input df"/"Profile_Library" remain unused;
# the input is always the raw/unnormalized file and the library is
# always derived from Input_Type (see library_name above). ---

# --- Load target.txt (Y/N sample labels + optional Patient column) ---
target_yes        <- character(0)
sample_to_patient <- NULL
target_file <- file.path(script_dir, "target.txt")
if (file.exists(target_file)) {
  tdf <- read.delim(target_file, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  if (all(c("Sample","Target") %in% names(tdf))) {
    tdf$Target <- toupper(trimws(tdf$Target))
    target_yes <- tdf$Sample[tdf$Target == "Y"]
    message("Target (Y)   : ", length(target_yes), " samples")
  }
  if ("Patient" %in% names(tdf)) {
    tdf$Patient <- as.character(trimws(tdf$Patient))
    sample_to_patient <- setNames(tdf$Patient, tdf$Sample)
    message("Patient IDs  : ", length(unique(tdf$Patient)), " unique patients found")
  }
} else {
  message("target.txt not found — clade detection and TREND will be skipped.")
}

# --- Load sampletree_control.txt (clade thresholds) ---
thr_min_targeted <- 2L
thr_max_clade    <- 50L
thr_max_others   <- 3L
thr_assign_policy <- "best"
ctrl_file <- file.path(script_dir, "sampletree_control.txt")
if (file.exists(ctrl_file)) {
  cfg_c <- read_methods(ctrl_file)
  if (!is.null(cfg_c[["min_targeted"]]))  thr_min_targeted  <- as.integer(cfg_c[["min_targeted"]])
  if (!is.null(cfg_c[["max_clade_size"]])) thr_max_clade    <- as.integer(cfg_c[["max_clade_size"]])
  if (!is.null(cfg_c[["max_others"]]))     thr_max_others   <- as.integer(cfg_c[["max_others"]])
  if (!is.null(cfg_c[["assign_policy"]])) {
    ap <- tolower(trimws(cfg_c[["assign_policy"]]))
    if (ap %in% c("best","first","largest","smallest")) thr_assign_policy <- ap
    else warning("Unknown assign_policy '", ap, "' — using 'best'")
  }
  message("Thresholds   : min_targeted=", thr_min_targeted,
          ", max_clade=", thr_max_clade, ", max_others=", thr_max_others,
          ", assign_policy=", thr_assign_policy)
}

# --- Load input CSV (from Inputs/ subfolder) ---
csv_path <- file.path(script_dir, "Inputs", filename)
if (!file.exists(csv_path)) stop("Input CSV not found: ", csv_path)
raw <- as.matrix(read.csv(csv_path, row.names=1,
                           check.names=FALSE, stringsAsFactors=FALSE))
raw[] <- as.numeric(raw)   # convert in-place to preserve dimnames
message("Data loaded  : ", nrow(raw), " samples x ", ncol(raw), " features")

# --- Remove all-zero columns (universal: applies to every metric) ---
zero_cols <- colSums(raw, na.rm=TRUE) == 0
if (any(zero_cols)) {
  message("Removed ", sum(zero_cols), " all-zero column(s): ",
          paste(colnames(raw)[zero_cols], collapse=", "))
  raw <- raw[, !zero_cols, drop=FALSE]
  message("Remaining    : ", nrow(raw), " samples x ", ncol(raw), " features")
}

# --- Differential input = normalized data (same transformation as tree) ---
otu_diff <- apply_normalization(raw, norm)
message("Diff matrix  : ", nrow(otu_diff), " x ", ncol(otu_diff),
        "  [", norm, "-normalized]")

# --- Output root folder: Outputs/{FileBase}/ ---
file_base <- tools::file_path_sans_ext(basename(filename))
out_base  <- file.path(script_dir, "Outputs", file_base)
message("Output root  : ", out_base)

# --- Phylo tree ---
use_phy  <- FALSE
ref_tree <- NULL
coph_mat <- NULL

phylo_val <- trimws(phylo_field %||% "None")
# Guard: if value is not a recognised "None" token and doesn't end in .nwk,
# it is almost certainly a copy-paste mistake (e.g. PhyloTree = Bray).
if (!phylo_val %in% c("None","none","0","NA","") &&
    !grepl("\\.nwk$", phylo_val, ignore.case=TRUE)) {
  warning("PhyloTree value '", phylo_val,
          "' does not look like a .nwk file — treating as None. ",
          "Set PhyloTree = None or PhyloTree = yourfile.nwk in methods.txt.")
  phylo_val <- "None"
}
if (!phylo_val %in% c("None","none","0","NA","")) {
  tree_path <- file.path(script_dir, "Inputs", phylo_val)
  if (!file.exists(tree_path)) stop("PhyloTree file not found: ", tree_path)
  ref_tree <- ape::read.tree(tree_path)
  if (dist_method %in% phylo_metrics) {
    aligned  <- align_to_tree(raw, ref_tree, label=filename)
    raw      <- aligned$mat
    ref_tree <- aligned$tree
    coph_mat <- cophenetic(ref_tree)
    use_phy  <- TRUE
    message("Phylo tree   : ", tree_path,
            " (", length(ref_tree$tip.label), " taxa after alignment)")
  }
}

if (dist_method %in% phylo_metrics && !use_phy)
  stop("Distance '", dist_method,
       "' requires a tree. Set PhyloTree = <file>.nwk in methods.txt")

# =============================================================
# Control run (SD = 0)
# =============================================================
# Phylo methods: skip the Normalization/SD0 nesting entirely -- there's
# no SD sweep (support is a fixed raw-count bootstrap), and only one
# normalization runs per job, so that nesting is just clutter. Output
# goes directly to Outputs/{FileBase}/{Distance}/. Non-phylo methods
# keep the full nesting unchanged, since those ARE swept across
# multiple Normalization/Distance/SD combinations for the same file.
ctrl_outdir <- if (.is_phylo_dist) {
  file.path(out_base, dist_method)
} else {
  file.path(out_base, norm, dist_method, "SD0")
}
message("\n--- Running Control (SD=0, n=", n_total, ") ---")

ctrl_res  <- run_condition(
  mat=raw, sd=0, outdir=ctrl_outdir,
  norm=norm, dist_method=dist_method,
  use_phy=use_phy, ref_tree=ref_tree, coph=coph_mat,
  input_type=input_type, target_yes=target_yes,
  sample_to_patient=sample_to_patient)
message("Control done. Mean support = ", round(ctrl_res$mean_sup, 1), "%")

ctrl_coph   <- cophenetic(ctrl_res$main_tr)
run_results <- list(ctrl_res)

# --- Per-patient clustering (TrueFalse.csv) -- REMOVED for this pipeline.
# This SampleBioShift copy answers one question only: do target.txt's
# Y-targeted samples cluster together, real-world-blind to which patient
# each one came from? Patient identity is not part of that question, so
# build_patient_truefalse_table() (and its TrueFalse.csv/patient_clusters/
# output) no longer runs here. The target_yes-gated block right below --
# find_clades_global() + plot_highlighted_tree() + run_trend_for_job() --
# is the entire clustering/differential-analysis path now.


# --- Highlighted tree + Clades + TREND (on SD0 main tree, run once) --
# Gated on BOTH a real target_yes set (global Y/N clade detection needs
# one to mean anything) AND methods.txt's "Differential analysis" field:
# "No" means the plain sample tree (support_tree.nwk/pdf/png above) is
# ALL that gets built for this block -- no clades/, no
# support_tree_highlighted.jpeg, no trend_outputs/. "Yes" (or the field
# absent) runs the full path down to the sig-element output BioShift
# consumes.
if (length(target_yes) > 0 && run_diff_analysis) {
  message("\n--- Clade detection ---")
  clade_info <- find_clades_global(
    tree=ctrl_res$main_tr, tag=file_base, method=dist_method,
    outdir=ctrl_outdir,
    min_targeted=thr_min_targeted, max_clade=thr_max_clade,
    max_others=thr_max_others, targeted_set=target_yes,
    assign_policy=thr_assign_policy)

  plot_highlighted_tree(
    tree=ctrl_res$main_tr,
    targeted_set=target_yes,
    contamination_set=clade_info$contam_tips,
    outfile=file.path(ctrl_outdir, "support_tree_highlighted.jpeg"),
    tag=file_base, method=dist_method,
    clade_member_targets=unique(unlist(lapply(clade_info$kept, function(ci) ci$tgt))))
  message("  Highlighted tree saved.")

  if (!is.null(raw)) {
    message("\n--- TREND analysis (", library_name, ", original/unnormalized input) ---")
    run_trend_for_job(
      otu=raw,
      clade_info=clade_info,
      outdir=ctrl_outdir,
      tag=file_base, method=dist_method,
      matrix_label=file_base,
      library_name=library_name,
      target_yes=target_yes)
    # finalize_all_merges()/collect_merged_inputs() (Input_MERGED_*.csv,
    # Merged_Group_Inputs/) -- REMOVED. SampleBioShift.py's own
    # merge_input_csvs() rebuilds the same merge directly from each
    # group_*/Input_*.csv (explicitly skipping any Input_MERGED_* it
    # finds), so this R-side copy was never actually consumed.
    message("  TREND complete.")
  }
} else if (length(target_yes) > 0 && !run_diff_analysis) {
  message("\n--- Differential analysis = No -- stopping after the plain sample tree. ",
          "(target.txt has a real Y-target, but clade detection/highlighting/TREND are skipped.) ---")
}

# Write control log.txt
write_run_log(
  outdir  = ctrl_outdir,
  params  = list(filename=filename, input_type=input_type,
                 dist_method=dist_method, norm=norm,
                 sd=if (.is_phylo_dist) "NA (phylo bootstrap)" else 0,
                 noise_type=if (.is_phylo_dist) "None (raw-count bootstrap)" else noise_type,
                 library_name=library_name,
                 n_total=n_total, n_samples=nrow(raw),
                 n_features=ncol(raw), n_target=length(target_yes)),
  res         = ctrl_res,
  clade_info  = if (exists("clade_info")) clade_info else NULL)
message("  log.txt written -> ", ctrl_outdir)

# =============================================================
# Test run (SD > 0) -- non-phylo methods only. Phylo methods already
# got their full support estimate from the raw-count bootstrap inside
# the Control run above; there's no separate "noise test" condition
# for them since SD doesn't apply.
# =============================================================
if (.is_phylo_dist) {
  message("\n(No SD Test run -- not applicable for phylogenetic Distance_Metric.)")
} else if (sd_val > 0) {
  sd_label    <- formatC(sd_val, format="g")   # e.g. "0.01"
  test_outdir <- file.path(out_base, norm, dist_method,
                            paste0("SD", sd_label))
  message("\n--- Running Test (SD=", sd_val, ", n=", n_total, ") ---")

  test_res <- run_condition(
    mat=raw, sd=sd_val, outdir=test_outdir,
    norm=norm, dist_method=dist_method,
    use_phy=use_phy, ref_tree=ref_tree, coph=coph_mat,
    input_type=input_type, target_yes=target_yes,
    sample_to_patient=sample_to_patient)
  message("Test done. Mean support = ", round(test_res$mean_sup, 1), "%")

  # RF distance: test replicates vs control tree
  reps       <- test_res$rep_trees_kept
  rf_vs_ctrl <- if (length(reps) > 0) {
    all_tr <- structure(c(list(ctrl_res$main_tr), reps), class="multiPhylo")
    rf2    <- suppressWarnings(as.matrix(ape::dist.topo(all_tr, method="PH85")))
    mean(rf2[1, -1])
  } else NA

  # Cophenetic correlation: test replicates vs control tree
  cc <- if (length(reps) > 0) {
    cors <- sapply(reps, function(t2) {
      tc2 <- tryCatch(cophenetic(t2), error=function(e) NULL)
      if (is.null(tc2)) return(NA)
      cmn <- intersect(rownames(ctrl_coph), rownames(tc2))
      cor(as.vector(ctrl_coph[cmn, cmn]),
          as.vector(tc2[cmn, cmn]),
          use="pairwise.complete.obs")
    })
    mean(cors, na.rm=TRUE)
  } else NA

  # --- Per-patient clustering (TrueFalse.csv, Test) -- REMOVED, same
  # reasoning as the Control block above: patient identity isn't part of
  # this pipeline's question anymore. ---

  # --- Highlighted tree + Clades + TREND (Test main tree) -- same gate
  # as the Control block above (target_yes AND run_diff_analysis). ---
  if (length(target_yes) > 0 && run_diff_analysis) {
    message("\n--- Clade detection (Test) ---")
    clade_info_test <- find_clades_global(
      tree=test_res$main_tr, tag=file_base, method=dist_method,
      outdir=test_outdir,
      min_targeted=thr_min_targeted, max_clade=thr_max_clade,
      max_others=thr_max_others, targeted_set=target_yes,
      assign_policy=thr_assign_policy)

    plot_highlighted_tree(
      tree=test_res$main_tr,
      targeted_set=target_yes,
      contamination_set=clade_info_test$contam_tips,
      outfile=file.path(test_outdir, "support_tree_highlighted.jpeg"),
      tag=file_base, method=dist_method,
      clade_member_targets=unique(unlist(lapply(clade_info_test$kept, function(ci) ci$tgt))))
    message("  Highlighted tree saved (Test).")

    if (!is.null(raw)) {
      message("\n--- TREND analysis (Test, ", library_name, ", original/unnormalized input) ---")
      run_trend_for_job(
        otu=raw,
        clade_info=clade_info_test,
        outdir=test_outdir,
        tag=file_base, method=dist_method,
        matrix_label=file_base,
        library_name=library_name,
        target_yes=target_yes)
      # finalize_all_merges()/collect_merged_inputs() -- REMOVED, same
      # reasoning as the Control block above.
      message("  TREND complete (Test).")
    }
  } else if (length(target_yes) > 0 && !run_diff_analysis) {
    message("\n--- Differential analysis = No -- stopping after the plain sample tree (Test). ---")
  }

  # Write test log.txt
  write_run_log(
    outdir  = test_outdir,
    params  = list(filename=filename, input_type=input_type,
                   dist_method=dist_method, norm=norm, sd=sd_val,
                   noise_type=noise_type, library_name=library_name,
                   n_total=n_total, n_samples=nrow(raw),
                   n_features=ncol(raw), n_target=length(target_yes)),
    res          = test_res,
    clade_info   = if (exists("clade_info_test")) clade_info_test else NULL,
    rf_vs_ctrl   = rf_vs_ctrl,
    coph_cor     = cc)
  message("  log.txt written -> ", test_outdir)

  run_results <- c(run_results, list(test_res))
}

message("\n=== Complete ===")
message("Noise type : ", if (.is_phylo_dist) "None (raw-count bootstrap)" else noise_type)
message("Library    : ", library_name)
message("Outputs in : ", out_base)
