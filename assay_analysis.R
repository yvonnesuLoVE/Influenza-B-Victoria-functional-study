#!/usr/bin/env Rscript
# Combined pH, thermal stability and receptor-binding workflow.
# Usage: Rscript scripts/assay_analysis.R config/analysis_config.R ph|temperature|receptor|all

fail <- function(...) stop(..., call. = FALSE)
need_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) fail("Install required R packages: ", paste(missing, collapse = ", "))
}
read_assay <- function(path, replicate = TRUE) {
  if (!file.exists(path)) fail("Input not found: ", path)
  dat <- openxlsx::read.xlsx(path, rowNames = TRUE, check.names = FALSE)
  if (!"Type" %in% names(dat) && "tYPE" %in% names(dat)) names(dat)[names(dat) == "tYPE"] <- "Type"
  required <- c("Type", if (replicate) "Replicate")
  if (!all(required %in% names(dat))) fail(path, " requires columns: ", paste(required, collapse = ", "))
  if (anyNA(dat$Type) || any(!nzchar(trimws(as.character(dat$Type))))) fail("Missing Type in ", path)
  if (replicate && (anyNA(dat$Replicate) || any(!nzchar(trimws(as.character(dat$Replicate)))))) fail("Missing Replicate in ", path)
  dat$Type <- as.character(dat$Type)
  dat
}
check_numeric <- function(dat, cols) {
  if (!length(cols) || !all(cols %in% names(dat))) fail("Required measurement columns are absent.")
  if (any(!vapply(dat[cols], is.numeric, logical(1)))) fail("Measurement columns must be numeric.")
  if (any(vapply(dat[cols], function(x) any(is.infinite(x)), logical(1)))) fail("Infinite input measurements.")
  if (any(vapply(dat[cols], function(x) any(x < 0, na.rm = TRUE), logical(1)))) fail("Negative input measurements.")
}
transform_values <- function(x, method) {
  switch(method, identity = x, log2 = suppressWarnings(log2(x)),
         log10 = suppressWarnings(log10(x)), fail("Unknown transformation: ", method))
}
finite_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")
save_plot <- function(p, path, width = 12, height = 8) {
  ggplot2::ggsave(path, p, width = width, height = height, units = "in", limitsize = FALSE)
}
metadata_for <- function(types, cfg) {
  if (is.null(cfg$virus_metadata)) {
    return(data.frame(Type = types, Clade = "All", Colour = grDevices::hcl.colors(length(types), "Dark 3"), Order = seq_along(types)))
  }
  m <- utils::read.csv(cfg$virus_metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Type", "Clade", "Colour", "Order") %in% names(m))) fail("Metadata requires Type, Clade, Colour, Order.")
  if (anyDuplicated(m$Type)) fail("Duplicate Type in virus metadata.")
  if (anyNA(match(types, m$Type))) fail("Virus metadata is missing: ", paste(setdiff(types, m$Type), collapse = ", "))
  m <- m[match(types, m$Type), ]
  if (anyNA(m) || any(!nzchar(m$Clade)) || !is.numeric(m$Order)) fail("Incomplete virus metadata.")
  tryCatch(grDevices::col2rgb(m$Colour), error = function(e) fail("Invalid metadata Colour."))
  m[order(m$Order), ]
}
nonzero_range <- function(x, include_zero = FALSE) {
  x <- x[is.finite(x)]
  if (include_zero) x <- c(0, x)
  if (!length(x)) return(c(0, 1))
  r <- range(x)
  if (diff(r) == 0) r <- r + c(-1, 1) * max(0.01, abs(r[1]) * 0.05)
  r
}

# The model is identical for both stability assays except for the condition labels.
# Type:Replicate makes replicate identifiers unique within each virus.
fit_stability <- function(long, out, prefix) {
  d <- long[is.finite(long$Value), c("Type", "Replicate", "Condition", "Value")]
  d$Type <- factor(d$Type); d$Replicate <- factor(d$Replicate)
  d$Condition <- droplevels(d$Condition)
  if (nlevels(d$Type) < 2 || nlevels(d$Condition) < 2) fail("Mixed model needs at least two viruses and conditions.")
  model <- lmerTest::lmer(Value ~ Type * Condition + (1 | Type:Replicate), data = d)
  saveRDS(model, file.path(out, paste0(prefix, "_model.rds")))
  capture.output(summary(model), file = file.path(out, paste0(prefix, "_model_summary.txt")))
  a <- as.data.frame(anova(model)); a$Term <- rownames(a)
  write_csv(a, file.path(out, paste0(prefix, "_anova.csv")))
  em_type <- emmeans::emmeans(model, ~ Type | Condition)
  em_cond <- emmeans::emmeans(model, ~ Condition | Type)
  ct <- as.data.frame(emmeans::contrast(em_type, "pairwise", adjust = "tukey"))
  cc <- as.data.frame(emmeans::contrast(em_cond, "pairwise", adjust = "tukey"))
  write_csv(ct, file.path(out, paste0(prefix, "_types_within_condition.csv")))
  write_csv(cc, file.path(out, paste0(prefix, "_conditions_within_type.csv")))
  # Build matrices from emmeans pair ordering, without splitting virus labels.
  # Virus names may themselves contain hyphens or spaces.
  matrices <- function(emm, pairs, slices, prefix) {
    labels <- as.data.frame(emm)
    result <- list()
    for (slice in unique(labels[[slices]])) {
      d_slice <- labels[as.character(labels[[slices]]) == as.character(slice), , drop = FALSE]
      levels <- as.character(d_slice[[if (slices == "Condition") "Type" else "Condition"]])
      # emmeans pairwise ordering is combn(levels, 2), in factor order.
      z <- pairs[as.character(pairs[[slices]]) == as.character(slice), , drop = FALSE]
      if (length(levels) < 2) next
      ij <- utils::combn(seq_along(levels), 2)
      if (nrow(z) != ncol(ij)) fail("Unexpected emmeans contrast ordering.")
      mat <- matrix(NA_real_, length(levels), length(levels), dimnames = list(levels, levels))
      mat[t(ij)] <- z$p.value; mat[t(ij[2:1, , drop = FALSE])] <- z$p.value
      result[[paste0(prefix, slice)]] <- data.frame(Level = levels, mat, check.names = FALSE)
    }
    result
  }
  sheets <- c(matrices(em_type, ct, "Condition", "Types@"), matrices(em_cond, cc, "Type", "Conditions@"))
  # Prefix sheets with an index so truncation/sanitization cannot create duplicate names.
  names(sheets) <- vapply(seq_along(sheets), function(i) {
    nm <- names(sheets)[i]
    for (bad in c("[", "]", ":", "*", "?", "/", "\\", "'")) nm <- gsub(bad, "_", nm, fixed = TRUE)
    paste0(sprintf("%03d_", i), substr(nm, 1, 27))
  }, character(1))
  openxlsx::write.xlsx(sheets, file.path(out, paste0(prefix, "_pairwise_pvalues.xlsx")), overwrite = TRUE)
  diagnostics <- data.frame(n_observations = stats::nobs(model), singular = lme4::isSingular(model),
                            dropped_nonfinite = nrow(long) - nrow(d))
  write_csv(diagnostics, file.path(out, paste0(prefix, "_model_diagnostics.csv")))
  if (lme4::isSingular(model)) warning(prefix, ": singular model; inspect diagnostics before interpretation.", call. = FALSE)
}

panel_definitions <- function(M, assay, meta, cfg) {
  if (!is.null(cfg$panel_membership)) {
    p <- utils::read.csv(cfg$panel_membership, stringsAsFactors = FALSE)
    required <- c("Assay", "Panel", "Type", "Order", "AbsoluteDelta")
    if (!all(required %in% names(p))) fail("Panel membership requires: ", paste(required, collapse = ", "))
    p <- p[p$Assay == assay, required]
    p$AbsoluteDelta <- as.logical(p$AbsoluteDelta)
    if (anyNA(p) || any(!p$Type %in% rownames(M))) fail("Invalid panel membership for ", assay)
    if (anyDuplicated(p[c("Panel", "Type")])) fail("Duplicate virus within a panel.")
    return(p)
  }
  if (isTRUE(cfg$use_legacy_panels)) {
    if (nrow(M) < 19) fail("Legacy panel indices require at least 19 alphabetically ordered Types.")
    indices <- list(`Pair 2-3` = c(18, 10), `Pair 2-6` = c(18, 5),
                    `Pair 7-13` = c(1, 2), `Pair 2-13` = c(18, 2),
                    `Clade V1A.1` = c(17,18,10,3,5,1,4,11,13,15,2,12,14),
                    `Clade V1A.3` = if (assay == "ph") c(17,19,9,16,6,8,7) else c(19,9,16,6,8,7))
    return(do.call(rbind, lapply(seq_along(indices), function(i) {
      ids <- indices[[i]]
      data.frame(Assay = assay, Panel = names(indices)[i], Type = rownames(M)[ids],
                 Order = seq_along(ids), AbsoluteDelta = i > 4)
    })))
  }
  data.frame(Assay = assay, Panel = meta$Clade, Type = meta$Type, Order = meta$Order, AbsoluteDelta = TRUE)
}

plot_stability <- function(long, assay, meta, cfg, out) {
  # Explicit alphabetical Type order matches dplyr's original grouped summaries.
  types <- sort(unique(as.character(long$Type)))
  conditions <- levels(long$Condition)
  summary <- long %>% dplyr::group_by(Type, Condition, NumericCondition) %>%
    dplyr::summarise(Mean = finite_mean(Value), .groups = "drop")
  M <- matrix(NA_real_, length(types), length(conditions), dimnames = list(types, conditions))
  M[cbind(match(summary$Type, types), match(summary$Condition, conditions))] <- summary$Mean
  if (!any(is.finite(M))) fail("No finite values to plot for ", assay)
  write_csv(data.frame(Type = rownames(M), M, check.names = FALSE), file.path(out, paste0(assay, "_heatmap_values.csv")))
  targets <- if (assay == "ph") cfg$ph_targets else cfg$temperature_targets
  nums <- vapply(conditions, function(x) unique(long$NumericCondition[as.character(long$Condition) == x])[1], numeric(1))
  idx <- vapply(targets, function(x) which.min(abs(nums - x)), integer(1))
  if (length(unique(idx)) != 2) fail("Target conditions resolve to the same measurement column.")
  lo <- idx[which.min(nums[idx])]; hi <- idx[which.max(nums[idx])]
  write_csv(data.frame(Requested = targets, Used = nums[idx]), file.path(out, paste0(assay, "_delta_targets.csv")))
  delta <- M[, hi] - M[, lo]
  panels <- panel_definitions(M, assay, meta, cfg)
  write_csv(panels, file.path(out, paste0(assay, "_resolved_panels.csv")))
  method <- if (assay == "ph") cfg$ph_plot_transform else cfg$temperature_plot_transform
  label <- if (method == "identity") "Mean TCID" else paste0("Mean ", method, "(TCID)")
  palette <- if (assay == "ph") rev(c("#081d58", "#225ea8", "#41b6c4", "#a1dab4", "#ffffcc")) else rev(viridisLite::magma(8, begin = .15, end = 1))
  colour <- circlize::colorRamp2(seq(nonzero_range(M)[1], nonzero_range(M)[2], length.out = length(palette)), palette)
  hs <- lapply(unique(panels$Panel), function(panel) {
    p <- panels[panels$Panel == panel, ]; p <- p[order(p$Order), ]
    if (length(unique(p$AbsoluteDelta)) != 1) fail("Panel AbsoluteDelta must be constant within a panel.")
    block <- M[p$Type, , drop = FALSE]
    values <- delta[p$Type]
    absolute <- p$AbsoluteDelta[1]
    if (absolute) values <- abs(values)
    delta_label <- paste0(if (absolute) "|" else "", nums[hi], " - ", nums[lo], if (absolute) "|" else "")
    right <- ComplexHeatmap::rowAnnotation(Delta = ComplexHeatmap::anno_barplot(
      values, ylim = nonzero_range(values, TRUE), gp = grid::gpar(fill = "#4C72B0", col = NA),
      border = FALSE), annotation_label = delta_label, width = grid::unit(26, "mm"))
    top <- bottom <- NULL
    if (nrow(block) == 2) {
      top <- ComplexHeatmap::HeatmapAnnotation(Lines = ComplexHeatmap::anno_lines(
        t(block), ylim = nonzero_range(block), gp = grid::gpar(col = c("#1F77B4", "#D62728")),
        add_points = TRUE), height = grid::unit(16, "mm"))
      bottom <- ComplexHeatmap::HeatmapAnnotation(`Row difference` = ComplexHeatmap::anno_barplot(
        as.numeric(block[1, ] - block[2, ]), ylim = nonzero_range(block[1, ] - block[2, ], TRUE),
        gp = grid::gpar(fill = "#4C72B0", col = NA), border = FALSE), height = grid::unit(16, "mm"))
    }
    ComplexHeatmap::Heatmap(block, name = paste0(assay, "_", match(panel, unique(panels$Panel))),
      heatmap_legend_param = list(title = label), col = colour,
      cluster_rows = FALSE, cluster_columns = FALSE, row_title = panel,
      right_annotation = right, top_annotation = top, bottom_annotation = bottom,
      row_names_side = "left", column_names_rot = 0, na_col = "white",
      rect_gp = grid::gpar(col = "black", lwd = .5),
      width = grid::unit(ncol(block) * 5, "mm"), height = grid::unit(nrow(block) * 5, "mm"))
  })
  if (!length(hs)) fail("No panels defined for ", assay)
  combined <- Reduce(function(a, b) ComplexHeatmap::`%v%`(a, b), hs)
  grDevices::pdf(file.path(out, paste0(assay, "_heatmaps.pdf")), width = 12,
                 height = max(7, nrow(panels) * 5 / 25.4 + length(hs) * 1.8))
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(combined, ht_gap = grid::unit(6, "mm"), merge_legend = TRUE)
}

run_stability <- function(assay, cfg) {
  need_packages(c("ComplexHeatmap", "circlize", "viridisLite", "lmerTest", "lme4", "emmeans"))
  dat <- read_assay(if (assay == "ph") cfg$ph_file else cfg$temperature_file)
  if (assay == "ph") {
    cols <- c("pH5", "pH5.6", "pH6", "pH6.5", "pH7")
    nums <- as.numeric(sub("^pH", "", cols))
  } else {
    cols <- grep("^\\s*[0-9]+(?:\\.[0-9]+)?\\s*°?C\\s*$", names(dat), value = TRUE, perl = TRUE)
    nums <- as.numeric(gsub("[^0-9.]", "", cols))
    o <- order(nums); cols <- cols[o]; nums <- nums[o]
  }
  check_numeric(dat, cols)
  if (anyDuplicated(nums)) fail("Duplicate numeric conditions.")
  # pH heatmaps in the final source were descending; models retain ascending order.
  long <- tidyr::pivot_longer(dat, dplyr::all_of(cols), names_to = "Condition", values_to = "Raw")
  long$NumericCondition <- nums[match(long$Condition, cols)]
  long$Condition <- factor(long$Condition, levels = cols)
  long$Value <- transform_values(long$Raw, "log2")
  override <- if (assay == "temperature") long$NumericCondition %in% cfg$temperature_force_zero else rep(FALSE, nrow(long))
  if (assay == "temperature" && any(!cfg$temperature_force_zero %in% nums)) fail("Forced-zero temperature absent from input.")
  long$Value[override] <- 0
  long$Value[!is.finite(long$Value)] <- NA_real_
  dir.create(out <- file.path(cfg$output_dir, assay), recursive = TRUE, showWarnings = FALSE)
  write_csv(long, file.path(out, paste0(assay, "_model_input.csv")))
  fit_stability(long, out, assay)
  plot_long <- long
  method <- if (assay == "ph") cfg$ph_plot_transform else cfg$temperature_plot_transform
  plot_long$Value <- transform_values(plot_long$Raw, method)
  plot_long$Value[override] <- 0
  # Retain the final source's thermal plotting convention: log2(0) -> 0.
  # This is exported separately from the model input, where non-finite -> NA.
  if (assay == "temperature") plot_long$Value[is.infinite(plot_long$Value) & plot_long$Value < 0] <- 0
  plot_long$Value[!is.finite(plot_long$Value)] <- NA_real_
  if (assay == "ph") plot_long$Condition <- factor(plot_long$Condition, levels = rev(cols))
  write_csv(plot_long, file.path(out, paste0(assay, "_plot_input.csv")))
  plot_stability(plot_long, assay, metadata_for(sort(unique(dat$Type)), cfg), cfg, out)
}

run_receptor <- function(cfg) {
  need_packages("patchwork")
  dat <- read_assay(cfg$receptor_file, replicate = FALSE)
  # The first spreadsheet column was row names in the source. Each row is one
  # replicate; no pairing between distinct row IDs is invented here.
  dat$Replicate_ID <- rownames(dat)
  cols <- grep("^[36]SLN_", names(dat), value = TRUE)
  check_numeric(dat, cols)
  long <- tidyr::pivot_longer(dat, dplyr::all_of(cols), names_to = "Measurement", values_to = "Absorbance")
  long$Ligand <- sub("_.*$", "", long$Measurement)
  # Strip the ligand prefix BEFORE parsing numbers (3SLN_64 is 64, not 3).
  long$HAU <- suppressWarnings(as.numeric(sub("^[36]SLN_", "", long$Measurement)))
  if (any(!is.finite(long$HAU)) || any(long$HAU <= 0)) fail("Receptor columns must be 3SLN_<HAU> or 6SLN_<HAU>.")
  measured <- unique(long[c("Ligand", "HAU")])
  if (!setequal(unique(long$Ligand), c("3SLN", "6SLN"))) fail("Both receptor ligands are required.")
  if (anyDuplicated(data.frame(Ligand = sub("_.*$", "", cols), HAU = as.numeric(sub("^[36]SLN_", "", cols))))) fail("Duplicate receptor dose columns.")
  for (ligand in c("3SLN", "6SLN")) {
    if (!all(cfg$receptor_doses %in% measured$HAU[measured$Ligand == ligand])) fail("Requested doses absent for ", ligand)
  }
  meta <- metadata_for(unique(dat$Type), cfg)
  long$Type <- factor(long$Type, levels = meta$Type)
  long$Clade <- meta$Clade[match(long$Type, meta$Type)]
  dir.create(out <- file.path(cfg$output_dir, "receptor"), recursive = TRUE, showWarnings = FALSE)
  write_csv(long, file.path(out, "receptor_long.csv"))
  # Historical 'AUC' was an unweighted SUM, not numerical integration.
  # Keep that calculation explicitly named and report missing counts.
  sums <- long %>% dplyr::group_by(Type, Replicate_ID, Ligand, Clade) %>%
    dplyr::summarise(SummedAbsorbance = if (all(is.na(Absorbance))) NA_real_ else sum(Absorbance, na.rm = TRUE),
                     n_measured = sum(!is.na(Absorbance)), n_missing = sum(is.na(Absorbance)), .groups = "drop")
  write_csv(sums, file.path(out, "receptor_summed_absorbance_all_doses.csv"))
  selected <- long[long$HAU %in% cfg$receptor_doses, ]
  selected$Dose <- factor(selected$HAU, levels = cfg$receptor_doses)
  sm <- selected %>% dplyr::group_by(Type, Clade, Ligand, Dose) %>%
    dplyr::summarise(Mean = finite_mean(Absorbance), SD = stats::sd(Absorbance, na.rm = TRUE),
                     n = sum(!is.na(Absorbance)), .groups = "drop")
  write_csv(sm, file.path(out, "receptor_summary.csv"))
  palette <- stats::setNames(meta$Colour, meta$Type)
  dodge <- ggplot2::position_dodge(width = .8)
  bars <- ggplot2::ggplot(sm, ggplot2::aes(Dose, Mean, fill = Type)) +
    ggplot2::geom_col(position = dodge, width = cfg$bar_width) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SD, ymax = Mean + SD), position = dodge, width = .2) +
    ggplot2::geom_hline(yintercept = cfg$receptor_llod, linetype = "dotted") +
    ggplot2::facet_grid(Clade ~ Ligand) + ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(x = "HAU", y = "Absorbance (450 nm)", fill = "Virus") + ggplot2::theme_classic()
  save_plot(bars, file.path(out, "receptor_grouped_bars.pdf"), 14, max(6, 4 * length(unique(meta$Clade))))
  sm$MirrorMean <- ifelse(sm$Ligand == "3SLN", -sm$Mean, sm$Mean)
  sm$MirrorMin <- ifelse(sm$Ligand == "3SLN", -(sm$Mean + sm$SD), sm$Mean - sm$SD)
  sm$MirrorMax <- ifelse(sm$Ligand == "3SLN", -(sm$Mean - sm$SD), sm$Mean + sm$SD)
  mirror <- ggplot2::ggplot(sm, ggplot2::aes(Dose, MirrorMean, fill = Type)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_col(position = dodge, width = cfg$bar_width) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = MirrorMin, ymax = MirrorMax), position = dodge, width = .2) +
    ggplot2::facet_wrap(~ Clade, ncol = 1) + ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = abs) +
    ggplot2::labs(x = "HAU", y = "Absorbance (450 nm); 6SLN up / 3SLN down", fill = "Virus") + ggplot2::theme_minimal()
  box <- ggplot2::ggplot(sums, ggplot2::aes(Type, SummedAbsorbance, fill = Type)) +
    ggplot2::geom_boxplot(width = .6) +
    ggplot2::geom_point(ggplot2::aes(shape = Ligand), position = ggplot2::position_jitter(width = .1, height = 0, seed = cfg$seed), alpha = .7, size = 1.8) +
    ggplot2::scale_fill_manual(values = palette) + ggplot2::scale_shape_manual(values = c(`3SLN` = 17, `6SLN` = 16)) +
    ggplot2::facet_wrap(~ Clade, ncol = 1, scales = "free_y") + ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Summed absorbance across all measured doses", fill = "Virus") + ggplot2::theme_minimal()
  save_plot(mirror, file.path(out, "receptor_mirror.pdf"), 12, 8)
  save_plot(box, file.path(out, "receptor_sum_boxplots.pdf"), 12, max(8, nrow(meta) * .3))
  combined <- patchwork::wrap_plots(mirror, box, ncol = 2, guides = "collect")
  save_plot(combined, file.path(out, "receptor_combined.pdf"), 18, max(8, nrow(meta) * .3))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2 || !args[2] %in% c("ph", "temperature", "receptor", "all")) {
    fail("Usage: Rscript scripts/assay_analysis.R config/analysis_config.R ph|temperature|receptor|all")
  }
  e <- new.env(parent = baseenv()); sys.source(args[1], envir = e)
  if (!exists("config", e, inherits = FALSE) || !is.list(e$config)) fail("Configuration must define a list called config.")
  cfg <- e$config
  need_packages(c("openxlsx", "dplyr", "tidyr", "ggplot2"))
  suppressPackageStartupMessages(library(dplyr))
  set.seed(cfg$seed)
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(args[1], file.path(cfg$output_dir, "analysis_config_used.R"), overwrite = TRUE)
  on.exit(capture.output(sessionInfo(), file = file.path(cfg$output_dir, "sessionInfo.txt")), add = TRUE)
  modes <- if (args[2] == "all") c("ph", "temperature", "receptor") else args[2]
  for (mode in modes) {
    message("Running ", mode)
    if (mode == "receptor") run_receptor(cfg) else run_stability(mode, cfg)
  }
}
if (sys.nframe() == 0L) main()
