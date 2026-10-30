#!/usr/bin/env Rscript
# Run explicitly; analysis scripts never install packages automatically.
cran <- c("openxlsx", "dplyr", "tidyr", "ggplot2", "patchwork", "circlize",
          "viridisLite", "lme4", "lmerTest", "emmeans", "ape", "BiocManager")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
bio <- c("ComplexHeatmap", "ggtree")
missing <- bio[!vapply(bio, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE)
