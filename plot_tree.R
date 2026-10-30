#!/usr/bin/env Rscript
# Optional generic replacement for the mixed-project exploratory tree script.
# Usage: Rscript scripts/plot_tree.R tree.nwk metadata.csv output.pdf [drop_tips.txt]
# metadata.csv: label (exact tree tip), Clade, Group. Group controls point shape.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3, 4)) stop("Usage: Rscript scripts/plot_tree.R tree.nwk metadata.csv output.pdf [drop_tips.txt]", call. = FALSE)
for (p in c("ape", "ggtree", "ggplot2")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Required package: ", p, call. = FALSE)
}
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggplot2))
tree <- ape::read.tree(args[1])
if (!inherits(tree, "phylo")) stop("Provide exactly one Newick tree.")
if (anyDuplicated(tree$tip.label)) stop("Tree contains duplicate tip labels.")
if (length(args) == 4) {
  drops <- readLines(args[4], warn = FALSE)
  drops <- drops[nzchar(drops) & !startsWith(drops, "#")]
  missing <- setdiff(drops, tree$tip.label)
  if (length(missing)) stop("Drop list contains unknown tips: ", paste(missing, collapse = ", "))
  if (length(drops)) tree <- ape::drop.tip(tree, drops)
}
metadata <- utils::read.csv(args[2], stringsAsFactors = FALSE)
if (!all(c("label", "Clade", "Group") %in% names(metadata))) stop("Metadata requires label, Clade, Group.")
if (anyDuplicated(metadata$label) || anyNA(match(tree$tip.label, metadata$label))) stop("Metadata tip labels are duplicated or incomplete.")
metadata <- metadata[match(tree$tip.label, metadata$label), ]
if (anyNA(metadata[c("label", "Clade", "Group")])) stop("Missing tree metadata values.")
p <- ggtree(tree) %<+% metadata + geom_tippoint(aes(colour = Clade, shape = Group), size = 2) +
  geom_tiplab(size = 2, offset = max(ape::node.depth.edgelength(tree)) * .01) +
  geom_treescale() + theme(legend.position = "right")
dir.create(dirname(args[3]), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(args[3], p, width = 12, height = max(6, length(tree$tip.label) * .15), limitsize = FALSE)
utils::write.csv(data.frame(label = tree$tip.label), paste0(tools::file_path_sans_ext(args[3]), "_tips.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = paste0(tools::file_path_sans_ext(args[3]), "_sessionInfo.txt"))
