#!/usr/bin/env Rscript

library("argparse")
library(network)
library(ggnetwork)

parser <- ArgumentParser(description = "Plot Dcifer result")
parser$add_argument(
  "--metadata",
  type = "character",
  required = TRUE,
  help = "Path to the metadata file (tsv, must have ID and Group columns)"
)
parser$add_argument(
  "--outdir",
  type = "character",
  required = TRUE,
  help = "Output directory to save the plots"
)
parser$add_argument(
    "--indir",
    type = "character",
    required = TRUE,
    help = "Input directory containing the between-infection relatedness estimate file (outdir of run_Dcifer)"
)

args <- parser$parse_args()
## Change FIXME lines to the appropriate parameters
if (!dir.exists(args$outdir)) {
  dir.create(args$outdir, recursive = TRUE)
}

# FIXME: change to path to metadata
metadata.file <- args$metadata

# FIXME: adjust based on file format
metadata <- read.delim(metadata.file)


# FIXME: change to path of between-infection relatedness estimate
mall.estimate.file <- file.path(args$indir, "mhap_between_relatedness_estimate.tsv")

mall.estimate <- read.delim(mall.estimate.file)

## STOP and make sure the metadata has the ID in the first column
## the ID must be the ones used in Dcifer analysis


# FIXME: change according sample and group column in metadata
metadata.sample.column <- "ID"
metadata.group.column <- "Group"

mall.estimate.meta <-
  merge(
    mall.estimate,
    metadata,
    by.x = "sample_id1",
    by.y = metadata.sample.column
  )
mall.estimate.meta <-
  merge(
    mall.estimate.meta,
    metadata,
    by.x = "sample_id2",
    by.y = metadata.sample.column
  )


mall.estimate.meta.within <-
  mall.estimate.meta[
    mall.estimate.meta[[paste0(metadata.group.column, ".x")]] ==
      mall.estimate.meta[[paste0(metadata.group.column, ".y")]],
    
  ]

mall.estimate.meta.within[[metadata.group.column]] <-
  mall.estimate.meta.within[[paste0(metadata.group.column, ".x")]]

# FIXME: adjust name and size (in inches) of between-infection relatedness plot
between.relatedness.plot.file <- file.path(args$outdir, "mhap_between_relatedness.pdf")

pdf(file = between.relatedness.plot.file, width = 4, height = 4)

ggplot(
  mall.estimate.meta.within,
  aes(x = as.factor(.data[[metadata.group.column]]), y = relatedness)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.35, height = 0, alpha = 0.5) +
  ylim(0, 1) +
  xlab(metadata.group.column) +
  ylab("between relatedness") +
  theme_classic()

dev.off()


z <- as.factor(metadata[[metadata.group.column]])

# FIXME: select the colour palette according to the number of groups
# The "Okabe-Ito" colour palette can accommodate up to 9 groups (default)
# The "Polychrome 36" colour palette can accommodate up to 36 groups
# You can also specify your own colour palette to use
nlevels(z)

# exclude black from Okabe-Ito palette
cols <- palette.colors(palette = "Okabe-Ito")[-1]

palette <- setNames(cols[1:nlevels(z)], levels(z))

# FIXME: determine path to save plots
prefix.filename <- args$outdir

# FIXME: adjust the IBD thresholds as needed
IBD.thresholds <- c(1, 1/2, 1/4, 1/8, 1/16) * 0.95

# make sure the IBD thresholds are in decreasing order
IBD.thresholds <- sort(IBD.thresholds, decreasing = TRUE)

g <-
  network(
    mall.estimate.meta[, c("sample_id1", "sample_id2", "relatedness")],
    directed = FALSE
  )

for (i in seq_along(IBD.thresholds)) {
  IBD.threshold <- IBD.thresholds[i]
  
  net <- network.copy(g)
  delete.edges(net, which(get.edge.value(net, "relatedness") < IBD.threshold))
  
  # make plots reproducible
  set.seed(1)
  net <- ggnetwork(net)
  
  net <-
    merge(net, metadata, by.x = "vertex.names", by.y = metadata.sample.column, sort = FALSE)
  
  # sort by relatedness to emphasize close relationships
  net <- net[order(net[["relatedness"]]), ]
  
  threshold.index <- i - 1

  # FIXME: adjust size (in inches) of network plot
  pdf( file.path(prefix.filename, paste0("mhap_network_", IBD.threshold, ".pdf")), width = 9, height = 8)

  # FIXME: relatedness network framework, adjust metadata as necessary
  print(
    ggplot(net, aes(x, y, xend = xend, yend = yend, fill = as.factor(.data[[metadata.group.column]]))) +
      geom_edges(aes(colour = relatedness, linewidth = relatedness)) +
      geom_nodes(size = 5, shape = 21) +
      theme_blank() +
      labs(title = paste0(IBD.threshold * 100, "%"), fill = metadata.group.column) +
      guides(linewidth = "none") +
      scale_fill_manual(breaks = levels(z), values = palette) +
      # using low and high colours of Greys palette from ColorBrewer
      scale_colour_binned(
        name = "IBD",
        breaks = rev(IBD.thresholds[1:threshold.index]),
        labels = rev(paste0(IBD.thresholds[1:threshold.index] * 100, "%")),
        limits = c(0, 1),
        low = "#F0F0F0",
        high = "#252525"
      ) +
      scale_linewidth_binned(
        name = "IBD",
        breaks = rev(IBD.thresholds[1:threshold.index]),
        labels = rev(paste0(IBD.thresholds[1:threshold.index] * 100, "%")),
        limits = c(0, 1),
        range = c(0.5, 4)
      )
  )
  
  dev.off()
}
