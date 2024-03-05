suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggplot2))

# Function to plot a Newick tree and save it as a JPEG image
plot_tree_from_newick <- function(newick_file, output_file) {
  # Read the Newick tree
  tree <- read.tree(newick_file)

  # Create a ggtree object from the tree
  tree_gg <- ggtree(tree,
                    layout='rectangular') +
                geom_tiplab(align = TRUE) +
                geom_tippoint(shape=21, fill="white")

  # Create the plot
  p <- tree_gg +
    geom_treescale() +
    coord_cartesian(clip = 'off')

  # Save the plot as a JPEG image
  ggsave(filename = output_file, plot = p, device = "jpeg")
}


plot_tree_from_newick(snakemake@input[[1]], snakemake@output[[1]])
