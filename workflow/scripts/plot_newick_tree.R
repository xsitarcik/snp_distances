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

# Check if the correct number of command-line arguments are provided

if (length(commandArgs(trailingOnly = TRUE)) != 2) {
  cat("Usage: Rscript tree_plott.R input_newick_file output_jpeg_file\n")
  quit(status = 1)
}

# Get input and output file names from command-line arguments
input_newick_file <- commandArgs(trailingOnly = TRUE)[1]
output_jpeg_file <- commandArgs(trailingOnly = TRUE)[2]

# Call the function to plot the tree and save it as a JPEG image
plot_tree_from_newick(input_newick_file, output_jpeg_file)
