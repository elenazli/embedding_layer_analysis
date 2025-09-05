#!/usr/bin/env Rscript
# =============================================================================
# Embedding Layer Variant Analysis: Layer 14 vs Layer 28 Comparison
# =============================================================================
# Author: Elena Li (@elenazli)
# Description: Compares neural network embeddings across layers 14 and 28
#              to understand how genetic variants affect protein representations
# 
# Usage: Rscript 14v28cos_sim.R <SUBJECT_NAME>
# Example: Rscript 14v28cos_sim.R BRCA1
# =============================================================================

# Load required libraries
suppressMessages({
  library(reticulate)
  library(stringr)
  library(ggplot2)
})

# Import numpy for array operations
np <- import("numpy")

# =============================================================================
# Utility Functions
# =============================================================================

#' Read codon table and create mapping
#' @param codon_table_file Path to tab-delimited codon table
#' @return Named vector mapping codons to amino acids
read_codon_table <- function(codon_table_file) {
  # Read the tab-delimited file into a data frame
  codon_df <- read.delim(codon_table_file, stringsAsFactors = FALSE)
  
  # Create a named vector: names are codons, values are amino acids
  codon_map <- setNames(codon_df$aa, codon_df$codon)
  
  return(codon_map)
}

#' Calculate cosine similarity between two vectors
#' @param vec1 First vector
#' @param vec2 Second vector
#' @return Cosine similarity value
cosine_similarity <- function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

# =============================================================================
# Main Analysis Function
# =============================================================================

main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Rscript 14v28cos_sim.R <SUBJECT_NAME>")
  }
  
  subject <- args[1]
  cat(sprintf("Starting analysis for subject: %s\n", subject))
  
  # Set up file paths and variables
  lwr <- tolower(subject)
  lwr_14 <- sprintf("%s_14", lwr)
  lwr_28 <- sprintf("%s_28", lwr)
  
  # Load codon mapping
  codons <- read_codon_table("data/input/codon_table")
  
  # Get list of variant files
  files_dir <- list.files(
    path = sprintf("jobs/variant-embeddings-%s/output", lwr_14), 
    pattern = "\\.npy$"
  )
  
  if (length(files_dir) == 0) {
    stop(sprintf("No variant files found for subject %s", subject))
  }
  
  # Load subject metadata
  subject_table <- read.csv("data/example/subject_summary_table_test.csv", header = TRUE)
  index <- which(subject_table$aid == subject)
  
  if (length(index) == 0) {
    stop(sprintf("Subject %s not found in metadata table", subject))
  }
  
  gene <- subject_table$gene[index]
  position <- as.numeric(subject_table$position[index]) 
  start <- as.numeric(subject_table$start[index])
  end <- as.numeric(subject_table$end[index])
  source_codon <- subject_table$src_codon[index]
  dist_from_gene_end <- (end - start) - position
  
  cat(sprintf("Gene: %s, Position: %d, Source codon: %s\n", gene, position, source_codon))
  
  # Load source embeddings for both layers
  cat("Loading source embeddings...\n")
  df_14_src <- np$load(sprintf(
    "jobs/subject-layer-14/output/input_%s_%s_source_embeddings_blocks_14_mlp_l3.npy", 
    subject, gene
  ))
  df_28_src <- np$load(sprintf(
    "jobs/subject-layer-28/output/input_%s_%s_source_embeddings_blocks_28_mlp_l3.npy", 
    subject, gene
  ))
  
  # Convert to R matrices
  mat_14_src <- py_to_r(df_14_src)[1, , ]
  mat_28_src <- py_to_r(df_28_src)[1, , ]
  
  # Create output directory
  output_dir <- sprintf("figures/%s/14v28", subject)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", output_dir))
  }
  
  # Initialize summary dataframe
  variant_summary <- data.frame(
    Variant = character(),
    Codon = character(),
    Amino_Acid = character(),
    Min_Diff = numeric(),
    Max_Diff = numeric(),
    Range = numeric(),
    Mean_Diff = numeric(),
    Median_Diff = numeric(),
    Std_Diff = numeric(),
    Abs_Mean_Diff = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each variant
  cat(sprintf("Processing %d variants...\n", length(files_dir)))
  
  for (i in seq_along(files_dir)) {
    file <- files_dir[i]
    file_28 <- sub("14", "28", file)
    codon <- substr(as.character(file), 11, 13)
    label <- sprintf("%s_%s", codons[[codon]], codon)
    
    cat(sprintf("Processing variant %d/%d: %s\n", i, length(files_dir), label))
    
    # Load variant embeddings
    df_14_var <- np$load(sprintf("jobs/variant-embeddings-%s/output/%s", lwr_14, file))
    mat_14_var <- py_to_r(df_14_var)[1, , ]
    df_28_var <- np$load(sprintf("jobs/variant-embeddings-%s/output/%s", lwr_28, file_28))
    mat_28_var <- py_to_r(df_28_var)[1, , ]
    
    # Compute cosine similarities for both layers
    cosine_similarities_14 <- sapply(1:nrow(mat_14_src), function(i) {
      cosine_similarity(mat_14_src[i, ], mat_14_var[i, ])
    })
    
    cosine_similarities_28 <- sapply(1:nrow(mat_28_src), function(i) {
      cosine_similarity(mat_28_src[i, ], mat_28_var[i, ])
    })
    
    # Calculate layer difference
    cosine_similarities_diff <- cosine_similarities_28 - cosine_similarities_14
    
    # Calculate summary statistics
    min_diff <- min(cosine_similarities_diff)
    max_diff <- max(cosine_similarities_diff)
    range_diff <- max_diff - min_diff
    mean_diff <- mean(cosine_similarities_diff)
    median_diff <- median(cosine_similarities_diff)
    std_diff <- sd(cosine_similarities_diff)
    abs_mean_diff <- mean(abs(cosine_similarities_diff))
    
    # Add to summary dataframe
    variant_summary <- rbind(variant_summary, data.frame(
      Variant = label,
      Codon = codon,
      Amino_Acid = codons[[codon]],
      Min_Diff = min_diff,
      Max_Diff = max_diff,
      Range = range_diff,
      Mean_Diff = mean_diff,
      Median_Diff = median_diff,
      Std_Diff = std_diff,
      Abs_Mean_Diff = abs_mean_diff,
      stringsAsFactors = FALSE
    ))
    
    # Create individual variant plot
    cosine_sim <- data.frame(
      Position = 1:nrow(mat_14_src),
      Diff_Cos_Sim = cosine_similarities_diff
    )
    
    normal_pos <- ggplot(cosine_sim, aes(x = Position, y = Diff_Cos_Sim)) +
      geom_point(size = 0.5, alpha = 0.7) +  
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 10)
      ) +
      ggtitle(sprintf("%s: %s (Layer 28 - Layer 14)", subject, label)) +
      ylab("Cosine Similarity Difference") +
      xlab("Sequence Position")
    
    # Save individual plot
    ofn <- sprintf("figures/%s/14v28/%s.pdf", subject, label)
    ggsave(ofn, plot = normal_pos, width = 8, height = 6)
  }
  
  # Sort summary by absolute mean difference (most impactful variants first)
  variant_summary <- variant_summary[order(-variant_summary$Abs_Mean_Diff), ]
  
  # Save summary to CSV
  summary_file <- sprintf("figures/%s/14v28/variant_summary_%s.csv", subject, subject)
  write.csv(variant_summary, summary_file, row.names = FALSE)
  cat(sprintf("Saved summary statistics to: %s\n", summary_file))
  
  # Print top 10 most impactful variants
  cat("\n=== Top 10 Most Impactful Variants ===\n")
  print(head(variant_summary, 10))
  
  # Create summary plot
  range_plot <- ggplot(variant_summary, aes(x = reorder(Variant, Abs_Mean_Diff), y = Abs_Mean_Diff)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    ) +
    ggtitle(sprintf("%s: Variant Impact Ranking", subject)) +
    ylab("Absolute Mean Difference (Layer 28 - Layer 14)") +
    xlab("Variant")
  
  range_plot_file <- sprintf("figures/%s/14v28/variant_summary_plot_%s.pdf", subject, subject)
  ggsave(range_plot_file, plot = range_plot, width = 12, height = max(8, nrow(variant_summary) * 0.3))
  cat(sprintf("Saved summary plot to: %s\n", range_plot_file))
  
  cat(sprintf("\nAnalysis complete for subject %s!\n", subject))
  cat(sprintf("Results saved in: figures/%s/14v28/\n", subject))
}

# =============================================================================
# Execute main function
# =============================================================================
if (!interactive()) {
  main()
}
