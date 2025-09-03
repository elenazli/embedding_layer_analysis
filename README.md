# Embedding Layer Variant Analysis

**Author:** Elena Li ([@elenazli](https://github.com/elenazli))

A comprehensive analysis tool for comparing neural network embeddings across different layers to understand how genetic variants affect protein representations in the Evo2 deep learning model.

## Overview

This project analyzes how genetic variants (single nucleotide polymorphisms) affect neural network embeddings at different layers (14 vs 28) of a protein language model. By comparing cosine similarities between source and variant embeddings across layers, we can understand how the model's "understanding" of genetic variation evolves through the network.

## Methodology

### Key Concepts
- **Layer Comparison**: Analyzes embeddings from layers 14 and 28 of a neural network
- **Cosine Similarity**: Measures how similar variant embeddings are to source embeddings
- **Position-Specific Analysis**: Examines how variants affect different sequence positions
- **Statistical Summarization**: Provides comprehensive metrics for variant impact assessment

### Analysis Pipeline
1. **Data Loading**: Loads source and variant embeddings from numpy arrays
2. **Similarity Calculation**: Computes cosine similarities for each position
3. **Layer Comparison**: Calculates differences between layer 14 and 28 similarities
4. **Statistical Analysis**: Generates summary statistics for each variant
5. **Visualization**: Creates individual plots and summary visualizations
6. **Output Generation**: Saves results as CSV and PDF files

## Key Metrics

- **Mean Difference**: Average difference in cosine similarity between layers
- **Range**: Spread of differences across sequence positions
- **Standard Deviation**: Consistency of effects across positions
- **Absolute Mean Difference**: Overall magnitude of layer-specific effects

## Start

### Prerequisites
- R (>= 4.0)
- Python (>= 3.7) with numpy
- Required R packages: `reticulate`, `stringr`, `ggplot2`

### Installation
```bash
# Clone the repository
git clone https://github.com/elenazli/embedding_layer_variant_analysis.git
cd embedding_layer_variant_analysis

# Install R dependencies
Rscript -e "install.packages(c('reticulate', 'stringr', 'ggplot2'))"

# Install Python dependencies
pip install -r requirements.txt
```

### Usage
```bash
# Run analysis for a specific subject
Rscript src/14v28cos_sim.R <SUBJECT_NAME>

# Example
Rscript src/14v28cos_sim.R BAA
```

### Expected Outputs
- **Individual plots**: `figures/{SUBJECT}/14v28/{VARIANT}.pdf`
- **Summary plot**: `figures/{SUBJECT}/14v28/variant_summary_plot_{SUBJECT}.pdf`
- **CSV results**: `figures/{SUBJECT}/14v28/variant_summary_{SUBJECT}.csv`
- **Console output**: Top 10 most impactful variants

## Project Structure

```
embedding_layer_variant_analysis/
├── README.md                 # This file
├── requirements.txt          # Python dependencies
├── environment.yml          # Conda environment
├── src/
│   ├── 14v28cos_sim.R      # Main analysis script
│   └── utils/               # Utility functions
├── data/
│   ├── input/               # Input data files
│   └── example/             # Example data
├── scripts/
│   └── run_analysis.sh     # Execution script
├── docs/
│   └── methodology.md      # Detailed methodology
└── examples/
    └── sample_output/      # Example results
```

## Data Requirements

### Input Files
- **Embedding files**: `.npy` files containing neural network embeddings
- **Subject table**: CSV with subject metadata
- **Codon table**: Tab-delimited file mapping codons to amino acids

### File Structure
```
jobs/
├── subject-layer-14/output/
│   └── input_{SUBJECT}_{GENE}_source_embeddings_blocks_14_mlp_l3.npy
├── subject-layer-28/output/
│   └── input_{SUBJECT}_{GENE}_source_embeddings_blocks_28_mlp_l3.npy
└── variant-embeddings-{SUBJECT}_14/output/
    └── {VARIANT_FILES}.npy
```

## Interpretation Guide

### High Impact Variants
- **High Absolute Mean Difference**: Variants with large overall effects
- **High Range**: Variants with position-dependent effects
- **High Standard Deviation**: Variants with inconsistent effects

### Biological Significance
- **Positive differences**: Variants become more similar to source at later layers
- **Negative differences**: Variants become less similar to source at later layers
- **Position-specific patterns**: May indicate functional domains or structural regions

## Example Results

The analysis generates several types of outputs:

1. **Individual Variant Plots**: Show position-specific differences for each variant
2. **Summary Statistics**: CSV file with comprehensive metrics for all variants
3. **Ranking Plot**: Bar chart showing variants ranked by impact
4. **Console Summary**: Top 10 most impactful variants printed to terminal

---

*This tool helps researchers understand how deep learning models interpret genetic variation and can provide insights into protein function and disease mechanisms.*
