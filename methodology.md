# Methodology: Embedding Layer Variant Analysis

## Overview

This document describes the detailed methodology for analyzing how genetic variants affect neural network embeddings across different layers of a protein language model.

## Background

### Neural Network Embeddings
- **Layer 14**: Earlier layer capturing basic sequence patterns and local features
- **Layer 28**: Later layer capturing higher-order abstractions and complex relationships
- **Embedding Evolution**: How the model's "understanding" changes as information flows through the network

### Genetic Variants
- **Single Nucleotide Polymorphisms (SNPs)**: Single base pair changes in DNA
- **Codon Variants**: All 64 possible codons at a specific position
- **Amino Acid Changes**: Resulting protein sequence modifications

## Analysis Pipeline

### 1. Data Loading
```r
# Load source embeddings for both layers
df_14_src <- np$load("jobs/subject-layer-14/output/input_{SUBJECT}_{GENE}_source_embeddings_blocks_14_mlp_l3.npy")
df_28_src <- np$load("jobs/subject-layer-28/output/input_{SUBJECT}_{GENE}_source_embeddings_blocks_28_mlp_l3.npy")

# Load variant embeddings
df_14_var <- np$load("jobs/variant-embeddings-{SUBJECT}_14/output/{VARIANT_FILE}.npy")
df_28_var <- np$load("jobs/variant-embeddings-{SUBJECT}_28/output/{VARIANT_FILE}.npy")
```

### 2. Cosine Similarity Calculation
For each position `i` in the sequence:

```r
cosine_similarity_14[i] = (src_14[i] · var_14[i]) / (||src_14[i]|| × ||var_14[i]||)
cosine_similarity_28[i] = (src_28[i] · var_28[i]) / (||src_28[i]|| × ||var_28[i]||)
```

Where:
- `src_14[i]` and `src_28[i]` are source embeddings at position `i` for layers 14 and 28
- `var_14[i]` and `var_28[i]` are variant embeddings at position `i` for layers 14 and 28
- `·` denotes dot product
- `||·||` denotes vector magnitude

### 3. Layer Difference Calculation
```r
layer_difference[i] = cosine_similarity_28[i] - cosine_similarity_14[i]
```

### 4. Statistical Analysis
For each variant, calculate:

- **Min_Diff**: Minimum layer difference across all positions
- **Max_Diff**: Maximum layer difference across all positions
- **Range**: Max_Diff - Min_Diff (spread of differences)
- **Mean_Diff**: Average layer difference
- **Median_Diff**: Median layer difference
- **Std_Diff**: Standard deviation of layer differences
- **Abs_Mean_Diff**: Mean of absolute layer differences

## Interpretation

### Positive Differences (Layer 28 > Layer 14)
- Variant becomes **more similar** to source at later layers
- Suggests the model learns to "normalize" or "understand" the variant better
- Could indicate functionally similar variants despite sequence differences

### Negative Differences (Layer 28 < Layer 14)
- Variant becomes **less similar** to source at later layers
- Suggests the model learns to distinguish the variant as functionally different
- Could indicate pathogenic or functionally disruptive variants

### High Range Values
- **Position-dependent effects**: Variant affects different sequence positions very differently
- **Context sensitivity**: Variant's effect depends on surrounding sequence context
- **Functional significance**: May indicate variants affecting specific protein domains

### High Absolute Mean Difference
- **Overall impact**: Large magnitude of layer-specific effects
- **Functional importance**: Variants with significant biological consequences
- **Model sensitivity**: Positions where the model shows high sensitivity to changes

## Biological Context

### Protein Structure
- **Active sites**: Positions critical for protein function
- **Binding domains**: Regions involved in protein-protein interactions
- **Structural regions**: Areas important for protein folding and stability

### Variant Effects
- **Conservative changes**: Amino acid substitutions with similar properties
- **Radical changes**: Amino acid substitutions with very different properties
- **Silent mutations**: Changes that don't alter the amino acid sequence

### Functional Implications
- **Pathogenic variants**: Changes associated with disease
- **Benign variants**: Changes with no functional consequences
- **Regulatory variants**: Changes affecting gene expression or regulation

## Statistical Considerations

### Normalization
- Cosine similarity is naturally normalized (range: -1 to 1)
- Layer differences provide relative comparisons
- Absolute values focus on magnitude regardless of direction

### Multiple Testing
- 64 variants per subject (all possible codons)
- Multiple positions per sequence
- Consider correction for multiple comparisons in downstream analysis

### Reproducibility
- Fixed random seeds for consistent results
- Standardized input preprocessing
- Version-controlled analysis pipeline

## Limitations

### Model Dependencies
- Results specific to the trained protein language model
- May not generalize to other architectures or training data
- Limited by model's training distribution

### Biological Validation
- Computational predictions require experimental validation
- Context-dependent effects may be model artifacts
- Functional significance needs biological interpretation

### Data Requirements
- Requires high-quality embedding data
- Sensitive to input preprocessing
- Computational resources for large-scale analysis

## Future Directions

### Enhanced Analysis
- Multi-layer comparisons (not just 14 vs 28)
- Attention weight integration
- Structural context incorporation

### Biological Integration
- Functional annotation databases
- Protein structure mapping
- Disease variant databases

### Methodological Improvements
- Statistical significance testing
- Effect size quantification
- Uncertainty estimation
