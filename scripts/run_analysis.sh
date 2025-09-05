#!/bin/bash
# =============================================================================
# Embedding Layer Variant Analysis - Execution Script
# =============================================================================
# Author: Elena Li (@elenazli)
# Description: Shell script to run the embedding layer variant analysis
# 
# Usage: ./run_analysis.sh <SUBJECT_NAME>
# Example: ./run_analysis.sh BAA
# =============================================================================

set -e  # Exit on any error

# Check if subject name is provided
if [ $# -eq 0 ]; then
    echo "Error: No subject name provided"
    echo "Usage: $0 <SUBJECT_NAME>"
    echo "Example: $0 BAA"
    exit 1
fi

SUBJECT=$1
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "============================================================================="
echo "Embedding Layer Variant Analysis"
echo "Author: Elena Li (@elenazli)"
echo "============================================================================="
echo "Subject: $SUBJECT"
echo "Project Root: $PROJECT_ROOT"
echo ""

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    echo "Error: R is not installed or not in PATH"
    echo "Please install R and ensure Rscript is available"
    exit 1
fi

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH"
    echo "Please install Python 3"
    exit 1
fi

# Check if required R packages are installed
echo "Checking R dependencies..."
Rscript -e "
required_packages <- c('reticulate', 'stringr', 'ggplot2')
missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if(length(missing_packages) > 0) {
    cat('Missing R packages:', paste(missing_packages, collapse=', '), '\n')
    cat('Installing missing packages...\n')
    install.packages(missing_packages, repos='https://cran.rstudio.com/')
} else {
    cat('All required R packages are installed.\n')
}
"

# Check if required Python packages are installed
echo "Checking Python dependencies..."
python3 -c "
import sys
required_packages = ['numpy']
missing_packages = []
for package in required_packages:
    try:
        __import__(package)
    except ImportError:
        missing_packages.append(package)

if missing_packages:
    print(f'Missing Python packages: {missing_packages}')
    print('Please install missing packages using: pip install ' + ' '.join(missing_packages))
    sys.exit(1)
else:
    print('All required Python packages are installed.')
"

# Check if data files exist
echo "Checking data files..."
if [ ! -f "$PROJECT_ROOT/data/input/codon_table" ]; then
    echo "Error: Codon table not found at $PROJECT_ROOT/data/input/codon_table"
    exit 1
fi

if [ ! -f "$PROJECT_ROOT/data/example/subject_summary_table_test.csv" ]; then
    echo "Error: Subject summary table not found at $PROJECT_ROOT/data/example/subject_summary_table_test.csv"
    exit 1
fi

# Check if embedding files exist (approximate check)
if [ ! -d "jobs" ]; then
    echo "Warning: 'jobs' directory not found in current location"
    echo "Make sure you're running this script from the directory containing your embedding data"
    echo "Expected structure: jobs/variant-embeddings-{subject}_14/output/"
fi

echo ""
echo "Starting analysis..."
echo "============================================================================="

# Change to project root directory
cd "$PROJECT_ROOT"

# Run the R script
Rscript src/14v28cos_sim.R "$SUBJECT"

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "============================================================================="
    echo "Analysis completed successfully!"
    echo "Results saved in: figures/$SUBJECT/14v28/"
    echo ""
    echo "Generated files:"
    echo "  - Individual variant plots: figures/$SUBJECT/14v28/*.pdf"
    echo "  - Summary statistics: figures/$SUBJECT/14v28/variant_summary_$SUBJECT.csv"
    echo "  - Summary plot: figures/$SUBJECT/14v28/variant_summary_plot_$SUBJECT.pdf"
    echo "============================================================================="
else
    echo ""
    echo "============================================================================="
    echo "Analysis failed! Please check the error messages above."
    echo "============================================================================="
    exit 1
fi
