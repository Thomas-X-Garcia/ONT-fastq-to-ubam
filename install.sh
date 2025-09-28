#!/bin/bash
# Installation script for ONT-fastq-to-ubam

set -e

echo "ONT-fastq-to-ubam Installation"
echo "=============================="
echo

# Check Python version
echo "Checking Python version..."
if ! python3 --version >/dev/null 2>&1; then
    echo "Error: Python 3 is not installed. Please install Python 3.8 or higher."
    exit 1
fi

python_version=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
echo "Found Python $python_version"

# Install pysam
echo
echo "Installing pysam..."
pip3 install --user pysam || pip install --user pysam || {
    echo "Error: Failed to install pysam. Try:"
    echo "  sudo pip3 install pysam"
    echo "  or"
    echo "  conda install -c bioconda pysam"
    exit 1
}

# Make scripts executable
echo
echo "Setting permissions..."
chmod +x ONT-fastq-to-ubam.py
chmod +x ONT-fastq-to-ubam_shard_wrapper.sh

# Check samtools for parallel wrapper
echo
echo "Checking for samtools (required for parallel processing)..."
if command -v samtools >/dev/null 2>&1; then
    samtools_version=$(samtools --version | head -1)
    echo "Found $samtools_version"
else
    echo "Warning: samtools not found. The parallel wrapper requires samtools."
    echo "To install samtools:"
    echo "  Ubuntu/Debian: sudo apt-get install samtools"
    echo "  CentOS/RHEL:   sudo yum install samtools"
    echo "  macOS:         brew install samtools"
fi

# Check for optional tools
echo
echo "Checking optional dependencies..."
if command -v pigz >/dev/null 2>&1; then
    echo "✓ pigz found (faster .gz decompression)"
else
    echo "ℹ pigz not found (optional, for faster .gz handling)"
fi

echo
echo "Installation complete!"
echo
echo "Usage:"
echo "  ./ONT-fastq-to-ubam.py input.fastq output.ubam --sample-name SAMPLE"
echo
echo "For parallel processing:"
echo "  ./ONT-fastq-to-ubam_shard_wrapper.sh -i input.fq.gz -o output.ubam -s SAMPLE -n 8 -p ./ONT-fastq-to-ubam.py"
echo
echo "For help:"
echo "  ./ONT-fastq-to-ubam.py --help"
echo "  ./ONT-fastq-to-ubam_shard_wrapper.sh -h"