#!/bin/bash
# Example usage of ONT-fastq-to-ubam

echo "ONT-fastq-to-ubam Example Usage"
echo "================================"
echo

# Example 1: Basic conversion
echo "Example 1: Basic single-threaded conversion"
echo "-------------------------------------------"
echo "./ONT-fastq-to-ubam.py \\"
echo "  input_reads.fastq.gz \\"
echo "  output.ubam \\"
echo "  --sample-name SAMPLE001"
echo

# Example 2: Parallel processing for large files
echo "Example 2: High-speed parallel conversion (8 cores)"
echo "---------------------------------------------------"
echo "./ONT-fastq-to-ubam_shard_wrapper.sh \\"
echo "  -i large_dataset.fastq.gz \\"
echo "  -o large_dataset.ubam \\"
echo "  -s SAMPLE001 \\"
echo "  -n 8 \\"
echo "  -p ./ONT-fastq-to-ubam.py"
echo

# Example 3: Maximum speed with all available cores
echo "Example 3: Maximum speed using all CPU cores"
echo "--------------------------------------------"
echo "./ONT-fastq-to-ubam_shard_wrapper.sh \\"
echo "  -i huge_dataset.fastq.gz \\"
echo "  -o huge_dataset.ubam \\"
echo "  -s SAMPLE001 \\"
echo "  -n \$(nproc) \\"
echo "  -p ./ONT-fastq-to-ubam.py \\"
echo "  -- --no-individual-tags  # Skip individual tags for speed"
echo

# Example 4: Preparing for wf-human-variation
echo "Example 4: Preparing data for epi2me-labs/wf-human-variation"
echo "------------------------------------------------------------"
echo "# Step 1: Convert to uBAM"
echo "./ONT-fastq-to-ubam.py ont_reads.fq.gz prepared.ubam --sample-name HG002"
echo
echo "# Step 2: Run the workflow"
echo "nextflow run epi2me-labs/wf-human-variation \\"
echo "  --bam prepared.ubam \\"
echo "  --ref reference.fasta \\"
echo "  --snp --sv --mod --cnv --str"
echo

# Example 5: Verify output
echo "Example 5: Verify the output"
echo "----------------------------"
echo "# Check header and read group"
echo "samtools view -H output.ubam | grep '^@RG'"
echo
echo "# Check preserved metadata tags"
echo "samtools view output.ubam | head -1 | tr '\\t' '\\n' | grep '^X'"
echo
echo "# Count reads"
echo "samtools view -c output.ubam"