#!/bin/bash

# Define the directory containing the VCF files
vcf_dir="/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/HPRC/trgt-v0.8.0-9bd9f00/"

# Define the output file
output_file="/uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/Git/STRchive_manuscript/data/combined_output.txt"

# Remove the output file if it exists
rm -f "$output_file"

# Loop through each VCF file in the directory
for vcf_file in "$vcf_dir"/*.vcf.gz; do
    # Get the filename without the directory path
    filename=$(basename "$vcf_file")

    # Run the Python script for the current VCF file and append the output to the output file
    python scripts/tr-solve2TRGT.py "$vcf_file" >> "$output_file" 2>&1

done
