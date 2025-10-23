#!/usr/bin/env python3
"""
Fast BAM to SNV Matrix Extractor
Optimized for speed with minimal arguments

Usage: python fast_bam_to_snv.py input.bam output_prefix
"""

import pysam
import numpy as np
import sys
import time
from collections import Counter

def extract_snv_matrix_fast(bam_file, output_prefix):
    """
    Fast extraction of SNV matrix from BAM file
    Uses optimized approach with minimal processing
    """
    print(f"Processing BAM file: {bam_file}")
    start_time = time.time()
    
    # Get basic info from BAM header
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        chromosome = bam.references[0]  # Use first chromosome
        chr_length = bam.lengths[0]
        print(f"Processing chromosome: {chromosome} ({chr_length:,} bp)")
    
    # Step 1: Collect all read sequences quickly
    print("Step 1: Reading alignments...")
    read_data = []
    total_reads = 0
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(chromosome):
            total_reads += 1
            
            # Skip unmapped, low quality, secondary alignments
            if (read.is_unmapped or read.is_secondary or 
                read.is_supplementary or read.mapping_quality < 20):
                continue
            
            # Only process simple alignments (no complex CIGAR)
            if len(read.cigartuples) != 1 or read.cigartuples[0][0] != 0:
                continue
            
            # Extract sequence as integers: A=1, C=2, G=3, T=4
            sequence = []
            for nucleotide in read.query_sequence:
                if nucleotide == 'A': sequence.append(1)
                elif nucleotide == 'C': sequence.append(2)
                elif nucleotide == 'G': sequence.append(3)
                elif nucleotide == 'T': sequence.append(4)
                else: sequence.append(0)
            
            read_data.append({
                'seq': sequence,
                'start': read.reference_start,
                'end': read.reference_end
            })
    
    print(f"Collected {len(read_data)} valid reads from {total_reads} total")
    
    if len(read_data) == 0:
        print("No valid reads found!")
        return
    
    # Step 2: Find coverage region automatically
    print("Step 2: Determining coverage region...")
    min_start = min(r['start'] for r in read_data)
    max_end = max(r['end'] for r in read_data)
    region_length = max_end - min_start
    
    print(f"Coverage region: {min_start:,} - {max_end:,} ({region_length:,} bp)")
    
    # Step 3: Build position matrix efficiently
    print("Step 3: Building position matrix...")
    position_matrix = {}
    
    for read in read_data:
        for i, nucleotide in enumerate(read['seq']):
            if nucleotide == 0:  # Skip gaps/unknown
                continue
            
            pos = read['start'] + i
            if pos < min_start or pos >= max_end:
                continue
                
            if pos not in position_matrix:
                position_matrix[pos] = []
            position_matrix[pos].append(nucleotide)
    
    # Step 4: Find SNV positions quickly
    print("Step 4: Calling SNVs...")
    snv_positions = []
    
    for pos in sorted(position_matrix.keys()):
        nucleotides = position_matrix[pos]
        if len(nucleotides) < 3:  # Need minimum coverage
            continue
            
        counts = Counter(nucleotides)
        total = len(nucleotides)
        
        # Check if position has multiple alleles above 10% frequency
        significant_alleles = sum(1 for count in counts.values() if count/total >= 0.1)
        
        if significant_alleles > 1:
            snv_positions.append(pos)
    
    print(f"Found {len(snv_positions)} SNV positions")
    
    if len(snv_positions) == 0:
        print("No SNVs detected!")
        return
    
    # Step 5: Create SNV matrix
    print("Step 5: Building SNV matrix...")
    snv_matrix = []
    
    for read in read_data:
        snv_row = []
        
        for snv_pos in snv_positions:
            # Check if read covers this SNV position
            if read['start'] <= snv_pos < read['end']:
                offset = snv_pos - read['start']
                if offset < len(read['seq']):
                    snv_row.append(read['seq'][offset])
                else:
                    snv_row.append(0)
            else:
                snv_row.append(0)  # No coverage
        
        # Only keep reads that cover at least one SNV
        if any(x != 0 for x in snv_row):
            snv_matrix.append(snv_row)
    
    snv_matrix = np.array(snv_matrix)
    
    # Step 6: Save outputs
    print("Step 6: Saving outputs...")
    
    # Save SNV matrix
    np.savetxt(f"{output_prefix}_SNV_matrix.txt", snv_matrix, fmt='%d')
    
    # Save SNV positions
    with open(f"{output_prefix}_SNV_pos.txt", 'w') as f:
        f.write(' '.join(str(pos) for pos in snv_positions))
    
    # Save basic info
    with open(f"{output_prefix}_info.txt", 'w') as f:
        f.write(f"Chromosome: {chromosome}\n")
        f.write(f"Total reads processed: {len(read_data)}\n")
        f.write(f"SNV positions: {len(snv_positions)}\n")
        f.write(f"Final matrix: {snv_matrix.shape[0]} reads × {snv_matrix.shape[1]} SNVs\n")
        f.write(f"Coverage region: {min_start}-{max_end}\n")
    
    end_time = time.time()
    print(f"Complete! Processing time: {end_time - start_time:.2f} seconds")
    print(f"Final matrix: {snv_matrix.shape[0]} reads × {snv_matrix.shape[1]} SNVs")
    print(f"Files saved: {output_prefix}_SNV_matrix.txt, {output_prefix}_SNV_pos.txt, {output_prefix}_info.txt")

def main():
    if len(sys.argv) != 3:
        print("Usage: python fast_bam_to_snv.py input.bam output_prefix")
        print("Example: python fast_bam_to_snv.py sample.bam results")
        return 1
    
    bam_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    if not os.path.exists(bam_file):
        print(f"Error: BAM file {bam_file} not found")
        return 1
    
    extract_snv_matrix_fast(bam_file, output_prefix)
    return 0

if __name__ == "__main__":
    import os
    exit(main())