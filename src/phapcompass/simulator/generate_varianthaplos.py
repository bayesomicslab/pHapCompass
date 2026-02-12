#!/usr/bin/env python3
"""
Generate _varianthaplos.txt files by comparing FASTA files in directories.

This script processes directories containing multiple FASTA files (haplotypes)
and generates varianthaplos.txt files in the same format as haplogenerator3.

Usage:
    python generate_varianthaplos.py <directory> [<directory2> ...]
    
Each directory should contain FASTA files representing different haplotypes.
The script will compare all FASTA files within each directory and output
a _varianthaplos.txt file listing all variant positions.
"""

import os
import sys
import gzip
import argparse
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO


def read_fasta_file(filepath):
    """
    Read a FASTA file (supports .gz files) and return sequences as dict.
    
    Args:
        filepath: Path to FASTA file
        
    Returns:
        dict: {contig_id: sequence_string}
    """
    sequences = {}
    
    if str(filepath).endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    try:
        with opener(filepath, mode) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences[record.id] = str(record.seq)
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return {}
    
    return sequences


def find_variants(haplotype_sequences, contig_id):
    """
    Find variant positions by comparing sequences across haplotypes.
    
    Args:
        haplotype_sequences: list of sequence strings for the same contig
        contig_id: name of the contig
        
    Returns:
        list: [(position, ref_allele, alt_alleles, haplotype_calls)]
    """
    if not haplotype_sequences or len(set(len(seq) for seq in haplotype_sequences)) > 1:
        print(f"Warning: Sequences for {contig_id} have different lengths or are empty", file=sys.stderr)
        return []
    
    variants = []
    seq_length = len(haplotype_sequences[0])
    num_haplotypes = len(haplotype_sequences)
    
    for pos in range(seq_length):
        # Get all alleles at this position
        alleles_at_pos = [seq[pos].upper() for seq in haplotype_sequences]
        unique_alleles = list(set(alleles_at_pos))
        
        # Skip if all alleles are the same (no variant)
        if len(unique_alleles) <= 1:
            continue
            
        # Skip if any allele is 'N' or other ambiguous character
        if any(allele not in 'ATCG' for allele in unique_alleles):
            continue
        
        # Determine reference and alternate alleles
        # Use the most common allele as reference
        allele_counts = {allele: alleles_at_pos.count(allele) for allele in unique_alleles}
        ref_allele = max(allele_counts.keys(), key=lambda x: allele_counts[x])
        alt_alleles = [allele for allele in unique_alleles if allele != ref_allele]
        
        # For now, handle only bi-allelic sites (most common case)
        if len(alt_alleles) == 1:
            alt_allele = alt_alleles[0]
            
            # Create haplotype calls (1 = alt allele, 0 = ref allele)
            hap_calls = [1 if allele == alt_allele else 0 for allele in alleles_at_pos]
            
            # 1-based position for output
            variants.append((pos + 1, ref_allele, alt_allele, hap_calls))
    
    return variants


def process_directory(directory_path, output_filename=None):
    """
    Process a directory containing FASTA files and generate varianthaplos.txt.
    
    Args:
        directory_path: Path to directory containing FASTA files
        output_filename: Optional custom output filename
    """
    directory_path = Path(directory_path)
    
    if not directory_path.exists() or not directory_path.is_dir():
        print(f"Error: {directory_path} is not a valid directory", file=sys.stderr)
        return
    
    # Find all FASTA files in the directory
    fasta_patterns = ['*.fa', '*.fasta', '*.fa.gz', '*.fasta.gz']
    fasta_files = []
    for pattern in fasta_patterns:
        fasta_files.extend(directory_path.glob(pattern))
    
    if not fasta_files:
        print(f"No FASTA files found in {directory_path}", file=sys.stderr)
        return
    
    # Sort files for consistent ordering
    fasta_files = sorted(fasta_files)
    
    print(f"Processing {len(fasta_files)} FASTA files in {directory_path}")
    print(f"Files: {[f.name for f in fasta_files]}")
    
    # Read all sequences from all files
    all_sequences = {}  # {file_index: {contig_id: sequence}}
    contig_ids = set()
    
    for i, fasta_file in enumerate(fasta_files):
        sequences = read_fasta_file(fasta_file)
        if sequences:
            all_sequences[i] = sequences
            contig_ids.update(sequences.keys())
    
    if not all_sequences:
        print(f"No sequences could be read from {directory_path}", file=sys.stderr)
        return
    
    # Ensure all files have the same contigs
    for file_idx, sequences in all_sequences.items():
        missing_contigs = contig_ids - set(sequences.keys())
        if missing_contigs:
            print(f"Warning: File {fasta_files[file_idx].name} missing contigs: {missing_contigs}", 
                  file=sys.stderr)
    
    # Generate output filename
    if output_filename is None:
        output_filename = directory_path / f"{directory_path.name}_varianthaplos.txt"
    else:
        output_filename = Path(output_filename)
    
    # Process each contig
    total_variants = 0
    total_length = 0
    
    with open(output_filename, 'w') as outfile:
        # Write header
        num_haplotypes = len(fasta_files)
        hap_columns = '\t'.join([f'hap_{i+1}' for i in range(num_haplotypes)])
        
        # Calculate total sequence length for header
        for contig_id in sorted(contig_ids):
            if 0 in all_sequences and contig_id in all_sequences[0]:
                total_length += len(all_sequences[0][contig_id])
        
        header = f"Block length  {total_length}\tvar_id\tcontig\tvarpos\tref_allele\talt_allele\t{hap_columns}\n"
        outfile.write(header)
        
        var_id = 1
        
        # Process each contig
        for contig_id in sorted(contig_ids):
            # Get sequences for this contig from all files
            haplotype_sequences = []
            for file_idx in range(len(fasta_files)):
                if file_idx in all_sequences and contig_id in all_sequences[file_idx]:
                    haplotype_sequences.append(all_sequences[file_idx][contig_id])
                else:
                    print(f"Warning: Missing {contig_id} in file {fasta_files[file_idx].name}", 
                          file=sys.stderr)
                    # Use empty sequence as placeholder
                    haplotype_sequences.append("")
            
            # Find variants for this contig
            if haplotype_sequences and all(haplotype_sequences):
                variants = find_variants(haplotype_sequences, contig_id)
                
                # Write variants to file
                for pos, ref_allele, alt_allele, hap_calls in variants:
                    hap_calls_str = '\t'.join(map(str, hap_calls))
                    line = f"{var_id}\t{contig_id}\t{pos}\t{ref_allele}\t{alt_allele}\t{hap_calls_str}\n"
                    outfile.write(line)
                    var_id += 1
                    total_variants += 1
    
    print(f"Generated {output_filename}")
    print(f"Total variants found: {total_variants}")
    print(f"Total sequence length: {total_length}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate _varianthaplos.txt files from FASTA directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process single directory
    python generate_varianthaplos.py /path/to/haplotypes/
    
    # Process multiple directories
    python generate_varianthaplos.py dir1/ dir2/ dir3/
    
    # Process with custom output name
    python generate_varianthaplos.py /path/to/haplotypes/ -o custom_variants.txt

Each directory should contain FASTA files representing different haplotypes
of the same genome/region. The script will identify variant positions by
comparing sequences across all FASTA files in each directory.
        """
    )
    
    parser.add_argument(
        'directories',
        nargs='+',
        help='Directories containing FASTA files to process'
    )
    
    parser.add_argument(
        '-o', '--output',
        help='Custom output filename (only for single directory)'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    if len(args.directories) > 1 and args.output:
        print("Error: Custom output filename can only be used with a single directory", file=sys.stderr)
        sys.exit(1)
    
    # Process each directory
    for directory in args.directories:
        try:
            if args.verbose:
                print(f"\n=== Processing {directory} ===")
            
            output_filename = args.output if len(args.directories) == 1 else None
            process_directory(directory, output_filename)
            
        except Exception as e:
            print(f"Error processing {directory}: {e}", file=sys.stderr)
            if args.verbose:
                import traceback
                traceback.print_exc()


if __name__ == "__main__":
    main()