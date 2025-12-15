#!/usr/bin/env python3
"""
Batch process allopolyploid directories to generate varianthaplos.txt files.

This script processes the output directories from run_haplogenerator_3_allo.sh
and generates varianthaplos.txt files for each sample directory.

Usage:
    python batch_varianthaplos.py <allo_root_directory>
    
The script will find all sample directories (containing multiple .fa.gz files)
and generate a varianthaplos.txt file for each one.
"""

import os
import sys
import argparse
from pathlib import Path
import subprocess
from multiprocessing import Pool
from functools import partial
import time


def find_sample_directories(root_path):
    """
    Find all directories that contain multiple FASTA files (samples).
    
    Args:
        root_path: Root directory to search
        
    Returns:
        list: List of Path objects for sample directories
    """
    root_path = Path(root_path)
    sample_dirs = []
    
    for dirpath, dirnames, filenames in os.walk(root_path):
        dirpath = Path(dirpath)
        
        # Count FASTA files in this directory
        fasta_files = []
        for pattern in ['*.fa', '*.fasta', '*.fa.gz', '*.fasta.gz']:
            fasta_files.extend(dirpath.glob(pattern))
        
        # If directory has 2 or more FASTA files, it's likely a sample directory
        if len(fasta_files) >= 2:
            sample_dirs.append(dirpath)
    
    return sorted(sample_dirs)


def process_single_directory(sample_dir, force_overwrite=False, verbose=False):
    """
    Process a single sample directory to generate varianthaplos.txt file.
    
    Args:
        sample_dir: Path to the sample directory
        force_overwrite: Whether to overwrite existing varianthaplos.txt files
        verbose: Enable verbose output
        
    Returns:
        tuple: (status, sample_dir, message) where status is 'success', 'skip', or 'error'
    """
    try:
        sample_dir = Path(sample_dir)
        
        # Check if varianthaplos.txt already exists
        existing_file = sample_dir / f"{sample_dir.name}_varianthaplos.txt"
        if existing_file.exists() and not force_overwrite:
            return ('skip', sample_dir, 'output exists')
        
        # Call the generate_varianthaplos.py script
        script_path = Path(__file__).parent / "generate_varianthaplos.py"
        cmd = [sys.executable, str(script_path), str(sample_dir)]
        
        if verbose:
            cmd.append("--verbose")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            return ('success', sample_dir, 'processed successfully')
        else:
            error_msg = result.stderr.strip() or result.stdout.strip() or 'unknown error'
            return ('error', sample_dir, error_msg)
    
    except Exception as e:
        return ('error', sample_dir, str(e))


def process_allo_directories(root_path, force_overwrite=False, verbose=False):
    """
    Process all sample directories in an allo output structure using multiprocessing.
    
    Args:
        root_path: Root directory containing allo output
        force_overwrite: Whether to overwrite existing varianthaplos.txt files
        verbose: Enable verbose output
    """
    root_path = Path(root_path)
    
    if not root_path.exists():
        print(f"Error: {root_path} does not exist", file=sys.stderr)
        return
    
    print(f"Searching for sample directories in {root_path}")
    sample_dirs = find_sample_directories(root_path)
    
    if not sample_dirs:
        print(f"No sample directories found in {root_path}", file=sys.stderr)
        return
    
    print(f"Found {len(sample_dirs)} sample directories")
    print(f"Processing with 10 parallel workers...")
    
    success_count = 0
    skip_count = 0
    error_count = 0
    total_dirs = len(sample_dirs)
    
    # Create a partial function with the fixed arguments
    worker_func = partial(process_single_directory, 
                         force_overwrite=force_overwrite, 
                         verbose=verbose)
    
    start_time = time.time()
    
    # Use multiprocessing pool with 10 workers
    with Pool(processes=10) as pool:
        # Use imap for progress tracking
        results = pool.imap(worker_func, sample_dirs)
        
        # Process results as they come in
        for i, (status, sample_dir, message) in enumerate(results, 1):
            if status == 'success':
                success_count += 1
                if verbose:
                    print(f"✓ Successfully processed {sample_dir}")
            elif status == 'skip':
                skip_count += 1
                if verbose:
                    print(f"Skipping {sample_dir} ({message})")
            elif status == 'error':
                error_count += 1
                print(f"✗ Error processing {sample_dir}: {message}", file=sys.stderr)
            
            # Show progress every 50 directories (or at the end) if not verbose
            if not verbose and (i % 50 == 0 or i == total_dirs):
                elapsed = time.time() - start_time
                rate = i / elapsed if elapsed > 0 else 0
                eta = (total_dirs - i) / rate if rate > 0 else 0
                print(f"Progress: {i}/{total_dirs} directories processed "
                      f"({i/total_dirs*100:.1f}%) - "
                      f"Rate: {rate:.1f}/sec - "
                      f"ETA: {eta/60:.1f}min")
    
    elapsed = time.time() - start_time
    print(f"\nCompleted in {elapsed/60:.1f} minutes")
    print(f"Summary:")
    print(f"  Successfully processed: {success_count}")
    print(f"  Skipped (existing): {skip_count}")
    print(f"  Errors: {error_count}")
    print(f"  Total directories: {total_dirs}")


def main():
    parser = argparse.ArgumentParser(
        description="Batch generate varianthaplos.txt files for allopolyploid directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python batch_varianthaplos.py /path/to/allo_output/
    
This will process all sample directories under the allo output directory
and generate varianthaplos.txt files for each one.
        """
    )
    
    parser.add_argument(
        'root_directory',
        help='Root directory containing allopolyploid sample directories'
    )
    
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite existing varianthaplos.txt files'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    try:
        process_allo_directories(args.root_directory, args.force, args.verbose)
    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()