#!/usr/bin/env python3
import argparse
import os
import sys
import logging
from .run_pHapCompass_short import run_pHapCompass_short
from .run_pHapCompass_long import run_pHapCompass_long
from .read_input import *
from pathlib import Path
THIS_DIR = Path(__file__).resolve().parent          
PACKAGE_ROOT = THIS_DIR.parent                      
REPO_ROOT = PACKAGE_ROOT.parent   

EXTRACT_HAIRS_BIN = (REPO_ROOT / "third_party" / "extract_poly" / "build" / "extractHAIRS").resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="pHapCompass main entry point for short- and long-read haplotype assembly.")

    # Inputs
    parser.add_argument("--bam-path", type=str,  default=None, help="Path to BAM file. If provided with --vcf-path and no --fragment-file-path, extractHAIRS will be run.")
    parser.add_argument("--frag-path", type=str, default=None, help="Path to fragment file (.frag). If provided, BAM is optional.")
    parser.add_argument("--vcf-path", type=str, required=True, help="Path to VCF file (required: used for genotypes and/or ploidy inference).")
    parser.add_argument("--ploidy", type=int,  default=None, help="Ploidy (e.g., 2, 3, 4). If not provided, inferred from VCF.")
    parser.add_argument("--result-path", type=str, required=True, help="Output path (file path or prefix) for saving final haplotypes/results.")
    parser.add_argument("--epsilon", type=float, default=0.00001, help="Sequencing error rate ε (used by both models). Default: 0.00001")
    parser.add_argument("--data-type", type=str, choices=["short", "long"], required=True, help="Type of data: 'short' for short reads, 'long' for long reads.")
    # Uncertainty (optional int, default 3 if flag used without value)
    parser.add_argument("--uncertainty", nargs="?", const=3, type=int, default=None, help=("Enable uncertainty quantification. Optional value specifies the number of haplotype samples (default: 3 if flag is provided without a value)."))


    # Hyperparameters for short-read model
    parser.add_argument("--mw", type=float, default=10.0, help="MEC weight for short-read model.")
    parser.add_argument("--lw", type=float, default=1.0, help="Likelihood weight for short-read model.")
    parser.add_argument("--sw", type=float, default=1.0, help="Samples weight for short-read model.")

    # Hyperparameters for long-read model
    parser.add_argument("--delta", type=float, default=5, help="δ hyperparameter for long-read model.")
    parser.add_argument("--learning-rate", type=float, default=0.02,help="Learning rate for long-read model.")

    # extractHAIRS options
    parser.add_argument("--mbq", type=int, default=13, help="Minimum base quality for extractHAIRS (default: 13).")
    parser.add_argument("--extracthairs-bin", type=str, default="extractHAIRS", help="Path to extractHAIRS binary (default: 'extractHAIRS' on PATH).")


    parser.add_argument("--verbose", action="store_true", help="Enable verbose mode (equivalent to --log-level DEBUG).")
    parser.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level.",)

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    log_level = "DEBUG" if args.verbose else args.log_level

    logging.basicConfig(level=getattr(logging, log_level), format="%(asctime)s | %(levelname)s | %(message)s")

    if not os.path.exists(args.vcf_path):
        logging.error("VCF file not found: %s", args.vcf_path)
        sys.exit(1)

    if args.frag_path is None and args.bam_path is None:
        logging.error("You must provide either --frag-path or --bam-path (with --vcf-path).")
        sys.exit(1)

    # Ensure output directory exists (if any)
    out_dir = os.path.dirname(args.result_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Infer ploidy if needed
    if args.ploidy is None:
        logging.info("Inferring ploidy from VCF: %s", args.vcf_path)
        inferred_ploidy = infer_ploidy_from_vcf(args.vcf_path)
        logging.info("Inferred ploidy: %d", inferred_ploidy)
        args.ploidy = inferred_ploidy

    # Generate genotype file (same folder as result-path, or current dir if none)
    # geno_dir = out_dir if out_dir else "."
    # args.genotype_path = os.path.join(geno_dir, "genotype.csv")
    # logging.info("Writing genotype file to: %s", args.genotype_path)

    args.genotype = vcf_gt_to_csv(args.vcf_path, ploidy=args.ploidy)


    # Prepare fragment file
    if args.frag_path is None:
        if args.bam_path is None:
            logging.error("Need --bam-path to generate fragment file with extractHAIRS.")
            sys.exit(1)

        bam_root = bam_root_name(args.bam_path)  # X from X.bam
        # Folder = same folder as result-path; if none, use current directory
        frag_dir = out_dir if out_dir else "."
        if args.mbq == 13:
            frag_filename = f"{bam_root}.frag"
        else:
            frag_filename = f"{bam_root}_mbq{args.mbq}.frag"
        frag_path = os.path.join(frag_dir, frag_filename)

        logging.info("No fragment file provided. Will generate using extractHAIRS at: %s", frag_path)
        args.extracthairs_bin = str(EXTRACT_HAIRS_BIN)
        run_extract_hairs(
            bam_path=args.bam_path,
            vcf_path=args.vcf_path,
            frag_path=frag_path,
            ploidy=args.ploidy,
            epsilon=args.epsilon,
            mbq=args.mbq,
            extracthairs_bin=args.extracthairs_bin,
        )
        args.frag_path = frag_path
    else:
        if not os.path.exists(args.frag_path):
            logging.error("Fragment file not found: %s", args.frag_path)
            sys.exit(1)
        logging.info("Using existing fragment file: %s", args.frag_path)

    logging.info("Data type: %s", args.data_type)
    logging.info("Result path: %s", args.result_path)

    if args.data_type == "short":
        logging.info("Running short-read model with mw=%.4f, lw=%.4f, sw=%.4f", args.mw, args.lw, args.sw)
        run_pHapCompass_short(args)

    elif args.data_type == "long":
        logging.info("Running long-read model with delta=%.4f, learning_rate=%.6f", args.delta, args.learning_rate)
        run_pHapCompass_long(args)

    if args.uncertainty is not None:
        logging.info("Uncertainty quantification enabled with %d samples.", args.uncertainty)
    else:
        logging.info("Uncertainty quantification disabled.")


    logging.info("Finished pHapCompass run.")


if __name__ == "__main__":
    main()
