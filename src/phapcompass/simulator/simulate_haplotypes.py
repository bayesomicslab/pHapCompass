#!/usr/bin/env python3
"""
Simulate haplotypes for auto- and allopolyploid configurations.

This version assumes haplogenerator_3.py lives inside the package at:
    phapcompass.simulator.haplogenerator_3

So we invoke it via:
    python -m phapcompass.simulator.haplogenerator_3 ...

Key changes vs your old script:
- Removed haplogenerator_path arg (no filesystem paths).
- Uses sys.executable + -m module invocation.
- Adds check=True to fail fast.
- Fixes sg_rate scoping in allopolyploidy_simulator.
"""

import argparse
import subprocess
import itertools
import os
import sys


HAPLOGEN_MODULE = "phapcompass.simulator.haplogenerator_3"
BATCH_VARIANTHAPLOS_MODULE = "phapcompass.simulator.batch_varianthaplos"


def simulate_one(reference_path, output_dir, model, mutation_rate, mut_map, ploidy, shifted=False):
    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        sys.executable, "-m", HAPLOGEN_MODULE,
        "-f", reference_path,
        "-o", output_dir,
        "--model", model,
        "--s", f"[{mutation_rate}, 0, 0]",
        "-m", str(mut_map),
        "-p", str(ploidy),
    ]

    if shifted:
        cmd += ["--region-start", "500000", "--region-end", "1000000"]

    subprocess.run(cmd, check=True)
    return


def autopolyploidy_simulator(args):
    if args.mutation_rates is None:
        args.mutation_rates = [0.01, 0.001, 0.005]
    if args.ploidies is None:
        args.ploidies = [2, 3, 4, 6]

    for ploidy, mutation_rate, sample in itertools.product(args.ploidies, args.mutation_rates, range(args.num_samples)):
        sample_index = str(sample)
        output_dir = os.path.join(args.output_dir, f"auto/ploidy_{ploidy}/mutation_{mutation_rate}/{sample_index}/{sample_index}")
        simulate_one(args.reference_path, output_dir, model="poisson", mutation_rate=mutation_rate, mut_map='{"A":"C","G":"A","C":"T","T":"G"}', ploidy=ploidy, shifted=True)
    return


def allopolyploidy_simulator(args):
    if args.mutation_rates is None:
        args.mutation_rates = [0.00005, 0.0001]
    if args.sg_rates is None:
        args.sg_rates = [0.0005, 0.0001]

    # NOTE: You hard-coded supported allopolyploid ploidies in the original script
    # (3, 4, 6) with fixed subgenome configurations. Keep that behavior.
    args.ploidies = [3, 4, 6]

    # For each subgenome mutation rate, generate subgenomes and then derive haplotypes
    for sg_rate in args.sg_rates:
        # 1) simulate subgenomes A/B/C
        for sg in ["A", "B", "C"]:
            output_dir = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/{sg}/{sg}")
            simulate_one(args.reference_path, output_dir, model="poisson", mutation_rate=sg_rate, mut_map='{"A":"C","G":"A","C":"T","T":"G"}', ploidy=1, shifted=True)

        # 2) paths to subgenome references (now sg_rate is defined correctly)
        ref_A_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/A/A_hap1.fa")
        ref_B_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/B/B_hap1.fa")
        ref_C_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/C/C_hap1.fa")

        # ploidy 3 = A(2) + B(1)
        for mutation_rate, sample in itertools.product(args.mutation_rates, range(args.num_samples)):
            sample_index = str(sample)
            base_out = os.path.join(args.output_dir, f"allo/ploidy_3/sg_mut_{sg_rate}/mutation_{mutation_rate}/{sample_index}/{sample_index}")
            simulate_one(ref_A_path, base_out + "A", "poisson", mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, shifted=False)
            simulate_one(ref_B_path, base_out + "B", "poisson", mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 1, shifted=False)

        # ploidy 4 = A(2) + B(2)
        for mutation_rate, sample in itertools.product(args.mutation_rates, range(args.num_samples)):
            sample_index = str(sample)
            base_out = os.path.join(args.output_dir, f"allo/ploidy_4/sg_mut_{sg_rate}/mutation_{mutation_rate}/{sample_index}/{sample_index}")
            simulate_one(ref_A_path, base_out + "A", "poisson", mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, shifted=False)
            simulate_one(ref_B_path, base_out + "B", "poisson", mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, shifted=False)

        # ploidy 6 = A(2) + B(2) + C(2)
        for mutation_rate, sample in itertools.product(args.mutation_rates, range(args.num_samples)):
            sample_index = str(sample)
            base_out = os.path.join(args.output_dir, f"allo/ploidy_6/sg_mut_{sg_rate}/mutation_{mutation_rate}/{sample_index}/{sample_index}")
            simulate_one(ref_A_path, base_out + "A", "poisson", mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, shifted=False)
            simulate_one(ref_B_path, base_out + "B", "poisson", mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, shifted=False)
            simulate_one(ref_C_path, base_out + "C", "poisson", mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, shifted=False)

    # 3) post-process (kept from your original script, but invoked as a module)
    subprocess.run([sys.executable, "-m", BATCH_VARIANTHAPLOS_MODULE, os.path.join(args.output_dir, "allo")], check=True)
    return


def build_parser():
    parser = argparse.ArgumentParser(description="Simulate haplotypes")

    parser.add_argument("--reference_path", type=str, required=True, help="Path to the reference genome FASTA")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the simulated haplotypes")
    parser.add_argument("--structure", type=str, required=True, choices=["autopolyploidy", "allopolyploidy"], help="Structure of the haplotypes")
    parser.add_argument("--num_samples", type=int, required=True, help="Number of samples to simulate for each configuration")
    parser.add_argument("--mutation_rates", type=float, nargs="+", default=None, required=False, help="List of mutation rates to simulate")
    parser.add_argument("--sg_rates", type=float, nargs="+", default=None, required=False, help="List of subgenome mutation rates to simulate in allopolyploidy")
    parser.add_argument("--ploidies", type=int, nargs="+", default=None, required=False, help="List of ploidy levels to simulate (used for autopolyploidy only)")

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.structure == "autopolyploidy":
        autopolyploidy_simulator(args)
    elif args.structure == "allopolyploidy":
        allopolyploidy_simulator(args)
    else:
        raise ValueError(f"Unknown structure: {args.structure}")


if __name__ == "__main__":
    main()
