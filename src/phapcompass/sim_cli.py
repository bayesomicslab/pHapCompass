#!/usr/bin/env python3
import argparse
from phapcompass.simulator import simulate_haplotypes


def build_parser():
    p = argparse.ArgumentParser(
        prog="phapcompass simulation",
        description="Simulation utilities for phapcompass"
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    hap = sub.add_parser("haplotypes", help="Simulate haplotype references")
    hap.add_argument("--reference_path", required=True, type=str)
    hap.add_argument("--output_dir", required=True, type=str)
    hap.add_argument("--structure", required=True, choices=["autopolyploidy", "allopolyploidy"])
    hap.add_argument("--num_samples", required=True, type=int)
    hap.add_argument("--mutation_rates", nargs="+", type=float, default=None)
    hap.add_argument("--sg_rates", nargs="+", type=float, default=None)
    hap.add_argument("--ploidies", nargs="+", type=int, default=None)

    reads = sub.add_parser("reads", help="Simulate reads (not wired yet)")
    reads.add_argument("--not_implemented", action="store_true", help=argparse.SUPPRESS)

    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "haplotypes":
        mod_argv = [
            "--reference_path", args.reference_path,
            "--output_dir", args.output_dir,
            "--structure", args.structure,
            "--num_samples", str(args.num_samples),
        ]
        if args.mutation_rates is not None:
            mod_argv += ["--mutation_rates", *[str(x) for x in args.mutation_rates]]
        if args.sg_rates is not None:
            mod_argv += ["--sg_rates", *[str(x) for x in args.sg_rates]]
        if args.ploidies is not None:
            mod_argv += ["--ploidies", *[str(x) for x in args.ploidies]]

        simulate_haplotypes.main(mod_argv)
        return 0

    if args.cmd == "reads":
        raise SystemExit("reads simulation is not wired yet. Use `phapcompass-sim haplotypes ...` for now.")

    raise SystemExit(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
