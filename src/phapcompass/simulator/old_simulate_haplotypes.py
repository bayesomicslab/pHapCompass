import argparse
import subprocess
import itertools
import os

def simulate_one(reference_path, output_dir, model, mutation_rate, mut_map, ploidy, haplogenerator_path, shifted=False):
    # strip off final number in output_dir if present
    output_file = output_dir
    output_dir = os.path.dirname(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if shifted:
        subprocess.run([
            "python", haplogenerator_path,
            "-f", reference_path,
            "-o", output_file,
            "--model", model,
            "--subfreq", f'[{mutation_rate}, 0, 0]',
            "-m", str(mut_map),
            "-p", str(ploidy),
            "--region-start", "500000",
            "--region-end", "1000000"
        ])
    else:
        subprocess.run([
            "python", haplogenerator_path,
            "-f", reference_path,
            "-o", output_file,
            "--model", model,
            "--subfreq", f'[{mutation_rate}, 0, 0]',
            "-m", str(mut_map),
            "-p", str(ploidy)
        ])
    return

def autopolyploidy_simulator(args):
    if args.mutation_rates is None:
        args.mutation_rates = [0.01, 0.001, 0.005]
    if args.ploidies is None:
        args.ploidies = [2, 3, 4, 6]
    for item in itertools.product(args.ploidies, args.mutation_rates, range(args.num_samples)):
        ploidy = item[0]
        mutation_rate = item[1]
        sample_index = str(item[2])
        output_dir = os.path.join(args.output_dir, f"auto/ploidy_{ploidy}/mutation_{mutation_rate}/{sample_index}/{sample_index}")
        simulate_one(args.reference_path, output_dir, 'poisson', mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', ploidy, args.haplogenerator_path, shifted=args.shifted)
    return

def allopolyploidy_simulator(args):
    if args.mutation_rates is None:
        args.mutation_rates = [0.00005, 0.0001]
    if args.sg_rates is None:
        args.sg_rates = [0.0005, 0.0001]
    args.ploidies = [3, 4, 6] # figure out how to handle inputting specific subgenome configurations later

    # first simulate subgenomes
    for mut in args.sg_rates:
        for sg in ['A', 'B', 'C']:
            output_dir = os.path.join(args.output_dir, f"allo/subgenome/mutation_{mut}/{sg}/{sg}")
            simulate_one(args.reference_path, output_dir, 'poisson', mut, '{"A":"C","G":"A","C":"T","T":"G"}', 1, args.haplogenerator_path, shifted=args.shifted)

    # ploidy 3
    for item in itertools.product(args.sg_rates, args.mutation_rates, range(args.num_samples)):
        sg_rate = item[0]
        mutation_rate = item[1]
        sample_index = str(item[2])
        ref_A_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/A/A_hap1.fa")
        ref_B_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/B/B_hap1.fa")
        ref_C_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/C/C_hap1.fa")
        output_dir = os.path.join(args.output_dir, f"allo/ploidy_3/sg_mut_{sg_rate}/mutation_{mutation_rate}/{sample_index}/{sample_index}")
        simulate_one(ref_A_path, output_dir+'A', 'poisson', mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, args.haplogenerator_path)
        simulate_one(ref_B_path, output_dir+'B', 'poisson', mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 1, args.haplogenerator_path)

    # ploidy 4
    for item in itertools.product(args.sg_rates, args.mutation_rates, range(args.num_samples)):
        sg_rate = item[0]
        mutation_rate = item[1]
        sample_index = str(item[2])
        ref_A_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/A/A_hap1.fa")
        ref_B_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/B/B_hap1.fa")
        ref_C_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/C/C_hap1.fa")
        output_dir = os.path.join(args.output_dir, f"allo/ploidy_4/sg_mut_{sg_rate}/mutation_{mutation_rate}/{sample_index}/{sample_index}")
        simulate_one(ref_A_path, output_dir+'A', 'poisson', mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, args.haplogenerator_path)
        simulate_one(ref_B_path, output_dir+'B', 'poisson', mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, args.haplogenerator_path)

    for item in itertools.product(args.sg_rates, args.mutation_rates, range(args.num_samples)):
        sg_rate = item[0]
        mutation_rate = item[1]
        sample_index = str(item[2])
        ref_A_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/A/A_hap1.fa")
        ref_B_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/B/B_hap1.fa")
        ref_C_path = os.path.join(args.output_dir, f"allo/subgenome/mutation_{sg_rate}/C/C_hap1.fa")
        output_dir = os.path.join(args.output_dir, f"allo/ploidy_6/sg_mut_{sg_rate}/mutation_{mutation_rate}/{sample_index}/{sample_index}")
        simulate_one(ref_A_path, output_dir+'A', 'poisson', mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, args.haplogenerator_path)
        simulate_one(ref_B_path, output_dir+'B', 'poisson', mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, args.haplogenerator_path)
        simulate_one(ref_C_path, output_dir+'C', 'poisson', mutation_rate, '{"A":"C","G":"A","C":"T","T":"G"}', 2, args.haplogenerator_path)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    batch_script = os.path.join(script_dir, "batch_varianthaplos.py")
    subprocess.run(["python", batch_script, os.path.join(args.output_dir, f"allo")])
    return

def main():
    parser = argparse.ArgumentParser(description="Simulate haplotypes")
    
    parser.add_argument('--reference_path', type=str, required=True, help='Path to the reference genome')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the simulated haplotypes')
    parser.add_argument('--haplogenerator_path', type=str, required=False, default='haplogenerator/haplogenerator_3.py', help='Path to the haplogenerator tool')
    parser.add_argument('--structure', type=str, required=True, help='Structure of the haplotypes (autopolyploidy or allopolyploidy)')
    parser.add_argument('--num_samples', type=int, required=True, help='Number of samples to simulate for each configuration')
    parser.add_argument('--mutation_rates', type=float, nargs='+', default=None, required=False, help='List of mutation rates to simulate')
    parser.add_argument('--shifted', type=bool, default=False, required=False, help='When true, simulate mutations in a shifted region of the genome (true for data in RECOMB 2026 paper)')
    parser.add_argument('--sg_rates', type=float, nargs='+', default=None, required=False, help='List of subgenome mutation rates to simulate in allopolyploidy')
    parser.add_argument('--ploidies', type=int, nargs='+', default=None, required=False, help='List of ploidy levels to simulate')

    args = parser.parse_args()

    if args.structure == "autopolyploidy":
        autopolyploidy_simulator(args)
    elif args.structure == "allopolyploidy":
        allopolyploidy_simulator(args)

if __name__ == "__main__":
    main()