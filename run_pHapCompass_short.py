from pHapCompass_short.combinatorials import *
from pHapCompass_short.FFBS import *
from pHapCompass_short.viterbi import *
from pHapCompass_short.predict_hap import *
from pHapCompass_short.matching import *
from evaluations import *
from utils import *
from read_input import *
import pandas as pd
import numpy as np



def run_pHapCompass_short(args):
    """
    End-to-end short-reads pipeline with canonical naming:
      - Nodes: '{i}-{j}' with 1-based SNP positions, i<j
      - Edges: '{i}-{j}--{j}-{k}' (topological)
    Uses vectorized emissions/transitions, forward pass, and FFBS sampler.
    """

    frag = args.fragment_file_path
    vcf  = args.vcf_path
    ploidy   = args.ploidy
    # args.error_rate  = args.epsilon
    args.error_rate = 0.001
    result_path  = args.result_path

    print(f'Working on {frag} ...')

    gen_df = pd.read_csv(args.genotype_path)
    n_snps = len(gen_df)

    # ------------------- 2) Read fragments (.frag) into sparse CSR ------------------------
    cfg  = InputConfigSparse(data_path=args.data_path, genotype_path=args.genotype_path, ploidy=args.ploidy)
    frag = SparseFragment(cfg, positions_from_genotype=list(range(n_snps)))         # NOTE: _ingest_block must use START_IS_ONE_BASED=True
    data_matrix = frag.csr                    # (reads × SNPs) with {0,1(REF),2(ALT)}

    # SNP-index mappings:
    col_to_snp = np.array(frag.col_to_snp, dtype=np.int64)


    g = gen_df.sum(axis=1).to_numpy(dtype=np.int16)   # length = #SNP positions
    g_per_snp = g.astype(np.int16)

    # ------------------- 4) Build pair layer (nodes, counts4, edges, triple index) -------
    pair_layer = build_pair_layer(data_matrix, min_cocov=1)     # uses O^T O etc.

    # genotype per node, aligned by POSITIONS (not raw columns!)
    nodes_cols = pair_layer.nodes.astype(int)         # (P,2) column indices (0-based)
    gi = g_per_snp[col_to_snp[nodes_cols[:, 0]]]      # g at SNP position for col_i
    gj = g_per_snp[col_to_snp[nodes_cols[:, 1]]]
    genotype_pairs = list(zip(gi.tolist(), gj.tolist()))

    # ------------------- 5) Bank + node emissions (vectorized) ---------------------------
    bank = build_node_phasing_bank(args.ploidy, args.error_rate, genotype_pairs)
    state_names = build_state_names_from_bank(pair_layer, gi, gj, bank, args.ploidy, True, frag.col_to_snp)

    transitions_dict, transitions_extra = compute_transitions_optimized(data_matrix, pair_layer, state_names, args.ploidy, args.error_rate, col_to_snp, frag.snp_to_col, return_extra=True)
    emissions = build_pair_emissions_from_state_names(state_names=state_names, ploidy=args.ploidy, error_rate=args.error_rate)

    # Build nodes/edges lists from canonical dictionaries
    nodes  = list(emissions.keys())
    edges  = [(u, v) for (u, v) in (k.split("--") for k in transitions_dict.keys())]

    # ------------------- 8) Slices (DAG layering) ----------------------------------------
    slices, _ = assign_slices_and_interfaces(nodes, edges)

    # evidence indices per node/edge from M using position→column map
    assignment_dict = assign_evidence_to_states_and_transitions_from_M(nodes, edges, data_matrix, frag.snp_to_col)

    if args.uncertainty is not None:
        predicted_haplotypes = []
        block_ids = []
        likelihoods = []
        forward_messages = compute_forward_messages_from_M(slices,edges, assignment_dict, emissions, transitions_dict, data_matrix, frag.snp_to_col)
        for _ in range(args.uncertainty):
            samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
            predicted_haplotype, block_id = predict_haplotypes_mec_based(nodes, edges, samples, ploidy, args.genotype_path, data_matrix, frag.snp_to_col, transitions_extra, args, 
            priority="combined", likelihood_weight= args.lw, mec_weight=args.mw, ffbs_weight=args.sw, allow_ffbs_override=True, verbose=False)
            predicted_haplotype = predicted_haplotype.apply(pd.to_numeric, errors='coerce') 
            likelihood = compute_likelihood()  # ... > I should complete

            predicted_haplotypes.append(predicted_haplotype)
            block_ids.append(block_id)
            likelihoods.append(likelihood)

    else:

        best_states = viterbi_decode(slices, edges, assignment_dict, emissions, transitions_dict, data_matrix, frag.snp_to_col)
        predicted_haplotype, block_ids = predict_haplotypes_mec_based(nodes, edges, best_states, ploidy, args.genotype_path, data_matrix, frag.snp_to_col, transitions_extra, args, 
        priority="combined", likelihood_weight= args.lw, mec_weight=args.mw, ffbs_weight=args.sw, allow_ffbs_override=True, verbose=False)
        predicted_haplotype = predicted_haplotype.apply(pd.to_numeric, errors='coerce') 
        likelihood = compute_likelihood() #... > I should complete
    write_phased_vcf(input_vcf_path=vcf, output_vcf_path=result_path, predicted_haplotypes=predicted_haplotypes if args.uncertainty else predicted_haplotype,
        block_ids=block_ids if args.uncertainty else block_ids, likelihoods=likelihoods if args.uncertainty else likelihood, ploidy=ploidy)





