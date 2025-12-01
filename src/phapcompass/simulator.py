import os, re, gzip, random, subprocess, shutil, gzip
import numpy as np
import pandas as pd
from collections import namedtuple


def extract_chromosome(fasta_path, chrom, out_path=None):
    """
    Extract one chromosome/contig from a plain-text FASTA by matching the
    first token after '>' (e.g., '>Chr8 something' -> 'Chr8').

    If out_path is given, writes the FASTA block there.
    Otherwise returns the FASTA-formatted string.
    """
    out_lines = []
    found = False
    write_block = False

    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                # End current block if we were writing and we hit a new header
                if found and not write_block:
                    break
                # Decide whether this header matches the requested chrom
                key = line[1:].split()[0] if len(line) > 1 else ""
                write_block = (key == chrom)
                if write_block:
                    # out_lines.append(line)
                    out_lines.append('>'+key+'\n')
                    found = True
            else:
                if write_block:
                    out_lines.append(line)

    if not found:
        raise ValueError(f"Contig '{chrom}' not found in {fasta_path}")

    if out_path:
        with open(out_path, "w", encoding="utf-8") as out:
            out.writelines(out_lines)
    else:
        return "".join(out_lines)


def old_tsv_chrom_to_vcf(tsv_path, chromosome, out_vcf_path, sample_name="SAMPLE"):
    """
    Convert a whitespace-delimited variants table to a VCF for one chromosome.
    Expected leading columns:
      var_id contig varpos ref_allele alt_allele hap_1 hap_2 [hap_3 hap_4 ...]
    Extra columns are ignored. Lines with too few columns are skipped.

    Writes a phased GT built from all hap_* columns present.
    """
    strip_prefix = re.compile(r'^\s*Block\s+length\s+\d+\s+')

    # --- minimal gzip-aware I/O ---
    def _open_in(p):
        return gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, "r", encoding="utf-8")
    def _open_out(p):
        return gzip.open(p, "wt", encoding="utf-8") if p.endswith(".gz") else open(p, "w", encoding="utf-8")

    with _open_out(out_vcf_path) as out, _open_in(tsv_path) as fh:

        # Minimal header
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##contig=<ID={chromosome}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype from hap_* columns (0=REF,1=ALT)">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        first_line = True
        for line in fh:
            if first_line:
                line = strip_prefix.sub("", line, count=1)  # remove "Block length <num> "
                first_line = False

            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()  # whitespace-split

            if parts[0].lower() == "var_id":
                continue  # skip header row
            if len(parts) < 7:           # supports diploid+
                continue

            var_id, contig, pos, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            if contig != chromosome:
                continue

            haps = parts[5:]             # take all hap_* columns present
            gt_fields = []
            for x in haps:
                if x in {"0", "1"}:
                    gt_fields.append(x)
                elif x == ".":
                    gt_fields.append(".")
                else:
                    gt_fields.append(".")  # unknown -> missing
            gt = "|".join(gt_fields)

            out.write(f"{contig}\t{pos}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")


def tsv_chrom_to_vcf(tsv_path, chromosome, out_vcf_path, sample_name="SAMPLE"):
    """
    Convert a whitespace-delimited variants table to a VCF for one chromosome.
    NOW ALSO FILTERS OUT HOMOZYGOUS SITES (all 0 or all 1)
    """
    strip_prefix = re.compile(r'^\s*Block\s+length\s+\d+\s+')

    def _open_in(p):
        return gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, "r", encoding="utf-8")
    def _open_out(p):
        return gzip.open(p, "wt", encoding="utf-8") if p.endswith(".gz") else open(p, "w", encoding="utf-8")

    with _open_out(out_vcf_path) as out, _open_in(tsv_path) as fh:
        # Minimal header
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##contig=<ID={chromosome}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype from hap_* columns (0=REF,1=ALT)">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        first_line = True
        num_written = 0
        num_skipped = 0
        
        for line in fh:
            if first_line:
                line = strip_prefix.sub("", line, count=1)
                first_line = False

            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()

            if parts[0].lower() == "var_id":
                continue
            if len(parts) < 7:
                continue

            var_id, contig, pos, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            if contig != chromosome:
                continue

            haps = parts[5:]
            gt_fields = []
            for x in haps:
                if x in {"0", "1"}:
                    gt_fields.append(x)
                elif x == ".":
                    gt_fields.append(".")
                else:
                    gt_fields.append(".")
            
            # SKIP HOMOZYGOUS SITES (all same allele)
            unique_gts = set(gt_fields)
            if len(unique_gts) == 1:
                num_skipped += 1
                continue
            
            gt = "|".join(gt_fields)
            out.write(f"{contig}\t{pos}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
            num_written += 1
        
        print(f"  VCF: wrote {num_written} variants, skipped {num_skipped} homozygous sites")


def generate_fasta_vcfs_from_tsv_new_auto():
    """
    Complete regeneration pipeline from varianthaplos files
    """
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    ploidies = [2, 3, 4, 6, 8]
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    chrom = 'Chr1'
    
    for ploidy in ploidies:
        for mr in mut_rates:
            for sample in samples:
                print(f'Ploidy: {ploidy}, mr: {mr}, sample: {sample}')
                haplotype_path = os.path.join(hap_path, f'ploidy_{ploidy}', f'mut_{mr}', str(sample).zfill(2))
                
                # Remove old files
                files_to_remove =  ['HPOPG_Chr1.vcf',  'Chr1.fa.gz',  'Chr1.vcf.gz', f'Chr1_ploidy_{ploidy}_sample_{str(sample).zfill(2)}.fa.gz']
                for ff in files_to_remove:
                    file_path = os.path.join(haplotype_path, ff)
                    if os.path.exists(file_path): 
                        os.remove(file_path)
                        print(f'  Removed: {ff}')

                # Input files
                varianthaplos_file = os.path.join(haplotype_path, f'{str(sample).zfill(2)}_varianthaplos.txt.gz')
                haplotype_names = [os.path.join(haplotype_path, f'{str(sample).zfill(2)}_hap{i+1}.fa.gz') for i in range(ploidy)]
                
                # Output files
                out_vcf_path = os.path.join(haplotype_path, f'{chrom}.vcf')
                fasta_file_name = os.path.join(haplotype_path, f'{chrom}_ploidy_{ploidy}_sample_{str(sample).zfill(2)}.fa')
                chrom_fasta_path = os.path.join(haplotype_path, f'{chrom}.fa')
                vcf_unphased = os.path.join(haplotype_path, chrom + '_unphased.vcf')
                genotype_path = os.path.join(haplotype_path, 'genotype.csv')
                
                # STEP 1: Filter homozygous sites from varianthaplos
                filtered_varianthaplos = filter_homozygous_sites(varianthaplos_file, ploidy)
                
                # STEP 2: Create VCF from filtered varianthaplos
                tsv_chrom_to_vcf(filtered_varianthaplos, chrom, out_vcf_path, sample_name=str(sample).zfill(2))
                
                # STEP 3: Create reference FASTA from haplotype 1 + VCF
                create_reference_from_haplotype_and_vcf(haplotype_names[0], out_vcf_path, chrom, chrom_fasta_path)
                
                # STEP 4: Create concatenated FASTA for ART
                write_concat_chromosome_simple(haplotype_names, chrom, fasta_file_name)
                
                # Clean up temp file
                if filtered_varianthaplos != varianthaplos_file:
                    os.remove(filtered_varianthaplos)

                # STEP 5: generate genotype csv file
                vcf_to_genotype(out_vcf_path, genotype_path, sample_name=str(sample).zfill(2), ploidy=ploidy)
                
                # STEP 6: generate unphased vcf to be used in other methods
                unphase_and_sort_vcf(out_vcf_path, vcf_unphased, sample_column=str(sample).zfill(2))             
                
                print(f'  ✓ Created VCF: {out_vcf_path}')
                print(f'  ✓ Created reference: {chrom_fasta_path}')
                print(f'  ✓ Created concat FASTA: {fasta_file_name}')


def generate_fasta_vcfs_from_tsv_new_allo():
    """
    Complete regeneration pipeline from varianthaplos files
    """
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    subgenome_configs = [0.01, 0.05]
    ploidies = [3, 4, 6, 8]
    samples = range(20)
    chrom = 'Chr1'
    within_mutation_rates = [0.001, 0.005, 0.0075, 0.01, 0.05]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}
    samples = range(20)
    chrom = 'Chr1'
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for sample in samples:
                    print('Ploidy:', ploidy, 'Subgenome_mr:', sgc, 'within_mr:', mr, 'sample:', sample)
                    # stop
                    haplotype_path = os.path.join(hap_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr) , str(sample).zfill(2))
                    # variation_txt_path = [os.path.join(haplotype_path, str(sample).zfill(2) + '_' + sbg + '_varianthaplos.txt.gz') for sbg in subgenomes]
                    # variation_txt_name = os.path.join(haplotype_path, str(sample).zfill(2) + '_varianthaplos.txt.gz')
                    
                    # Remove old files
                    files_to_remove =  ['HPOPG_Chr1.vcf',  'Chr1.fa.gz',  'Chr1.vcf.gz', f'Chr1_ploidy_{ploidy}_sample_{str(sample).zfill(2)}.fa.gz']
                    for ff in files_to_remove:
                        file_path = os.path.join(haplotype_path, ff)
                        if os.path.exists(file_path): 
                            os.remove(file_path)
                            print(f'  Removed: {ff}')
                        
                    
                    # Input files
                    fasta_files = [f for f in os.listdir(haplotype_path) if 'fa.gz' in f]
                    fasta_paths = [os.path.join(haplotype_path, f) for f in fasta_files]
                    fasta_names = ['haplotype_' + f.split('.fa.gz')[0].split('_')[2].split('hap')[1] + '_subgenome_' + f.split('.fa.gz')[0].split('_')[1] for f in fasta_files]
                    varianthaplos_file = os.path.join(haplotype_path, f'{str(sample).zfill(2)}_varianthaplos.txt.gz')
                    # haplotype_names = [os.path.join(haplotype_path, f'{str(sample).zfill(2)}_hap{i+1}.fa.gz') for i in range(ploidy)]
                
                    # Output files
                    out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                    fasta_file_name = os.path.join(haplotype_path, f'{chrom}_ploidy_{ploidy}_sample_{str(sample).zfill(2)}.fa')
                    chrom_fasta_path = os.path.join(haplotype_path, f'{chrom}.fa')
                    vcf_unphased = os.path.join(haplotype_path, chrom + '_unphased.vcf')
                    genotype_path = os.path.join(haplotype_path, 'genotype.csv')
                    
                    # STEP 1: Filter homozygous sites from varianthaplos
                    filtered_varianthaplos = filter_homozygous_sites(varianthaplos_file, ploidy)
                    
                    # STEP 2: Create VCF from filtered varianthaplos
                    tsv_chrom_to_vcf(filtered_varianthaplos, chrom, out_vcf_path, sample_name=str(sample).zfill(2))
                    
                    # STEP 3: Create reference FASTA from haplotype 1 + VCF
                    create_reference_from_haplotype_and_vcf(fasta_paths[0], out_vcf_path, chrom, chrom_fasta_path)
                    
                    # STEP 4: Create concatenated FASTA for ART
                    write_concat_chromosome_simple(fasta_paths, chrom, fasta_file_name, names_list=fasta_names)

                    # Clean up temp file
                    if filtered_varianthaplos != varianthaplos_file:
                        os.remove(filtered_varianthaplos)

                    # STEP 5: generate genotype csv file
                    vcf_to_genotype(out_vcf_path, genotype_path, sample_name=str(sample).zfill(2), ploidy=ploidy)
                    
                    # STEP 6: generate unphased vcf to be used in other methods
                    unphase_and_sort_vcf(out_vcf_path, vcf_unphased, sample_column=str(sample).zfill(2))             
                    
                    print(f'  ✓ Created VCF: {out_vcf_path}')
                    print(f'  ✓ Created reference: {chrom_fasta_path}')
                    print(f'  ✓ Created concat FASTA: {fasta_file_name}')


def filter_homozygous_sites(varianthaplos_path, ploidy):
    """
    Remove rows where all haplotypes have the same allele (all 0s or all 1s)
    Returns path to filtered file
    """
    import tempfile
    
    strip_prefix = re.compile(r'^\s*Block\s+length\s+\d+\s+')
    
    def _open_in(p):
        return gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, "r", encoding="utf-8")
    
    # Create temp file for filtered output
    temp_fd, temp_path = tempfile.mkstemp(suffix='.txt.gz')
    os.close(temp_fd)
    
    lines_kept = 0
    lines_removed = 0
    
    with gzip.open(temp_path, 'wt', encoding='utf-8') as out, _open_in(varianthaplos_path) as fh:
        first_line = True
        for line in fh:
            if first_line:
                # Keep the header with Block length
                out.write(line)
                first_line = False
                continue
            
            s = line.strip()
            if not s or s.startswith("#"):
                out.write(line)
                continue
            
            parts = s.split()
            
            # Header row
            if parts[0].lower() == "var_id":
                out.write(line)
                continue
            
            # Data rows: var_id contig varpos ref_allele alt_allele hap_1 hap_2 [hap_3 ...]
            if len(parts) < 5 + ploidy:
                out.write(line)  # Keep malformed lines
                continue
            
            # Extract haplotype alleles
            hap_alleles = parts[5:5+ploidy]
            
            # Check if all same (homozygous)
            unique_alleles = set(hap_alleles)
            
            # Remove if all 0 or all 1
            if len(unique_alleles) == 1 and list(unique_alleles)[0] in ['0', '1']:
                lines_removed += 1
                continue
            
            # Keep heterozygous sites
            lines_kept += 1
            out.write(line)
    
    print(f'  Filtered varianthaplos: kept {lines_kept}, removed {lines_removed} homozygous sites')
    return temp_path


def create_reference_from_haplotype_and_vcf(haplotype_path, vcf_path, chrom, out_ref_path):
    """
    Create reference FASTA by:
    1. Taking haplotype sequence as template
    2. Forcing REF allele at each VCF position
    
    This ensures reference matches what VCF calls "REF"
    """
    def _open_in(p):
        return gzip.open(p, "rt", encoding="utf-8") if p.endswith(".gz") else open(p, "r", encoding="utf-8")
    
    # Load haplotype sequence
    seq_parts = []
    take = False
    with _open_in(haplotype_path) as f:
        for line in f:
            if line.startswith('>'):
                header_chrom = line[1:].split()[0]
                take = (header_chrom == chrom)
                continue
            if take:
                seq_parts.append(''.join(line.strip().split()))
    
    if not seq_parts:
        raise ValueError(f"Chromosome {chrom} not found in {haplotype_path}")
    
    seq = bytearray(''.join(seq_parts).upper().encode('ascii'))
    
    # Apply REF alleles from VCF
    num_applied = 0
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            contig, pos, var_id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            
            if contig != chrom:
                continue
            
            # VCF is 1-based, convert to 0-based
            pos_0based = int(pos) - 1
            
            # Only apply single-nucleotide variants
            if len(ref) == 1 and ref in 'ACGT':
                if 0 <= pos_0based < len(seq):
                    seq[pos_0based] = ref.encode('ascii')[0]
                    num_applied += 1
    
    # Write reference
    with open(out_ref_path, 'w') as f:
        f.write(f'>{chrom}\n')
        f.write(seq.decode('ascii') + '\n')
    
    print(f'  Applied {num_applied} REF alleles to reference')


def merge_tsvs_to_vcf(tsv_paths, chromosome, out_vcf_path, sample_name="SAMPLE"):
    """
    Merge multiple whitespace-delimited variant TSVs for one chromosome into a single VCF.
    Each TSV has leading columns:
      var_id contig varpos ref_allele alt_allele hap_1 hap_2 [hap_3 hap_4 ...]
    Lines before the header like "Block length <num> ..." are ignored (first line only).
    - Keeps only rows where contig == `chromosome`
    - Builds GT from all hap_* columns present (phased: 0|1|...).
    - Sorts by POS (numeric).
    - If the same POS appears more than once across any files, raises ValueError and writes nothing.

    Input TSVs may be .gz or plain. Output VCF may be .gz or plain (based on out_vcf_path suffix).
    """

    strip_prefix = re.compile(r'^\s*Block\s+length\s+\d+\s+')

    def _open_in(p):
        return gzip.open(p, "rt", encoding="utf-8") if str(p).endswith(".gz") \
               else open(p, "r", encoding="utf-8")

    def _open_out(p):
        return gzip.open(p, "wt", encoding="utf-8") if str(p).endswith(".gz") \
               else open(p, "w", encoding="utf-8")

    records = []           # will hold tuples: (pos_int, contig, var_id, ref, alt, gt_string)
    seen_pos = set()       # to detect duplicates across files

    for tsv_path in tsv_paths:
        first_line = True
        with _open_in(tsv_path) as fh:
            for line in fh:
                if first_line:
                    line = strip_prefix.sub("", line, count=1)
                    first_line = False

                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()

                if parts[0].lower() == "var_id":
                    continue
                if len(parts) < 7:
                    continue

                var_id, contig, pos, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                if contig != chromosome:
                    continue

                try:
                    p_int = int(pos)
                except ValueError:
                    continue

                # duplicate-position check across all inputs
                if p_int in seen_pos:
                    raise ValueError(f"Duplicate position detected at {chromosome}:{p_int} across input files. Aborting without writing.")
                seen_pos.add(p_int)

                haps = parts[5:]  # all hap_* present
                gt_fields = []
                for x in haps:
                    if x in {"0", "1"}:
                        gt_fields.append(x)
                    elif x == ".":
                        gt_fields.append(".")
                    else:
                        gt_fields.append(".")
                gt = "|".join(gt_fields)

                records.append((p_int, contig, var_id, ref, alt, gt))

    # sort by numeric position
    records.sort(key=lambda r: r[0])

    # write output VCF
    with _open_out(out_vcf_path) as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##contig=<ID={chromosome}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype from hap_* columns (0=REF,1=ALT)">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        for p_int, contig, var_id, ref, alt, gt in records:
            out.write(f"{contig}\t{p_int}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")


def print_fasta_lengths(fasta_path, header_prefix=None):
    """
    Print 'header<TAB>length' for each FASTA record in fasta_path.
    If header_prefix is given (e.g., 'haplotype_1'), only print records whose
    header (first token after '>') starts with that prefix.
    """
    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                print(line, end="")           # print header line itself
            else:
                print(len(line.rstrip("\r\n")))  # length of sequence line (no newline)


def write_concat_chromosome_simple(fasta_paths, chrom, out_path, names_list=None):
    """
    From each FASTA in `fasta_paths` (plain or .gz), extract the block whose first header
    token == `chrom` and write a single multi-FASTA to `out_path` (gz if out_path endswith .gz).
    If `names_list` is provided, use those as headers in order; otherwise infer 'haplotype_{N}'
    from the filename (matching hap1/hap_2/haplotype-3, case-insensitive). Sequences are one line.
    """
    def _open_in(path):
        return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, "r", encoding="utf-8")

    def _open_out(path):
        return gzip.open(path, "wt", encoding="utf-8") if path.endswith(".gz") else open(path, "w", encoding="utf-8")

    fasta_paths = [os.path.expanduser(p.strip()) for p in fasta_paths]
    if names_list is not None:
        if len(names_list) != len(fasta_paths):
            raise ValueError("names_list length must match fasta_paths length")
        headers = [str(n).strip() for n in names_list]
    else:
        # infer hap number from filename
        headers = []
        for fp in fasta_paths:
            fname = os.path.basename(fp)
            m = re.search(r'hap(?:lotype[_-]?)?(\d+)', fname, re.IGNORECASE)
            if not m:
                raise ValueError(f"Cannot determine hap number from filename: {fname} (provide names_list)")
            headers.append(f"haplotype_{int(m.group(1))}")

    with _open_out(out_path) as out:
        for fp, header in zip(fasta_paths, headers):
            seq_chunks, writing, found = [], False, False
            with _open_in(fp) as fh:
                for line in fh:
                    if line.startswith(">"):
                        key = line[1:].split()[0] if len(line) > 1 else ""
                        writing = (key == chrom)
                        if writing:
                            found = True
                        continue
                    if writing:
                        seq_chunks.append("".join(line.split()))
            if not found:
                raise ValueError(f"Contig '{chrom}' not found in file: {fp}")
            seq = "".join(seq_chunks).upper()
            out.write(f">{header}\n{seq}\n")


def merge_tsvs_to_vcf_allow_same_refalt(tsv_paths, chromosome, out_vcf_path, sample_name="SAMPLE"):
    """
    Merge multiple whitespace-delimited variant TSVs for one chromosome into a single VCF.
    Each TSV has columns:
      var_id contig varpos ref_allele alt_allele hap_1 hap_2 [hap_3 ...]
    - Keeps contig == `chromosome`
    - GT is built from all hap_* columns present (phased '|' join)
    - Sorts by POS (numeric)
    - If the same POS appears across inputs with DIFFERENT REF/ALT, raise ValueError.
      If same POS has SAME REF/ALT, keep the FIRST occurrence and skip the rest.

    Input TSVs and output VCF can be .gz or plain (based on filename suffix).
    """
    strip_prefix = re.compile(r'^\s*Block\s+length\s+\d+\s+')

    def _open_in(p):
        return gzip.open(p, "rt", encoding="utf-8") if str(p).endswith(".gz") else open(p, "r", encoding="utf-8")

    def _open_out(p):
        return gzip.open(p, "wt", encoding="utf-8") if str(p).endswith(".gz") else open(p, "w", encoding="utf-8")

    records = []                 # (pos_int, contig, var_id, ref, alt, gt)
    pos_to_refalt = {}           # pos_int -> (ref, alt) of the FIRST accepted record

    for tsv_path in tsv_paths:
        first_line = True
        with _open_in(tsv_path) as fh:
            for line in fh:
                if first_line:
                    line = strip_prefix.sub("", line, count=1)  # strip leading "Block length <num> " if present
                    first_line = False

                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                parts = s.split()
                if parts[0].lower() == "var_id":
                    continue
                if len(parts) < 7:
                    continue

                var_id, contig, pos, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                if contig != chromosome:
                    continue
                try:
                    p_int = int(pos)
                except ValueError:
                    continue

                # REF/ALT comparison logic
                refU, altU = ref.upper(), alt.upper()
                if p_int in pos_to_refalt:
                    ref0, alt0 = pos_to_refalt[p_int]
                    if refU != ref0 or altU != alt0:
                        raise ValueError(f"Conflict at {chromosome}:{p_int} -> ({ref0},{alt0}) vs ({refU},{altU})")
                    # same REF/ALT at the same position: skip this duplicate record
                    continue
                else:
                    pos_to_refalt[p_int] = (refU, altU)

                # Build GT from all hap_* columns present
                haps = parts[5:]
                gt_fields = []
                for x in haps:
                    if x in {"0", "1"}:
                        gt_fields.append(x)
                    elif x == ".":
                        gt_fields.append(".")
                    else:
                        gt_fields.append(".")
                gt = "|".join(gt_fields)
                records.append((p_int, contig, var_id, ref, alt, gt))

    # sort by position
    records.sort(key=lambda r: r[0])

    # write VCF
    with _open_out(out_vcf_path) as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##contig=<ID={chromosome}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype from hap_* columns (0=REF,1=ALT)">\n')
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        for p_int, contig, var_id, ref, alt, gt in records:
            out.write(f"{contig}\t{p_int}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")


def fasta_from_vcf_ref(hap_fasta_path, vcf_path, chrom, out_fasta_path):
    """
    Build a single-FASTA for `chrom` by taking the haplotype FASTA as a template,
    forcing the REF allele at each VCF POS (VCF is 1-based). Writes header `>{chrom}`
    and the sequence on ONE line. Supports .fa/.fa.gz for input and output.
    Assumes SNVs; non-SNV records are ignored.
    """

    def _open_text(path, mode):
        # mode: "r" or "w"
        if path.endswith(".gz"):
            return gzip.open(path, mode + "t", encoding="utf-8")
        return open(path, mode, encoding="utf-8")

    # 1) Load template sequence for `chrom`
    seq_parts, take = [], False
    with _open_text(hap_fasta_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                take = (line[1:].split()[0] == chrom)
                continue
            if take:
                seq_parts.append("".join(line.split()))
    if not seq_parts:
        raise ValueError(f"Chromosome '{chrom}' not found in {hap_fasta_path}.")
    seq = bytearray("".join(seq_parts).upper().encode("ascii"))

    # 2) Apply REF at variant sites from VCF
    with _open_text(vcf_path, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            c, pos, _id, REF, ALT = parts[0], parts[1], parts[2], parts[3], parts[4]
            if c != chrom or len(REF) != 1:
                continue  # SNVs only
            try:
                p = int(pos) - 1  # VCF POS is 1-based
            except ValueError:
                continue
            if 0 <= p < len(seq):
                seq[p:p+1] = REF.upper().encode("ascii")

    # 3) Write FASTA: single header + single sequence line (NO WRAPPING)
    with _open_text(out_fasta_path, "w") as out:
        out.write(f">{chrom}\n")
        out.write(seq.decode("ascii") + "\n")


def generate_fasta_vcfs_auto():
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    ploidies = [2, 3, 4, 6, 8]
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    chrom = 'Chr1'
    for ploidy in ploidies:
        for mr in mut_rates:
            for sample in samples:
                # stop
                print('Ploidy:', ploidy, 'mr:', mr, 'sample:', sample)
                haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                haplotype_names = [os.path.join(haplotype_path, str(sample).zfill(2) + '_' + 'hap' + str(i + 1) + '.fa.gz') for i in range(ploidy)]
                # variation_txt_path = os.path.join(haplotype_path, str(sample).zfill(2) + '_varianthaplos.txt.gz')
                out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(2) +'.fa')
                chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
               
                write_concat_chromosome_simple(haplotype_names, chrom, fasta_file_name)
                # # tsv_chrom_to_vcf(variation_txt_path, chrom, out_vcf_path, sample_name="SAMPLE")
                fasta_from_vcf_ref(haplotype_names[0], out_vcf_path, chrom, chrom_fasta_path)


def generate_fasta_vcfs_allo():
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    subgenome_configs = [0.01, 0.05]
    # subgenomes = ['A', 'B']
    # subgenomes = {3: ['A', 'B'], 4: ['A', 'B'], 6: ['A', 'B', 'C'], 8: ['A', 'B']}
    ploidies = [3, 4, 6, 8]
    within_mutation_rates = [0.001, 0.005, 0.0075, 0.01, 0.05]
    # ploidy_folder_name = {4: '4ploid_2+2', 8: '8ploid_4+4'}
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}
    samples = range(20)
    chrom = 'Chr1'
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for sample in samples:
                    print('Ploidy:', ploidy, 'Subgenome_mr:', sgc, 'within_mr:', mr, 'sample:', sample)
                    # stop
                    haplotype_path = os.path.join(hap_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr) , str(sample).zfill(2))
                    # variation_txt_path = [os.path.join(haplotype_path, str(sample).zfill(2) + '_' + sbg + '_varianthaplos.txt.gz') for sbg in subgenomes]
                    # variation_txt_name = os.path.join(haplotype_path, str(sample).zfill(2) + '_varianthaplos.txt.gz')
                    chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                    # fasta_paths = [os.path.join(haplotype_path, str(sample).zfill(2) + '_' + sbg + '_hap' + str(i + 1) + '.fa.gz') for i in range(int(ploidy/2)) for sbg in subgenomes]
                    fasta_files = [f for f in os.listdir(haplotype_path) if 'fa.gz' in f]
                    fasta_paths = [os.path.join(haplotype_path, f) for f in fasta_files]
                    fasta_names = ['haplotype_' + f.split('.fa.gz')[0].split('_')[2].split('hap')[1] + '_subgenome_' + f.split('.fa.gz')[0].split('_')[1] for f in fasta_files]
                    
                    # fasta_names = ['haplotype_' + str(i + 1) + '_sungenome_' + sbg for i in range(int(ploidy/2)) for sbg in subgenomes[ploidy]]
                    fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(2) +'.fa')
                    out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                    # ready_vcf_file = os.path.join(haplotype_path, str(sample).zfill(2) + '_variants.vcf.gz')
                    write_concat_chromosome_simple(fasta_paths, chrom, fasta_file_name, names_list=fasta_names)
                    # tsv_chrom_to_vcf(variation_txt_name, chrom, out_vcf_path, sample_name=str(sample).zfill(2))
                    fasta_from_vcf_ref(fasta_paths[0], out_vcf_path, chrom, chrom_fasta_path)
                    # compare_two_vcfs_first_mismatch_print_reason(out_vcf_path, ready_vcf_file, variation_txt_name, sample1=str(sample).zfill(2), sample2=str(sample).zfill(2))


def index_references_auto():
    commands = ''
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    ploidies = [2, 3, 4, 6, 8]
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    chrom = 'Chr1'
    for ploidy in ploidies:
        for mr in mut_rates:
            for sample in samples:
                print('Ploidy:', ploidy, 'mr:', mr, 'sample:', sample)
                haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                chrom_fasta_path_mmi = os.path.join(haplotype_path, chrom + '.mmi')
                idx_comm = f'bwa index {chrom_fasta_path}\n'
                idx_comm2 = f'minimap2 -d {chrom_fasta_path_mmi} {chrom_fasta_path}\n'
                commands += idx_comm
                commands += idx_comm2
    # Save to shell script
    sh_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/index_ref_auto.sh'
    with open(sh_path, 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write(commands)


def index_references_allo():
    sh_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/index_ref_allo.sh'
    commands = ''
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    within_mut_rates  = [0.001, 0.005, 0.0075, 0.01, 0.05]
    subgenome_mut_rates = [0.01, 0.05]
    ploidies = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}
    samples = range(20)
    chrom = 'Chr1'
    for sgr in subgenome_mut_rates:
        for ploidy in ploidies.keys():
            for mr in within_mut_rates:
                for sample in samples:
                    print(f"Ploidy: {ploidy}, subgenome μ={sgr:g}, within μ={mr:g}, sample {sample:02d}")
                    haplotype_path = os.path.join(hap_path, f"subgenome_config_mut{sgr}", f"{ploidies[ploidy]}" ,f"within_mut_{mr}", f"{sample:02d}")                   
                    chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                    chrom_fasta_path_mmi = os.path.join(haplotype_path, chrom + '.mmi')
                    idx_comm = f'bwa index {chrom_fasta_path}\n'
                    idx_comm2 = f'minimap2 -d {chrom_fasta_path_mmi} {chrom_fasta_path}\n'
                    commands += idx_comm
                    commands += idx_comm2
    # Save to shell script
    with open(sh_path, 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write(commands)
                    

def unphase_vcf_files_auto():
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    ploidies = [2, 3, 4, 6, 8]
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    chrom = 'Chr1'
    for ploidy in ploidies:
        for mr in mut_rates:
            for sample in samples:
                print('Ploidy:', ploidy, 'mr:', mr, 'sample:', sample)
                haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                vcf_phased = os.path.join(haplotype_path, chrom + '.vcf')
                vcf_unphased = os.path.join(haplotype_path, chrom + '_unphased.vcf')
                genotype_path = os.path.join(haplotype_path, 'genotype.csv')
                
                remove_all_alt_variants(vcf_phased, sample_column='SAMPLE')
                vcf_to_genotype(vcf_phased, genotype_path, sample_name='SAMPLE', ploidy=ploidy)
                unphase_and_sort_vcf(vcf_phased, vcf_unphased, sample_column='SAMPLE')
                # incidents_all = np.sum(pd.read_csv(genotype_path).to_numpy().sum(axis=1) == ploidy)
                # incidents_zero = np.sum(pd.read_csv(genotype_path).to_numpy().sum(axis=1) == 0)
                # if incidents_all > 0: 
                #     print('====> Not heterozygous (all 1) # sites:', incidents_all, 'Ploidy:', ploidy, 'mr:', mr, 'sample:', sample)

                # if incidents_zero > 0: 
                #     print('====> Not heterozygous (all 0): # sites', incidents_zero, 'Ploidy:', ploidy, 'mr:', mr, 'sample:', sample)


def unphase_vcf_files_allo():
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    within_mut_rates  = [0.001, 0.005, 0.0075, 0.01, 0.05]
    subgenome_mut_rates = [0.01, 0.05]
    ploidies = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}


    samples = range(20)
    chrom = 'Chr1'

    for sgr in subgenome_mut_rates:
        for ploidy in ploidies.keys():
            for mr in within_mut_rates:
                for sample in samples:
                    # print('Ploidy:', ploidy, 'mr:', mr, 'sample:', sample)

                    haplotype_path = os.path.join(hap_path, f"subgenome_config_mut{sgr}", f"{ploidies[ploidy]}" ,f"within_mut_{mr}", f"{sample:02d}")                   
                    vcf_phased = os.path.join(haplotype_path, chrom + '.vcf')
                    vcf_unphased = os.path.join(haplotype_path, chrom + '_unphased.vcf')
                    genotype_path = os.path.join(haplotype_path, 'genotype.csv')
                    
                    remove_all_alt_variants(vcf_phased, sample_column=f"{sample:02d}")
                    vcf_to_genotype(vcf_phased, genotype_path, sample_name=f"{sample:02d}", ploidy=ploidy)
                    unphase_and_sort_vcf(vcf_phased, vcf_unphased, sample_column=f"{sample:02d}")
                    incidents_all = np.sum(pd.read_csv(genotype_path).to_numpy().sum(axis=1) == ploidy)
                    incidents_zero = np.sum(pd.read_csv(genotype_path).to_numpy().sum(axis=1) == 0)
                    if incidents_all > 0: 
                        print('====> Not heterozygous (all 1) # sites:', incidents_all, 'Ploidy:', ploidy, 'mr:', mr, 'sample:', sample)

                    if incidents_zero > 0: 
                        print('====> Not heterozygous (all 0): # sites', incidents_zero, 'Ploidy:', ploidy, 'mr:', mr, 'sample:', sample)


def vcf_to_genotype(vcf_path, output_path, sample_name, ploidy):
    """
    Reads a VCF file and creates a DataFrame where each column represents a haplotype.
    
    Args:
        vcf_path: Path to input VCF file
        output_path: Path to output CSV file
        sample_name: Name of the sample column to extract
        ploidy: Number of haplotypes (e.g., 2, 3, 4, 6, 8)
    """
    # Lists to store data
    haplotype_data = {f'haplotype_{i+1}': [] for i in range(ploidy)}
    
    with open(vcf_path, 'r') as infile:
        for line in infile:
            # Skip header lines
            if line.startswith('##'):
                continue
            
            # Get sample column index from header line
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                try:
                    sample_idx = headers.index(sample_name)
                except ValueError:
                    raise ValueError(f"Sample column '{sample_name}' not found in VCF header")
                continue
            
            # Process data lines
            fields = line.strip().split('\t')
            
            # Get the genotype field
            genotype = fields[sample_idx]
            
            # Split by ':' in case there are multiple FORMAT fields
            gt = genotype.split(':')[0]  # GT is always first
            
            # Split by either '|' (phased) or '/' (unphased)
            if '|' in gt:
                alleles = gt.split('|')
            elif '/' in gt:
                alleles = gt.split('/')
            else:
                raise ValueError(f"Unexpected genotype format: {gt}")
            
            # Verify ploidy matches
            if len(alleles) != ploidy:
                raise ValueError(f"Genotype {gt} has {len(alleles)} alleles, expected {ploidy}")
            
            # Add alleles to respective haplotype columns
            for i, allele in enumerate(alleles):
                haplotype_data[f'haplotype_{i+1}'].append(int(allele))
    
    # Create DataFrame
    df = pd.DataFrame(haplotype_data)
    
    # Save to CSV
    df.to_csv(output_path, index=False)


def remove_all_alt_variants(vcf_path, sample_column):
    """
    Removes variants where all alleles in the sample column are 1s (e.g., 1|1|1, 1/1).
    Overwrites the original VCF file and creates a gzipped version.
    
    Args:
        vcf_path: Path to VCF file (will be overwritten)
        sample_column: Name of the sample column to check (default: '00')
    """
    # Read all lines first
    lines_to_keep = []
    sample_idx = None
    
    with open(vcf_path, 'r') as infile:
        for line in infile:
            # Keep all header lines
            if line.startswith('##'):
                lines_to_keep.append(line)
                continue
            
            # Process column header line
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                try:
                    sample_idx = headers.index(sample_column)
                except ValueError:
                    raise ValueError(f"Sample column '{sample_column}' not found in VCF header")
                lines_to_keep.append(line)
                continue
            
            # Process data lines
            fields = line.strip().split('\t')
            genotype = fields[sample_idx]
            
            # Extract GT field (first field before ':')
            gt = genotype.split(':')[0]
            
            # Split by either '|' or '/'
            if '|' in gt:
                alleles = gt.split('|')
            elif '/' in gt:
                alleles = gt.split('/')
            else:
                # If no separator, keep the line (unexpected format)
                lines_to_keep.append(line)
                continue
            
            # Check if all alleles are '1'
            if all(allele == '1' for allele in alleles):
                # Skip this line (don't add to lines_to_keep)
                continue
            else:
                # Keep this line
                lines_to_keep.append(line)
    
    # Write the regular VCF file
    with open(vcf_path, 'w') as outfile:
        outfile.writelines(lines_to_keep)
    
    # Create gzipped version
    gz_path = vcf_path + '.gz'
    with open(vcf_path, 'rb') as f_in:
        with gzip.open(gz_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    # print(f"Created: {vcf_path}")
    # print(f"Created: {gz_path}")


def test_extractHAIRS_polyploid():
    """
    Test extractHAIRS on a small polyploid sample to verify it works correctly
    """
    # Test parameters - use smallest dataset
    test_ploidy = 3  # Start with triploid
    test_mr = 0.001
    test_coverage = 5
    test_sample = 0
    
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_auto_short'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    
    # Paths
    haplotype_path = os.path.join(hap_path, f'ploidy_{test_ploidy}', f'mut_{test_mr}', str(test_sample).zfill(2))
    vcf_file = os.path.join(haplotype_path, 'Chr1.vcf')
    bam_file = os.path.join(dataset_path, f'ploidy_{test_ploidy}', f'mut_{test_mr}', f'cov_{test_coverage}', f'{str(test_sample).zfill(2)}.bam')
    frag_file = os.path.join(dataset_path, f'ploidy_{test_ploidy}', f'mut_{test_mr}', f'cov_{test_coverage}', f'{str(test_sample).zfill(2)}_TEST.frag')
    
    # Run extractHAIRS
    cmd = f'{extract_location} --bam {bam_file} --VCF {vcf_file} --out {frag_file}'
    print(f"Testing command: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        print(f"Return code: {result.returncode}")
        
        # Check if fragment file was created and is not empty
        if os.path.exists(frag_file):
            file_size = os.path.getsize(frag_file)
            print(f"\n✓ Fragment file created: {frag_file}")
            print(f"  File size: {file_size} bytes")
            
            if file_size > 0:
                # Read and display first few lines
                with open(frag_file, 'r') as f:
                    lines = f.readlines()[:10]
                    print(f"\n  First {len(lines)} lines of fragment file:")
                    for line in lines:
                        print(f"    {line.rstrip()}")
                
                # Check format
                print("\n✓ Checking fragment format...")
                check_fragment_format(frag_file, test_ploidy)
            else:
                print("✗ Fragment file is empty!")
        else:
            print(f"✗ Fragment file not created: {frag_file}")
            
    except subprocess.TimeoutExpired:
        print("✗ Command timed out (>5 minutes)")
    except Exception as e:
        print(f"✗ Error: {e}")


def check_fragment_format(frag_file, expected_ploidy):
    """
    Check if fragment file format supports polyploid data
    Fragment format should have variant positions and allele calls
    """
    with open(frag_file, 'r') as f:
        for i, line in enumerate(f):
            # if i > 5:  # Check first few fragments
            #     break
            parts = line.strip().split()
            if len(parts) > 0:
                print(f"    Fragment {i}: {len(parts)} fields")
                # Fragment format: each line has info about read spanning variants
                # Should show allele calls (0 or 1) for each variant position
                
    print(f"  Expected ploidy: {expected_ploidy}")
    print("  ✓ Check manually: Do fragments show correct allele information?")


def simulate_long_read_dataset_auto():
    run_all_sh = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_auto/simulate_auto_long.sh'
    pbsim_location = '/home/mah19006/pHapCompass_v2/pbsim3/src/pbsim'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    model_path = '/home/mah19006/pHapCompass_v2/pbsim3/data/QSHMM-ONT-HQ.model' # QSHMM-ONT-HQ.model
    # model_path = '/home/mah19006/pHapCompass_v2/pbsim3/data/QSHMM-RSII.model'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_auto_long'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    commands_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_auto'
    if not os.path.exists(dataset_path):
        os.makedirs(dataset_path)

    if not os.path.exists(commands_path):
        os.makedirs(commands_path)
    
    coverages = [3, 5, 10, 20, 40]
    ploidies = [2, 3, 4, 6, 8]
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    chrom = 'Chr1'
    all_sh = ''
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                cov_path = os.path.join(dataset_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage))
                if not os.path.exists(cov_path):
                    os.makedirs(cov_path, exist_ok=True)
                # stop
                this_sample = ''
                sh_file_name = 'ploidy_' + str(ploidy) + '_mut_' + str(mr) + '_cov_' + str(coverage) + '.sh'
                for sample in samples:
                    # this_sample = ''
                    print('Ploidy:', ploidy, '- mr:', mr, '- Coverage:', coverage, '- sample:', sample)     
                    haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                    out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                    fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(2) +'.fa')
                    # chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                    mmi_path = os.path.join(haplotype_path, chrom + '.mmi')
                    comm0 = f'# Ploidy: {ploidy} - mutation rate: {mr} - Coverage: {coverage} - sample: {sample}\n'
                    comm1 = 'cd {}\n'.format(cov_path)
                    # pbsim_command = '{} --strategy wgs --method qshmm --qshmm {} --depth {} --genome {}  --accuracy-mean 0.999 0.0005 --prefix {}\n'.format(pbsim_location, model_path, coverage/ploidy, fasta_file_name, str(sample).zfill(2))
                    pbsim_command = '{} --strategy wgs --method qshmm --qshmm {} --depth {} --genome {} --prefix {}\n'.format(pbsim_location, model_path, coverage/ploidy, fasta_file_name, str(sample).zfill(2))
                    align_command = f'minimap2 -t 16 -x map-ont -a {mmi_path} {cov_path}/{str(sample).zfill(2)}_00*.fq.gz | samtools sort -@8 -o {cov_path}/{str(sample).zfill(2)}.bam\n'
                    index_command = f'samtools index {cov_path}/{str(sample).zfill(2)}.bam\n'
                    ext_h_command = f'{extract_location} --bam {cov_path}/{str(sample).zfill(2)}.bam --VCF {out_vcf_path} --out {cov_path}/{str(sample).zfill(2)}.frag\n\n\n'
                    this_sample += comm0 + comm1 + pbsim_command + align_command + index_command + ext_h_command
                with open(os.path.join(commands_path, sh_file_name), "a") as f:
                    f.write(this_sample)
                all_sh += 'bash ' + os.path.join(commands_path, sh_file_name) + '\n'
    with open(run_all_sh, "w") as f:
        f.write(all_sh)


def simulate_long_read_dataset_allo():
    run_all_sh = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_allo2/simulate_allo.sh'
    pbsim_location = '/home/mah19006/pHapCompass_v2/pbsim3/src/pbsim'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    model_path = '/home/mah19006/pHapCompass_v2/pbsim3/data/QSHMM-ONT-HQ.model' # QSHMM-ONT-HQ.model
    # model_path = '/home/mah19006/pHapCompass_v2/pbsim3/data/QSHMM-RSII.model'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_allo_long'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    commands_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_allo2' 
    
    if not os.path.exists(dataset_path):
        os.makedirs(dataset_path)

    if not os.path.exists(commands_path):
        os.makedirs(commands_path)
    
    subgenome_configs = [0.01, 0.05]
    ploidies = [3, 4, 6, 8]
    within_mutation_rates = [0.001, 0.005, 0.0075, 0.01, 0.05]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}
    samples = range(20)
    coverages = [5, 10, 20, 40]
    chrom = 'Chr1'
    all_sh = ''
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for coverage in coverages:
                    cov_path = os.path.join(dataset_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage))
                    if not os.path.exists(cov_path):
                        os.makedirs(cov_path, exist_ok=True)
                    # stop
                    done_sampels = [f for f in os.listdir(cov_path) if '.frag' in f]
                    if len(done_sampels) < 20:
                        this_sample = ''
                        sh_file_name = 'ploidy_' + str(ploidy) + '_subgenome_' + str(sgc) + '_mut_' + str(mr) + '_cov_' + str(coverage) + '.sh'
                        for sample in samples:

                            print('Ploidy:', ploidy, '- subgenome:', sgc, '- mr:', mr, '- Coverage:', coverage, '- sample:', sample)     
                            # if not os.path.exists(os.path.join(cov_path, str(sample).zfill(2) + '.frag')):
                            #     print('not done.')
                                # stop
                            # haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                            haplotype_path = os.path.join(hap_path, f"subgenome_config_mut{sgc}", f"{ploidy_folder_name[ploidy]}" ,f"within_mut_{mr}", f"{sample:02d}")
                            out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                            fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(2) +'.fa')
                            # chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                            mmi_path = os.path.join(haplotype_path, chrom + '.mmi')
                            comm0 = f'# Ploidy: {ploidy} - mutation rate: {mr} - Coverage: {coverage} - sample: {sample}\n'
                            comm1 = 'cd {}\n'.format(cov_path)
                            # pbsim_command = '{} --strategy wgs --method qshmm --qshmm {} --depth {} --genome {} --accuracy-mean 0.999 0.0005 --prefix {}\n'.format(pbsim_location, model_path, coverage/ploidy, fasta_file_name, str(sample).zfill(2))
                            pbsim_command = '{} --strategy wgs --method qshmm --qshmm {} --depth {} --genome {} --prefix {}\n'.format(pbsim_location, model_path, coverage/ploidy, fasta_file_name, str(sample).zfill(2))
                            align_command = f'minimap2 -t 16 -x map-ont -a {mmi_path} {cov_path}/{str(sample).zfill(2)}_00*.fq.gz | samtools sort -@8 -o {cov_path}/{str(sample).zfill(2)}.bam\n'
                            index_command = f'samtools index {cov_path}/{str(sample).zfill(2)}.bam\n'
                            ext_h_command = f'{extract_location} --bam {cov_path}/{str(sample).zfill(2)}.bam --VCF {out_vcf_path} --out {cov_path}/{str(sample).zfill(2)}.frag\n\n\n'
                            this_sample += comm0 + comm1 + pbsim_command + align_command + index_command + ext_h_command
                        with open(os.path.join(commands_path, sh_file_name), "a") as f:
                            f.write(this_sample)
                        all_sh += 'bash ' + os.path.join(commands_path, sh_file_name) + '\n'
    with open(run_all_sh, "w") as f:
        f.write(all_sh)


def simulate_short_read_dataset_auto():
    run_all_sh = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/short/simulate_auto_short.sh'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_auto_short'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    commands_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/short' 
    art_path = 'art_illumina'
    coverages = [3, 5, 10, 20, 40]
    ploidies = [2, 3, 4, 6, 8]
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    chrom = 'Chr1'
    read_length = 125
    mean_insert_length = 350
    std_insert_length = 50
    all_sh = ''
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                # stop
                cov_path = os.path.join(dataset_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage))
                if not os.path.exists(cov_path):
                    os.makedirs(cov_path, exist_ok=True)
                # stop
                this_sample = ''
                sh_file_name = 'ploidy_' + str(ploidy) + '_mut_' + str(mr) + '_cov_' + str(coverage) + '_short.sh'
                for sample in samples:
                    # this_sample = ''
                    rs = random.randint(1, 2**32)
                    print('Ploidy:', ploidy, '- mr:', mr, '- Coverage:', coverage, '- sample:', sample)     
                    haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                    out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                    fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(2) +'.fa')
                    chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                    # mmi_path = os.path.join(haplotype_path, chrom + '.mmi')
                    fastq_file_name = os.path.join(cov_path, str(sample).zfill(2))
                    comm0 = f'# Ploidy: {ploidy} - mutation rate: {mr} - Coverage: {coverage} - sample: {sample}\n'
                    art_command = f'{art_path} -ss HS25 -rs {rs} -i {fasta_file_name} -p -na -l {read_length} -f {coverage/ploidy} -m {mean_insert_length} -s {std_insert_length} -o {fastq_file_name};'
                    align_command = f'bwa mem {chrom_fasta_path} {fastq_file_name}1.fq {fastq_file_name}2.fq > {fastq_file_name}.sam; samtools view -Sb {fastq_file_name}.sam | samtools sort -o {fastq_file_name}.bam; samtools index {fastq_file_name}.bam; rm {fastq_file_name}.sam;'
                    ext_h_command = f'{extract_location} --bam {fastq_file_name}.bam --VCF {out_vcf_path} --out {fastq_file_name}.frag\n\n\n'
                    this_sample += comm0 + art_command + align_command + ext_h_command
                    # check_extractHAIRS_systematic(f'{fastq_file_name}.bam', out_vcf_path, f'{fastq_file_name}.frag')
                with open(os.path.join(commands_path, sh_file_name), "a") as f:
                    f.write(this_sample)
                all_sh += 'bash ' + os.path.join(commands_path, sh_file_name) + '\n'
    with open(run_all_sh, "a") as f:
        f.write(all_sh)


def simulate_short_read_dataset_allo():

    run_all_sh = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/short_allo2/simulate_short_allo.sh'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_allo_short'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    commands_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/short_allo2' 
    art_path = 'art_illumina'
    
    if not os.path.exists(dataset_path):
        os.makedirs(dataset_path)

    if not os.path.exists(commands_path):
        os.makedirs(commands_path)
    
    subgenome_configs = [0.01, 0.05]
    ploidies = [3, 4, 6, 8]
    within_mutation_rates = [0.001, 0.005, 0.0075, 0.01, 0.05]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}
    samples = range(20)
    coverages = [3, 5, 10, 20, 40]
    read_length = 125
    mean_insert_length = 350
    std_insert_length = 50
    chrom = 'Chr1'
    all_sh = ''
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for coverage in coverages:
                    cov_path = os.path.join(dataset_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage))
                    if not os.path.exists(cov_path):
                        os.makedirs(cov_path, exist_ok=True)
                    # stop
                    done_sampels = [f for f in os.listdir(cov_path) if '.frag' in f]
                    if len(done_sampels) < 20:
                        this_sample = ''
                        sh_file_name = 'ploidy_' + str(ploidy) + '_subgenome_' + str(sgc) + '_mut_' + str(mr) + '_cov_' + str(coverage) + '.sh'
                        for sample in samples:

                            print('Ploidy:', ploidy, '- subgenome:', sgc, '- mr:', mr, '- Coverage:', coverage, '- sample:', sample)     
                            
                            rs = random.randint(1, 2**32)
                            haplotype_path = os.path.join(hap_path, f"subgenome_config_mut{sgc}", f"{ploidy_folder_name[ploidy]}" ,f"within_mut_{mr}", f"{sample:02d}")
                            out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                            fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(2) +'.fa')
                            chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                            mmi_path = os.path.join(haplotype_path, chrom + '.mmi')
                            fastq_file_name = os.path.join(cov_path, str(sample).zfill(2))
                            comm0 = f'# Ploidy: {ploidy} - subgenome mutation rate: {sgc} - within mutation rate: {mr} - Coverage: {coverage} - sample: {sample}\n'
                            art_command = f'{art_path} -ss HS25 -rs {rs} -i {fasta_file_name} -p -na -l {read_length} -f {coverage/ploidy} -m {mean_insert_length} -s {std_insert_length} -o {fastq_file_name};'
                            align_command = f'bwa mem {chrom_fasta_path} {fastq_file_name}1.fq {fastq_file_name}2.fq > {fastq_file_name}.sam; samtools view -Sb {fastq_file_name}.sam | samtools sort -o {fastq_file_name}.bam; samtools index {fastq_file_name}.bam; rm {fastq_file_name}.sam;'
                            ext_h_command = f'{extract_location} --bam {fastq_file_name}.bam --VCF {out_vcf_path} --out {fastq_file_name}.frag\n\n\n'
                            this_sample += comm0 + art_command + align_command + ext_h_command
                        with open(os.path.join(commands_path, sh_file_name), "a") as f:
                            f.write(this_sample)
                        all_sh += 'bash ' + os.path.join(commands_path, sh_file_name) + '\n'
    with open(run_all_sh, "w") as f:
        f.write(all_sh)


def simulate_auto_short_read_validation():
    run_all_sh = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/grid/simulate_auto_short.sh'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/evaluation'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    commands_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/grid' 
    art_path = 'art_illumina'
    coverages = [5, 10, 20]
    ploidies = [2, 3, 4, 6]
    mut_rates = [0.001]
    samples = range(5)
    chrom = 'Chr1'
    read_length = 125
    mean_insert_length = 350
    std_insert_length = 50
    all_sh = ''
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                # stop
                cov_path = os.path.join(dataset_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage))
                if not os.path.exists(cov_path):
                    os.makedirs(cov_path, exist_ok=True)
                # stop
                this_sample = ''
                sh_file_name = 'ploidy_' + str(ploidy) + '_mut_' + str(mr) + '_cov_' + str(coverage) + '_short.sh'
                for sample in samples:
                    # this_sample = ''
                    rs = random.randint(1, 2**32)
                    print('Ploidy:', ploidy, '- mr:', mr, '- Coverage:', coverage, '- sample:', sample)     
                    haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                    out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                    fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(2) +'.fa')
                    chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                    # mmi_path = os.path.join(haplotype_path, chrom + '.mmi')
                    fastq_file_name = os.path.join(cov_path, str(sample).zfill(2))
                    comm0 = f'# Ploidy: {ploidy} - mutation rate: {mr} - Coverage: {coverage} - sample: {sample}\n'
                    art_command = f'{art_path} -ss HS25 -rs {rs} -i {fasta_file_name} -p -na -l {read_length} -f {coverage/ploidy} -m {mean_insert_length} -s {std_insert_length} -o {fastq_file_name};'
                    align_command = f'bwa mem {chrom_fasta_path} {fastq_file_name}1.fq {fastq_file_name}2.fq > {fastq_file_name}.sam; samtools view -Sb {fastq_file_name}.sam | samtools sort -o {fastq_file_name}.bam; samtools index {fastq_file_name}.bam; rm {fastq_file_name}.sam;'
                    ext_h_command = f'{extract_location} --bam {fastq_file_name}.bam --VCF {out_vcf_path} --out {fastq_file_name}.frag\n\n\n'
                    this_sample += comm0 + art_command + align_command + ext_h_command
                    # check_extractHAIRS_systematic(f'{fastq_file_name}.bam', out_vcf_path, f'{fastq_file_name}.frag')
                with open(os.path.join(commands_path, sh_file_name), "a") as f:
                    f.write(this_sample)
                all_sh += 'bash ' + os.path.join(commands_path, sh_file_name) + '\n'
    with open(run_all_sh, "a") as f:
        f.write(all_sh)


def check_extractHAIRS_systematic(bam_file, vcf_file, frag_file):
    """
    Systematically verify extractHAIRS is working correctly for polyploid data
    """
    print("=" * 80)
    print("SYSTEMATIC extractHAIRS VERIFICATION")
    print("=" * 80)
    
    # 1. Check fragment file statistics
    print("\n1. FRAGMENT FILE STATISTICS")
    print("-" * 40)
    
    with open(frag_file, 'r') as f:
        lines = f.readlines()
    
    num_fragments = len(lines)
    print(f"Total fragments: {num_fragments}")
    
    # Count variants covered
    variant_ids = set()
    block_counts = {}
    allele_lengths = []
    
    for line in lines:
        parts = line.strip().split()
        num_blocks = int(parts[0])
        block_counts[num_blocks] = block_counts.get(num_blocks, 0) + 1
        
        # Extract variant IDs (every other field starting from index 2)
        for i in range(2, len(parts)):
            if parts[i].isdigit() and len(parts[i]) >= 4:
                variant_ids.add(int(parts[i]))
            elif i > 2 and parts[i-1].isdigit() and len(parts[i-1]) >= 4:
                # This is allele string after variant ID
                allele_lengths.append(len(parts[i]))
    
    print(f"Unique variants covered: {len(variant_ids)}")
    print(f"Variant ID range: {min(variant_ids)} - {max(variant_ids)}")
    print(f"\nBlock distribution:")
    for blocks, count in sorted(block_counts.items()):
        print(f"  {blocks} block(s): {count} fragments ({100*count/num_fragments:.1f}%)")
    
    if allele_lengths:
        from collections import Counter
        allele_dist = Counter(allele_lengths)
        print(f"\nAlleles per fragment (SNPs covered):")
        for length, count in sorted(allele_dist.items())[:10]:
            print(f"  {length} SNPs: {count} fragments ({100*count/num_fragments:.1f}%)")
    
    # 2. Check VCF statistics
    print("\n2. VCF FILE STATISTICS")
    print("-" * 40)
    
    vcf_variants = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            var_id = int(parts[2])
            gt = parts[9].split(':')[0]
            vcf_variants.append((var_id, gt))
    
    print(f"Total variants in VCF: {len(vcf_variants)}")
    
    # Count heterozygous variants
    het_variants = [v for v in vcf_variants if '0' in v[1] and '1' in v[1]]
    hom_ref = [v for v in vcf_variants if v[1].replace('|', '').replace('/', '') == '0' * v[1].count('0')]
    hom_alt = [v for v in vcf_variants if v[1].replace('|', '').replace('/', '') == '1' * v[1].count('1')]
    
    print(f"Heterozygous variants: {len(het_variants)} ({100*len(het_variants)/len(vcf_variants):.1f}%)")
    print(f"Homozygous REF (all 0): {len(hom_ref)} ({100*len(hom_ref)/len(vcf_variants):.1f}%)")
    print(f"Homozygous ALT (all 1): {len(hom_alt)} ({100*len(hom_alt)/len(vcf_variants):.1f}%)")
    
    # 3. Coverage of heterozygous sites
    print("\n3. HETEROZYGOUS VARIANT COVERAGE")
    print("-" * 40)
    
    het_var_ids = set(v[0] for v in het_variants)
    covered_het = variant_ids.intersection(het_var_ids)
    
    print(f"Heterozygous variants in fragments: {len(covered_het)} / {len(het_var_ids)} ({100*len(covered_het)/len(het_var_ids):.1f}%)")
    
    # Find uncovered heterozygous variants
    uncovered = het_var_ids - variant_ids
    if uncovered and len(uncovered) <= 20:
        print(f"\nUncovered heterozygous variants (first 20): {sorted(list(uncovered))[:20]}")
    
    # 4. Manual verification of specific fragments
    print("\n4. MANUAL VERIFICATION - Sample Fragments")
    print("-" * 40)
    
    # Pick a few fragments to manually verify
    sample_fragments = lines[:5]
    
    for i, frag_line in enumerate(sample_fragments):
        print(f"\nFragment {i+1}:")
        parts = frag_line.strip().split()
        read_name = parts[1]
        print(f"  Read: {read_name}")
        
        # Extract variant info
        variant_info = []
        j = 2
        while j < len(parts):
            if parts[j].isdigit() and len(parts[j]) >= 4:
                var_id = parts[j]
                alleles = parts[j+1] if j+1 < len(parts) else "?"
                variant_info.append((var_id, alleles))
                j += 2
            else:
                j += 1
        
        print(f"  Variants: {variant_info}")
        
        # Get the read from BAM
        cmd = f"samtools view {bam_file} | grep '{read_name}' | head -2"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.stdout:
            bam_lines = result.stdout.strip().split('\n')
            print(f"  BAM entries: {len(bam_lines)}")
            for bam_line in bam_lines:
                bam_parts = bam_line.split('\t')
                pos = bam_parts[3]
                flag = int(bam_parts[1])
                is_reverse = (flag & 16) != 0
                is_first = (flag & 64) != 0
                mate_info = "R1" if is_first else "R2"
                strand = "rev" if is_reverse else "fwd"
                print(f"    {mate_info} at pos {pos} ({strand})")
            
            # Check which variants this read should cover
            print(f"  Expected to cover variants: {variant_info}")
    
    # 5. Check for systematic issues
    print("\n5. SYSTEMATIC CHECKS")
    print("-" * 40)
    
    # Check if fragments only use certain alleles
    all_alleles = []
    for line in lines[:1000]:  # Sample first 1000
        parts = line.strip().split()
        for i in range(2, len(parts)):
            if not parts[i].isdigit() and len(parts[i]) <= 10 and all(c in '01' for c in parts[i]):
                all_alleles.extend(list(parts[i]))
    
    if all_alleles:
        from collections import Counter
        allele_counts = Counter(all_alleles)
        print(f"Allele distribution in fragments:")
        print(f"  0 (REF): {allele_counts.get('0', 0)} ({100*allele_counts.get('0', 0)/len(all_alleles):.1f}%)")
        print(f"  1 (ALT): {allele_counts.get('1', 0)} ({100*allele_counts.get('1', 0)/len(all_alleles):.1f}%)")
        
        # For ploidy 3 with GT 0|0|1, we expect ~2/3 REF (0) and ~1/3 ALT (1)
        # For GT 0|1|0, same thing
        # For GT 1|0|0, same thing
        print(f"\n  Expected for ploidy 3: ~67% REF, ~33% ALT (for heterozygous 0|0|1 sites)")
    
    print("\n" + "=" * 80)
    print("VERIFICATION COMPLETE")
    print("=" * 80)


def snp_count(vcf_path: str) -> int:
    """Return the number of variant lines in a biallelic, single-sample VCF."""
    n = 0
    op = gzip.open if vcf_path.endswith(".gz") else open
    with op(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): 
                continue
            n += 1
    return n


def hamming_on_snps(fasta_path: str, vcf_path: str):
    """
    Compute pairwise Hamming distances (fraction of SNP mismatches)
    among biallelic haplotypes. Assumes all sequences are equal length
    and only SNP positions differ (A/B alleles).
    Returns a pandas DataFrame (values in [0,1]).
    """
  

    # --- SNP positions (1-based -> 0-based) ---
    positions = []
    _open_v = gzip.open if vcf_path.endswith(".gz") else open
    with _open_v(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                positions.append(int(line.split("\t", 3)[1]) - 1)
    if not positions:
        raise ValueError(f"No SNPs found in {vcf_path}")

    idx = np.array(positions, dtype=np.int64)
    m = len(idx)  # number of SNPs

    # --- haplotypes ---
    labels, seqs = [], []
    _open_f = gzip.open if fasta_path.endswith(".gz") else open
    with _open_f(fasta_path, "rt") as fh:
        seq, label = [], None
        for line in fh:
            if line.startswith(">"):
                if label is not None:
                    seqs.append("".join(seq).strip())
                    seq = []
                label = line[1:].strip().split()[0]
                labels.append(label or f"hap_{len(labels)+1}")
            else:
                seq.append(line.strip())
        if label is not None:
            seqs.append("".join(seq).strip())

    n = len(seqs)
    if n < 2:
        raise ValueError("Need at least two haplotypes to compare")

    # --- subset to SNPs ---
    arr = np.vstack([np.frombuffer(s.encode("ascii"), dtype="S1")[idx] for s in seqs])

    # --- compute fraction of mismatches explicitly ---
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            mism = np.count_nonzero(arr[i] != arr[j])
            D[i, j] = D[j, i] = mism / m

    return pd.DataFrame(D, index=labels, columns=labels)


def intersnp_distances(vcf_path: str):
    """
    Return the list of distances between consecutive SNPs (in bp), sorted by position.
    For 0 or 1 SNP, returns [].
    """
    pos = []
    op = gzip.open if vcf_path.endswith(".gz") else open
    with op(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): 
                continue
            pos.append(int(line.split("\t", 3)[1]))
    pos.sort()
    if len(pos) < 2:
        return []
    diffs = [pos[i] - pos[i-1] for i in range(1, len(pos))]
    return diffs


def diversity_stats(vcf_path: str, fasta_path: str):
    """
    # Diversity stats from VCF (π, Watterson’s θ, heterozygosity) ------------
    Compute per-sample nucleotide diversity π and Watterson’s θ using the phased GT field
    across a biallelic VCF, plus basic counts.
    Interpretation note: here 'n' is the number of haplotypes (ploidy).
    Returns dict with: S, n_hap, seq_len, pi, theta_w, het_mean
    """
    # sequence length from FASTA (assumes all haplotypes same length)
    L = 0
    opf = gzip.open if fasta_path.endswith(".gz") else open
    with opf(fasta_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"): 
                continue
            L += len(line.strip())

    S = 0
    n_hap = None
    het_vals = []   # per-site heterozygosity = 2p(1-p) for biallelic
    opv = gzip.open if vcf_path.endswith(".gz") else open
    with opv(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): 
                continue
            S += 1
            f = line.rstrip("\n").split("\t")
            gt = f[9]  # single-sample VCF
            # FORMAT is GT only per your example; handle "0|1|0..."
            alleles = gt.split(":")[0].replace("/", "|").split("|")
            if n_hap is None:
                n_hap = len(alleles)
            alt_count = sum(1 for a in alleles if a == "1")
            p = alt_count / n_hap
            het_vals.append(2 * p * (1 - p))  # biallelic expected heterozygosity == per-site π

    if n_hap is None:
        # empty VCF; define safe outputs
        return {"S": 0, "n_hap": 0, "seq_len": L, "pi": 0.0, "theta_w": 0.0, "het_mean": 0.0}

    # π as mean per-site pairwise differences over variable sites, then scaled by L if desired.
    # Here we report π per site over the *whole* sequence length (common in genomics).
    pi_per_variable_site = float(np.mean(het_vals)) if het_vals else 0.0
    pi = pi_per_variable_site * (S / max(L, 1))  # scale by fraction of variable sites over length

    # Watterson’s θ = S / a_n, where a_n = sum_{k=1}^{n-1} 1/k, with n = number of sequences (haplotypes)
    a_n = sum(1.0 / k for k in range(1, max(n_hap, 1)))
    theta_w = (S / a_n) / max(L, 1) if a_n > 0 else 0.0

    return {
        "S": S,
        "n_hap": n_hap,
        "seq_len": L,
        "pi": pi,
        "theta_w": theta_w,
        "het_mean": float(np.mean(het_vals)) if het_vals else 0.0
    }


def shared_segment_lengths(fasta_path: str):
    """
    # Haplotype “block” / shared identical segment lengths (pairwise) --------
    For every haplotype pair, compute run-lengths of *identical* bases (contiguous segments with no differences).
    Returns dict: {'pairs': {(i,j): [run_len1, run_len2, ...]}, 'pooled': [all run lengths]}
    """
    # read haplotypes
    seqs = []
    opf = gzip.open if fasta_path.endswith(".gz") else open
    with opf(fasta_path, "rt") as fh:
        cur = []
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).strip())
                    cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur).strip())

    if len(seqs) < 2:
        return {"pairs": {}, "pooled": []}

    L = len(seqs[0])
    arr = np.vstack([np.frombuffer(s.encode("ascii"), dtype="S1") for s in seqs])

    out = {}
    pooled = []
    n = len(seqs)
    for i in range(n):
        for j in range(i+1, n):
            eq = (arr[i] == arr[j]).astype(np.int8)  # 1 for equal, 0 for diff
            runs = []
            run = 0
            for v in eq:
                if v == 1:
                    run += 1
                else:
                    if run > 0:
                        runs.append(run)
                        run = 0
            if run > 0:
                runs.append(run)
            out[(i, j)] = runs
            pooled.extend(runs)
    return {"pairs": out, "pooled": pooled}


def allele_count_distribution(vcf_path: str):
    """
    # Allele count / frequency distribution from VCF -------------------------
    For each SNP, count #ALT haplotypes from GT like '0|1|...'.
    Returns dict with:
      - 'ploidy'
      - 'counts': list of alt-count per site
      - 'hist': {k: num_sites_with_k_alts} for k=0..ploidy
    """
    counts = []
    ploidy = None
    op = gzip.open if vcf_path.endswith(".gz") else open
    with op(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"): 
                continue
            f = line.rstrip("\n").split("\t")
            gt = f[9].split(":")[0]
            alleles = gt.replace("/", "|").split("|")
            if ploidy is None:
                ploidy = len(alleles)
            counts.append(sum(1 for a in alleles if a == "1"))
    if ploidy is None:
        ploidy = 0
    hist = {k: 0 for k in range(ploidy + 1)}
    for c in counts:
        if c in hist:
            hist[c] += 1
    return {"ploidy": ploidy, "counts": counts, "hist": hist}


def distance_matrix_and_mds(fasta_path: str, use_snps_only: bool = False, vcf_path: str | None = None, n_components: int = 2):
    """
    # Distance matrix + classical MDS (no sklearn) ---------------------------
    Build normalized Hamming distance matrix among haplotypes and run classical MDS
    to get low-dim embedding. Returns dict with 'D' (matrix) and 'coords' (N x n_components).
    """
    # read haplotypes
    seqs = []
    opf = gzip.open if fasta_path.endswith(".gz") else open
    with opf(fasta_path, "rt") as fh:
        cur = []
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur).strip())
                    cur = []
            else:
                cur.append(line.strip())
        if cur:
            seqs.append("".join(cur).strip())
    n = len(seqs)
    if n == 0:
        return {"D": [], "coords": []}
    if n == 1:
        return {"D": [[0.0]], "coords": [[0.0]*n_components]}

    L = len(seqs[0])
    arr = np.vstack([np.frombuffer(s.encode("ascii"), dtype="S1") for s in seqs])

    # optional SNP-only restriction
    if use_snps_only:
        if vcf_path is None:
            raise ValueError("use_snps_only=True requires vcf_path")
        pos = []
        opv = gzip.open if vcf_path.endswith(".gz") else open
        with opv(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"): 
                    continue
                pos.append(int(line.split("\t", 3)[1]) - 1)
        if not pos:
            D = np.zeros((n, n), dtype=float)
            return {"D": D.tolist(), "coords": np.zeros((n, n_components)).tolist()}
        arr = arr[:, np.array(pos, dtype=np.int64)]
        eff_len = arr.shape[1]
    else:
        eff_len = L

    # distance matrix
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i+1, n):
            mism = np.sum(arr[i] != arr[j])
            d = mism / eff_len
            D[i, j] = D[j, i] = d

    # classical MDS via double-centering
    J = np.eye(n) - np.ones((n, n))/n
    B = -0.5 * J @ (D**2) @ J
    # eigen-decomposition
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    # keep top components with non-negative eigenvalues
    keep = min(n_components, np.sum(eigvals > 0))
    if keep == 0:
        coords = np.zeros((n, n_components))
    else:
        Lmb = np.diag(np.sqrt(eigvals[:keep]))
        coords_k = eigvecs[:, :keep] @ Lmb
        # pad if keep < n_components
        if keep < n_components:
            coords = np.zeros((n, n_components))
            coords[:, :keep] = coords_k
        else:
            coords = coords_k
    return {"D": D.tolist(), "coords": coords.tolist()}


def compute_reference_stats_auto():
    reference_path =  '/home/mah19006/projects/HaplOrbit/reference/simulated_haplotypes/auto'
    save_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/info_auto'
    ploidies = [2, 3, 4, 6, 8]
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    snp_count_df = pd.DataFrame(columns=['ploidy', 'mutation_rate', 'sample', 'snp_count'], index=range(len(ploidies)* len(mut_rates)* len(samples)))
    counter = 0
    hamming_info = {p: {mr: {s: None for s in samples} for mr in mut_rates} for p in ploidies}
    snp_density_info  = {p: {mr: {s: None for s in samples} for mr in mut_rates} for p in ploidies}
    for ploidy in ploidies:
        for mr in mut_rates:
            for sample in samples:
                this_ref_path = os.path.join(reference_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr) , str(sample).zfill(2))
                vcf_file_path = os.path.join(this_ref_path, 'Chr1.vcf')
                fasta_file_path = os.path.join(this_ref_path, 'Chr1_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(2) + '.fa')
                
                snp_count_df.loc[counter, :] = ploidy, mr, sample, snp_count(vcf_file_path)
                counter += 1
                hamming_info[ploidy][mr][sample] = hamming_on_snps(fasta_file_path, vcf_file_path)
                snp_density_info[ploidy][mr][sample] = intersnp_distances(vcf_file_path)

    # save snp count
    snp_count_df.to_csv(os.path.join(save_path, 'snp_counts.csv'))
    # save hamming_info
    pd.to_pickle(hamming_info, os.path.join(save_path, "hamming_info.pkl"))
    # save snp_density_info
    pd.to_pickle(snp_density_info, os.path.join(save_path, "snp_density_info.pkl"))


def compute_reference_stats_allo():
    reference_path =  '/home/mah19006/projects/HaplOrbit/reference/simulated_haplotypes/allo'
    save_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/info_allo'
    # subgenome_configs = [] # [0.01, 0.05]
    subgenome_configs = ['0.0005', '0.0001']
    within_mutation_rates = ['0.00005', '0.0001']
    
    ploidies = [3, 4, 6] #, 8]
    # within_mutation_rates = [0.001, 0.005, 0.0075, 0.01, 0.05]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}
    samples = range(10)
    snp_count_df = pd.DataFrame(columns=['ploidy', 'subgenome_configs', 'within_mutation_rates','sample', 'snp_count'], 
    index=range(len(ploidies) * len(subgenome_configs) * len(within_mutation_rates) * len(samples)))
    counter = 0
    # hamming_info = {p: {mr: {s: None for s in samples} for mr in mut_rates} for p in ploidies}
    # snp_density_info  = {p: {mr: {s: None for s in samples} for mr in mut_rates} for p in ploidies}


    hamming_info = {
        sgc: {p: {mr: {s: None for s in samples} for mr in within_mutation_rates} for p in ploidies}
        for sgc in subgenome_configs
    }

    snp_density_info = {
        sgc: {p: {mr: {s: None for s in samples} for mr in within_mutation_rates} for p in ploidies}
        for sgc in subgenome_configs
    }

    for sgc in subgenome_configs:
        for ploidy in ploidies:
            for mr in within_mutation_rates:
                for sample in samples:
                    
                    this_ref_path = os.path.join(reference_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr) , str(sample).zfill(1))
                    vcf_file_path = os.path.join(this_ref_path, 'Chr1.vcf')
                    fasta_file_path = os.path.join(this_ref_path, 'Chr1_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(1) + '.fa')
                    
                    snp_count_df.loc[counter, :] = ploidy, sgc, mr, sample, snp_count(vcf_file_path)
                    counter += 1
                    hamming_info[sgc][ploidy][mr][sample] = hamming_on_snps(fasta_file_path, vcf_file_path)
                    snp_density_info[sgc][ploidy][mr][sample] = intersnp_distances(vcf_file_path)

    # save snp count
    snp_count_df.to_csv(os.path.join(save_path, 'snp_counts.csv'))
    # save hamming_info
    pd.to_pickle(hamming_info, os.path.join(save_path, "hamming_info.pkl"))
    # save snp_density_info
    pd.to_pickle(snp_density_info, os.path.join(save_path, "snp_density_info.pkl"))


def compare_two_vcfs_first_mismatch_print_reason(vcf1, vcf2, tsv_path, sample1=None, sample2=None):
    """
    Prints one concise line on the first difference:
      - '<vcf1_basename> | other <reason>'
      - '<vcf2_basename> | other <reason>'
      - 'both | <reason>'
      - 'none | <reason>'
    Prints nothing if the VCFs are identical. Returns True if it printed, else False.
    """

    def _open(path):
        return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, "r", encoding="utf-8")

    def _parse_vcf(path, sample_name=None):
        """Yield (chrom,pos, ref, alt, gt_list) in file order for the chosen sample (first by default)."""
        with _open(path) as f:
            sample_col = None
            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 10:
                        raise ValueError(f"No sample column in {os.path.basename(path)}")
                    sample_col = 9 if sample_name is None else cols.index(sample_name)
                    continue
                parts = line.rstrip("\n").split("\t")
                if not parts or len(parts) < (sample_col or 9) + 1:
                    continue
                chrom, pos = parts[0], int(parts[1])
                ref, alt = parts[3], parts[4]
                gt_field = parts[sample_col].split(":", 1)[0]
                gt = [] if gt_field in (".", "./.", ".|.") else gt_field.replace("|","/").split("/")
                yield (chrom, pos, ref, alt, gt)

    _hdr_strip = re.compile(r'^\s*Block\s+length\s+\d+\s+')

    def _load_tsv(tsv_path):
        """Load TSV into dict[(chrom,pos)] = (ref, alt, gt_list)."""
        recs = {}
        opener = gzip.open if tsv_path.endswith(".gz") else open
        mode   = "rt" if tsv_path.endswith(".gz") else "r"
        with opener(tsv_path, mode, encoding="utf-8") as f:
            for line in f:
                s = _hdr_strip.sub('', line.strip())
                if not s or s.startswith(("var_id", "#")):
                    continue
                parts = s.split()
                if len(parts) < 6:
                    continue
                _, chrom, pos, ref, alt, *haps = parts
                recs[(chrom, int(pos))] = (ref, alt, [str(int(x)) for x in haps])
        return recs

    def _matches_tsv(rec, tsv_rec):
        """Strict match: REF/ALT equal and GT vector exactly equal."""
        ref, alt, gt = rec
        rT, aT, gT = tsv_rec
        return (ref, alt) == (rT, aT) and gt == gT



    name1 = os.path.basename(vcf1)
    name2 = os.path.basename(vcf2)
    TSV = _load_tsv(tsv_path)

    dict1 = { (c,p):(r,a,g) for (c,p,r,a,g) in _parse_vcf(vcf1, sample1) }
    dict2 = { (c,p):(r,a,g) for (c,p,r,a,g) in _parse_vcf(vcf2, sample2) }

    keys1, keys2 = set(dict1), set(dict2)

    # Helper to emit one-line verdict
    def _emit(verdict, reason):
        print(f"{verdict} | {reason}")

    # 1) Missing in vcf2
    for k in sorted(keys1 - keys2):
        chrom, pos = k
        if k in TSV:
            m1 = _matches_tsv(dict1[k], TSV[k])
            if m1:
                _emit(name1, f"other missing at {chrom}:{pos}")
            else:
                _emit("none", f"TSV mismatch at {chrom}:{pos}; other missing")
        else:
            _emit("none", f"TSV lacks {chrom}:{pos}; other missing")
        return True

    # 2) Missing in vcf1
    for k in sorted(keys2 - keys1):
        chrom, pos = k
        if k in TSV:
            m2 = _matches_tsv(dict2[k], TSV[k])
            if m2:
                _emit(name2, f"other missing at {chrom}:{pos}")
            else:
                _emit("none", f"TSV mismatch at {chrom}:{pos}; other missing")
        else:
            _emit("none", f"TSV lacks {chrom}:{pos}; other missing")
        return True

    # 3) Shared: check content and decide
    for k in sorted(keys1 & keys2):
        chrom, pos = k
        r1,a1,g1 = dict1[k]
        r2,a2,g2 = dict2[k]

        # classify difference
        diff_type = None
        short = ""
        if (r1,a1) != (r2,a2):
            diff_type = "refalt"
            short = f"REF/ALT {r1}>{a1} vs {r2}>{a2} at {chrom}:{pos}"
        elif len(g1) != len(g2):
            diff_type = "ploidy"
            short = f"ploidy {len(g1)} vs {len(g2)} at {chrom}:{pos}"
        elif g1 != g2:
            diff_type = "gt"
            short = f"GT {'/'.join(g1)} vs {'/'.join(g2)} at {chrom}:{pos}"

        if diff_type:
            if k not in TSV:
                _emit("none", f"TSV lacks {chrom}:{pos}; {short}")
                return True
            t = TSV[k]
            m1 = _matches_tsv((r1,a1,g1), t)
            m2 = _matches_tsv((r2,a2,g2), t)
            if m1 and m2:
                _emit("both", short)
            elif m1:
                _emit(name1, f"other {short}")
            elif m2:
                _emit(name2, f"other {short}")
            else:
                _emit("none", f"TSV mismatch at {chrom}:{pos}; {short}")
            return True

    # Identical: print nothing
    return False


def unphase_and_sort_vcf(vcf_path, output_path, sample_column='Sample'):
    """
    Reads a VCF file, unphases the genotypes in the specified sample column,
    sorts the alleles while preserving their counts, and sets QUAL to 75.
    
    Args:
        vcf_path: Path to input VCF file
        output_path: Path to output VCF file
        sample_column: Name of the sample column to modify (default: '00')
    """
    with open(vcf_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            # Write header lines as-is
            if line.startswith('#'):
                outfile.write(line)
                # Get sample column index from header line
                if line.startswith('#CHROM'):
                    headers = line.strip().split('\t')
                    try:
                        sample_idx = headers.index(sample_column)
                    except ValueError:
                        raise ValueError(f"Sample column '{sample_column}' not found in VCF header")
                continue
            
            # Process data lines
            fields = line.strip().split('\t')
            
            # Set QUAL column (index 5) to 75
            fields[5] = '75'
            
            # Get the genotype field (should be in sample column)
            genotype = fields[sample_idx]
            
            # Split by ':' in case there are multiple FORMAT fields
            gt_fields = genotype.split(':')
            gt = gt_fields[0]  # GT is always first
            
            # Replace phase separator '|' with unphased separator '/'
            # and sort the alleles
            if '|' in gt:
                alleles = gt.split('|')
                # Sort alleles (as strings, which works fine for single digits)
                alleles_sorted = sorted(alleles)
                # Join with unphased separator
                gt_unphased = '/'.join(alleles_sorted)
                
                # Replace the GT field
                gt_fields[0] = gt_unphased
                fields[sample_idx] = ':'.join(gt_fields)
            
            # Write modified line
            outfile.write('\t'.join(fields) + '\n')


def extractHAIRS_mbq4_long_auto():
    run_all_sh = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_auto_extractHAIRS2/extractHAIRS_mbq7.sh'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_auto_long'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'
    commands_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_auto_extractHAIRS2'

    if not os.path.exists(commands_path):
        os.makedirs(commands_path)
    
    coverages = [3, 5, 10, 20, 40]
    ploidies = [2, 3, 4, 6, 8]
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    chrom = 'Chr1'
    all_sh = ''
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                cov_path = os.path.join(dataset_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage))
                this_sample = ''
                sh_file_name = 'ploidy_' + str(ploidy) + '_mut_' + str(mr) + '_cov_' + str(coverage) + '.sh'
                for sample in samples:
                    # this_sample = ''
                    print('Ploidy:', ploidy, '- mr:', mr, '- Coverage:', coverage, '- sample:', sample)     
                    haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                    out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                    ext_h_command = f'{extract_location} --bam {cov_path}/{str(sample).zfill(2)}.bam --VCF {out_vcf_path} --mbq 7 --out {cov_path}/{str(sample).zfill(2)}_mbq7.frag\n\n\n'
                    this_sample += ext_h_command
                with open(os.path.join(commands_path, sh_file_name), "a") as f:
                    f.write(this_sample)
                all_sh += 'bash ' + os.path.join(commands_path, sh_file_name) + '\n'
    with open(run_all_sh, "w") as f:
        f.write(all_sh)


def extractHAIRS_mbq4_long_allo():

    run_all_sh = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_allo_extractHAIRS2/extractHAIRS_mbq7.sh'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_allo_long'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    commands_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_allo_extractHAIRS2' 
    
    if not os.path.exists(commands_path):
        os.makedirs(commands_path)
    subgenome_configs = [0.01, 0.05]
    ploidies = [3, 4, 6, 8]
    within_mutation_rates = [0.001, 0.005, 0.0075, 0.01, 0.05]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}
    samples = range(20)
    coverages = [5, 10, 20, 40]
    chrom = 'Chr1'
    all_sh = ''
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for coverage in coverages:
                    cov_path = os.path.join(dataset_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage))
                    # stop
                    this_sample = ''
                    sh_file_name = 'ploidy_' + str(ploidy) + '_subgenome_' + str(sgc) + '_mut_' + str(mr) + '_cov_' + str(coverage) + '.sh'
                    for sample in samples:
                        print('Ploidy:', ploidy, '- subgenome:', sgc, '- mr:', mr, '- Coverage:', coverage, '- sample:', sample)     
                        haplotype_path = os.path.join(hap_path, f"subgenome_config_mut{sgc}", f"{ploidy_folder_name[ploidy]}" ,f"within_mut_{mr}", f"{sample:02d}")
                        out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                        ext_h_command = f'{extract_location} --bam {cov_path}/{str(sample).zfill(2)}.bam --VCF {out_vcf_path} --mbq 7 --out {cov_path}/{str(sample).zfill(2)}_mbq7.frag\n\n\n'
                        this_sample += ext_h_command
                    with open(os.path.join(commands_path, sh_file_name), "a") as f:
                        f.write(this_sample)
                    all_sh += 'bash ' + os.path.join(commands_path, sh_file_name) + '\n'
    with open(run_all_sh, "w") as f:
        f.write(all_sh)


def generate_fasta_vcfs_from_tsv_new_allo_new():
    """
    Complete regeneration pipeline from varianthaplos files
    """
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    subgenome_configs = ['0.0005', '0.0001']
    ploidies = [3, 4, 6]
    samples = range(10)
    chrom = 'Chr1'
    within_mutation_rates = ['0.00005', '0.0001']
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C'}
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for sample in samples:
                    print('Ploidy:', ploidy, 'Subgenome_mr:', sgc, 'within_mr:', mr, 'sample:', sample)
                    # stop
                    haplotype_path = os.path.join(hap_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr) , str(sample).zfill(1))
                    # variation_txt_path = [os.path.join(haplotype_path, str(sample).zfill(2) + '_' + sbg + '_varianthaplos.txt.gz') for sbg in subgenomes]
                    # variation_txt_name = os.path.join(haplotype_path, str(sample).zfill(2) + '_varianthaplos.txt.gz')
                    
                    # Remove old files
                    files_to_remove =  ['HPOPG_Chr1.vcf',  'Chr1.fa.gz',  'Chr1.vcf.gz', f'Chr1_ploidy_{ploidy}_sample_{str(sample).zfill(2)}.fa.gz']
                    for ff in files_to_remove:
                        file_path = os.path.join(haplotype_path, ff)
                        if os.path.exists(file_path): 
                            os.remove(file_path)
                            print(f'  Removed: {ff}')
                        
                    
                    # Input files
                    fasta_files = [f for f in os.listdir(haplotype_path) if 'fa.gz' in f]
                    fasta_paths = [os.path.join(haplotype_path, f) for f in fasta_files]
                    fasta_names = ['haplotype_' + f.split('.fa.gz')[0].split('_')[2].split('hap')[1] + '_subgenome_' + f.split('.fa.gz')[0].split('_')[1] for f in fasta_files]
                    varianthaplos_file = os.path.join(haplotype_path, f'{str(sample).zfill(1)}_varianthaplos.txt.gz')
                    # haplotype_names = [os.path.join(haplotype_path, f'{str(sample).zfill(2)}_hap{i+1}.fa.gz') for i in range(ploidy)]
                
                    # Output files
                    out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                    fasta_file_name = os.path.join(haplotype_path, f'{chrom}_ploidy_{ploidy}_sample_{str(sample).zfill(1)}.fa')
                    chrom_fasta_path = os.path.join(haplotype_path, f'{chrom}.fa')
                    vcf_unphased = os.path.join(haplotype_path, chrom + '_unphased.vcf')
                    genotype_path = os.path.join(haplotype_path, 'genotype.csv')
                    
                    # STEP 1: Filter homozygous sites from varianthaplos
                    filtered_varianthaplos = filter_homozygous_sites(varianthaplos_file, ploidy)
                    
                    # STEP 2: Create VCF from filtered varianthaplos
                    tsv_chrom_to_vcf(filtered_varianthaplos, chrom, out_vcf_path, sample_name=str(sample).zfill(1))
                    
                    # STEP 3: Create reference FASTA from haplotype 1 + VCF
                    create_reference_from_haplotype_and_vcf(fasta_paths[0], out_vcf_path, chrom, chrom_fasta_path)
                    
                    # STEP 4: Create concatenated FASTA for ART
                    write_concat_chromosome_simple(fasta_paths, chrom, fasta_file_name, names_list=fasta_names)

                    # Clean up temp file
                    if filtered_varianthaplos != varianthaplos_file:
                        os.remove(filtered_varianthaplos)

                    # STEP 5: generate genotype csv file
                    vcf_to_genotype(out_vcf_path, genotype_path, sample_name=str(sample).zfill(1), ploidy=ploidy)
                    
                    # STEP 6: generate unphased vcf to be used in other methods
                    unphase_and_sort_vcf(out_vcf_path, vcf_unphased, sample_column=str(sample).zfill(1))             
                    
                    print(f'  ✓ Created VCF: {out_vcf_path}')
                    print(f'  ✓ Created reference: {chrom_fasta_path}')
                    print(f'  ✓ Created concat FASTA: {fasta_file_name}')


def index_references_allo_new():
    sh_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/index_ref_allo_new.sh'
    commands = ''
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    within_mut_rates  = ['0.00005', '0.0001']     
    subgenome_mut_rates = ['0.0005', '0.0001']
    ploidies = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C'}
    samples = range(10)
    chrom = 'Chr1'
    for sgr in subgenome_mut_rates:
        for ploidy in ploidies.keys():
            for mr in within_mut_rates:
                for sample in samples:
                    print(f"Ploidy: {ploidy}, subgenome μ={sgr}, within μ={mr}, sample {sample:01d}")
                    haplotype_path = os.path.join(hap_path, f"subgenome_config_mut{sgr}", f"{ploidies[ploidy]}" ,f"within_mut_{mr}", f"{sample:01d}")                   
                    chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                    chrom_fasta_path_mmi = os.path.join(haplotype_path, chrom + '.mmi')
                    idx_comm = f'bwa index {chrom_fasta_path}\n'
                    idx_comm2 = f'minimap2 -d {chrom_fasta_path_mmi} {chrom_fasta_path}\n'
                    commands += idx_comm
                    commands += idx_comm2
    # Save to shell script
    with open(sh_path, 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write(commands)


def simulate_short_read_dataset_allo_new():

    run_all_sh = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/short_allo_new/simulate_short_allo.sh'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_allo_short_new'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    commands_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/short_allo_new' 
    art_path = 'art_illumina'
    
    if not os.path.exists(dataset_path):
        os.makedirs(dataset_path)

    if not os.path.exists(commands_path):
        os.makedirs(commands_path)
    
    # subgenome_configs = [0.01, 0.05]
    subgenome_configs = ['0.0005', '0.0001']
    within_mutation_rates = ['0.00005', '0.0001']
    ploidies = [3, 4, 6]
    # within_mutation_rates = [0.001, 0.005, 0.0075, 0.01, 0.05]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C'}
    samples = range(10)
    coverages = [5, 10, 20, 40]
    read_length = 125
    mean_insert_length = 350
    std_insert_length = 50
    chrom = 'Chr1'
    all_sh = ''
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for coverage in coverages:
                    cov_path = os.path.join(dataset_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage))
                    if not os.path.exists(cov_path):
                        os.makedirs(cov_path, exist_ok=True)
                    # stop
                    done_sampels = [f for f in os.listdir(cov_path) if '.frag' in f]
                    if len(done_sampels) < 10:
                        this_sample = ''
                        sh_file_name = 'ploidy_' + str(ploidy) + '_subgenome_' + str(sgc) + '_mut_' + str(mr) + '_cov_' + str(coverage) + '.sh'
                        for sample in samples:

                            print('Ploidy:', ploidy, '- subgenome:', sgc, '- mr:', mr, '- Coverage:', coverage, '- sample:', sample)     
                            
                            rs = random.randint(1, 2**32)
                            haplotype_path = os.path.join(hap_path, f"subgenome_config_mut{sgc}", f"{ploidy_folder_name[ploidy]}" ,f"within_mut_{mr}", f"{sample:01d}")
                            out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                            fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(1) +'.fa')
                            chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                            mmi_path = os.path.join(haplotype_path, chrom + '.mmi')
                            fastq_file_name = os.path.join(cov_path, str(sample).zfill(1))
                            comm0 = f'# Ploidy: {ploidy} - subgenome mutation rate: {sgc} - within mutation rate: {mr} - Coverage: {coverage} - sample: {sample}\n'
                            art_command = f'{art_path} -ss HS25 -rs {rs} -i {fasta_file_name} -p -na -l {read_length} -f {coverage/ploidy} -m {mean_insert_length} -s {std_insert_length} -o {fastq_file_name};'
                            align_command = f'bwa mem {chrom_fasta_path} {fastq_file_name}1.fq {fastq_file_name}2.fq > {fastq_file_name}.sam; samtools view -Sb {fastq_file_name}.sam | samtools sort -o {fastq_file_name}.bam; samtools index {fastq_file_name}.bam; rm {fastq_file_name}.sam;'
                            ext_h_command = f'{extract_location} --bam {fastq_file_name}.bam --VCF {out_vcf_path} --out {fastq_file_name}.frag\n\n\n'
                            this_sample += comm0 + art_command + align_command + ext_h_command
                        with open(os.path.join(commands_path, sh_file_name), "a") as f:
                            f.write(this_sample)
                        all_sh += 'bash ' + os.path.join(commands_path, sh_file_name) + '\n'
    with open(run_all_sh, "w") as f:
        f.write(all_sh)


def simulate_long_read_dataset_allo_new():
    run_all_sh = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_allo_new/simulate_long_allo.sh'
    pbsim_location = '/home/mah19006/pHapCompass_v2/pbsim3/src/pbsim'
    extract_location = '~/pHapCompass_v2/Hap10/extract_poly/build/extractHAIRS'
    model_path = '/home/mah19006/pHapCompass_v2/pbsim3/data/QSHMM-ONT-HQ.model' # QSHMM-ONT-HQ.model
    # model_path = '/home/mah19006/pHapCompass_v2/pbsim3/data/QSHMM-RSII.model'
    dataset_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_allo_long_new'
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/allo'
    commands_path = '/mnt/research/aguiarlab/proj/HaplOrbit/scripts/simulation/long_allo_new' 
    
    if not os.path.exists(dataset_path):
        os.makedirs(dataset_path)

    if not os.path.exists(commands_path):
        os.makedirs(commands_path)
    
    subgenome_configs = ['0.0005', '0.0001']
    within_mutation_rates = ['0.00005', '0.0001']
    ploidies = [3, 4, 6]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C'}
    samples = range(10)
    coverages = [5, 10, 20, 40]
    chrom = 'Chr1'
    all_sh = ''
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for coverage in coverages:
                    cov_path = os.path.join(dataset_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage))
                    if not os.path.exists(cov_path):
                        os.makedirs(cov_path, exist_ok=True)
                    # stop
                    done_sampels = [f for f in os.listdir(cov_path) if '.frag' in f]
                    if len(done_sampels) < 10:
                        this_sample = ''
                        sh_file_name = 'ploidy_' + str(ploidy) + '_subgenome_' + str(sgc) + '_mut_' + str(mr) + '_cov_' + str(coverage) + '.sh'
                        for sample in samples:

                            print('Ploidy:', ploidy, '- subgenome:', sgc, '- mr:', mr, '- Coverage:', coverage, '- sample:', sample)     
                            # if not os.path.exists(os.path.join(cov_path, str(sample).zfill(2) + '.frag')):
                            #     print('not done.')
                                # stop
                            # haplotype_path = os.path.join(hap_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2))
                            haplotype_path = os.path.join(hap_path, f"subgenome_config_mut{sgc}", f"{ploidy_folder_name[ploidy]}" ,f"within_mut_{mr}", f"{sample:01d}")
                            out_vcf_path = os.path.join(haplotype_path, chrom + '.vcf')
                            fasta_file_name = os.path.join(haplotype_path, chrom + '_ploidy_' + str(ploidy) + '_sample_' + str(sample).zfill(1) +'.fa')
                            # chrom_fasta_path = os.path.join(haplotype_path, chrom + '.fa')
                            mmi_path = os.path.join(haplotype_path, chrom + '.mmi')
                            comm0 = f'# Ploidy: {ploidy} - mutation rate: {mr} - Coverage: {coverage} - sample: {sample}\n'
                            comm1 = 'cd {}\n'.format(cov_path)
                            # pbsim_command = '{} --strategy wgs --method qshmm --qshmm {} --depth {} --genome {} --accuracy-mean 0.999 0.0005 --prefix {}\n'.format(pbsim_location, model_path, coverage/ploidy, fasta_file_name, str(sample).zfill(2))
                            pbsim_command = '{} --strategy wgs --method qshmm --qshmm {} --depth {} --genome {} --prefix {}\n'.format(pbsim_location, model_path, coverage/ploidy, fasta_file_name, str(sample).zfill(1))
                            align_command = f'minimap2 -t 16 -x map-ont -a {mmi_path} {cov_path}/{str(sample).zfill(1)}_00*.fq.gz | samtools sort -@8 -o {cov_path}/{str(sample).zfill(1)}.bam\n'
                            index_command = f'samtools index {cov_path}/{str(sample).zfill(1)}.bam\n'
                            ext_h_command = f'{extract_location} --bam {cov_path}/{str(sample).zfill(1)}.bam --VCF {out_vcf_path} --mbq 4 --out {cov_path}/{str(sample).zfill(1)}_mbq4.frag\n\n\n'
                            this_sample += comm0 + comm1 + pbsim_command + align_command + index_command + ext_h_command
                        with open(os.path.join(commands_path, sh_file_name), "a") as f:
                            f.write(this_sample)
                        all_sh += 'bash ' + os.path.join(commands_path, sh_file_name) + '\n'
    with open(run_all_sh, "w") as f:
        f.write(all_sh)
