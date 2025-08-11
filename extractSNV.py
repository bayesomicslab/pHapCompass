#!/usr/bin/env python3
import pysam
from cyvcf2 import VCF

# ── Config ─────────────────────────────────────────────────────────────────────
vcf_path     = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878/contig_100/ploidy_3/NA12878_ploidy3_contig100.vcf"
bam_path     = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878/contig_100/ploidy_3/cov_10/bam/00.bam"
out_prefix   = "quickSNV"
min_mapq     = 20
min_snvs_per_read = 2

# ── 1) load SNV positions ──────────────────────────────────────────────────────
vcf = VCF(vcf_path)
# read positions into a sorted list of 0-based coords
snv_pos = [rec.POS - 1 for rec in vcf]  
snv_pos.sort()
nSNV = len(snv_pos)

# write positions file
with open(f"{out_prefix}_SNV_pos.txt","w") as f:
    f.write(" ".join(map(str, snv_pos))+"\n")

# ── 2) open BAM for streaming ───────────────────────────────────────────────────
bam = pysam.AlignmentFile(bam_path, "rb")

# ── 3) emit fragment file ──────────────────────────────────────────────────────
with open(f"{out_prefix}_SNV_matrix.txt","w") as out:
    out.write(f"{nSNV}\n")           # first line = #SNVs
    for read in bam.fetch():         # one read at a time
        if read.is_unmapped or read.mapping_quality < min_mapq:
            continue

        # get list of (read_idx, ref_pos) from the CIGAR
        # read.get_aligned_pairs() returns list of (read_pos, ref_pos)
        # we want ref_pos -> read_pos
        pos2q = {rp: qp for qp, rp in read.get_aligned_pairs(matches_only=True)}

        pairs = []
        # for each SNV, check if this read covers it
        for j, ref in enumerate(snv_pos):
            if ref in pos2q:
                base = read.query_sequence[pos2q[ref]]
                code = {"A":1,"C":2,"G":3,"T":4}.get(base.upper(), 0)
                if code:
                    pairs.append(f"{j},{code}")

        if len(pairs) >= min_snvs_per_read:
            out.write(" ".join(pairs) + "\n")
