from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple, Iterable, Optional, Union
from scipy.sparse import csr_matrix, csc_matrix
import numpy as np
import pandas as pd
import logging
import subprocess
import csv
import gzip
import io
import os



@dataclass
class InputConfigSparse:
    data_path: str                 # .frag file (plain text)
    genotype_path: Optional[str]   # CSV of haplotypes (optional for this step)
    ploidy: int                    # K
    alleles: Tuple[int, int] = (0, 1)

    # Encoding convention: 0 = no-call, 1 = REF, 2 = ALT
    refalt_encoding: Dict[str, int] = None

    def __post_init__(self):
        if self.refalt_encoding is None:
            self.refalt_encoding = {'0': 1, '1': 2, '-': 0, '.': 0}


class SparseFragment:
    """
    Parse extractHAIRS/extractHAIRS_poly style .frag files into a CSR matrix.

    Rows   = reads (fragments)
    Cols   = SNP positions (integer variant indices)
    Values = {0: no-call, 1: REF, 2: ALT}

    We support records that contain multiple (start, allele-string) "blocks"
    per read line (as in extract-poly output). For each block, we expand the
    allele string across consecutive SNP indices starting at 'start'.

    Example .frag line (conceptual):
        <B> <READ_NAME> <...unused...> <start1> <allele_string1> <start2> <allele_string2> ...

        where B is the number of blocks that follow (each block = start + string)

    Notes:
    - Duplicate assignments to the same (read, SNP) cell are checked; if conflicting,
      we retain the first and count a conflict (exposed in .n_conflicts).
    - Unknown characters map to 0 (no-call) by default; you can extend `refalt_encoding`.
    """

    def __init__(self,
                 cfg: InputConfigSparse,
                 positions_from_genotype: Optional[Iterable[int]] = None):
        self.cfg = cfg

        # Row/col mappings
        self.read_ids: List[str] = []
        self.col_to_snp: List[int] = []       # col index -> SNP position (int)
        self.snp_to_col: Dict[int, int] = {}  # SNP position (int) -> col index

        # Sparse matrix (built via ._rows/. _cols/. _data then finalized)
        self._rows: List[int] = []
        self._cols: List[int] = []
        self._data: List[int] = []
        self._finalized: bool = False

        # Stats
        self.n_lines: int = 0
        self.n_blocks: int = 0
        self.n_assignments: int = 0
        self.n_conflicts: int = 0

        # Optionally preseed SNP columns from genotype header
        if positions_from_genotype is not None:
            for pos in positions_from_genotype:
                self._ensure_col(pos)

        # Parse immediately
        self._parse_frag(cfg.data_path)
        self._finalize()

    # ---------- Public API ----------

    @property
    def csr(self) -> csr_matrix:
        if not self._finalized:
            raise RuntimeError("Sparse matrix not finalized yet.")
        return self._csr

    @property
    def csc(self) -> csc_matrix:
        if not self._finalized:
            raise RuntimeError("Sparse matrix not finalized yet.")
        # Convert lazily
        if not hasattr(self, "_csc"):
            self._csc = self._csr.tocsc()
        return self._csc

    @property
    def shape(self) -> Tuple[int, int]:
        return self.csr.shape

    def co_covered_pairs(self, min_reads: int = 1) -> np.ndarray:
        """
        Return SNP pairs (as 2-column int array of *column indices*) that share at least `min_reads` co-coverage.
        Useful later for building the quotient/pair node set without any graph library.
        """
        occ = (self.csr > 0).astype(np.int32)
        cov = (occ.T @ occ)           # denseish for small windows; if huge, we can block
        cov.setdiag(0)
        ii, jj = cov.nonzero()
        keep = cov.data >= min_reads
        return np.vstack([ii[keep], jj[keep]]).T

    def rows_covering_all(self, cols: List[int]) -> np.ndarray:
        """
        Return row indices of reads that cover ALL columns in `cols`.
        """
        sub = self.csr[:, cols]
        return np.where(np.array(sub.getnnz(axis=1)).ravel() == len(cols))[0]

    def pair_counts_4(self, col_i: int, col_j: int) -> Tuple[int, int, int, int]:
        """
        Fast 4-pattern counts for a SNP pair (by *column* indices):
            n00, n01, n10, n11   (0/1 means REF/ALT; here REF→1, ALT→2 in the matrix)

        We ignore no-calls (0) by intersecting nonzero rows and then tallying.
        """
        Mc = self.csc  # ensure CSC once
        ci = Mc[:, col_i]
        cj = Mc[:, col_j]

        ri = set(ci.indices.tolist())
        rj = set(cj.indices.tolist())
        rows = np.fromiter(ri & rj, dtype=np.int64)  # rows with both nonzero

        if rows.size == 0:
            return 0, 0, 0, 0

        vi = self.csr[rows, col_i].A.ravel()
        vj = self.csr[rows, col_j].A.ravel()

        # Map {1,2}→{0,1} (REF→0, ALT→1) for counting pattern bins
        bi = (vi == 2).astype(np.int8)
        bj = (vj == 2).astype(np.int8)
        patt = (bi << 1) | bj  # 2-bit code: 0=00,1=01,2=10,3=11

        cnt = np.bincount(patt, minlength=4)
        return int(cnt[0]), int(cnt[1]), int(cnt[2]), int(cnt[3])

    # ---------- Internal helpers ----------

    def _parse_frag(self, path: str):
        def _open_any(p: str):
            if p.endswith(".gz"):
                return io.TextIOWrapper(gzip.open(p, "rb"))
            return open(p, "r")

        with _open_any(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                self.n_lines += 1

                parts = line.split()
                # We accept both “with read name” and “without read name” styles.
                # Common extract-poly format seen in your older code:
                # parts[0] = B (num blocks)
                # parts[1] = read name (sometimes)
                # then repeating pairs: start_i, allele_string_i
                try:
                    n_blocks = int(parts[0])
                except ValueError:
                    raise ValueError(f"Malformed .frag line {self.n_lines}: expected integer in column 1")

                # Decide if a read name is present in col 2
                start_col = 2 if self._looks_like_readname(parts[1]) else 1
                read_name = parts[1] if start_col == 2 else f"read_{self.n_lines}"
                row = self._ensure_row(read_name)

                # Parse blocks
                idx = start_col
                for b in range(n_blocks):
                    if idx + 1 >= len(parts):
                        raise ValueError(f"Line {self.n_lines}: expected start+allele_string for block {b+1}")
                    try:
                        start_pos = int(parts[idx])
                    except ValueError:
                        raise ValueError(f"Line {self.n_lines}: start position not int at token {idx+1}")
                    allele_str = parts[idx + 1]
                    self._ingest_block(row, start_pos, allele_str)
                    self.n_blocks += 1
                    idx += 2

    def _looks_like_readname(self, token: str) -> bool:
        # Heuristic: non-integer tokens treated as read names.
        try:
            int(token)
            return False
        except ValueError:
            return True

    def _ingest_block(self, row: int, start_pos: int, allele_str: str):
        """
        Expand a (start_pos, allele_string) run into per-SNP updates.
        start_pos refers to the *SNP index* used throughout your draft/code (0-based or 1-based).
        We will treat it as 0-based here; if your .frag is 1-based, set START_IS_ONE_BASED=True.
        """
        START_IS_ONE_BASED = True  # flip here if needed to match your .frag
        pos0 = start_pos - 1 if START_IS_ONE_BASED else start_pos


        for k, ch in enumerate(allele_str):
            snp = pos0 + k
            val = self.cfg.refalt_encoding.get(ch, 0)  # default to no-call if unknown
            if val == 0:
                continue
            col = self._ensure_col(snp)
            self._assign(row, col, val)

    def _ensure_row(self, read_name: str) -> int:
        """Get/create row index for a read."""
        try:
            return self._readname_to_row[read_name]
        except AttributeError:
            self._readname_to_row = {}
        except KeyError:
            pass

        row = len(self.read_ids)
        self.read_ids.append(read_name)
        self._readname_to_row[read_name] = row
        return row

    def _ensure_col(self, snp_pos: int) -> int:
        """Get/create column index for a SNP position (integer index)."""
        if snp_pos in self.snp_to_col:
            return self.snp_to_col[snp_pos]
        col = len(self.col_to_snp)
        self.col_to_snp.append(snp_pos)
        self.snp_to_col[snp_pos] = col
        return col

    def _assign(self, row: int, col: int, val: int):
        """
        Append an assignment. If (row,col) already assigned, enforce consistency:
        - If same value, ignore extra (duplicate) write.
        - If conflicting, keep the first and count a conflict.
        """
        self.n_assignments += 1

        # Conflict detection: look back only for the very last assignment to this cell to avoid O(nnz) scan.
        # For correctness we use a map when density of duplicates is >0.
        if not hasattr(self, "_cell_seen"):
            self._cell_seen = {}

        key = (row, col)
        if key in self._cell_seen:
            prev_val_idx = self._cell_seen[key]
            prev_val = self._data[prev_val_idx]
            if prev_val != val:
                self.n_conflicts += 1
            return  # keep the first value
        else:
            self._rows.append(row)
            self._cols.append(col)
            self._data.append(val)
            self._cell_seen[key] = len(self._data) - 1

    def _finalize(self):
        n_rows = len(self.read_ids)
        n_cols = len(self.col_to_snp)
        if n_rows == 0 or n_cols == 0:
            self._csr = csr_matrix((0, 0), dtype=np.uint8)
        else:
            self._csr = csr_matrix(
                (np.asarray(self._data, dtype=np.uint8),
                 (np.asarray(self._rows, dtype=np.int64),
                  np.asarray(self._cols, dtype=np.int64))),
                shape=(n_rows, n_cols),
            )
        self._finalized = True

    # ---------- Convenience ----------

    def summary(self) -> str:
        return (
            f"SparseFragment: reads={len(self.read_ids)}, snps={len(self.col_to_snp)}, "
            f"nnz={self._csr.nnz}, lines={self.n_lines}, blocks={self.n_blocks}, "
            f"assignments={self.n_assignments}, conflicts={self.n_conflicts}"
        )


def vcf_gt_to_csv(vcf_path: str,
                  csv_path: str,
                  ploidy: int,
                  sample_name: str = None,
                  contig_filter: str | None = None) -> int:
    """
    Read a VCF and write a CSV with columns haplotype_0..haplotype_{ploidy-1}
    containing the genotype digits from the sample's GT field (e.g., 0/1/1/1).

    Parameters
    ----------
    vcf_path : str
        Path to the VCF (.vcf or .vcf.gz).
    csv_path : str
        Output CSV path.
    sample_name : str, optional
        If provided, select this sample column; otherwise use the only sample present.
    ploidy : int
        Expected number of alleles in GT (default 4).
    contig_filter : str | None
        If provided, keep only rows with CHROM == contig_filter.

    Returns
    -------
    int
        Number of rows written.
    """
    # opener for gz/plain
    opener = gzip.open if vcf_path.endswith(".gz") else open

    # headers for CSV
    fieldnames = [f"haplotype_{i}" for i in range(ploidy)]

    n_written = 0
    with opener(vcf_path, "rt") as fin, open(csv_path, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()

        sample_idx = None
        gt_index_in_format = None

        for line in fin:
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                # Parse header to identify sample column
                header = line.rstrip("\n").split("\t")
                # VCF fixed columns + FORMAT + samples...
                #   0     1   2   3   4   5     6      7      8      9..
                # #CHROM POS ID  REF ALT QUAL FILTER INFO   FORMAT  sample(s)
                if len(header) < 10:
                    raise ValueError("VCF has no sample columns.")
                samples = header[9:]
                if sample_name is None:
                    if len(samples) != 1:
                        raise ValueError(
                            f"VCF contains {len(samples)} samples; specify sample_name."
                        )
                    sample_idx = 9  # first and only sample
                else:
                    try:
                        sample_idx = header.index(sample_name)
                    except ValueError:
                        raise ValueError(f"sample_name '{sample_name}' not found in VCF header.")
                continue

            # Data lines
            parts = line.rstrip("\n").split("\t")
            if len(parts) < sample_idx + 1:
                continue

            chrom = parts[0]
            if contig_filter is not None and chrom != contig_filter:
                continue

            fmt = parts[8]         # FORMAT string like "GT:DP:AD:RO:QR:AO:QA"
            sample_field = parts[sample_idx]

            # Locate GT once per distinct FORMAT
            # (FORMAT is usually stable, but we can re-scan if it changes)
            fmt_keys = fmt.split(":")
            try:
                gt_index_in_format = fmt_keys.index("GT")
            except ValueError:
                # No GT => skip
                continue

            sample_values = sample_field.split(":")
            if gt_index_in_format >= len(sample_values):
                continue

            gt = sample_values[gt_index_in_format]
            # Accept phased or unphased separators
            if "|" in gt:
                alleles = gt.split("|")
            else:
                alleles = gt.split("/")

            if len(alleles) != ploidy:
                # Skip if ploidy mismatch (e.g., missing or wrong model)
                continue
            if any(a == "." or a == "" for a in alleles):
                continue
            # Keep only 0/1; if multi-allelic snuck in, map non-'0' to '1' (optional)
            try:
                nums = [int(a) for a in alleles]
            except ValueError:
                # Non-integer allele codes; skip
                continue
            # For biallelic VCF this is already 0/1; if a>1 appears, treat as 1
            nums = [0 if x == 0 else 1 for x in nums]

            row = {f"haplotype_{i}": nums[i] for i in range(ploidy)}
            writer.writerow(row)
            n_written += 1


def infer_ploidy_from_vcf(vcf_path: str, sample_index: int = 0) -> int:
    """
    Infer ploidy from the first non-header variant line of a VCF.
    Assumes GT is present and same ploidy across samples.
    """
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    with open_fn(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise ValueError(f"VCF {vcf_path} has no genotype columns.")
            sample_fmt = fields[8].split(":")
            try:
                gt_index = sample_fmt.index("GT")
            except ValueError:
                raise ValueError(f"No GT field in FORMAT column of {vcf_path}.")
            sample_field = fields[9 + sample_index]
            sample_vals = sample_field.split(":")
            gt = sample_vals[gt_index]
            if gt in (".", "./.", ".|."):
                raise ValueError(f"Missing GT in first variant for sample {sample_index}.")
            sep = "/" if "/" in gt else "|"
            alleles = gt.split(sep)
            return len(alleles)
    raise ValueError(f"Could not infer ploidy from VCF {vcf_path}: no variant lines found.")


def run_extract_hairs(bam_path: str, vcf_path: str, frag_path: str, ploidy: Optional[int] = None, epsilon: float = 0.01, mbq: int = 13,  extracthairs_bin: str = "extractHAIRS") -> None:
    """
    Call extractHAIRS to generate a fragment file from BAM + VCF.
    Adjust options if needed to match your setup.
    """
    cmd = [
        extracthairs_bin,
        "--bam", bam_path,
        "--vcf", vcf_path,
        "--out", frag_path,
        "--error", str(epsilon),
    ]
    if ploidy is not None:
        cmd.extend(["--ploidy", str(ploidy)])
    if mbq != 13:
        cmd.extend(["--mbq", str(mbq)])

    logging.info("Running extractHAIRS: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error("extractHAIRS failed with return code %s", e.returncode)
        raise


def bam_root_name(bam_path: str) -> str:
    """
    Get a clean base name X from X.bam or X.bam.gz etc.
    """
    base = os.path.basename(bam_path)
    # handle .gz, .bgz etc.
    root, ext = os.path.splitext(base)
    if ext in [".gz", ".bgz"]:
        root, ext2 = os.path.splitext(root)
        base = root + ext2
    else:
        base = root + ext
    for suffix in (".bam", ".cram", ".sam"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def write_phased_vcf(
    input_vcf_path: str,
    output_vcf_path: str,
    predicted_haplotypes: Union[pd.DataFrame, List[pd.DataFrame]],
    block_ids: Union[pd.Series, List[pd.Series]],
    likelihoods: Optional[Union[float, List[float]]] = None,
    ploidy: int = None,
    sample_name: str = "SAMPLE"
):
    """
    Write phased haplotypes to VCF format.
    
    Parameters
    ----------
    input_vcf_path : str
        Path to input VCF file (can be .vcf or .vcf.gz)
    output_vcf_path : str
        Path to output VCF file (.vcf.gz)
    predicted_haplotypes : pd.DataFrame or List[pd.DataFrame]
        Single DataFrame or list of DataFrames with columns [hap_0, hap_1, ..., hap_k]
        where k = ploidy-1. Values are 0 (REF) or 1 (ALT).
    block_ids : pd.Series or List[pd.Series]
        Single Series or list of Series with phase set IDs (or NaN for unphased)
    likelihoods : float or List[float], optional
        Single likelihood or list of likelihoods for each sample
    ploidy : int, optional
        Ploidy level (inferred from haplotype columns if not provided)
    sample_name : str
        Sample name for VCF column (default: "SAMPLE")
    """
    
    # Normalize inputs to lists
    if isinstance(predicted_haplotypes, pd.DataFrame):
        predicted_haplotypes = [predicted_haplotypes]
    if isinstance(block_ids, pd.Series):
        block_ids = [block_ids]
    if likelihoods is not None and not isinstance(likelihoods, list):
        likelihoods = [likelihoods]
    
    n_samples = len(predicted_haplotypes)
    
    # Infer ploidy if not provided
    if ploidy is None:
        ploidy = len([col for col in predicted_haplotypes[0].columns if col.startswith('hap_')])
    
    # Validate inputs
    assert len(block_ids) == n_samples, "Number of block_ids must match predicted_haplotypes"
    if likelihoods is not None:
        assert len(likelihoods) == n_samples, "Number of likelihoods must match predicted_haplotypes"
    
    # Open input VCF
    open_func = gzip.open if input_vcf_path.endswith('.gz') else open
    
    # Read input VCF header and data
    header_lines = []
    data_lines = []
    
    with open_func(input_vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line.rstrip('\n'))
            else:
                data_lines.append(line.rstrip('\n'))
    
    # Prepare output
    output_lines = []
    
    # Process header
    for line in header_lines:
        if line.startswith('##FORMAT=<ID=GT'):
            # Keep GT definition but update description
            output_lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">')
        elif line.startswith('##FORMAT=<ID=PS'):
            # Keep PS definition
            output_lines.append(line)
        elif line.startswith('#CHROM'):
            # Add our custom FORMAT fields before the #CHROM line
            if likelihoods is not None:
                output_lines.append('##FORMAT=<ID=LK,Number=1,Type=Float,Description="Phasing likelihood/confidence">')
            
            # Add PS if not already present
            has_ps = any('##FORMAT=<ID=PS' in l for l in header_lines)
            if not has_ps:
                output_lines.append('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">')
            
            # Modify column header to use our sample name
            cols = line.split('\t')
            if len(cols) > 9:
                cols[9] = sample_name
            output_lines.append('\t'.join(cols))
        else:
            output_lines.append(line)
    
    # Process data lines
    for idx, data_line in enumerate(data_lines):
        fields = data_line.split('\t')
        
        # Parse existing fields
        chrom = fields[0]
        pos = fields[1]
        snp_id = fields[2]
        ref = fields[3]
        alt = fields[4]
        qual = fields[5]
        filt = fields[6]
        info = fields[7]
        format_field = fields[8] if len(fields) > 8 else 'GT'
        
        # Build GT and PS for all samples
        gt_values = []
        ps_values = []
        lk_values = []
        
        for sample_idx in range(n_samples):
            hap_df = predicted_haplotypes[sample_idx]
            block_series = block_ids[sample_idx]
            
            # Check if this position is phased
            if idx >= len(block_series) or pd.isna(block_series.iloc[idx]):
                # Unphased - use unphased genotype
                # Get genotype from haplotype columns
                hap_cols = [col for col in hap_df.columns if col.startswith('hap_')]
                if idx < len(hap_df):
                    alleles = [str(int(hap_df[col].iloc[idx])) if pd.notna(hap_df[col].iloc[idx]) else '.' 
                              for col in hap_cols]
                    gt_values.append('/'.join(alleles))
                else:
                    gt_values.append('/'.join(['.'] * ploidy))
                ps_values.append('.')
            else:
                # Phased
                hap_cols = [col for col in hap_df.columns if col.startswith('hap_')]
                if idx < len(hap_df):
                    alleles = [str(int(hap_df[col].iloc[idx])) if pd.notna(hap_df[col].iloc[idx]) else '.' 
                              for col in hap_cols]
                    gt_values.append('|'.join(alleles))
                    ps_values.append(str(int(block_series.iloc[idx])))
                else:
                    gt_values.append('|'.join(['.'] * ploidy))
                    ps_values.append('.')
            
            # Add likelihood if available
            if likelihoods is not None:
                lk_values.append(f"{likelihoods[sample_idx]:.6f}")
        
        # Combine multiple samples with ':'
        combined_gt = ':'.join(gt_values)
        combined_ps = ':'.join(ps_values)
        
        # Build FORMAT string
        if n_samples == 1:
            # Single sample
            if ps_values[0] == '.':
                # Unphased
                format_str = 'GT'
                if likelihoods is not None:
                    format_str += ':LK'
                    sample_str = f"{gt_values[0]}:{lk_values[0]}"
                else:
                    sample_str = gt_values[0]
            else:
                # Phased
                format_str = 'GT:PS'
                if likelihoods is not None:
                    format_str += ':LK'
                    sample_str = f"{gt_values[0]}:{ps_values[0]}:{lk_values[0]}"
                else:
                    sample_str = f"{gt_values[0]}:{ps_values[0]}"
        else:
            # Multiple samples (uncertainty mode)
            format_str = 'GT:PS'
            if likelihoods is not None:
                format_str += ':LK'
                sample_str = f"{combined_gt}:{combined_ps}:{':'.join(lk_values)}"
            else:
                sample_str = f"{combined_gt}:{combined_ps}"
        
        # Reconstruct line
        output_line = '\t'.join([chrom, pos, snp_id, ref, alt, qual, filt, info, format_str, sample_str])
        output_lines.append(output_line)
    
    # Write output
    with gzip.open(output_vcf_path, 'wt') as f:
        for line in output_lines:
            f.write(line + '\n')
    
    print(f"Phased VCF written to: {output_vcf_path}")




