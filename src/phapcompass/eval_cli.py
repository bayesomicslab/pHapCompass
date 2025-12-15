#!/usr/bin/env python3
import argparse
import logging
import tempfile
from pathlib import Path
from .evaluations import *
import numpy as np
import pysam
import logging
import subprocess
import shutil
from typing import Optional, Tuple
from .read_input import InputConfigSparse, SparseFragment


# -----------------------------
# Utility: ploidy inference safe
# -----------------------------
def infer_ploidy_from_vcf_safe(vcf_path: str) -> Optional[int]:
    """
    Best-effort ploidy inference from VCF. Returns None on failure.
    """
    try:
        from .read_input import infer_ploidy_from_vcf
        return infer_ploidy_from_vcf(vcf_path)
    except Exception:
        return None


# -----------------------------
# VCF -> haplotype matrix parsers
# -----------------------------
def parse_pred_vcf_to_H_and_blocks(
    pred_vcf: str,
    ploidy: int,
    sample_name: str = "SAMPLE",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse predicted VCF into:
      H_pred: (K x L) float array (NaN for unphased or missing sites)
      block_ids: (L,) float array (NaN for unphased sites)

    Phased sites are defined by presence of PS.
    """
    vcf = pysam.VariantFile(pred_vcf)

    samples = list(vcf.header.samples)
    if not samples:
        raise ValueError(f"Predicted VCF has no samples: {pred_vcf}")
    if sample_name not in samples:
        # Fall back to first sample if needed
        sample_name = samples[0]

    gt_rows = [[] for _ in range(ploidy)]
    block_ids = []

    for rec in vcf.fetch():
        srec = rec.samples[sample_name]
        gt = srec.get("GT", None)  # tuple like (0,0,1) or contains None
        ps = srec.get("PS", None)

        # "phased" defined by PS presence
        phased = ps is not None

        if gt is None or not isinstance(gt, tuple) or any(a is None for a in gt) or len(gt) != ploidy or not phased:
            for k in range(ploidy):
                gt_rows[k].append(np.nan)
            block_ids.append(np.nan)
            continue

        # Map allele codes to 0/1 (biallelic; multi-allelic treated as ALT)
        bits = [0 if a == 0 else 1 for a in gt]

        for k in range(ploidy):
            gt_rows[k].append(float(bits[k]))
        block_ids.append(float(ps))

    H_pred = np.asarray(gt_rows, dtype=float)
    block_ids = np.asarray(block_ids, dtype=float)
    return H_pred, block_ids


def parse_truth_vcf_to_H(
    truth_vcf: str,
    ploidy: int,
    sample_name: Optional[str] = None) -> np.ndarray:
    """
    Parse ground-truth VCF into H_true (K x L).

    Truth VCF may not have PS. We parse GT for all sites and return:
      H_true: (K x L) float array (NaN for missing/invalid GT)
    """
    vcf = pysam.VariantFile(truth_vcf)

    samples = list(vcf.header.samples)
    if not samples:
        raise ValueError(f"Truth VCF has no samples: {truth_vcf}")

    if sample_name is None:
        sample_name = samples[0]
    elif sample_name not in samples:
        # Fall back
        sample_name = samples[0]

    gt_rows = [[] for _ in range(ploidy)]

    for rec in vcf.fetch():
        srec = rec.samples[sample_name]
        gt = srec.get("GT", None)

        if gt is None or not isinstance(gt, tuple) or any(a is None for a in gt) or len(gt) != ploidy:
            for k in range(ploidy):
                gt_rows[k].append(np.nan)
            continue

        bits = [0 if a == 0 else 1 for a in gt]
        for k in range(ploidy):
            gt_rows[k].append(float(bits[k]))

    H_true = np.asarray(gt_rows, dtype=float)
    return H_true


def valid_overlap_columns(H_true: np.ndarray, H_pred: np.ndarray) -> np.ndarray:
    """
    Return column indices where both H_true and H_pred are fully defined (no NaNs).
    """
    if H_true.shape != H_pred.shape:
        raise ValueError(f"Shape mismatch: H_true {H_true.shape} vs H_pred {H_pred.shape}")

    ok_true = ~np.any(np.isnan(H_true), axis=0)
    ok_pred = ~np.any(np.isnan(H_pred), axis=0)
    return np.where(ok_true & ok_pred)[0]


# -----------------------------
# VER evaluation
# -----------------------------
def compute_ver_from_vcfs(
    truth_vcf: str,
    pred_vcf: str,
    ploidy: int,
    truth_sample: Optional[str] = None,
    pred_sample: str = "SAMPLE") -> Tuple[float, int]:
    """
    Compute VER given truth VCF and predicted VCF.
    Blocks are taken from predicted VCF PS.
    Evaluates only on columns where both are defined.
    Returns (ver, evaluated_sites).
    """
    H_true = parse_truth_vcf_to_H(truth_vcf, ploidy, truth_sample)
    H_pred, block_ids = parse_pred_vcf_to_H_and_blocks(pred_vcf, ploidy, pred_sample)

    cols = valid_overlap_columns(H_true, H_pred)
    if cols.size == 0:
        raise ValueError("No overlapping phased sites between truth and prediction.")

    ver = vector_error_wrapper(H_true[:, cols], H_pred[:, cols], block_ids[cols])
    return float(ver), int(cols.size)


# -----------------------------
# Extract frag from BAM (for eval MEC)
# -----------------------------
def resolve_extracthairs_binary() -> str:
    """
    Resolve extractHAIRS binary path.
    Preference:
      1) bundled repo path (editable installs)
      2) PATH lookup
    """
    # Try bundled repo structure relative to this file (works if running from repo/editable)
    this_dir = Path(__file__).resolve().parent            # .../src/phapcompass
    repo_root_guess = this_dir.parent.parent              # .../pHapCompass
    bundled = repo_root_guess / "third_party" / "extract_poly" / "build" / "extractHAIRS"
    if bundled.is_file():
        return str(bundled)

    # Fall back to PATH
    resolved = shutil.which("extractHAIRS")
    if resolved:
        return resolved

    raise FileNotFoundError(
        "extractHAIRS binary not found.\n"
        f"Tried bundled path: {bundled}\n"
        "and PATH lookup for 'extractHAIRS'."
    )


def build_frag_from_bam_if_needed(
    bam_path: str,
    vcf_path: str,
    frag_path: str,
    ploidy: int,
    mbq: int = 13) -> None:
    """
    Generate a fragment file from BAM+VCF using extractHAIRS.
    Used only for MEC/geom-MEC evaluation when --bam is provided.
    """
    bin_path = resolve_extracthairs_binary()

    cmd = [
        bin_path,
        "--bam", bam_path,
        "--vcf", vcf_path,
        "--out", frag_path,
    ]
    if mbq != 13:
        cmd.extend(["--mbq", str(mbq)])

    logging.info("Running extractHAIRS: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


# -----------------------------
# Read list filtering for MEC on phased subset
# -----------------------------
def filter_read_list_to_columns(read_list: list, keep_cols: np.ndarray) -> list:
    """
    Filter and remap read_list to a subset of columns.

    read_list format:
      [pos0, alleles0, pos1, alleles1, ...]

    keep_cols: array of original column indices to keep.
    Returns a new read_list where positions are remapped to 0..len(keep_cols)-1.
    """
    keep_set = set(int(x) for x in keep_cols.tolist())
    col_map = {int(old): int(new) for new, old in enumerate(keep_cols.tolist())}

    new_list = []
    for i in range(0, len(read_list), 2):
        pos = np.asarray(read_list[i], dtype=int)
        alle = np.asarray(read_list[i + 1], dtype=int)

        mask = np.array([p in keep_set for p in pos], dtype=bool)
        if mask.sum() == 0:
            continue

        pos_f = pos[mask]
        alle_f = alle[mask]

        # remap positions
        pos_r = np.array([col_map[int(p)] for p in pos_f], dtype=int)

        # keep order by position
        order = np.argsort(pos_r)
        pos_r = pos_r[order]
        alle_f = alle_f[order]

        new_list.append(pos_r.tolist())
        new_list.append(alle_f.tolist())

    return new_list


# -----------------------------
# MEC evaluation (from pred VCF + frag)
# -----------------------------
def compute_mec_from_pred_vcf_and_frag(
    pred_vcf: str,
    frag_path: str,
    ploidy: int,
    pred_sample: str = "SAMPLE",
) -> Tuple[float, int]:
    """
    Compute MEC given predicted VCF and fragment file.
    MEC computed on columns where prediction is phased/defined.

    Returns (mec, evaluated_sites).
    """
    H_pred, block_ids = parse_pred_vcf_to_H_and_blocks(pred_vcf, ploidy, pred_sample)

    # Use only fully defined predicted columns (phased)
    ok_pred = ~np.any(np.isnan(H_pred), axis=0)
    cols = np.where(ok_pred)[0]
    if cols.size == 0:
        raise ValueError("No phased columns in predicted VCF for MEC evaluation.")

    # Build read matrix from frag
    cfg = InputConfigSparse(data_path=frag_path, genotype_path=np.empty((0, 0)), ploidy=ploidy)  # genotype not needed for parsing
    frag = SparseFragment(cfg)
    M = frag.csr
    read_list = build_read_list_from_M(M)

    # Filter/remap reads to the phased subset
    read_list_sub = filter_read_list_to_columns(read_list, cols)

    # Subset haplotypes and block ids
    H_sub = H_pred[:, cols]
    block_sub = block_ids[cols]

    mec = mec_full(H_sub, block_sub, read_list_sub, probabalistic=False)
    return float(mec), int(cols.size)


def compute_geom_mec_from_pred_vcf_and_frag(pred_vcf: str, frag_path: str, ploidy: int, pred_sample: str = "SAMPLE") -> Tuple[float, int]:
    """
    Compute geometric MEC given predicted VCF and fragment file.
    Returns (geom_mec, evaluated_sites).
    """
    H_pred, block_ids = parse_pred_vcf_to_H_and_blocks(pred_vcf, ploidy, pred_sample)

    ok_pred = ~np.any(np.isnan(H_pred), axis=0)
    cols = np.where(ok_pred)[0]
    if cols.size == 0:
        raise ValueError("No phased columns in predicted VCF for geometric MEC evaluation.")

    cfg = InputConfigSparse(data_path=frag_path, genotype_path=np.empty((0, 0)), ploidy=ploidy)
    frag = SparseFragment(cfg)
    data_matrix = frag.csr
    read_list = build_read_list_from_M(data_matrix)

    # read_list_sub = filter_read_list_to_columns(read_list, cols)

    # H_sub = H_pred[:, cols]
    # block_sub = block_ids[cols]

    gmec = mec_full_geometric_penalty(H_pred, block_ids, read_list, probabalistic=False)
    return float(gmec), int(cols.size)





def _setup_logging(verbose: bool, log_level: str) -> None:
    level = "DEBUG" if verbose else log_level
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser( prog="phapcompass-eval", description="Evaluation utilities for pHapCompass: VER, MEC, and geometric MEC.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")
    parser.add_argument("--log-level", type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level.")

    sub = parser.add_subparsers(dest="cmd", required=True)

    # ---------------- VER ----------------
    p_ver = sub.add_parser("ver", help="Compute VER between truth VCF and predicted VCF.")
    p_ver.add_argument("--truth-vcf", required=True, help="Ground-truth VCF (.vcf or .vcf.gz).")
    p_ver.add_argument("--pred-vcf", required=True, help="Predicted VCF (.vcf or .vcf.gz).")
    p_ver.add_argument("--ploidy", type=int, default=None, help="Optional ploidy override. If omitted, inferred from VCF.")
    p_ver.add_argument("--truth-sample", type=str, default=None, help="Sample name in truth VCF. If omitted, uses first sample.")
    p_ver.add_argument("--pred-sample", type=str, default="SAMPLE", help='Sample name in predicted VCF (default: "SAMPLE").')

    # ---------------- MEC ----------------
    p_mec = sub.add_parser("mec", help="Compute MEC for predicted VCF given reads (frag or bam).")
    p_mec.add_argument("--pred-vcf", required=True, help="Predicted VCF (.vcf or .vcf.gz).")
    p_mec.add_argument("--ploidy", type=int, default=None, help="Optional ploidy override. If omitted, inferred from predicted VCF.")
    p_mec.add_argument("--pred-sample", type=str, default="SAMPLE", help='Sample name in predicted VCF (default: "SAMPLE").')

    g_mec = p_mec.add_mutually_exclusive_group(required=True)
    g_mec.add_argument("--frag", type=str, default=None, help="Fragment file (.frag).")
    g_mec.add_argument("--bam", type=str, default=None, help="BAM file (if frag not provided).")

    # If BAM is used, we need the VCF used for extraction (unphased input VCF)
    p_mec.add_argument("--vcf", type=str, default=None, help="Input VCF used for extractHAIRS when --bam is provided (required with --bam).")
    p_mec.add_argument("--mbq", type=int, default=13, help="Min base quality for extractHAIRS (default: 13).")

    # ---------------- GEO-MEC ----------------
    p_gm = sub.add_parser("geom-mec", help="Compute geometric MEC for predicted VCF given reads (frag or bam).")
    p_gm.add_argument("--pred-vcf", required=True, help="Predicted VCF (.vcf or .vcf.gz).")
    p_gm.add_argument("--ploidy", type=int, default=None, help="Optional ploidy override. If omitted, inferred from predicted VCF.")
    p_gm.add_argument("--pred-sample", type=str, default="SAMPLE", help='Sample name in predicted VCF (default: "SAMPLE").')

    g_gm = p_gm.add_mutually_exclusive_group(required=True)
    g_gm.add_argument("--frag", type=str, default=None, help="Fragment file (.frag).")
    g_gm.add_argument("--bam", type=str, default=None, help="BAM file (if frag not provided).")

    p_gm.add_argument("--vcf", type=str, default=None, help="Input VCF used for extractHAIRS when --bam is provided (required with --bam).")
    p_gm.add_argument("--mbq", type=int, default=13, help="Min base quality for extractHAIRS (default: 13).")

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _setup_logging(args.verbose, args.log_level)

    if args.cmd == "ver":
        truth_vcf = args.truth_vcf
        pred_vcf = args.pred_vcf

        ploidy = args.ploidy
        if ploidy is None:
            # Prefer predicted VCF for inference (matches output)
            ploidy = infer_ploidy_from_vcf_safe(pred_vcf) or infer_ploidy_from_vcf_safe(truth_vcf)
            if ploidy is None:
                raise RuntimeError("Could not infer ploidy from VCF(s). Please provide --ploidy.")

        ver, n_eval = compute_ver_from_vcfs(
            truth_vcf=truth_vcf,
            pred_vcf=pred_vcf,
            ploidy=ploidy,
            truth_sample=args.truth_sample,
            pred_sample=args.pred_sample,
        )
        print(f"VER={ver:.6f}  (evaluated_sites={n_eval})")
        return

    if args.cmd in ("mec", "geom-mec"):
        pred_vcf = args.pred_vcf

        ploidy = args.ploidy
        if ploidy is None:
            ploidy = infer_ploidy_from_vcf_safe(pred_vcf)
            if ploidy is None:
                raise RuntimeError("Could not infer ploidy from predicted VCF. Please provide --ploidy.")

        # Determine frag path
        frag_path = None
        tmpdir = None

        if args.frag is not None:
            frag_path = args.frag
        else:
            # BAM mode => must have --vcf for extractHAIRS
            if args.vcf is None:
                raise ValueError("When using --bam, you must also provide --vcf (the VCF used for extractHAIRS).")

            tmpdir = tempfile.mkdtemp(prefix="phapcompass_eval_")
            frag_path = str(Path(tmpdir) / "eval_tmp.frag")

            build_frag_from_bam_if_needed(
                bam_path=args.bam,
                vcf_path=args.vcf,
                frag_path=frag_path,
                ploidy=ploidy,
                mbq=args.mbq,
            )

        if args.cmd == "mec":
            mec, n_eval = compute_mec_from_pred_vcf_and_frag(
                pred_vcf=pred_vcf,
                frag_path=frag_path,
                ploidy=ploidy,
                pred_sample=args.pred_sample,
            )
            print(f"MEC={mec:.6f}  (evaluated_sites={n_eval})")
        else:
            gmec, n_eval = compute_geom_mec_from_pred_vcf_and_frag(
                pred_vcf=pred_vcf,
                frag_path=frag_path,
                ploidy=ploidy,
                pred_sample=args.pred_sample,
            )
            print(f"GEOM_MEC={gmec:.6f}  (evaluated_sites={n_eval})")

        return


if __name__ == "__main__":
    main()
