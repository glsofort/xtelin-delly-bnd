#!/usr/bin/env python3
# ============================================================================
# extract_bnd.py — Delly v1.7.2 BND breakpoint coordinate extractor
# ============================================================================
# Extracts cross-chromosomal BND records from a Delly VCF/BCF and prints
# breakpoint coordinates in a structured, human-readable table.
#
# Target use case : SLC25A13 IVS16ins3kb breakpoint
# Expected BND    : chr7:95751080  <->  chr6:117996638  (GT=0/1, PASS)
# Tool            : Delly v1.7.2  (EMBL.DELLYv1.7.2)
# Reference       : hs37d5 bundle  (GRCh37 + decoys)
# ----------------------------------------------------------------------------
#
# Dependencies: python3 (stdlib only — no bcftools/pysam required)
# Usage:
#   python3 extract_bnd.py delly_out.vcf
#   python3 extract_bnd.py -t 6,7 delly_out.vcf          # only chr6<->chr7
#   python3 extract_bnd.py --target 7,95751080 delly.vcf  # match target by pos
#   python3 extract_bnd.py --pass-only delly.vcf
#   python3 extract_bnd.py delly_out.bcf                   # needs bcftools
#
# Outputs:
#   1. STDOUT — formatted summary table
#   2. JSON file  (.bnd.json) — machine-readable breakpoint list
#   3. VCF file   (.target_bnd.vcf) — filtered VCF (optional, --save-vcf)
# ============================================================================

import argparse
import csv
import json
import os
import re
import sys
from collections import OrderedDict

# ---------------------------------------------------------------------------
# Constants — known gene / variant lookup table
#   key: canonical pair_key "chr_lo:pos_lo-chr_hi:pos_hi"
#   ±2000bp tolerance used for matching
# ---------------------------------------------------------------------------
GENE_LOOKUP = [
    {
        "gene":       "SLC25A13",
        "variant":    "IVS16ins3kb",
        "chrA":       "7",   "posA":   95751080,
        "chrB":       "6",   "posB":  117996638,
        "omim_gene":  "603859",   # SLC25A13 gene
        "omim_disease": "605814", # Citrin deficiency
        "local_db":   "README_pipeline.md",
        "notes":      "NUS1 cDNA retrotransposal insertion (chr6→chr7, reverse-complement)",
        "tolerance":  2000,
    },
]

# Backward-compat alias for TARGET CHECK block
TARGET_SLC25A13 = GENE_LOOKUP[0]


def lookup_gene_info(rec):
    """Return dict with gene/omim/local_db for a BND record, or all '-' if unknown."""
    for entry in GENE_LOOKUP:
        tol = entry["tolerance"]
        # Match either orientation
        m1 = (rec.chrom == entry["chrA"] and rec.chr2 == entry["chrB"]
              and in_range(rec.pos,  entry["posA"], tol)
              and in_range(rec.pos2, entry["posB"], tol))
        m2 = (rec.chrom == entry["chrB"] and rec.chr2 == entry["chrA"]
              and in_range(rec.pos,  entry["posB"], tol)
              and in_range(rec.pos2, entry["posA"], tol))
        if m1 or m2:
            return {
                "gene":      entry["gene"],
                "omim":      entry["omim_gene"],
                "omim_full": f"{entry['omim_gene']} (gene) / {entry['omim_disease']} (disease)",
                "local_db":  entry["local_db"],
                "variant":   entry["variant"],
            }
    return {"gene": "-", "omim": "-", "omim_full": "-", "local_db": "-", "variant": "-"}

# Delly ALT bracket direction mapping
# C[X:Y[  ->  opening [ on right   =  5' of mate connects here (CT=5toX)
# C]X:Y]  ->  closing ] on right   =  3' of mate connects here (CT=3toX)
# [X:Y[C  ->  opening [ on left    =  5' of this connects to mate
# ]X:Y]C  ->  closing ] on left    =  3' of this connects to mate
ALT_PATTERN = re.compile(r"([\[\]])([^:\[\]]+):(\d+)([\[\]])")


# ============================================================================
# VCF record parser
# ============================================================================
class BNDRecord:
    """Lightweight parser for a single BND VCF line."""

    def __init__(self, line, header_fields):
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 10:
            raise ValueError(f"Too few columns ({len(cols)}) in line: {line[:80]}")

        self.chrom   = cols[0]
        self.pos     = int(cols[1])
        self.varid   = cols[2]
        self.ref     = cols[3]
        self.alt     = cols[4]
        self.qual    = cols[5]
        self.filter  = cols[6]
        self.info    = self._parse_info(cols[7])
        self.format  = cols[8]
        self.samples = self._parse_samples(cols[9:], self.format)

        # Delly cross-chr BND uses CHR2/POS2 in INFO
        self.chr2 = self.info.get("CHR2", None)
        self.pos2 = int(self.info.get("POS2", 0)) if self.info.get("POS2") else None
        self.ct   = self.info.get("CT", None)

        # Genotype (first sample)
        self.gt = self.samples[0].get("GT", "./.") if self.samples else "./."

    # ---- helpers ----
    @staticmethod
    def _parse_info(info_str):
        d = {}
        for field in info_str.split(";"):
            if "=" in field:
                k, v = field.split("=", 1)
                d[k] = v
            else:
                d[field] = True  # flag
        return d

    @staticmethod
    def _parse_samples(sample_cols, fmt_str):
        fmts = fmt_str.split(":")
        out = []
        for col in sample_cols:
            vals = col.split(":")
            d = {}
            for i, f in enumerate(fmts):
                d[f] = vals[i] if i < len(vals) else "."
            out.append(d)
        return out

    def is_cross_chr(self):
        return self.chr2 is not None and self.chr2 != self.chrom

    def pair_key(self):
        """Canonical (chr_lo, pos_lo, chr_hi, pos_hi) key for deduplication."""
        a = (self.chrom, self.pos)
        b = (self.chr2, self.pos2)
        if a > b:
            a, b = b, a
        return f"{a[0]}:{a[1]}-{b[0]}:{b[1]}"

    def passes(self):
        return self.filter == "PASS"

    # ---- pretty output ----
    def summary(self):
        pe  = self.info.get("PE", ".")
        sr  = self.info.get("SR", ".")
        mapq = self.info.get("MAPQ", ".")
        cip = self.info.get("CIPOS", ".")
        gq  = self.samples[0].get("GQ", ".") if self.samples else "."
        return (
            f"{self.chrom:>4}:{self.pos:<10}  <->  "
            f"{self.chr2:>4}:{self.pos2:<10}  "
            f"PE={pe:>3}  SR={sr:>2}  MAPQ={mapq:>3}  "
            f"CT={self.ct:<5}  GT={self.gt:<4}  GQ={gq:<6}  "
            f"FILTER={self.filter}"
        )


# ============================================================================
# VCF reader (supports plain .vcf; .bcf auto-converted via bcftools)
# ============================================================================
def iter_vcf_lines(path):
    """Yield non-header lines from a VCF/BCF file."""
    if not os.path.isfile(path):
        sys.exit(f"ERROR: input file not found: {path}")

    # BCF → pipe through bcftools
    if path.endswith(".bcf") or path.endswith(".bcf.gz"):
        import subprocess
        try:
            p = subprocess.Popen(
                ["bcftools", "view", path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, bufsize=1,
            )
        except FileNotFoundError:
            sys.exit("ERROR: reading .bcf requires 'bcftools' in PATH")
        for line in p.stdout:
            yield line.rstrip("\n")
        p.wait()
        if p.returncode != 0:
            err = p.stderr.read()
            sys.exit(f"ERROR bcftools view {path}: {err}")
        return

    # Plain VCF (optionally gzipped)
    if path.endswith(".gz"):
        import gzip
        fh = gzip.open(path, "rt")
    else:
        fh = open(path, "r")
    with fh:
        for line in fh:
            yield line.rstrip("\n")


# ============================================================================
# Filter helpers
# ============================================================================
def chr_pair_match(rec, pair_str):
    """True if BND connects the two chromosomes in pair_str (order-insensitive)."""
    c1, c2 = pair_str.split(",")
    c1, c2 = c1.strip(), c2.strip()
    return (rec.chrom == c1 and rec.chr2 == c2) or (rec.chrom == c2 and rec.chr2 == c1)


def in_range(pos, target, tolerance):
    return abs(int(pos) - int(target)) <= int(tolerance)


# ============================================================================
# Main
# ============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Extract Delly v1.7.2 BND breakpoint coordinates from VCF/BCF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s delly_out.vcf
  %(prog)s -t 6,7 delly_out.vcf
  %(prog)s -t 6,7 --pass-only delly_out.vcf
  %(prog)s --target 7,95751080 --tolerance 500 delly_out.vcf
  %(prog)s --save-vcf filtered.vcf delly_out.bcf
        """,
    )
    ap.add_argument("vcf", help="Delly output .vcf or .bcf")
    ap.add_argument(
        "-t", "--chr-pair", default=None,
        help="Only keep BNDs between these two chromosomes, e.g. '6,7'",
    )
    ap.add_argument(
        "--pass-only", action="store_true",
        help="Only keep records with FILTER=PASS",
    )
    ap.add_argument(
        "--target", default=None, metavar="CHR,POS",
        help="Match one breakpoint side to this coordinate (e.g. '7,95751080')",
    )
    ap.add_argument(
        "--tolerance", type=int, default=500,
        help="±bp tolerance for --target coordinate match (default: 500)",
    )
    ap.add_argument(
        "--all-bnd", action="store_true",
        help="Do NOT restrict to cross-chromosomal BNDs (keep intra-chr too)",
    )
    ap.add_argument(
        "--save-vcf", metavar="PATH", default=None,
        help="Write filtered VCF to this path",
    )
    ap.add_argument(
        "--sample-name", default=None,
        help="Override sample name in output (default: derived from VCF header)",
    )
    ap.add_argument(
        "--no-target-check", action="store_true",
        help="Skip SLC25A13 IVS16ins3kb expected-breakpoint sanity check",
    )
    args = ap.parse_args()

    # ------------------------------------------------------------------
    # 1. Parse VCF
    # ------------------------------------------------------------------
    header_fields = None
    header_lines  = []
    all_records   = []
    cross_records = []

    for line in iter_vcf_lines(args.vcf):
        if line.startswith("##"):
            header_lines.append(line)
            continue
        if line.startswith("#CHROM"):
            header_fields = line.split("\t")
            header_lines.append(line)
            continue
        if not line.strip():
            continue
        try:
            rec = BNDRecord(line, header_fields)
        except Exception as e:
            print(f"WARN: skipping malformed line: {e}", file=sys.stderr)
            continue
        all_records.append(rec)
        if rec.is_cross_chr():
            cross_records.append(rec)

    if not all_records:
        sys.exit(f"ERROR: no records found in {args.vcf}")

    # ------------------------------------------------------------------
    # 2. Apply filters
    # ------------------------------------------------------------------
    records = cross_records if not args.all_bnd else all_records

    if args.pass_only:
        records = [r for r in records if r.passes()]

    if args.chr_pair:
        records = [r for r in records if chr_pair_match(r, args.chr_pair)]

    if args.target:
        t_chr, t_pos = args.target.split(",")
        t_chr, t_pos = t_chr.strip(), int(t_pos.strip())
        records = [
            r for r in records
            if (r.chrom == t_chr and in_range(r.pos,  t_pos, args.tolerance))
            or (r.chr2  == t_chr and in_range(r.pos2, t_pos, args.tolerance))
        ]

    # ------------------------------------------------------------------
    # 3. Prepare output path
    # ------------------------------------------------------------------
    base = os.path.splitext(os.path.basename(args.vcf))[0]
    tsv_path = base + ".bnd.tsv"

    # ------------------------------------------------------------------
    # 4. Print pretty summary
    # ------------------------------------------------------------------
    sample_name = args.sample_name or base
    n_total    = len(all_records)
    n_cross    = len(cross_records)
    n_filtered = len(records)

    print("=" * 90)
    print(f"  Delly v1.7.2 BND Breakpoint Report")
    print(f"  Input      : {args.vcf}")
    print(f"  Sample     : {sample_name}")
    print(f"  Date       : {os.popen('date').read().strip()}")
    print("-" * 90)
    print(f"  Total BND records in VCF       : {n_total}")
    print(f"  Cross-chromosomal BNDs         : {n_cross}")
    print(f"  After user filters             : {n_filtered}")
    if args.chr_pair:
        print(f"  --chr-pair                     : {args.chr_pair}")
    if args.pass_only:
        print(f"  --pass-only                    : enabled")
    if args.target:
        print(f"  --target                       : {args.target}  (±{args.tolerance} bp)")
    print("=" * 90)

    if not records:
        print("\n  [NO MATCHING BND RECORDS]\n")
    else:
        print(f"\n  {'CHR_A':>5}  {'POS_A':>10}    {'CHR_B':>5}  {'POS_B':>10}  "
              f"{'PE':>4} {'SR':>3} {'MAPQ':>4} {'CT':>6} {'GT':>4} {'GQ':>5} FILTER")
        print("  " + "-" * 84)
        for r in records:
            pe   = r.info.get("PE", ".")
            sr   = r.info.get("SR", ".")
            mapq = r.info.get("MAPQ", ".")
            gq   = r.samples[0].get("GQ", ".") if r.samples else "."
            print(f"  {r.chrom:>5}  {r.pos:>10}  <->  "
                  f"{r.chr2:>5}  {r.pos2:>10}  "
                  f"{pe:>4} {sr:>3} {mapq:>4} {r.ct:>6} {r.gt:>4} {gq:>5} {r.filter}")

    # ------------------------------------------------------------------
    # 5. SLC25A13 target check
    # ------------------------------------------------------------------
    if not args.no_target_check:
        t = TARGET_SLC25A13
        print()
        print("-" * 90)
        print(f"  TARGET CHECK — {t['gene']} {t['variant']}")
        print(f"  Expected : {t['chrA']}:{t['posA']}  <->  {t['chrB']}:{t['posB']}")
        print(f"            (reverse-complement NUS1 cDNA insertion)")
        print("-" * 90)

        found_slc = []
        for r in cross_records:
            # Match in either orientation, ±2000 bp tolerance on pos
            orient1 = (r.chrom == t["chrA"] and r.chr2 == t["chrB"]
                       and in_range(r.pos,  t["posA"], 2000)
                       and in_range(r.pos2, t["posB"], 2000))
            orient2 = (r.chrom == t["chrB"] and r.chr2 == t["chrA"]
                       and in_range(r.pos,  t["posB"], 2000)
                       and in_range(r.pos2, t["posA"], 2000))
            if orient1 or orient2:
                found_slc.append(r)

        if found_slc:
            print(f"  ✓ FOUND {len(found_slc)} candidate(s)")
            for r in found_slc:
                match_gt   = (r.gt in ("0/1", "1/0"))
                match_pass = r.passes()
                flag = ""
                if not match_gt:   flag += " GT!=0/1"
                if not match_pass: flag += " FILTER!=PASS"
                tag = "  ★★★ TARGET" if match_gt and match_pass else f"  ⚠ WARNING:{flag}"
                print(f"    {r.summary()}{tag}")
        else:
            print("  ✗ NOT FOUND within ±2000 bp tolerance")
        print()

    # ------------------------------------------------------------------
    # 6. Write JSON  (6 keys: gene / mutation_position / gt / local_database / g_score / omim)
    # ------------------------------------------------------------------
    json_path = base + ".bnd.json"
    json_rows = []
    for r in records:
        info = lookup_gene_info(r)
        gq   = r.samples[0].get("GQ", "-") if r.samples else "-"
        g_score = gq if gq not in (".", "") else "-"
        json_rows.append({
            "gene":             info["gene"],
            "mutation_position": f"{r.chrom}:{r.pos}/{r.chr2}:{r.pos2}",
            "gt":               r.gt,
            "local_database":   info["local_db"],
            "g_score":          g_score,
            "omim":             info["omim"],
        })
    with open(json_path, "w") as fh:
        json.dump(json_rows, fh, indent=2, ensure_ascii=False)
    print(f"  [JSON written] {json_path}  ({len(records)} records)")

    # ------------------------------------------------------------------
    # 7. Save filtered VCF (optional)
    # ------------------------------------------------------------------
    if args.save_vcf:
        with open(args.save_vcf, "w") as fh:
            for hl in header_lines:
                fh.write(hl + "\n")
            for r in records:
                # Reconstruct original line from parsed record
                fh.write("\t".join([
                    r.chrom, str(r.pos), r.varid, r.ref, r.alt,
                    r.qual, r.filter,
                    ";".join(f"{k}={v}" if v is not True else k
                             for k, v in r.info.items()),
                    r.format,
                ] + [":".join(r.samples[0].get(f, ".") for f in r.format.split(":"))]) + "\n")
        print(f"  [VCF written]  {args.save_vcf}")

    print("=" * 90)


if __name__ == "__main__":
    main()
