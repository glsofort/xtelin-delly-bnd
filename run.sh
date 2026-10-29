#!/usr/bin/env bash
# ============================================================================
# xtelin Delly v1.7.2 BND detection Nextflow pipeline runner
# ============================================================================
# Usage:
#   bash run.sh                              # default sample 170TFW
#   bash run.sh /path/to/other.bam           # single sample
#   bash run.sh "/path/to/samples/*.bam"     # batch (glob)
#   RESUME=1 bash run.sh                     # resume failed run
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- defaults ---
DEFAULT_BAM="${SCRIPT_DIR}/../170TFW.deduped.bam"
DEFAULT_REF="/mnt/storage/bioinformatics-data/kian/hg19/hs37d5.fa"
EXTRACT_SCRIPT="${SCRIPT_DIR}/scripts/extract_bnd.py"
OUTDIR="${SCRIPT_DIR}/nf_results"

BAM_INPUT="${1:-$DEFAULT_BAM}"
REF="${DELLY_REF:-$DEFAULT_REF}"

# sample name: strip path + .deduped.bam / .bam suffixes
SAMPLE_NAME="$(basename "$BAM_INPUT" | sed 's/\.deduped\.bam$//; s/\.bam$//')"

echo "============================================"
echo " xtelin Delly v1.7.2 BND Nextflow pipeline"
echo "============================================"
echo " BAM         : $BAM_INPUT"
echo " FASTA       : $REF"
echo " extract.py  : $EXTRACT_SCRIPT"
echo " output      : $OUTDIR"
echo "--------------------------------------------"

# --- sanity checks ---
[[ -f "$REF" ]]          || { echo "ERROR: FASTA not found: $REF"; exit 1; }
[[ -f "$EXTRACT_SCRIPT" ]] || { echo "ERROR: extract_bnd.py not found: $EXTRACT_SCRIPT"; exit 1; }
[[ "$BAM_INPUT" == *.bam && ! -f "$BAM_INPUT" ]] \
  && { echo "ERROR: BAM not found: $BAM_INPUT"; exit 1; }

mkdir -p "$OUTDIR" nf_reports

# --- run ---
CMD=(
  nextflow run main.nf
    -profile docker
    --bam  "$BAM_INPUT"
    --ref  "$REF"
    --extract-script "$EXTRACT_SCRIPT"
    --outdir "$OUTDIR"
    -with-trace "nf_reports/${SAMPLE_NAME}.json"
)
[[ "${RESUME:-0}" == "1" ]] && CMD+=( -resume )

"${CMD[@]}"

echo ""
echo "============================================"
echo " Pipeline complete"
echo "--------------------------------------------"
echo " results : $OUTDIR"
echo "   delly/  - genome-wide BCF + VCF"
echo "   target/ - chr6<->chr7 BND  (.vcf + .bnd.json)"
echo "   ref/    - FASTA + .fai"
echo " reports : nf_reports/"
echo "============================================"
