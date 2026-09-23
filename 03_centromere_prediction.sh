#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# 5-kb window-based tandem repeat (TR) coverage analysis
# for chromosomes Chr1-Chr8
#
# Software:
#   TRF 4.10.0rc2
#   bedtools
#
# Input:
#   genome.fa
#
# Output:
#   TRF_Chr1-8_merged.bed
#   Chr1-8.lengths
#   windows_5kb.bed
#   windows_5kb_TR_raw.tsv
#   windows_5kb_TR_coverage.tsv
#   Chr1-8_5kb_TR_summary.tsv
# ============================================================


# -------------------------
# Input and software paths
# -------------------------

GENOME="genome.fa"

TRF="/opt/hdd/users/only/.conda/envs/quartet/bin/trf"
BEDTOOLS="/usr/bin/bedtools"

WINDOW_SIZE=5000


# -------------------------
# Check input and software
# -------------------------

if [[ ! -f "${GENOME}" ]]; then
    echo "ERROR: Genome file not found: ${GENOME}" >&2
    exit 1
fi

if [[ ! -x "${TRF}" ]]; then
    echo "ERROR: TRF executable not found: ${TRF}" >&2
    exit 1
fi

if [[ ! -x "${BEDTOOLS}" ]]; then
    echo "ERROR: bedtools executable not found: ${BEDTOOLS}" >&2
    exit 1
fi


# ============================================================
# Step 1. Identify tandem repeats using TRF
# ============================================================

echo "[1/5] Running TRF..."

"${TRF}" \
    "${GENOME}" \
    2 7 7 80 10 50 2000 \
    -d -h


# ============================================================
# Step 2. Generate TRF BED file
# ============================================================

echo "[2/5] Preparing TRF BED file..."

# This step assumes that the TRF results have been converted
# and merged into:
#     TRF_Chr1-8_merged.bed
#
# The BED file should contain:
#     chromosome  start  end
#
# Only Chr1-Chr8 are retained for downstream analysis.

if [[ ! -f "TRF_Chr1-8_merged.bed" ]]; then
    echo "ERROR: TRF_Chr1-8_merged.bed was not found." >&2
    echo "Please generate the merged TRF BED file before continuing." >&2
    exit 1
fi


# ============================================================
# Step 3. Generate non-overlapping 5-kb windows
# ============================================================

echo "[3/5] Generating 5-kb windows..."

awk 'BEGIN{OFS="\t"} $1 ~ /^Chr[1-8]$/ {print $1,$2}' \
    "${GENOME}.fai" \
    > Chr1-8.lengths

"${BEDTOOLS}" makewindows \
    -g Chr1-8.lengths \
    -w "${WINDOW_SIZE}" \
    > windows_5kb.bed

echo "Number of 5-kb windows:"
wc -l windows_5kb.bed


# ============================================================
# Step 4. Calculate TR coverage for each 5-kb window
# ============================================================

echo "[4/5] Calculating TR coverage..."

"${BEDTOOLS}" coverage \
    -a windows_5kb.bed \
    -b TRF_Chr1-8_merged.bed \
    > windows_5kb_TR_raw.tsv


# Format output
awk 'BEGIN{
        OFS="\t";
        print "Chr","Start","End","Window_bp","TR_covered_bp","TR_coverage"
     }
     {
        len=$3-$2;
        print $1,$2,$3,len,$6,$7
     }' \
    windows_5kb_TR_raw.tsv \
    > windows_5kb_TR_coverage.tsv


# ============================================================
# Step 5. Summarize TR coverage for each chromosome
# ============================================================

echo "[5/5] Generating chromosome-level summary..."

awk -F'\t' '
NR==1 {next}

{
    chr=$1
    cov=$6

    n[chr]++
    sum[chr]+=cov

    if(cov>=0.50) ge50[chr]++
    if(cov>=0.60) ge60[chr]++
    if(cov>=0.70) ge70[chr]++
    if(cov>=0.80) ge80[chr]++
    if(cov>=0.90) ge90[chr]++
    if(cov>=0.95) ge95[chr]++

    if(!(chr in min) || cov<min[chr])
        min[chr]=cov

    if(!(chr in max) || cov>max[chr])
        max[chr]=cov
}

END {
    print "Chr","Windows","Mean","Min","Max",
          ">=50%",">=60%",">=70%",">=80%",">=90%",">=95%"

    for(i=1;i<=8;i++) {

        chr="Chr" i

        printf "%s\t%d\t%.4f\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\n",
            chr,
            n[chr],
            sum[chr]/n[chr],
            min[chr],
            max[chr],
            ge50[chr],
            ge60[chr],
            ge70[chr],
            ge80[chr],
            ge90[chr],
            ge95[chr]
    }
}' \
    windows_5kb_TR_coverage.tsv \
    > Chr1-8_5kb_TR_summary.tsv


# ============================================================
# Final output
# ============================================================

echo
echo "Analysis completed successfully."
echo
echo "Output files:"
echo "  Chr1-8.lengths"
echo "  windows_5kb.bed"
echo "  windows_5kb_TR_raw.tsv"
echo "  windows_5kb_TR_coverage.tsv"
echo "  Chr1-8_5kb_TR_summary.tsv"
echo
echo "Chromosome-level summary:"
echo

cat Chr1-8_5kb_TR_summary.tsv
