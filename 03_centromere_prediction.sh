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

# Select windows with gene coverage = 0 and TR coverage >= 70%
awk -F'\t' '
BEGIN{OFS="\t"}

NR==FNR {
    if(NR>1)
        gene[$1 FS $2 FS $3]=$6
    next
}

NR>1 {
    key=$1 FS $2 FS $3
    if($6>=0.70 && gene[key]==0)
        print $1,$2,$3,$6,gene[key]
}
' \
windows_5kb_gene_coverage.tsv \
windows_5kb_TR_coverage.tsv \
> TR70_gene0_5kb.tsv

# Merge consecutive windows only when the end of one window exactly matches the start of the next
awk -F'\t' '
BEGIN{OFS="\t"}

NR==1 {
    chr=$1
    start=$2
    end=$3
    n=1
    tr_sum=$4
    next
}

{
    if($1==chr && $2==end) {
        end=$3
        n++
        tr_sum+=$4
    }
    else {
        print chr,start,end,end-start,n,tr_sum/n
        chr=$1
        start=$2
        end=$3
        n=1
        tr_sum=$4
    }
}

END {
    print chr,start,end,end-start,n,tr_sum/n
}
' TR70_gene0_5kb.tsv \
> TR70_gene0_5kb_merged.tsv

#View the start positions and lengths of TRs within all candidate regions in chromosome order.
only@admin1-Rack-Server:~/TMTM_10kb_centromere_analysis$ awk '$1=="ChrX" && $2 < xxxxxxxx && $3 > xxxxxxxx {
    print $1, $2, $3, $3-$2
}' OFS='\t' TRF_merged.bed


## Centromere prediction using quarTeT
set -euo pipefail

project_root=${1:-/opt/hdd/users/BH10160/tianmu_tiemu_centromere_quartet}
threads=${2:-8}

apptainer=${project_root}/../software/apptainer/bin/apptainer
if [[ ! -x "${apptainer}" ]]; then
    apptainer=/opt/hdd/users/BH10160/software/apptainer/bin/apptainer
fi

image=${project_root}/tools/quartet-bio-1.3.1.sif
input_dir=${project_root}/work/chromosome_inputs
genome=${input_dir}/tianmu_tiemu.chromosomes.fa
repeat_gff=${input_dir}/tianmu_tiemu.chromosomes.repeat.gff3
gene_gff=${input_dir}/tianmu_tiemu.chromosomes.gene.gff3
output_dir=${project_root}/results/centrominer_default
log_dir=${project_root}/logs
tmp_dir=${project_root}/work/tmp

mkdir -p "${output_dir}" "${log_dir}" "${tmp_dir}"

for required in "${apptainer}" "${image}" "${genome}" "${repeat_gff}" "${gene_gff}"; do
    if [[ ! -s "${required}" ]]; then
        echo "Missing or empty required file: ${required}" >&2
        exit 1
    fi
done

if find "${output_dir}" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
    echo "Output directory is not empty: ${output_dir}" >&2
    echo "Move or remove the existing run before starting a new one." >&2
    exit 1
fi

{
    echo "run_started=$(date -Iseconds)"
    echo "hostname=$(hostname)"
    echo "threads=${threads}"
    echo "image=${image}"
    sha256sum "${image}" "${genome}" "${repeat_gff}" "${gene_gff}"
    "${apptainer}" exec --cleanenv "${image}" quartet -h
    "${apptainer}" exec --cleanenv "${image}" blastn -version
    "${apptainer}" exec --cleanenv "${image}" Rscript --version
} > "${log_dir}/software_and_input_manifest.txt" 2>&1

cd "${output_dir}"

export APPTAINER_CACHEDIR=${project_root}/.apptainer/cache
export APPTAINER_TMPDIR=${project_root}/.apptainer/tmp
export APPTAINERENV_TMPDIR=${tmp_dir}

command=(
    "${apptainer}" exec --cleanenv --bind "${project_root}:${project_root}" "${image}"
    quartet CentroMiner
    -i "${genome}"
    --TE "${repeat_gff}"
    --gene "${gene_gff}"
    -t "${threads}"
    -p tianmu_tiemu
)

printf '%q ' "${command[@]}" > "${log_dir}/centrominer_command.txt"
printf '\n' >> "${log_dir}/centrominer_command.txt"

set +e
/usr/bin/time -v "${command[@]}" > >(tee "${log_dir}/centrominer.stdout.log") 2> >(tee "${log_dir}/centrominer.stderr.log" >&2)
status=$?
set -e

echo "${status}" > "${log_dir}/centrominer.exit_status"
date -Iseconds > "${log_dir}/run_finished.txt"

if [[ "${status}" -ne 0 ]]; then
    echo "CentroMiner failed with exit status ${status}." >&2
    exit "${status}"
fi

echo "CentroMiner completed successfully: ${output_dir}"

#
samtools=/opt/miniconda3/envs/bioenv/bin/samtools

$samtools view -F 2304 ont_vs_genome.sorted.bam | \
awk '
NR==FNR{
chr[++n]=$1
s[n]=$2
e[n]=$3
l[n]=$3-$2
next
}
{
for(i=1;i<=n;i++)
if($3==chr[i]){
rs=$4
re=$4+length($10)

ov=(re<e[i]?re:e[i])-(rs>s[i]?rs:s[i])

if(ov>0){
r=ov/l[i]
for(j=1;j<=10;j++)
if(r>=j/10)c[i,j]++
}
}
}
END{
print "Chr\t10%\t20%\t30%\t40%\t50%\t60%\t70%\t80%\t90%\t100%"
for(i=1;i<=n;i++){
printf "%s",chr[i]
for(j=1;j<=10;j++)
printf "\t%d",c[i,j]+0
print ""
}
}
' final_centromere_candidates.bed - \
> ONT_centromere_span_fraction.tsv
