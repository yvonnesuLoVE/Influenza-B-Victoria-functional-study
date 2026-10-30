#!/usr/bin/env bash
# Generic staged workflow distilled from the original sequencing notebook.
# Usage: bash scripts/sequencing.sh config/sequencing_config.sh STAGE
set -euo pipefail
usage() {
    echo "Usage: bash scripts/sequencing.sh CONFIG trim|index|map|consensus|qc|classify|headers|all"
}
[[ $# -eq 2 ]] || { usage >&2; exit 2; }
CONFIG=$1
STAGE=$2
case "$STAGE" in trim|index|map|consensus|qc|classify|headers|all) ;; *) usage >&2; exit 2 ;; esac
[[ -r "$CONFIG" ]] || { echo "Cannot read config: $CONFIG" >&2; exit 1; }
# A configuration is trusted shell code, not an input data file.
source "$CONFIG"
: "${SAMPLE_LIST:?}" "${REFERENCE_LIST:?}" "${RAW_DIR:?}" "${REF_DIR:?}" "${OUT_DIR:?}"
: "${THREADS:?}" "${TRIMMOMATIC:?}" "${BWA:?}" "${SAMTOOLS:?}" "${KRAKEN2:?}"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 1; }
require_tool() { command -v "$1" >/dev/null 2>&1 || { echo "Executable not found: $1" >&2; exit 1; }; }
require_file() { [[ -s "$1" ]] || { echo "Missing or empty input: $1" >&2; exit 1; }; }
# IDs are filenames, not paths. Accept CRLF lists; ignore blank/comment lines.
read_ids() {
    local line
    while IFS= read -r line || [[ -n "$line" ]]; do
        line=${line%$'\r'}
        [[ -z "$line" || "$line" == \#* ]] && continue
        [[ "$line" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo "Invalid ID: $line" >&2; return 1; }
        printf '%s\n' "$line"
    done < "$1"
}
require_file "$SAMPLE_LIST"; require_file "$REFERENCE_LIST"
mkdir -p "$OUT_DIR/logs" "$OUT_DIR/paired" "$OUT_DIR/unpaired" "$OUT_DIR/mapping" \
    "$OUT_DIR/consensus" "$OUT_DIR/consensus_named" "$OUT_DIR/qc" "$OUT_DIR/kraken2" "$OUT_DIR/references"
# Store normalized lists so failures in validation cannot disappear in a subshell.
read_ids "$SAMPLE_LIST" > "$OUT_DIR/logs/samples.txt"
read_ids "$REFERENCE_LIST" > "$OUT_DIR/logs/references.txt"
require_file "$OUT_DIR/logs/samples.txt"; require_file "$OUT_DIR/logs/references.txt"
for list in "$OUT_DIR/logs/samples.txt" "$OUT_DIR/logs/references.txt"; do
    [[ -z "$(sort "$list" | uniq -d)" ]] || { echo "Duplicate IDs in $list" >&2; exit 1; }
done
cp "$CONFIG" "$OUT_DIR/logs/config_used.sh"
LOG="$OUT_DIR/logs/run.log"
log() { printf '[%s] %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$*" | tee -a "$LOG"; }
trap 'echo "Failed at line $LINENO during stage $STAGE" >&2' ERR
trim_reads() {
    require_tool "$TRIMMOMATIC"
    "$TRIMMOMATIC" -version > "$OUT_DIR/logs/trimmomatic_version.txt" 2>&1
    while IFS= read -r sample; do
        r1="$RAW_DIR/$sample.R1.fastq.gz"; r2="$RAW_DIR/$sample.R2.fastq.gz"
        require_file "$r1"; require_file "$r2"
        log "Trim $sample"
        "$TRIMMOMATIC" PE -threads "$THREADS" "$r1" "$r2" \
            "$OUT_DIR/paired/$sample.trimmed.paired.R1.fastq.gz" "$OUT_DIR/unpaired/$sample.trimmed.unpaired.R1.fastq.gz" \
            "$OUT_DIR/paired/$sample.trimmed.paired.R2.fastq.gz" "$OUT_DIR/unpaired/$sample.trimmed.unpaired.R2.fastq.gz" \
            SLIDINGWINDOW:4:20 LEADING:30 TRAILING:30 MINLEN:70 \
            > "$OUT_DIR/logs/$sample.trim.log" 2>&1
    done < "$OUT_DIR/logs/samples.txt"
}
index_references() {
    require_tool "$BWA"
    # bwa-mem2's no-argument version/help command normally returns nonzero.
    "$BWA" > "$OUT_DIR/logs/bwa_version.txt" 2>&1 || true
    while IFS= read -r reference; do
        require_file "$REF_DIR/$reference.fasta"
        # Index a copy so input data are never modified.
        cp "$REF_DIR/$reference.fasta" "$OUT_DIR/references/$reference.fasta"
        "$BWA" index "$OUT_DIR/references/$reference.fasta" > "$OUT_DIR/logs/$reference.index.log" 2>&1
    done < "$OUT_DIR/logs/references.txt"
}
map_reads() {
    require_tool "$BWA"; require_tool "$SAMTOOLS"
    "$SAMTOOLS" --version > "$OUT_DIR/logs/samtools_version.txt"
    while IFS= read -r sample; do
        r1="$OUT_DIR/paired/$sample.trimmed.paired.R1.fastq.gz"
        r2="$OUT_DIR/paired/$sample.trimmed.paired.R2.fastq.gz"
        require_file "$r1"; require_file "$r2"
        while IFS= read -r reference; do
            fasta="$OUT_DIR/references/$reference.fasta"
            require_file "$fasta"
            prefix="$sample.$reference"
            log "Map $prefix"
            "$BWA" mem -t "$THREADS" "$fasta" "$r1" "$r2" \
                > "$OUT_DIR/mapping/$prefix.sam" 2> "$OUT_DIR/logs/$prefix.map.log"
            "$SAMTOOLS" sort -@ "$THREADS" -O BAM -o "$OUT_DIR/mapping/$prefix.sorted.bam" \
                "$OUT_DIR/mapping/$prefix.sam"
            "$SAMTOOLS" index "$OUT_DIR/mapping/$prefix.sorted.bam"
        done < "$OUT_DIR/logs/references.txt"
    done < "$OUT_DIR/logs/samples.txt"
}
call_consensus() {
    require_tool "$SAMTOOLS"
    while IFS= read -r sample; do
        while IFS= read -r reference; do
            prefix="$sample.$reference"; bam="$OUT_DIR/mapping/$prefix.sorted.bam"
            require_file "$bam"
            log "Consensus $prefix"
            # Historical simple-mode settings; do not replace with default Bayesian mode.
            "$SAMTOOLS" consensus -@ "$THREADS" --mode simple -f fasta -a -A -c 0.5 -H 0.5 \
                -o "$OUT_DIR/consensus/$prefix.consensus.fasta" "$bam"
        done < "$OUT_DIR/logs/references.txt"
    done < "$OUT_DIR/logs/samples.txt"
}
quality_control() {
    require_tool "$SAMTOOLS"
    while IFS= read -r sample; do
        while IFS= read -r reference; do
            prefix="$sample.$reference"; bam="$OUT_DIR/mapping/$prefix.sorted.bam"
            require_file "$bam"
            qc="$OUT_DIR/qc/$prefix"; mkdir -p "$qc"
            log "QC $prefix"
            "$SAMTOOLS" quickcheck -v "$bam" > "$qc/quickcheck.txt" 2>&1
            "$SAMTOOLS" index "$bam"
            "$SAMTOOLS" flagstat "$bam" > "$qc/flagstat.txt"
            "$SAMTOOLS" idxstats "$bam" > "$qc/idxstats.txt"
            "$SAMTOOLS" coverage "$bam" > "$qc/coverage.tsv"
            "$SAMTOOLS" depth -aa "$bam" > "$qc/depth.tsv"
            awk '{d=$3; n++; sum+=d; if(n==1||d<min)min=d; if(d>max)max=d;
                  if(d>=1)c1++; if(d>=10)c10++; if(d>=20)c20++; if(d>=50)c50++; if(d>=100)c100++}
                 END {if(n){printf "positions\t%d\nmean_depth\t%.2f\nmin_depth\t%d\nmax_depth\t%d\n",n,sum/n,min,max;
                      printf "breadth_1x\t%.2f%%\nbreadth_10x\t%.2f%%\nbreadth_20x\t%.2f%%\nbreadth_50x\t%.2f%%\nbreadth_100x\t%.2f%%\n",100*c1/n,100*c10/n,100*c20/n,100*c50/n,100*c100/n}}' \
                "$qc/depth.tsv" > "$qc/depth_summary.txt"
            "$SAMTOOLS" stats "$bam" > "$qc/stats.txt"
            awk '/^SN/' "$qc/stats.txt" > "$qc/stats_summary.txt"
            "$SAMTOOLS" view "$bam" | awk '{mq[$5]++} END{for(m in mq)print m "\t" mq[m]}' | sort -n > "$qc/mapq_distribution.tsv"
            awk '$3<10' "$qc/depth.tsv" > "$qc/depth_lt10.tsv"
        done < "$OUT_DIR/logs/references.txt"
    done < "$OUT_DIR/logs/samples.txt"
}
classify_reads() {
    require_tool "$KRAKEN2"
    : "${KRAKEN_DB:?Set KRAKEN_DB for classification}"
    [[ -d "$KRAKEN_DB" ]] || { echo "Database directory not found: $KRAKEN_DB" >&2; exit 1; }
    "$KRAKEN2" --version > "$OUT_DIR/logs/kraken2_version.txt" 2>&1
    while IFS= read -r sample; do
        r1="$OUT_DIR/paired/$sample.trimmed.paired.R1.fastq.gz"
        r2="$OUT_DIR/paired/$sample.trimmed.paired.R2.fastq.gz"
        require_file "$r1"; require_file "$r2"
        "$KRAKEN2" --threads "$THREADS" --db "$KRAKEN_DB" --paired --gzip-compressed \
            --report "$OUT_DIR/kraken2/$sample.report.txt" --output "$OUT_DIR/kraken2/$sample.output.txt" "$r1" "$r2"
    done < "$OUT_DIR/logs/samples.txt"
}
name_consensus() {
    while IFS= read -r sample; do
        while IFS= read -r reference; do
            prefix="$sample.$reference"
            cons="$OUT_DIR/consensus/$prefix.consensus.fasta"
            fasta="$REF_DIR/$reference.fasta"
            require_file "$cons"; require_file "$fasta"
            # Map every consensus record to its own reference record. Keep all
            # segments; never change only the first FASTA header.
            awk -v prefix="$prefix" '
                FNR==NR {if(/^>/){header=substr($0,2);split(header,a,/[ \t]/);id=a[1];gene="unknown";
                    if(index(header,"__")){sub(/^.*__/,"",header);split(header,b,/\|\|/);gene=b[1]}
                    genes[id]=gene}next}
                /^>/ {id=substr($0,2);sub(/[ \t].*$/,"",id);gene=(id in genes)?genes[id]:"unknown";
                      print ">" prefix "|" gene "|" id;next}
                {print}' "$fasta" "$cons" > "$OUT_DIR/consensus_named/$prefix.consensus.fasta"
        done < "$OUT_DIR/logs/references.txt"
    done < "$OUT_DIR/logs/samples.txt"
}
log "Start $STAGE"
case "$STAGE" in
    trim) trim_reads ;; index) index_references ;; map) map_reads ;;
    consensus) call_consensus ;; qc) quality_control ;; classify) classify_reads ;;
    headers) name_consensus ;;
    all) trim_reads; index_references; map_reads; call_consensus; quality_control; name_consensus ;;
esac
log "Finished $STAGE"
