#!/usr/bin/env bash
# Regenerate samtools_reference.csv (manuscript context row: raw decode+count
# throughput ceiling) for a given single BAM + multi BAM dir.
#   Usage: samtools_reference.sh SINGLE_BAM MULTI_DIR OUT_CSV [REPS]
# Rows: single @ threads 1/48/72 (samtools view -c -@ N-1), median of REPS;
#       multi  @ 72 cores as 24 concurrent `samtools view -c -@ 2`.
# Schema (matches manuscript loader): scenario,threads,median_s,mb_per_s
set -uo pipefail
SINGLE=$1; MULTI_DIR=$2; OUT=$3; REPS=${4:-5}

command -v samtools >/dev/null || { echo "samtools not found" >&2; exit 1; }
mb() { echo "scale=3; $(stat -c%s "$1") / 1048576" | bc; }

median() {  # newline-separated numbers on stdin
    sort -n | awk '{a[NR]=$1} END {print (NR%2 ? a[(NR+1)/2] : (a[NR/2]+a[NR/2+1])/2)}'
}

echo "scenario,threads,median_s,mb_per_s" > "$OUT"

SMB=$(mb "$SINGLE")
for T in 1 48 72; do
    AT=$((T - 1))
    times=""
    for r in $(seq "$REPS"); do
        s=$(date +%s.%N)
        samtools view -c -@ "$AT" "$SINGLE" > /dev/null
        e=$(date +%s.%N)
        times+="$(echo "$e - $s" | bc)"$'\n'
    done
    med=$(printf '%s' "$times" | median)
    mbs=$(echo "scale=2; $SMB / $med" | bc)
    echo "single,$T,$med,$mbs" >> "$OUT"
    echo "  single @$T: ${med}s (${mbs} MB/s)" >&2
done

# multi: 24 concurrent view -c -@2 (= 72 cores), median of REPS rounds
MMB=$(du -sb "$MULTI_DIR"/*.bam | awk '{s+=$1} END {print s/1048576}')
times=""
for r in $(seq "$REPS"); do
    s=$(date +%s.%N)
    for f in "$MULTI_DIR"/*.bam; do
        samtools view -c -@ 2 "$f" > /dev/null &
    done
    wait
    e=$(date +%s.%N)
    times+="$(echo "$e - $s" | bc)"$'\n'
done
med=$(printf '%s' "$times" | median)
mbs=$(echo "scale=2; $MMB / $med" | bc)
echo "multi,72,$med,$mbs" >> "$OUT"
echo "  multi @72: ${med}s (${mbs} MB/s)" >&2
echo "wrote $OUT" >&2
