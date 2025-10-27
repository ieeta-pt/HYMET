#!/bin/bash

# Accept arguments
INPUT_DIR="$1"
MASH_SCREEN="$2"
SCREEN_TAB="$3"
FILTERED_SCREEN="$4"
SORTED_SCREEN="$5"
TOP_HITS="$6"
SELECTED_GENOMES="$7"
INITIAL_THRESHOLD="$8"

# Step 1: Run mash screen and generate sorted_screen
THREAD_COUNT="${THREADS:-8}"

shopt -s nullglob
INPUT_FILES=()
for pattern in "$INPUT_DIR"/*.fna "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.fasta "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fq; do
  for file in $pattern; do
    [ -e "$file" ] || continue
    INPUT_FILES+=("$file")
  done
done
shopt -u nullglob

if [ "${#INPUT_FILES[@]}" -eq 0 ]; then
  echo "ERROR: no input sequences found under ${INPUT_DIR}" >&2
  exit 1
fi

mash screen -p "${THREAD_COUNT}" -v 0.9 "$MASH_SCREEN" "${INPUT_FILES[@]}" > "$SCREEN_TAB"
sort -u -k5,5 "$SCREEN_TAB" > "$FILTERED_SCREEN"
sort -gr "$FILTERED_SCREEN" > "$SORTED_SCREEN"

# Step 2: Adjust the threshold and select genomes
num_sequences=$(
python3 - "${INPUT_FILES[@]}" <<'PY'
import pathlib
import sys

files = sys.argv[1:]
total = 0
for file in files:
    path = pathlib.Path(file)
    suffix = path.suffix.lower()
    if suffix in {".fna", ".fa", ".fasta"}:
        with path.open() as fh:
            total += sum(1 for line in fh if line.startswith(">"))
    elif suffix in {".fastq", ".fq"}:
        with path.open() as fh:
            total += sum(1 for idx, line in enumerate(fh, 1) if idx % 4 == 1)
    else:
        with path.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    total += 1
print(total)
PY
)

if ! [[ "$num_sequences" =~ ^[0-9]+$ ]]; then
  echo "WARNING: falling back to sequence count of 1 (got '$num_sequences')" >&2
  num_sequences=1
fi

min_candidates=$(echo "$num_sequences * 3.25" | bc | awk '{printf("%d\n",$1 + 0.5)}')
min_candidates=$(( min_candidates < 5 ? 5 : min_candidates ))

best_threshold=$INITIAL_THRESHOLD
current_threshold=$INITIAL_THRESHOLD
threshold_found=0

echo "===================================="
echo "Number of input sequences: $num_sequences"
echo "Minimum expected candidates: $min_candidates"
echo "===================================="

while (( $(echo "$current_threshold >= 0.70" | bc -l) )); do
    count=$(awk -v t="$current_threshold" '$1 > t' "$SORTED_SCREEN" | wc -l)
    
    echo "Testing threshold: $current_threshold"
    echo "Candidates found: $count"

    if [ "$count" -ge "$min_candidates" ]; then
        best_threshold=$current_threshold
        threshold_found=1
        break
    fi
    
    current_threshold=$(echo "$current_threshold - 0.02" | bc -l)
done

if [ "$threshold_found" -eq 0 ]; then
    best_threshold=0.70
    count=$(awk -v t="$best_threshold" '$1 > t' "$SORTED_SCREEN" | wc -l)
    echo "No suitable threshold found. Using 0.70."
fi

# Filter with the best threshold found
awk -v threshold="$best_threshold" '$1 > threshold' "$SORTED_SCREEN" > "$TOP_HITS"
cut -f5 "$TOP_HITS" > "$SELECTED_GENOMES"

echo "===================================="
echo "Final threshold used: $best_threshold"
echo "Candidates found: $count"
echo "===================================="
