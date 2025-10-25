#!/bin/bash

# Accept arguments
INPUT_DIR="$1"          # This argument will no longer be used for .fna files directly.
REFERENCE_SET="$2"      # Path to combined_genomes.fasta.
NT_MMI="$3"             # Path to the minimap2 index (reference.mmi).
RESULTADOS_PAF="$4"     # Path to the alignment results (resultados.paf).

# Create the reference set index
if [ ! -s "$NT_MMI" ]; then
    echo "Creating index with minimap2..."
    minimap2 -I2g -d "$NT_MMI" "$REFERENCE_SET" # General index
    if [ $? -ne 0 ]; then
        echo "Error creating index with minimap2."
        exit 1
    fi
else
    echo "Using cached minimap2 index: $NT_MMI"
fi

INPUT_MODE="${INPUT_MODE:-contigs}"
THREAD_COUNT="${THREADS:-8}"

shopt -s nullglob
INPUT_FILES=()
case "$INPUT_MODE" in
  contigs)
    for pattern in "$INPUT_DIR"/*.fna "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.fasta; do
      for file in $pattern; do
        [ -e "$file" ] || continue
        INPUT_FILES+=("$file")
      done
    done
    ;;
  reads)
    for pattern in "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fq; do
      for file in $pattern; do
        [ -e "$file" ] || continue
        INPUT_FILES+=("$file")
      done
    done
    ;;
  *)
    echo "Unsupported INPUT_MODE '${INPUT_MODE}'" >&2
    exit 1
    ;;
esac
shopt -u nullglob

if [ "${#INPUT_FILES[@]}" -eq 0 ]; then
    echo "No input files found for mode ${INPUT_MODE} under ${INPUT_DIR}" >&2
    exit 1
fi

# Run alignment using minimap2 (for long reads)
echo "Running alignment with minimap2..."
if [ "$INPUT_MODE" = "reads" ]; then
    PRESET="${MINIMAP2_READS_PRESET:-sr}"
    minimap2 -I2g -t "${THREAD_COUNT}" -x "${PRESET}" "$NT_MMI" "${INPUT_FILES[@]}" >"$RESULTADOS_PAF"
else
    minimap2 -I2g -t "${THREAD_COUNT}" -x asm10 "$NT_MMI" "${INPUT_FILES[@]}" >"$RESULTADOS_PAF"
fi

# Check if the alignment was successful
if [ $? -ne 0 ]; then
    echo "Error running alignment with minimap2."
    exit 1
fi

echo "Alignment completed successfully! Results saved to $RESULTADOS_PAF."
