#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_FILE="$SCRIPT_DIR/../../solution/tests/test_kpar_fast.scan_kpara_1.root_1"
OUTPUT_FILE="$SCRIPT_DIR/example_kpar.png"

echo "Running: python3 $SCRIPT_DIR/quick_dispersion_plot.py $INPUT_FILE --x-axis kpar --x-scale log --y-scale linear --output $OUTPUT_FILE"
python3 "$SCRIPT_DIR/quick_dispersion_plot.py" \
    "$INPUT_FILE" \
    --x-axis kpar \
    --x-scale log \
    --y-scale linear \
    --output "$OUTPUT_FILE"
