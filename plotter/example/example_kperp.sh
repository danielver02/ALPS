#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_FILE="$SCRIPT_DIR/../../solution/tests/test_kperp.scan_kperp_1.root_1"
OUTPUT_FILE="$SCRIPT_DIR/example_kperp.png"

echo "Running: python3 $SCRIPT_DIR/quick_dispersion_plot.py $INPUT_FILE --x-axis kperp --x-scale log --y-scale log --output $OUTPUT_FILE"
python3 "$SCRIPT_DIR/quick_dispersion_plot.py" \
    "$INPUT_FILE" \
    --x-axis kperp \
    --x-scale log \
    --y-scale log \
    --output "$OUTPUT_FILE"
