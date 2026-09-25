#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status, 
# treat unset variables as an error, and fail pipelines on first error.
set -euo pipefail

# --- Directory Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
RAW_DIR="${DATA_DIR}/raw"
PROCESSED_DIR="${DATA_DIR}/processed"

mkdir -p "${RAW_DIR}" "${PROCESSED_DIR}"

echo "=================================================="
echo " Starting Bacterial Protein Preprocessing Pipeline "
echo "=================================================="

# -------------------------------------------------------------------
# STEP 1: Run Quality Filtering and Metadata Extraction Script First
# -------------------------------------------------------------------
echo ""
echo "[Step 1/4] Running getQualityBacteria.py to fetch PDB IDs & quality metadata..."

if [[ -f "getQualityBacteria.py" ]]; then
    python3 getQualityBacteria.py
else
    echo "Error: getQualityBacteria.py not found in current directory!" >&2
    exit 1
fi

echo "Step 1 complete! Bacterial metadata generated."



