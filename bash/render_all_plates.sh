#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Resolve repo root (parent of bash/)
# -------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RAW_DIR="$REPO_ROOT/data/raw"
RENDER_SCRIPT="$REPO_ROOT/run_render_one_plate.R"
HIST_DIR="$REPO_ROOT/outputs/history"

RESET_HISTORY=false

# -------------------------------
# Parse args
# -------------------------------
for arg in "$@"; do
  case $arg in
    --reset-history)
      RESET_HISTORY=true
      shift
      ;;
    *)
      ;;
  esac
done

# -------------------------------
# Sanity checks
# -------------------------------
if [[ ! -f "$RENDER_SCRIPT" ]]; then
  echo "ERROR: run_render_one_plate.R not found at $RENDER_SCRIPT"
  exit 1
fi

if [[ ! -d "$RAW_DIR" ]]; then
  echo "ERROR: raw data directory not found at $RAW_DIR"
  exit 1
fi

# -------------------------------
# Optional history reset
# -------------------------------
if [[ "$RESET_HISTORY" == true ]]; then
  echo "⚠️  Resetting QC history in $HIST_DIR"
  rm -f "$HIST_DIR"/*.rds
fi

FILES=$(ls "$RAW_DIR"/PlateRunResults_*.csv 2>/dev/null | sort)

if [[ -z "$FILES" ]]; then
  echo "No PlateRunResults CSV files found in $RAW_DIR"
  exit 0
fi

echo "Repo root: $REPO_ROOT"
echo "Rendering QC reports for plates:"
echo "$FILES"
echo "--------------------------------"

for f in $FILES; do
  echo "➡️  Rendering $(basename "$f")"
  Rscript "$RENDER_SCRIPT" "$f"
done

echo "✅ All QC reports rendered successfully."