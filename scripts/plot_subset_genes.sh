#!/usr/bin/env bash
# plot_subset_genes.sh — run plot_all.py on every data/<DS>/subset_genes/
# tree under <data_root>, writing the 3-panel comparison PDF into each
# subset_genes/ dir as <DS>_subset_genes_method_comparison.pdf.
#
# Usage:
#   bash plot_subset_genes.sh [data_root] [parallel]
#
# Defaults: data_root=./data parallel=12

set -u
data_root="${1:-data}"
PAR="${2:-12}"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd -P)
PY="$REPO_ROOT/.venv/bin/python3"
[ -x "$PY" ] || { echo "ERROR: .venv/bin/python3 missing; run 'uv sync' first" >&2; exit 1; }

N_DS=$(find "$data_root" -mindepth 2 -maxdepth 2 -type d -name subset_genes | wc -l)
echo "=== plot_subset_genes ==="
echo "data_root  = $data_root"
echo "parallel   = $PAR"
echo "datasets   = $N_DS"

t0=$(date +%s)
find "$data_root" -mindepth 2 -maxdepth 2 -type d -name subset_genes -print0 \
  | xargs -0 -n1 -P "$PAR" -I {} bash -c '
      ds_dir="$1"
      script="$2"
      ds=$(basename "$(dirname "$ds_dir")")
      "'"$PY"'" "$script" "${ds}_subset_genes" "$ds_dir" "$ds_dir" 2>&1 \
          | sed "s|^|[plot:$ds] |" | tail -3
    ' _ {} "$SCRIPT_DIR/plot_all.py"

echo "=== plot_subset_genes done in $(( $(date +%s) - t0 ))s ==="
