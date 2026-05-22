#!/usr/bin/env bash
# bench_subset_genes.sh — compute metrics for all methods on every
# data/<DS>/subset_genes/ tree found under <data_root>, in parallel.
#
# Usage:
#   bash bench_subset_genes.sh [data_root] [parallel_sparse] [parallel_dense]
#
# Defaults: data_root=./data parallel_sparse=32 parallel_dense=8
#
# Each <DS>/subset_genes/ directory is expected to contain:
#   raw.mtx.gz                             (required)
#   {pf,log,sqrt,pf_log,pf_log_pf,cpm_log,cp10k_log}.mtx.gz   (sparse methods)
#   {sctransform,cp10k_log_scale}.csv.gz   (dense methods, optional)
#   <DS>_subset_genes_metrics.json          (required; will be updated in-place)

set -u
data_root="${1:-data}"
PSPARSE="${2:-32}"
PDENSE="${3:-8}"

# Resolve script directory + .venv python so we don't pay uv startup per job.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd -P)
PY="$REPO_ROOT/.venv/bin/python3"
[ -x "$PY" ] || { echo "ERROR: .venv/bin/python3 missing; run 'uv sync' first" >&2; exit 1; }

# Functions invoked per-dataset by xargs.
run_metric() {
    local script="$1" raw="$2" mtx="$3" json="$4" ds_dir="$5" method="$6"
    [ -f "$mtx" ] || return 0
    [ -f "$json" ] || return 0
    "$PY" "$script" "$raw" "$mtx" "$json" 2>&1 \
        | sed "s|^|[$method:$(basename "$(dirname "$ds_dir")")] |"
}
export PY
export -f run_metric

dispatch() {
    local script="$1" method="$2" ext="$3" parallel="$4"
    local glob_dir
    find "$data_root" -mindepth 2 -maxdepth 2 -type d -name subset_genes -print0 \
      | xargs -0 -n1 -P "$parallel" -I {} bash -c '
          ds_dir="$1"; script="$2"; method="$3"; ext="$4"
          ds=$(basename "$(dirname "$ds_dir")")
          run_metric "$script" "$ds_dir/raw.mtx.gz" "$ds_dir/${method}.${ext}" \
                     "$ds_dir/${ds}_subset_genes_metrics.json" "$ds_dir" "$method"
        ' _ {} "$script" "$method" "$ext"
}

# Make run_metric available to subshells.
export PY
export -f run_metric

echo "=== bench_subset_genes ==="
echo "data_root  = $data_root"
echo "sparse jobs= $PSPARSE   dense jobs= $PDENSE"
N_DS=$(find "$data_root" -mindepth 2 -maxdepth 2 -type d -name subset_genes | wc -l)
echo "datasets   = $N_DS"

for method in raw pf log sqrt pf_log pf_log_pf cpm_log cp10k_log; do
    t0=$(date +%s)
    echo "[sparse $method] start"
    dispatch "$SCRIPT_DIR/metrics_methods_sparse.py" "$method" "mtx.gz" "$PSPARSE"
    echo "[sparse $method] done in $(( $(date +%s) - t0 ))s"
done

for method in sctransform cp10k_log_scale; do
    t0=$(date +%s)
    echo "[dense  $method] start"
    dispatch "$SCRIPT_DIR/metrics_methods_dense.py" "$method" "csv.gz" "$PDENSE"
    echo "[dense  $method] done in $(( $(date +%s) - t0 ))s"
done

echo "=== bench_subset_genes done ==="
