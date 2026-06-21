#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_DATA_DIR="${SCRIPT_DIR}/../../benchmark/output/clr_local/data/mouse_brain_10x"
DATA_DIR="${1:-${DEFAULT_DATA_DIR}}"

URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_filtered_feature_bc_matrix.h5"
FILENAME="neuron_1k_v3_filtered_feature_bc_matrix.h5"
SHA256="78166dda24a6103f6b690e0e5a392d89f9da429b84ee265ee4a0f4621f607695"

mkdir -p "${DATA_DIR}"
OUT="${DATA_DIR}/${FILENAME}"

if [[ ! -s "${OUT}" ]]; then
  curl -L --fail "${URL}" -o "${OUT}"
else
  echo "Using existing file: ${OUT}"
fi

if command -v shasum >/dev/null 2>&1; then
  printf "%s  %s\n" "${SHA256}" "${OUT}" | shasum -a 256 -c -
elif command -v sha256sum >/dev/null 2>&1; then
  printf "%s  %s\n" "${SHA256}" "${OUT}" | sha256sum -c -
else
  echo "Downloaded ${OUT}"
  echo "No SHA256 checker found; expected SHA256 is ${SHA256}"
fi
