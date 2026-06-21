#!/usr/bin/env python3
"""Sync current PFlog metrics into committed supplement flat JSON files.

The flat supplement files historically use the prefix ``PFlog1pPF`` for the
PFlog row.  For compatibility we keep that prefix, but the values are sourced
from the current ``PFlog (shift. CLR)`` entry produced by scclr:

    center(log1p(4 * alpha * x))

Only PFlog fields are updated.  Other normalization metrics are left untouched.
"""

from __future__ import annotations

import glob
import json
import os
import sys
from pathlib import Path


SOURCE_LABEL = "PFlog (shift. CLR)"
FLAT_PREFIX = "PFlog1pPF"


def sync_one(source_json: Path, flat_json: Path) -> bool:
    source = json.loads(source_json.read_text())
    entry = source.get(SOURCE_LABEL)
    if not isinstance(entry, dict):
        raise KeyError(f"{source_json} has no {SOURCE_LABEL!r} entry")

    flat = json.loads(flat_json.read_text())
    changed = False
    mapping = {
        "cov_gene": f"{FLAT_PREFIX}_cov_gene",
        "cov_cell": f"{FLAT_PREFIX}_cov_cell",
        "r2_depth": f"{FLAT_PREFIX}_r2_depth",
        "r_mono": f"{FLAT_PREFIX}_r_mono",
        "scclr_alpha": f"{FLAT_PREFIX}_scclr_alpha",
        "scclr_k": f"{FLAT_PREFIX}_scclr_k",
        "pflog_formula": f"{FLAT_PREFIX}_formula",
    }
    for source_key, flat_key in mapping.items():
        value = entry.get(source_key)
        if flat.get(flat_key) != value:
            flat[flat_key] = value
            changed = True

    impl = source.get("pflog_impl", "scclr PFlog: center(log1p(4*alpha*x))")
    if flat.get(f"{FLAT_PREFIX}_impl") != impl:
        flat[f"{FLAT_PREFIX}_impl"] = impl
        changed = True

    if changed:
        flat_json.write_text(json.dumps(flat, indent=2) + "\n")
    return changed


def main(data_root: str, supplement_metrics_dir: str) -> None:
    data_root_path = Path(data_root)
    metrics_dir = Path(supplement_metrics_dir)
    changed = missing = 0
    source_paths = sorted(
        Path(p)
        for p in glob.glob(str(data_root_path / "*" / "subset_genes" / "*_subset_genes_metrics.json"))
    )
    for source_json in source_paths:
        ds = source_json.name.removesuffix("_subset_genes_metrics.json")
        flat_json = metrics_dir / f"{ds}_all_genes_metrics_flat.json"
        if not flat_json.exists():
            missing += 1
            continue
        changed += int(sync_one(source_json, flat_json))
    print(f"synced {len(source_paths)} source metrics; changed {changed}; missing flat {missing}")


if __name__ == "__main__":
    data_root = sys.argv[1] if len(sys.argv) > 1 else "/home/sina/projects/synchromesh/data"
    metrics_dir = sys.argv[2] if len(sys.argv) > 2 else "analysis/supplement/metrics"
    main(data_root, metrics_dir)
