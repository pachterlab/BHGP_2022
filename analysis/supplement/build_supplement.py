#!/usr/bin/env python3
"""Generate the per-dataset supplement LaTeX (supplementary_per_dataset_methods.tex).

For each dataset (in data/datasets.txt order) it emits one numbered supplementary
figure (Supplementary Figure 1.X): the dataset's all-genes JSON metrics as a table,
the subset-genes method-comparison figure embedded via \\includegraphics, and a
caption carrying the accession.

Figures are embedded as 150-DPI PNGs (rasterized from the per-dataset PDFs) to
keep the compiled supplement small enough to upload to bioRxiv.

Inputs:
  - data/datasets.txt                              (dataset order)
  - analysis/supplement/metrics/<DS>_all_genes_metrics_flat.json   (dataset metrics)
  - analysis/supplement/figs_png/<DS>_subset_genes_method_comparison.png  (figure)

Usage: python build_supplement.py            (writes supplementary_per_dataset_methods.tex)
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = "/home/sina/projects/synchromesh/data"
DATASETS = os.path.join(REPO, "data", "datasets.txt")
OUT = os.path.join(HERE, "supplementary_per_dataset_methods.tex")

def esc(s):                                   # escape LaTeX specials in accessions
    return s.replace("\\", r"\textbackslash{}").replace("_", r"\_").replace("&", r"\&").replace("%", r"\%").replace("#", r"\#")

def comma(x):
    try: return f"{int(round(float(x))):,}"
    except Exception: return str(x)

def sig(x, n=3):
    try: return f"{float(x):.{n}g}"
    except Exception: return str(x)

# (json key, row label, formatter)
ROWS = [
    ("ncells", "cells", comma), ("ngenes", "genes", comma),
    ("nvals", "nonzero entries", comma), ("density", "density", sig),
    ("total_counts", "total counts", comma),
    ("avg_per_cell", "mean counts / cell", comma),
    ("avg_per_gene", "mean counts / gene", comma),
    ("min_cell", "min counts / cell", comma),
    ("max_cell", "max counts / cell", comma),
    ("overdispersion", "overdispersion", sig),
]

PREAMBLE = r"""\documentclass[10pt,a4paper]{article}
\usepackage[margin=1in]{geometry}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{booktabs}
\usepackage[hidelinks]{hyperref}
\usepackage{bookmark}
\hypersetup{
  bookmarksopen=true, bookmarksopenlevel=1,
  pdftitle={Supplementary Figures: per-dataset method comparisons},
  pdfauthor={Booeshaghi et al.},
}
% Number figures as 1.X (Supplementary Figure 1.1, 1.2, ...).
\renewcommand{\thefigure}{1.\arabic{figure}}
\captionsetup[figure]{name=Supplementary Figure, labelsep=period, font=small}
\setlength{\parindent}{0pt}
\begin{document}
\begin{center}
{\Large\bfseries Supplementary Figures: per-dataset method comparisons}\\[0.4em]
{\normalsize Booeshaghi et al.}\\[0.4em]
{\small 526 datasets, numbered Supplementary Figure 1.1--1.526. Use the PDF bookmarks (outline) to navigate by accession.}
\end{center}
\clearpage
"""

FIG = r"""\bookmark[level=1,page=\thepage]{{Fig 1.{idx}  {acc_bm}}}
\begin{{figure}}[p]
\centering
{{\small
\begin{{tabular}}{{@{{}}ll@{{}}}}
\toprule
\multicolumn{{2}}{{@{{}}l}}{{\textbf{{{acc}}}}}\\
\midrule
{rows}
\bottomrule
\end{{tabular}}}}\par\vspace{{0.7em}}
\includegraphics[width=\linewidth,height=0.62\textheight,keepaspectratio]{{{pdf}}}
\caption{{\textbf{{{acc}}}. Per-method comparison of variance stabilization (left), depth normalization (middle), and monotonicity (right) for each normalization on dataset {acc}.}}
\label{{suppfig:{idx}}}
\end{{figure}}
\clearpage
"""

def main(limit=None):
    datasets = [l.strip() for l in open(DATASETS) if l.strip()]
    if limit:
        datasets = datasets[:int(limit)]
    blocks, missing_fig, missing_metrics = [], [], []
    idx = 0
    for ds in datasets:
        pdf = os.path.join(HERE, "figs_png", f"{ds}_subset_genes_method_comparison.png")
        mfn = os.path.join(HERE, "metrics", f"{ds}_all_genes_metrics_flat.json")
        if not os.path.exists(pdf):
            missing_fig.append(ds); continue
        m = json.load(open(mfn)) if os.path.exists(mfn) else {}
        if not m: missing_metrics.append(ds)
        rows = "\n".join(
            rf"{lbl} & {fmt(m[k])} \\" for k, lbl, fmt in ROWS if k in m
        ) or r"\multicolumn{2}{@{}l}{(metrics unavailable)} \\"
        idx += 1
        blocks.append(FIG.format(acc=esc(ds), acc_bm=ds.replace("_", "-"), rows=rows, pdf=pdf, idx=idx))
    with open(OUT, "w") as f:
        f.write(PREAMBLE)
        f.write("\n".join(blocks))
        f.write("\n\\end{document}\n")
    print(f"wrote {OUT}: {idx} figures", file=sys.stderr)
    if missing_fig: print(f"  WARNING: {len(missing_fig)} datasets missing figure pdf: {missing_fig[:5]}...", file=sys.stderr)
    if missing_metrics: print(f"  WARNING: {len(missing_metrics)} datasets missing metrics json", file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
