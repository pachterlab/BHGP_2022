#!/usr/bin/env python3
"""Remake of Figure 1b panel 1 (variance-stabilization panel): per-method boxplot of
CV(q) across datasets, replacing the old 'CV of gene variances' (cov_gene) metric.

Functional transforms use the model-based delta-method q (scripts/compute_cvq_all.py).
sctransform has no functional q, so it is shown as its EMPIRICAL cv_q_emp = CV of the
full-vst Pearson-residual gene variances on the same filtered genes
(scripts/sct_cvq_one.py -> sct_cvq_filtered.csv). It is colored distinctly to flag the
empirical-vs-model distinction. Dataset filter matches analysis/summary.ipynb:
avg_per_cell >= angelidis_2019's, minus the NaN-metric datasets."""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG  = os.path.join(HERE, "..", "analysis", "figures")
CSV  = os.path.join(FIG, "cvq_per_method_subset_filtered.csv")
SCT  = os.path.join(FIG, "sct_cvq_filtered.csv")
OUT  = os.path.join(FIG, "fig1b_panel1_cvq")

matplotlib.rcParams["mathtext.fontset"] = "cm"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 15
GENE = matplotlib.colormaps["cividis"](0.5)
SCTC = "#C9527E"   # distinct color: sctransform is empirical, not model-based q

NAN = ['ERS4228665','GSM3576415','GSM3576416','GSM3576418','CRX102286','CRX102292',
       'GSM3576425','GSM3731527','GSM3396166','GSM3148578','GSM3576429','CRX118014',
       'ERX2756731','GSM3576398','GSM4037323','GSM3576410']

df = pd.read_csv(CSV)
OLD_PFLOG_NAME = "PFlog" + "PF"
if "PFlog" not in df.columns and OLD_PFLOG_NAME in df.columns:
    df = df.rename(columns={OLD_PFLOG_NAME: "PFlog"})
# Attach sctransform empirical cv_q_emp (filtered) if computed
have_sct = os.path.exists(SCT)
if have_sct:
    s = pd.read_csv(SCT)[["ds","cv_filt"]].rename(columns={"cv_filt":"sctransform"})
    df = df.merge(s, on="ds", how="left")

# order matches summary.ipynb `labels`; sctransform included only if available.
# sctransform gets an xlabel asterisk so the caption can flag it as an EMPIRICAL CoV
# (CV of full-vst residual gene variances) rather than the model-based delta-method q.
ORDER = ["raw","PF","sqrt","log1p","log1pCP10k","log1pCPM","scalelog1pCP10k"] \
        + (["sctransform"] if have_sct else []) + ["log1pPF","PFlog"]
DISP  = {"PFlog":"PFlog\n(shift. CLR)", "sctransform":"sctransform*"}

_amj = os.path.join(HERE, "..", "analysis", "supplement", "metrics", "angelidis_2019_all_genes_metrics_flat.json")
amin = float(json.load(open(_amj))["avg_per_cell"])
passing = df[(df.sbar >= amin) & (~df.ds.isin(NAN))].copy()
print(f"{len(passing)} of {len(df)} datasets pass (avg_per_cell >= {amin:.0f}, minus NaN list)")

data = [passing[m].replace([np.inf,-np.inf],np.nan).dropna().values for m in ORDER]
print("median:", {m: round(float(np.median(d)),3) for m,d in zip(ORDER,data) if len(d)})

fig, ax = plt.subplots(figsize=(11.5,5.5))
pos = np.arange(len(ORDER))
bp = ax.boxplot(data, positions=pos, widths=0.6, patch_artist=True, notch=True, showfliers=True,
                medianprops=dict(color="k", linewidth=2),
                whiskerprops=dict(color="k"), capprops=dict(color="k"))
# all boxes share the cividis "gene" color (as in the original panel); sctransform is
# distinguished only by the xlabel asterisk (+ caption), not by color.
for box, fl in zip(bp["boxes"], bp["fliers"]):
    box.set_facecolor(GENE); box.set_edgecolor("k")
    fl.set(marker="o", markersize=3, markerfacecolor=GENE, markeredgecolor="none", alpha=0.4)
ax.scatter(pos, [np.mean(d) for d in data], marker="D", s=45, facecolor="white",
           edgecolor="k", zorder=5, label="mean")

ax.set_yscale("log")
ax.set_xticks(pos)
ax.set_xticklabels([DISP.get(m,m) for m in ORDER], rotation=40, ha="right", fontsize=11)
ax.set_ylabel("Coefficient of variation gene variance")
ax.set_title(f"Variance stabilization — CV($q$) across {len(passing)} datasets", fontsize=14)
ax.set_ylim(0.03, None)
handles = [plt.Line2D([],[],marker="D",color="w",markerfacecolor="w",markeredgecolor="k",label="mean")]
if not have_sct:
    ax.axhline(0.04, color=SCTC, linestyle="--", linewidth=1.8, zorder=1)
    ax.text(len(ORDER)-0.5, 0.043, "sctransform: $q=0$ by construction",
            color=SCTC, ha="right", va="bottom", fontsize=11)
ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=11)
fig.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(f"{OUT}.{ext}", dpi=200, facecolor="white", bbox_inches="tight")
print(f"saved -> {OUT}.png / .pdf")
