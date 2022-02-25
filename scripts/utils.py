from sklearn.preprocessing import normalize, scale
import numpy as np
from scipy.sparse import csr_matrix
import anndata
from pysctransform import SCTransform
import pandas as pd


def read_str_list(fname, lst=list):
    with open(fname, "r") as f:
        for idx, line in enumerate(f.readlines()):
            lst.append(line.strip())


def read_int_list(fname, lst=list):
    with open(fname, "r") as f:
        for idx, line in enumerate(f.readlines()):
            lst.append(int(line.strip()))


def sanitize_mtx(mtx: np.ndarray):
    cell_count_mask = mtx.sum(1) > 0  # count for each cell
    gene_count_mask = mtx.sum(0) > 0  # count for each gene

    genes_detected_mask = (mtx > 0).sum(1) > 0  # n genes per cell
    cells_detected_mask = (mtx > 0).sum(0) > 0  # n cells per gene
    row_mask = np.logical_and(cell_count_mask, genes_detected_mask)
    col_mask = np.logical_and(gene_count_mask, cells_detected_mask)

    return (row_mask, col_mask)


def do_pf(mtx, target_sum=None):
    pf = mtx.sum(1).flatten()
    mtx_pf = mtx / (pf / pf.mean())[:, None]
    if target_sum:
        mtx_pf = normalize(mtx, norm="l1") * target_sum
    return mtx_pf


def do_log_pf(mtx, pc=0.5, iter=1):
    print(f"iter: {iter}")
    log = np.log(mtx + pc)
    pf = do_pf(log)

    iter -= 1
    if iter == 0:
        return pf
    pf_up = do_pf(np.exp(pf) - 1)
    return do_log_pf(pf_up, iter)


def norm(mtx, pc=0.5):
    d = {}
    rm, cm = sanitize_mtx(mtx)
    sanmtx = mtx[rm][:, cm]

    print("sctransform")
    genes = np.arange(sanmtx.shape[1])

    var = pd.DataFrame(genes, columns=["gids"])
    adata = anndata.AnnData(X=csr_matrix(sanmtx), var=var)
    adata.var_names = var["gids"].astype(str)

    residuals = SCTransform(adata, var_features_n=3000, vst_flavor="v2")
    columns = residuals.columns.values.astype(int)

    remap_genes = np.array([list(genes).index(i) for i in columns])

    # when we do the remap genes we have to drop some rows since then become zero
    remapmtx = sanmtx[:, remap_genes]
    rm, cm = sanitize_mtx(remapmtx)
    mtx = remapmtx[rm]

    residuals = residuals[rm]

    d["sctransform"] = residuals.values

    print("raw")
    d["raw"] = mtx

    print("pf")
    d["pf"] = do_pf(mtx)

    print("log")
    d["log"] = np.log(pc + mtx)

    print("pf -> log")
    d["pf_log"] = np.log(pc + do_pf(mtx))

    print("pf -> log -> pf")
    d["pf_log_pf"] = do_log_pf(do_pf(mtx), pc=pc)

    print("cp10k -> log")
    d["cp10k_log"] = np.log(pc + do_pf(mtx, target_sum=10_000))

    print("cp10k -> log -> scale")
    d["cp10k_log_scale"] = pd.DataFrame(
        scale(np.log(pc + do_pf(mtx, target_sum=10_000)))
    ).values

    print("cpm -> log")
    d["cpm_log"] = np.log(pc + do_pf(mtx, target_sum=1_000_000))

    return d
