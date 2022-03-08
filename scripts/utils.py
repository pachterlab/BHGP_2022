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

import anndata
from pysctransform import SCTransform

def norm(sanmtx, sanbcs = [], sangenes = []):
    nc, ng = sanmtx.shape
    pc = 1.0
    if len(sanbcs) == 0:
        sanbcs = np.arange(nc).astype(str)
    if len(sangenes) == 0:
        sangenes = np.arange(ng).astype(str)
    
    # assume you get a sanitized matrix
    data = {}

    # prepare anndata for sctransform
    var = pd.DataFrame(sangenes, columns=["gids"])
    obs = pd.DataFrame(sanbcs, columns=["bcs"])
    adata = anndata.AnnData(X=csr_matrix(sanmtx), var=var, obs=obs)
    adata.var_names = var["gids"].astype(str)
    adata.obs_names = obs["bcs"].astype(str)

    print("sctransform")    
    residuals = SCTransform(adata, var_features_n=ng, vst_flavor="v2")
    sctgenes = residuals.columns.values

    reorder_gidx = np.array([list(sangenes).index(i) for i in sctgenes])

    # when we do the remap genes we have to drop some rows since then become zero
    reorder_mtx = sanmtx[:, reorder_gidx]
    rm, cm = sanitize_mtx(reorder_mtx)

    # clean and ordered matrices (we dont drop genes)
    mtx = reorder_mtx[rm]
    bcs = sanbcs[rm]
    genes = sangenes[reorder_gidx]
    
    sct = residuals[rm].values
    
    # create data dict with all transformations
    data["sctransform"] = sct
    print("raw")
    data["raw"] = mtx
    print("pf")
    data["pf"] = do_pf(mtx)
    print("log")
    data["log"] = np.log(pc + mtx)
    print("pf_log")
    data["pf_log"] = np.log(pc + do_pf(mtx))
    print("pf_log_pf")
    data["pf_log_pf"] = do_log_pf(do_pf(mtx), pc=pc)
    print("cp10k_log")
    data["cp10k_log"] = np.log(pc + do_pf(mtx, target_sum=10_000))
    print("cp10k_log_scale")
    data["cp10k_log_scale"] = pd.DataFrame(
        scale(np.log(pc + do_pf(mtx, target_sum=10_000)))
    ).values
    print("cpm_log")
    data["cpm_log"] = np.log(pc + do_pf(mtx, target_sum=1_000_000))
    print("sqrt")
    data["sqrt"] = np.sqrt(mtx)
    return (data, mtx, bcs, genes)
