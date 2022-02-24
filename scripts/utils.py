from sklearn.preprocessing import normalize, scale
import numpy as np


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


# def norm(mtx):
#     d = {}
#     print("raw")
#     d["raw"] = mtx

#     print("pf")
#     d["pf"] = do_pf(mtx)

#     print("log1p")
#     d["log1p"] = np.log1p(mtx)

#     print("pf -> log1p")
#     d["pf -> log1p"] = np.log1p(do_pf(mtx))

#     print("pf -> log1p -> pf")
#     d["pf -> log1p -> pf"] = do_log1p_pf(do_pf(mtx))

#     print("cpm -> log1p -> scale")
#     d["cpm -> log1p -> scale"] = scale(np.log1p(do_pf(mtx, target_sum=1_000_000)))

#     return d
