#!/usr/bin/env python3

import sys
import os
from utils import do_pf, do_log_pf, norm, sanitize_mtx, notsct, sclcp
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
from sklearn.preprocessing import scale
import pandas as pd

import numpy as np


def main(in_matrix_fn, out_prefix):
    mtx = mmread(in_matrix_fn).toarray()
#    rm, cm = sanitize_mtx(mtx)

    # sanitize
#    sanmtx = mtx[rm][:, cm]
 
#    (data, _, _, _) = norm(sanmtx)
    # takes in raw.mtx.gz
    (data, _, _, _) = notsct(mtx)
    # takes in cp10_log.mtx
#    (data, _, _, _) = sclcp(mtx)

    titles = ["pf", "log", "pf_log", "pf_log_pf", "cpm_log", "cp10k_log", "sqrt"]
    for title in titles:
        print(f"saving {title}")
        out_fn = os.path.join(out_prefix, f"{title}.mtx")
        mmwrite(out_fn, csr_matrix(data[f"{title}"]))

    title = "cp10k_log_scale"
    print(f"saving {title}")
    out_fn = os.path.join(out_prefix, f"{title}.csv")
    pd.DataFrame(data[title]).to_csv(out_fn, index=False, header=False)

#    title = "sctransform"
#    print(f"saving {title}")
#    out_fn = os.path.join(out_prefix, f"{title}.csv")
#    pd.DataFrame(data[title]).to_csv(out_fn, index=False, header=False)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
