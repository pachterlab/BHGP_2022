#!/usr/bin/env python3

import sys
import os
from utils import sanitize_mtx, sct
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
from sklearn.preprocessing import scale
import pandas as pd

import numpy as np


def main(in_matrix_fn, out_prefix):
    mtx = mmread(in_matrix_fn).toarray()
    rm, cm = sanitize_mtx(mtx)

    # sanitize
    sanmtx = mtx[rm][:, cm]
 
    (data, _, _, _) = sct(sanmtx)

    titles = ["raw"] 
    for title in titles:
        print(f"saving {title}")
        out_fn = os.path.join(out_prefix, f"{title}.mtx")
        mmwrite(out_fn, csr_matrix(data[f"{title}"]))

    title = "sctransform"
    print(f"saving {title}")
    out_fn = os.path.join(out_prefix, f"{title}.csv.gz")
    pd.DataFrame(data[title]).to_csv(out_fn, index=False, header=False, compression="gzip")
    

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
