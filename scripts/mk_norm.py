#!/usr/bin/env python3

import sys
import os
from utils import do_pf, do_log1p_pf
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
from sklearn.preprocessing import scale
import pandas as pd

import numpy as np


def main(in_matrix_fn, out_prefix):

    mtx = mmread(in_matrix_fn).toarray()

    print("pf")
    out_fn = os.path.join(out_prefix, "pf.mtx")
    mmwrite(out_fn, csr_matrix(do_pf(mtx)))

    print("log1p")
    out_fn = os.path.join(out_prefix, "log1p.mtx")
    mmwrite(out_fn, csr_matrix(np.log1p(mtx)))

    print("pf -> log1p")
    out_fn = os.path.join(out_prefix, "pf_log1p.mtx")
    mmwrite(out_fn, csr_matrix(np.log1p(do_pf(mtx))))

    out_fn = os.path.join(out_prefix, "pf_log1p_pf.mtx")
    mmwrite(out_fn, csr_matrix(do_log1p_pf(do_pf(mtx))))

    print("cpm -> log1p -> scale")
    out_fn = os.path.join(out_prefix, "cpm_log1p_scale.csv")
    pd.DataFrame(scale(np.log1p(do_pf(mtx, target_sum=1_000_000)))).to_csv(
        out_fn, index=False, header=False
    )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
