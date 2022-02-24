#!/usr/bin/env python3

import sys
import os
from utils import do_pf, do_log_pf
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
from sklearn.preprocessing import scale
import pandas as pd

import numpy as np


def main(in_matrix_fn, out_prefix):
    # pseudocount
    pc = 0.5

    mtx = mmread(in_matrix_fn).toarray()

    print("pf")
    out_fn = os.path.join(out_prefix, "pf.mtx")
    mmwrite(out_fn, csr_matrix(do_pf(mtx)))

    print("log")
    out_fn = os.path.join(out_prefix, "log.mtx")
    mmwrite(out_fn, csr_matrix(np.log(pc + mtx)))

    print("pf -> log")
    out_fn = os.path.join(out_prefix, "pf_log.mtx")
    mmwrite(out_fn, csr_matrix(np.log(pc + do_pf(mtx))))

    print("pf -> log -> pf")
    out_fn = os.path.join(out_prefix, "pf_log_pf.mtx")
    mmwrite(out_fn, csr_matrix(do_log_pf(do_pf(mtx))))

    print("cp10k -> log")
    out_fn = os.path.join(out_prefix, "cp10k_log.mtx")
    mmwrite(out_fn, csr_matrix(np.log(pc + do_pf(mtx, target_sum=10_000))))

    print("cp10k -> log -> scale")
    out_fn = os.path.join(out_prefix, "cp10k_log_scale.csv")
    pd.DataFrame(scale(np.log(pc + do_pf(mtx, target_sum=10_000)))).to_csv(
        out_fn, index=False, header=False
    )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
