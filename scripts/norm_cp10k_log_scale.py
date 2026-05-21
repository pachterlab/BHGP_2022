#!/usr/bin/env python3

import sys
import os
from utils import cps
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
from sklearn.preprocessing import scale
import pandas as pd

import numpy as np


def main(in_matrix_fn, out_prefix):
    mtx = mmread(in_matrix_fn).toarray()
 
    (data, _, _, _) = cps(mtx)

    title = "cp10k_log_scale"
    print(f"saving {title}")
    out_fn = os.path.join(out_prefix, f"{title}.csv.gz")
    pd.DataFrame(data[title]).to_csv(out_fn, index=False, header=False, compression="gzip")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
