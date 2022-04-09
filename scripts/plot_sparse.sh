#!/usr/bin/env bash


mknorm(){
    echo $1 && ../scripts/plot_sparse.py $1 $1/sparse/ $1/
}

N=256
(
for d in *; do
    if [ -f $d/sparse/raw.mtx.gz ] && [ ! -f $d/${d}_full_method_comparison.pdf ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
)
