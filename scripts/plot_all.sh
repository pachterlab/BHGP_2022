#!/usr/bin/env bash


mknorm(){
    echo $1 && ../scripts/plot_all.py $1 $1/ $1/
}

N=5
(
for d in *; do
    if [ ! -f $d/${d}_method_comparison.pdf ] && [ -f $d/raw.mtx.gz ] && [ -f $d/sctransform.csv.gz ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
)
