#!/usr/bin/env bash

mknorm(){
    echo $1 && ../scripts/metrics_methods_dense.py ${1}/raw.mtx.gz ${1}/${2}.csv.gz ${1}/${1}_metrics.json
}

N=8
OBS=$(cat final | cut -f2 -d '/')
(
for MTX in sctransform cp10k_log_scale; do
    echo $MTX &&
    for d in $OBS; do
        if [ -f $d/$MTX.csv.gz ]; then
            ((i=i%N)); ((i++==0)) && wait
            mknorm $d $MTX &
        fi
    done
    wait
done
)
