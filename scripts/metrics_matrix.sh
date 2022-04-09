#!/usr/bin/env bash


mknorm(){
    echo $1 && ../scripts/metrics_matrix.py ${1}/raw.mtx.gz ${1}/${1}_metrics.json
}

N=256
OBS=$(cat final | cut -f2 -d '/')
(
for d in $OBS; do
    if [ -f $d/raw.mtx.gz ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
)
