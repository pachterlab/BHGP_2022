#!/usr/bin/env bash


mknorm(){
    echo $1 && ../scripts/norm_sparse.py $1/raw.mtx.gz $1/ && gzip $1/*.mtx
}

N=256

OBS=$(cat fixed | cut -f2 -d'/')

(
for d in $OBS; do
    if [ -f $d/raw.mtx.gz ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
)
