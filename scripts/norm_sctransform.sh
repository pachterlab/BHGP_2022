#!/usr/bin/env bash


mknorm(){
    echo $1 && ../scripts/norm_sctransform.py $1/matrix.mtx.gz $1/ && gzip $1/raw.mtx.gz
}

N=16
OBS=$(cat fixed | cut -f2 -d'/')
(
for d in $OBS; do
    if  [ ! -f $d/raw.mtx.gz ] && [ ! -f $d/sctransform.csv ] && [ -f $d/matrix.mtx.gz ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
wait
)
