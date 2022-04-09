#!/usr/bin/env bash


mknorm(){
    echo $1 && ../scripts/norm_cp10k_log_scale.py $1/cp10k_log.mtx.gz $1/ 
}

N=16
OBS=$(cat fixed | cut -f2 -d'/')
(
for d in $OBS; do
    if  [ -f $d/cp10k_log.mtx.gz ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
)
