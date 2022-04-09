#!/usr/bin/env bash


mknorm(){
    echo $1 && gzip $1/cp10k_log_scale.csv
}

N=256

OBS=$( ls -d * )

(
for d in $OBS; do
    if [ -f $d/cp10k_log_scale.csv ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
)
