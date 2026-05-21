#!/usr/bin/env bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

mknorm(){
    echo $1 && $parent_path/metrics_methods_dense.py ${1}/raw.mtx.gz ${1}/${2}.csv.gz ${1}/${1}_metrics.json
}

N=8

path=$(pwd)
OBS=$(find $path -mindepth 1 -maxdepth 1 -type d  \( ! -iname ".*" \) | sed 's|^\./||g')

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
