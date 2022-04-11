#!/usr/bin/env bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

mknorm(){
    echo $1 && $parent_path/metrics_methods_sparse.py ${1}/raw.mtx.gz ${1}/${2}.mtx.gz ${1}/${1}_metrics.json
}

N=256
path=$(pwd)
OBS=$(find $path -mindepth 1 -maxdepth 1 -type d  \( ! -iname ".*" \) | sed 's|^\./||g')
(
for MTX in raw pf log sqrt pf_log pf_log_pf cpm_log cp10k_log; do
    echo $MTX &&
    for d in $OBS; do
        if [ -f $d/$MTX.mtx.gz ]; then
            ((i=i%N)); ((i++==0)) && wait
            mknorm $d $MTX &
        fi
    done
    wait
done
)
