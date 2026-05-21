#!/usr/bin/env bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

mknorm(){
    echo $1 && $parent_path/norm_cp10k_log_scale.py $1/cp10k_log.mtx.gz $1/ 
}

N=16
path=$(pwd)
OBS=$(find $path -mindepth 1 -maxdepth 1 -type d  \( ! -iname ".*" \) | sed 's|^\./||g')
(
for d in $OBS; do
    if  [ -f $d/cp10k_log.mtx.gz ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
)
