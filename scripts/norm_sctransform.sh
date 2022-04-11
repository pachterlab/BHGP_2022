#!/usr/bin/env bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

mknorm(){
    echo $1 && $parent_path/norm_sctransform.py $1/matrix.mtx.gz $1/ && gzip $1/raw.mtx.gz
}

N=16
path=$(pwd)
OBS=$(find $path -mindepth 1 -maxdepth 1 -type d  \( ! -iname ".*" \) | sed 's|^\./||g')
(
for d in $OBS; do
    if  [ ! -f $d/raw.mtx.gz ] && [ ! -f $d/sctransform.csv ] && [ -f $d/matrix.mtx.gz ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
wait
)
