#!/usr/bin/env bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

mknorm(){
    echo $1 && $parent_path/plot_all.py $1 $1/ $1/
}

N=5
path=$(pwd)
OBS=$(find $path -mindepth 1 -maxdepth 1 -type d  \( ! -iname ".*" \) | sed 's|^\./||g')
(
for d in $OBS; do
    if [ ! -f $d/${d}_method_comparison.pdf ] && [ -f $d/raw.mtx.gz ] && [ -f $d/sctransform.csv.gz ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
)
