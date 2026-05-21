#!/usr/bin/env bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

mknorm(){
    echo $1 && $parent_path/plot_sparse.py $1 $1/sparse/ $1/
}

N=256
path=$(pwd)
OBS=$(find $path -mindepth 1 -maxdepth 1 -type d  \( ! -iname ".*" \) | sed 's|^\./||g')
(
for d in *; do
    if [ -f $d/sparse/raw.mtx.gz ] && [ ! -f $d/${d}_full_method_comparison.pdf ]; then
        ((i=i%N)); ((i++==0)) && wait
        mknorm $d &
    fi
done
)
