#!/bin/bash

# $1: path for result files
# $2: save path

files="$(ls $1/hydro_gas*)"

for f in $files
do
    /data1/wilhelm/py_envs/amuse_dev/amuse/amuse.sh src/plot_molecular_cloud_collapse_with_sinks.py -f $f -s $2
done
