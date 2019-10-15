#!/bin/bash

files="$(ls $1/hydro_gas*)"

for f in $files
do
    ../amuse/amuse.sh src/plot_molecular_cloud_collapse_with_sinks.py -f $f -s $2
done
