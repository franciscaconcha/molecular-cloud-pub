#!/bin/bash

# $1: Mcloud
# $2: Rcloud
# $3: results path
# $4: method for star formation, cluster or single

Nsph=(4000 8000 16000 32000)
Nruns=(1 2 3 4 5 6 7 8)

logfile="$3/log_M$1MSun""_R$2pc.txt"
echo "" > $logfile

for Ns in ${Nsph[@]}
do
    count=0
    for r in ${Nruns[@]}
    do
        if [ $(expr $Ns \* $r) -le 32000 ]
        then
            filepath="$3/M$1MSun""_R$2pc_N$Ns""/$r"
            start=`date "+%Y-%m-%d %H:%M:%S"`
            start_s=`date "+%s"`
            echo $filepath
            echo $filepath >> $logfile
            echo "START: $start" >> $logfile
            ../../amuse/amuse.sh cloud_collapse.py --Mcloud $1 --Rcloud $2 --Ncloud $Ns -s $filepath --tend 5 -m $4
            end=`date "+%Y-%m-%d %H:%M:%S"`
            end_s=`date "+%s"`
            echo "END: $end" >> $logfile
            echo "ELAPSED: $(($end_s - $start_s)) s" >> $logfile
            echo "" >> $logfile
        fi
    done
done
