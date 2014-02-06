#!/bin/bash

base=$1
input1=$2
input2=$3
output=$4
Ks=`echo $5 | sed s'/,/ /g'` #convert comma seperated to space separated
Kmerge=$6
transLength=$7

mkdir ${base}/oases ; cd ${base}/oases
for k in $Ks ; do
    echo "directory_${k}" >> log_${k}
    echo "running velveth" >> log_${k}
    velveth directory_${k} $k -fastq -separate ../../$input1 ../../$input2 >> log_${k}
    echo "running velvetg" >> log_${k}
    velvetg directory_$k -read_trkg yes >> log_${k}
    echo "running oases" >> log_${k}
    oases directory_$k >> log_${k} 
    echo "done" >> log_${k}
done

echo "running velveth" >> log
velveth mergedAssembly $Kmerge -long directory*/transcripts.fa >> log
echo "running velvetg" >> log
velvetg mergedAssembly -read_trkg yes -conserveLong yes >> log
echo "running oases" >> log
oases mergedAssembly -merge yes -min_trans_lgth $transLength >> log
echo "done" >> log
cd ../../

mv ${base}/oases/mergedAssembly/transcripts.fa ${output}
