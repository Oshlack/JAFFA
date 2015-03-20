#!/bin/bash

velveth=$1
velvetg=$2
oases=$3
base=$4
output=$5
Kall_original=$6
Kseq=`echo $Kall_original | awk 'BEGIN {FS=","} {print $1,$3,$2 }'`
Kall=`seq $Kseq` #convert to a list of numbers
Kmerge=$7
transLength=$8
#export OMP_THREAD_LIMIT=$9
#export OMP_NUM_THREADS=$9
inputs="../../${10}"
if [ $# -eq 11 ] ; then
   inputs="$inputs ../../${11}"
fi


# Lets start by checking that the input k-mer lengths are supported
Kmin=`$velveth | grep "MAXKMERLENGTH" | cut -d"=" -f2`
Ko=`$oases | grep "MAXKMERLENGTH" | cut -d"=" -f2`
if (( $Ko < $Kmin )) ; then Kmin=$Ko ; fi
for k in $Kall ; do
    if (( $k <= $Kmin )) ; then Ks="$Ks $k"
    else 
       echo "k=$k is above the minimum k-mer of your velvet/oases (${Kmin})."
       echo "Lets keep going without it..."
    fi
done

#now start the real stuff
mkdir ${base}/oases ; cd ${base}/oases

echo "$velveth directory $Kall_original -fastq -separate $inputs >> log_velveth"
$velveth directory $Kall_original -fastq -separate $inputs >> log_velveth

for k in $Ks ; do
    echo "Running Assembly for k="$k
    echo "directory_${k}" >> log_${k}
    echo "running $velvetg" >> log_${k}
    $velvetg directory_$k -read_trkg yes >> log_${k}
    echo "running $oases" >> log_${k}
    $oases directory_$k >> log_${k} 
    echo "done" >> log_${k}
done

echo "running velveth to merge with Kmerge = $Kmerge " > log_merge
$velveth mergedAssembly $Kmerge -long directory*/transcripts.fa >> log_merge
echo "running velvetg" >> log_merge
$velvetg mergedAssembly -read_trkg yes -conserveLong yes >> log_merge
echo "running oases" >> log_merge
$oases mergedAssembly -merge yes -min_trans_lgth $transLength >> log_merge
echo "done" >> log_merge
cd ../../

#check whether the jobs completed successfully.
for k in $Ks ; do
   if ! [ -s ${base}/oases/directory_${k}/transcripts.fa  ] ; then
      echo "Oases K-mer=${k} assembly failed" 
   else
      echo "Oases K-mer=${k} assembly succeeded" 
      rm -rf ${base}/oases/directory_${k}
   fi
done

if ! [ -s ${base}/oases/mergedAssembly/transcripts.fa  ] ; then
   echo "Oases mergeAssembly failed"
   exit 1
else
   echo "Oases mergeAssembly succeeded"
   #move the output if it was a success
   mv ${base}/oases/mergedAssembly/transcripts.fa ${output}
   rm -rf ${base}/oases/mergedAssembly
fi

