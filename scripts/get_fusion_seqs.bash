#!/bin/bash

if [[ -z ${@} ]]
then
  echo "USAGE: get_fusion_seqs <FUSION EVEN LINE> [OUTPUT NAME]"
  echo ""
  echo "Typically this script is used inside a while loop, e.g"
  echo ""
  echo "while read l;do get_fusion_seqs \$l;done < jaffa_results.csv"
  echo ""
  exit
fi

function get_fusion_seqs() {

  type=""

  IFS=","
  tokens=(${@})

  field1=$(echo ${tokens[0]} | tr -d '"')
  field2=$(echo ${tokens[1]} | tr -d '"')
  field3=$(echo ${tokens[14]} | tr -d '"')
  field4=$(echo ${tokens[15]} | tr -d '"')

  #DEBUG
  #echo "---"
  #echo "Processing ${field1}"
  #echo "Processing ${field2}"
  #echo "Processing ${field3}"
  #echo "Processing ${field4}"
  #echo "---"

  res=${tokens[18]}

  if [[ ${field1} =~ "sample" ]]
  then
    return
  fi

  if [[ -z ${res} ]]
  then
    res="jaffa_results"
  fi

  res="${res}.fasta"

  fusions_file=${field1}/${field1}${type}.fusions.fa

  new_id=">${field1}|${field2}|${field3}"
  echo ${new_id} >> ${res}
  brk_pnt=${field4}

  #NOTE not sure why grep ^> isn't working
  #sequence=$(grep -A1 "^>${field3}" ${fusions_file} | grep -v "^>")
  sequence=$(grep -A1 "^>${field3}" ${fusions_file} | grep -v ">")

  start=$(echo ${sequence} | cut -c 1-$((${brk_pnt}-1)))
  middle=$(echo ${sequence} | cut -c ${brk_pnt}-$((${brk_pnt}+1)) | tr '[:upper:]' '[:lower:]')
  string_length=$(echo ${#sequence})

  end=$(echo ${sequence} | cut -c $((${brk_pnt}+2))-${string_length})

  echo ${start}${middle}${end} >> ${res}
}

get_fusion_seqs ${@}
