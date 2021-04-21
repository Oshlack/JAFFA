#!/bin/bash

## This simple script will match fusions to cell barcodes 
## It requires the single cell long read data to first be processed with FLAMES
## method for barcode extraction, ./match_cell_barcode
## JAFFAL should then be run on the output fastq.gz reads file.
## A list of fusions of interest (identified from JAFFAL or otherwise)
## should be made and then this script run with the intermediate output file
## <sample>/<sample.txt>

if [[ -z ${@} ]]
then
  echo ""
  echo "USAGE: ./get_cell_barcodes_by_fusion.bash <list of fusions> <sample/sample.txt>"
  echo ""
  echo "Where <list of fusions> is a list of fusions with colon separated gene names. e.g. BCR:ABL1"
  echo "and <sample> is a long read single cell sample which has been processed by"
  echo "FLAMES method, match_cell_barcode, as well as JAFFAL" 
  exit
fi

cat $1 | while read fusion ; do
   grep $fusion $2 | cut -f1 -d "_" | sort -u | while read barcode ; do
       echo -e "$fusion\t$barcode"
   done
done
