#!/bin/bash

## If you don't have the reference files already downloaded, this will download
## a minimum set and create all the required files.
## Please adjust the variable as required.

JAFFA_DIR=./

## For hg38 and gencode49
GENOME_NAME=hg38
TRANS_NAME=gencode49
wget --content-disposition https://figshare.com/ndownloader/files/61624573

## For T2T and liftOverGenes uncomment the three lines below
#GENOME_NAME=hs1
#TRANS_NAME=catLiftOffGenesV1
#wget --content-disposition https://figshare.com/ndownloader/files/61624642

## For mm39 and gencode uncomment the three lines below
#GENOME_NAME=mm39
#TRANS_NAME=gencodeVM38
#wget --content-disposition https://figshare.com/ndownloader/files/61624690

tar -xvf JAFFA_REFERENCE_FILES_*.small.tar.gz

echo "Preparing transcriptome reference files"
#The steps below is uncessary for the offical references files, but is included custom references.
#$JAFFA_DIR/tools/bin/reformat fastawrap=0 in=${TRANSCRIPTOME_FASTA} out=stdout.fa | sed 's/ /__/g' > ${GENOME_NAME}_${TRANS_NAME}.fa

echo "Make the bowtie index for the transcriptome sequences"
$JAFFA_DIR/tools/bin/bowtie2-build ${GENOME_NAME}_${TRANS_NAME}.fa ${GENOME_NAME}_${TRANS_NAME}
echo "Make the blast index for the transcriptome sequences"
$JAFFA_DIR/tools/bin/makeblastdb -in ${GENOME_NAME}_${TRANS_NAME}.fa -dbtype nucl -out ${GENOME_NAME}_${TRANS_NAME}_blast

echo "Masking version of genome with genes masked"
$JAFFA_DIR/tools/bin/bedtools maskfasta -fi ${GENOME_NAME}.fa -fo Masked_${GENOME_NAME}.fa -bed ${GENOME_NAME}_${TRANS_NAME}.bed
echo "Make the bowtie index for the masked genome sequences"
$JAFFA_DIR/tools/bin/bowtie2-build Masked_${GENOME_NAME}.fa Masked_${GENOME_NAME}
