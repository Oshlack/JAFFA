#!/bin/bash

## Specifiy the paths to the required files
GENOME_FASTA=hg38.fa
TRANSCRIPTOME_FASTA=hg38_gencode49.fa
TRANSCRIPTOME_TABLE=hg38_gencode49.tab
TRANSCRIPTOME_BED=hg38_gencode49.bed

GENOME_NAME=hg38
TRANS_NAME=gencode49
JAFFA_DIR=./

echo "Preparing transcriptome reference files"

#Steps below are uncessary for the offical references files, but may need to be altered for custom references.
#rename files
#echo "Renaming files"
#cp $TRANSCRIPTOME_BED ${GENOME_NAME}_${TRANS_NAME}.bed
#cp $TRANSCRIPTOME_TABLE ${GENOME_NAME}_${TRANS_NAME}.tab
#echo "Reformatting the transcriptome fasta to one line per transcript"
#$JAFFA_DIR/tools/bin/reformat fastawrap=0 in=${TRANSCRIPTOME_FASTA} out=stdout.fa | sed 's/ /__/g' > ${GENOME_NAME}_${TRANS_NAME}.fa

echo "Make the bowtie index for the transcriptome sequences"
$JAFFA_DIR/tools/bin/bowtie2-build ${GENOME_NAME}_${TRANS_NAME}.fa ${GENOME_NAME}_${TRANS_NAME}
echo "Make the blast index for the transcriptome sequences"
$JAFFA_DIR/tools/bin/makeblastdb -in $TRANSCRIPTOME_FASTA -dbtype nucl -out ${GENOME_NAME}_${TRANS_NAME}_blast

echo "Masking version of genome with genes masked"
$JAFFA_DIR/tools/bin/bedtools maskfasta -fi $GENOME_FASTA -fo Masked_${GENOME_NAME}.fa -bed ${GENOME_NAME}_${TRANS_NAME}.bed
echo "Make the bowtie index for the masked genome sequences"
$JAFFA_DIR/tools/bin/bowtie2-build Masked_${GENOME_NAME}.fa Masked_${GENOME_NAME}
