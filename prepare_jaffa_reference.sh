#!/bin/bash
set -euo pipefail

# prepare_jaffa_reference.sh
#
# Single, self-contained JAFFA reference builder. One script, one build,
# usable by every pipeline variant:
#   JAFFA_direct.groovy / JAFFA_assembly.groovy / JAFFA_hybrid.groovy
#                                       (short-read: bowtie2 + BLAST + masked)
#   JAFFAL.groovy                       (long-read: minimap2)
#
# This builds a reference for most organism with reasonable annotation quality.
# The shared helper functions and the full build flow all
# live in this one file. It builds the transcriptome FASTA, .tab coordinate
# table, exon BED, masked genome, bowtie2 indexes, BLAST db, and the
# minimap2 genome index, so a single reference serves both the short-read
# and long-read pipelines and run_check can verify all relevant files
# unconditionally.
#
# Self-contained: sanitizes only unsafe sequence names (pass-through by
# default), then builds the reference directly, following the official
# Oshlack/JAFFA wiki recipe for non-supported genomes. All classification /
# sanitization / validation logic runs through the compiled
# prepare_ref_helper binary.
#
# Unsafe conditions checked per sequence name:
#   - contains "|"   (breaks BLAST/samtools accession-style parsing)
#   - contains ":"   (breaks make_final_table's chrom1:base1:chrom2:base2
#                      breakpoint grouping key)
#   - contains "/"   (breaks sed-delimiter reverse-mapping commands)
#   - duplicate name (breaks unique lookup by chromosome)
#
# Usage:
#   bash prepare_jaffa_reference.sh \
#       <JAFFA_PATH> <REF_DIR> <GENOME_NAME> <ANNOTATION_NAME> [THREADS]
#
# The genome FASTA is expected at <REF_DIR>/<GENOME_NAME>.fa
#
# THREADS is optional (default 1) and only affects bowtie2-build -- the
# other tools here don't meaningfully parallelize at this scale.
#
# Requires <REF_DIR>/<GENOME_NAME>_<ANNOTATION_NAME>_raw.gtf (or _raw.gff3)
# to already exist. Requires tools/bin/prepare_ref_helper to be compiled
# first:
#   g++ -std=c++11 -O2 -o tools/bin/prepare_ref_helper src/prepare_ref_helper.c++
#
# NOTE: since transcript headers are plain ">transcriptID" with no
# __range= suffix, JAFFA_stages.groovy's anno_prefix must be:
#   anno_prefix="'(.*)'"

usage() {
  echo "Usage: $0 <JAFFA_PATH> <REF_DIR> <GENOME_NAME> <ANNOTATION_NAME> [THREADS]"
  exit 1
}

[ $# -eq 4 ] || [ $# -eq 5 ] || usage

JAFFA_PATH=$1
REF_DIR=$2
GENOME_NAME=$3
ANNOTATION_NAME=$4
THREADS=${5:-1}

GENOME_FA=$REF_DIR/${GENOME_NAME}.fa

TOOLS=$JAFFA_PATH/tools/bin
HELPER=$TOOLS/prepare_ref_helper

PREFIX=${GENOME_NAME}_${ANNOTATION_NAME}
RAW_GTF=$REF_DIR/${PREFIX}_raw.gtf
RAW_GFF=$REF_DIR/${PREFIX}_raw.gff3
SANITIZED_FA=$REF_DIR/${GENOME_NAME}_sanitized.fa
SANITIZED_ANNO=$REF_DIR/${PREFIX}_sanitized.gtf
RENAME_LOG=$REF_DIR/${GENOME_NAME}_sanitize_log.tsv

############################################################################
# Helper functions
############################################################################

check_helper() {
  if [ ! -x "$HELPER" ]; then
    echo "ERROR: $HELPER not found or not executable."
    echo "Compile it first:"
    echo "  g++ -std=c++11 -O2 -o $HELPER $JAFFA_PATH/src/prepare_ref_helper.c++"
    exit 1
  fi
}

resolve_annotation() {
  if [[ -f "$RAW_GTF" ]]; then
    echo "$RAW_GTF"
  elif [[ -f "$RAW_GFF" ]]; then
    echo "  GFF3 detected, converting to GTF first..." >&2
    "$TOOLS/gffread_bin" "$RAW_GFF" -T -o "${RAW_GFF%.gff3}.gtf"
    echo "${RAW_GFF%.gff3}.gtf"
  else
    echo "ERROR: no raw GTF or GFF3 found at $RAW_GTF or $RAW_GFF" >&2
    exit 1
  fi
}

# Structural validation for GFF3 input: warns on non-unique IDs / missing
# Parent= on some children, hard-fails if Parent= is absent on all children
# (gffread cannot associate exons with transcripts without it). No-op for
# GTF input.
check_gff3_structure() {
  local GFF3=$1

  case "$GFF3" in
    *.gff3|*.gff) ;;
    *) return 0 ;;
  esac

  echo "=== Checking GFF3 structural integrity (ID/Parent linkage) ===" >&2

  local total_ids dup_ids parent_lines child_lines
  total_ids=$(awk -F'\t' '$3=="mRNA" || $3=="transcript"' "$GFF3" | grep -oP 'ID=\K[^;]+' | wc -l)
  dup_ids=$(awk -F'\t' '$3=="mRNA" || $3=="transcript"' "$GFF3" | grep -oP 'ID=\K[^;]+' | sort | uniq -d | wc -l)
  child_lines=$(awk -F'\t' '$3=="CDS" || $3=="exon"' "$GFF3" | wc -l)
  parent_lines=$(awk -F'\t' '($3=="CDS" || $3=="exon") && /Parent=/' "$GFF3" | wc -l)

  local problems=0

  if [ "$total_ids" -gt 0 ] && [ "$dup_ids" -gt 0 ]; then
    echo "  WARNING: $dup_ids of $total_ids ID= values are non-unique." >&2
    echo "           gffread may mis-group features or emit colliding" >&2
    echo "           transcript IDs downstream." >&2
    problems=1
  fi

  if [ "$child_lines" -gt 0 ] && [ "$parent_lines" -eq 0 ]; then
    echo "  ERROR: no Parent= attribute found on any of $child_lines" >&2
    echo "         CDS/exon feature(s). gffread cannot associate exons" >&2
    echo "         with transcripts without Parent= linkage -- this is" >&2
    echo "         not a standard-conformant GFF3 and cannot be built" >&2
    echo "         from directly." >&2
    problems=2
  elif [ "$child_lines" -gt 0 ]; then
    local missing=$((child_lines - parent_lines))
    if [ "$missing" -gt 0 ]; then
      echo "  WARNING: $missing of $child_lines CDS/exon feature(s) are" >&2
      echo "           missing Parent=. Those records will likely be" >&2
      echo "           dropped or mis-grouped by gffread." >&2
      problems=1
    fi
  fi

  if [ "$problems" -eq 2 ]; then
    echo "" >&2
    echo "This GFF3 needs a manual, organism-specific repair before it can" >&2
    echo "be used -- e.g. deriving a unique transcript key from another" >&2
    echo "attribute (locus_tag, protein_id) if one is present and" >&2
    echo "correlates 1:1 with transcripts, and/or sourcing transcript" >&2
    echo "sequences from a separately-provided CDS/transcript FASTA rather" >&2
    echo "than deriving them via gffread from this file." >&2
    exit 1
  elif [ "$problems" -eq 1 ]; then
    echo "  Proceeding, but expect gffread/validate to surface problems." >&2
  else
    echo "  OK: ID= values unique, Parent= present on all child features." >&2
  fi
}

# Sets GFA_OUT and GTF_OUT to the genome/annotation to actually build from
# (original, or sanitized copies if unsafe characters were found).
classify_and_sanitize() {
  echo "=== Checking for unsafe characters in sequence names ==="
  local NEEDS_SANITIZE
  NEEDS_SANITIZE=$("$HELPER" classify "$GENOME_FA")
  echo "  Sanitization needed: $NEEDS_SANITIZE"

  local INPUT_ANNO
  INPUT_ANNO=$(resolve_annotation)
  check_gff3_structure "$INPUT_ANNO"

  if [ "$NEEDS_SANITIZE" == "NO" ]; then
    echo "=== Genome names are already safe -- using original files directly ==="
    GFA_OUT="$GENOME_FA"
    GTF_OUT="$INPUT_ANNO"
    SANITIZED=0
    return 0
  fi

  echo "=== Sanitizing only unsafe names (minimal touch) ==="
  if [ ! -f "$SANITIZED_FA" ]; then
    "$HELPER" sanitize "$GENOME_FA" "$SANITIZED_FA" "$RENAME_LOG"
  fi
  if [ ! -f "$SANITIZED_ANNO" ]; then
    "$HELPER" rename-anno "$INPUT_ANNO" "$RENAME_LOG" "$SANITIZED_ANNO"
  fi
  GFA_OUT="$SANITIZED_FA"
  GTF_OUT="$SANITIZED_ANNO"
  SANITIZED=1
}

# Extract transcript sequences, reformat headers, validate uniqueness.
# Sets TRANS_FA to the resulting plain-header transcript FASTA.
build_transcript_fasta() {
  local GFA=$1
  local GTF=$2
  local TRANS_ANNO_FA=$REF_DIR/${PREFIX}_anno.fasta
  TRANS_FA=$REF_DIR/${PREFIX}.fa

  echo "=== Extracting transcript sequences with gffread ==="
  if [ ! -f "$TRANS_ANNO_FA" ]; then
    "$TOOLS/gffread_bin" "$GTF" -g "$GFA" -w "$TRANS_ANNO_FA"
  fi

  echo "=== Reformatting transcript FASTA headers (ID only, single line) ==="
  if [ ! -f "$TRANS_FA" ]; then
    cat "$TRANS_ANNO_FA" \
      | awk '{if(/^>/){print $1}else{print}}' \
      | "$TOOLS/reformat" fastawrap=0 in=stdin.fa out=stdout.fa > "$TRANS_FA"
  fi

  echo "=== Validating transcript IDs are unique ==="
  "$HELPER" validate "$TRANS_FA"
}

# Sets TAB_FILE to the .tab gene coordinate table.
build_tab_file() {
  local GTF=$1
  local GPD_FILE=$REF_DIR/${PREFIX}.gpd
  TAB_FILE=$REF_DIR/${PREFIX}.tab

  echo "=== Building .tab gene coordinate table ==="
  if [ ! -f "$TAB_FILE" ]; then
    # Ensembl GTFs omit gene_name for genes without a symbol (mostly lncRNAs), which
    # leaves name2 empty and JAFFA reports those partners as chrom:pos. Fall back to
    # gene_id, as GENCODE GTFs already do, so they are reported as ENSG... instead.
    local NAMED_GTF=$REF_DIR/${PREFIX}.gene_name_filled.gtf
    awk 'BEGIN{FS=OFS="\t"} !/^#/ && $9 !~ /gene_name "/ && match($9, /gene_id "[^"]+"/) {$9 = $9 " gene_name " substr($9, RSTART+8, RLENGTH-8) ";"} {print}' "$GTF" > "$NAMED_GTF"
    GTF=$NAMED_GTF
    "$TOOLS/gtfToGenePred" "$GTF" "$GPD_FILE" -genePredExt -geneNameAsName2
    rm -f "$NAMED_GTF"
    awk 'BEGIN{FS=OFS="\t"; print "#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames"} {print "1",$0}' "$GPD_FILE" > "$TAB_FILE"
  fi
}

# Sets BED_FILE to the sorted exon BED. Used as masking input (short-read
# pipeline) or as --junc-bed input (JAFFAL).
build_exon_bed() {
  local GTF=$1
  BED_FILE=$REF_DIR/${PREFIX}.bed

  echo "=== Building exon BED file ==="
  if [ ! -f "$BED_FILE" ]; then
    grep -P "\texon\t" "$GTF" \
      | awk '{print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' \
      | "$TOOLS/bedtools" sort -i - > "$BED_FILE"
  fi
}

write_placeholder_known_fusion_tables() {
  for f in known_mitelman known_cosmic known_cosmic_tier known_gtex; do
    local OUT=$REF_DIR/${f}.txt
    [ -f "$OUT" ] || touch "$OUT"
  done
}

############################################################################
# Main build flow
############################################################################

check_helper

if [ ! -x "$TOOLS/minimap2" ]; then
  echo "ERROR: $TOOLS/minimap2 not found or not executable."
  echo "This reference builds a minimap2 index for the long-read (JAFFAL)"
  echo "pipeline -- see install_linux64.sh."
  exit 1
fi

mkdir -p "$REF_DIR"

classify_and_sanitize   # sets GFA_OUT, GTF_OUT, SANITIZED

build_transcript_fasta "$GFA_OUT" "$GTF_OUT"   # sets TRANS_FA
build_tab_file "$GTF_OUT"                       # sets TAB_FILE
build_exon_bed "$GTF_OUT"                       # sets BED_FILE -- masking
                                                # input (short-read) and
                                                # --junc-bed input (JAFFAL)

MASKED_FA=$REF_DIR/Masked_${GENOME_NAME}.fa
BLAST_PREFIX=$REF_DIR/${PREFIX}_blast

echo "=== Masking genome and building bowtie2 indexes (threads=$THREADS) ==="
if [ ! -f "$MASKED_FA" ]; then
  "$TOOLS/bedtools" maskfasta -fi "$GFA_OUT" -fo "$MASKED_FA" -bed "$BED_FILE"
fi

if [ ! -f "${TRANS_FA%.fa}.1.bt2" ]; then
  "$TOOLS/bowtie2-build" --threads "$THREADS" "$TRANS_FA" "${TRANS_FA%.fa}"
fi

if [ ! -f "${MASKED_FA%.fa}.1.bt2" ]; then
  "$TOOLS/bowtie2-build" --threads "$THREADS" "$MASKED_FA" "${MASKED_FA%.fa}"
fi

echo "=== Building BLAST db ==="
if [ ! -f "${BLAST_PREFIX}.nsq" ]; then
  "$TOOLS/makeblastdb" -in "$TRANS_FA" -dbtype nucl -out "$BLAST_PREFIX"
fi

echo "=== Building minimap2 genome index (for JAFFAL) ==="
GENOME_MMI=$REF_DIR/${GENOME_NAME}.mmi
if [ ! -f "$GENOME_MMI" ]; then
  "$TOOLS/minimap2" -x splice -d "$GENOME_MMI" "$GFA_OUT"
fi

# JAFFA_stages.groovy hardcodes genomeFasta=<fastaBase>/<genome>.fa; make
# sure a ".fa"-named genome exists even if the source was ".fasta"/other.
GENOME_FA_LINK=$REF_DIR/${GENOME_NAME}.fa
if [ ! -e "$GENOME_FA_LINK" ]; then
  ln -s "$(cd "$(dirname "$GFA_OUT")" && pwd)/$(basename "$GFA_OUT")" "$GENOME_FA_LINK"
fi

write_placeholder_known_fusion_tables

echo "Done. JAFFA reference written to $REF_DIR with prefix '$PREFIX'"
echo "  Serves both short-read (direct/assembly/hybrid) and long-read (JAFFAL) pipelines."
echo "  Transcriptome FASTA:   $TRANS_FA"
echo "  Masked genome:         $MASKED_FA (+ bowtie2 indexes)"
echo "  BLAST db:              $BLAST_PREFIX"
echo "  minimap2 genome index: $GENOME_MMI"
echo "  Junction BED:          $BED_FILE"
echo "anno_prefix should be: anno_prefix=\"'(.*)'\""
if [ "$SANITIZED" -eq 1 ]; then
  echo "Names with unsafe characters were changed. See $RENAME_LOG for the (small) list."
  echo "To reverse in output files:"
  echo "  awk -F'\t' 'NR==FNR{m[\$2]=\$1;next}{for(k in m) gsub(k,m[k]);print}' $RENAME_LOG jaffa_results.csv > jaffa_results_original_names.csv"
else
  echo "No renaming was applied -- output uses original chromosome names directly."
fi
