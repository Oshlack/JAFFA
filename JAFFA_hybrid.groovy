/***********************************************************
 ** This is the JAFFA pipeline file for fusion detection
 ** which combines assembly with raw reads. ie. it's a 
 ** fusion (sorry :) ) between the JAFFA.groovy and 
 ** JAFF_no_assemble.groovy pipeline.
 ** Run like so:
 **    bpipe run <path_to_this_file> <path_to_fastq/fasta_files>
 ** See our website for details	on running options:
 ** https://code.google.com/p/jaffa-project/.
 **
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 12th Feb 2014
 *********************************************************/

commands="trimmomatic oases velveth velvetg R bowtie2 blat fasta_formatter samtools fastx_collapser"
load "JAFFA_stages.groovy"

run{ run_check + fastq_input_format * [ 
   		      make_dir_using_fastq_names +
      		      prepare_reads + 
		      run_assembly +
		      cat_reads + 
		      remove_dup + 
		      get_unmapped_from_assembly +
		      merge_assembly_with_unmapped +
		      align_transcripts_to_annotation + 
		      filter_transcripts + 
		      extract_fusion_sequences + 
		      map_reads + 
		      get_spanning_reads +
		      align_transcripts_to_genome + 
		      get_final_list 
		      ] + compile_all_results 
}
