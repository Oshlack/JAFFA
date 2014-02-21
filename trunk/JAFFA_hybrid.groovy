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
load "/vlsci/VR0193/shared/nadiad/JAFFA/JAFFA/JAFFA_stages.groovy"

run{ run_check + fastq_input_format * [ 
   		      make_dir_using_fastq_names +
      		      prepare_reads + 
		      run_assembly + //start the first part - assembly
		      align_transcripts_to_annotation + 
		      filter_transcripts + 
		      extract_fusion_sequences + 
		      map_reads + 
		      get_spanning_reads +
		      cat_reads + // start the second part - unmapped reads
		      remove_dup + 
		      get_assembly_unmapped +
		      align_transcripts_to_annotation + 
		      filter_transcripts + 
		      extract_fusion_sequences + 
		      make_simple_reads_table +
		      merge_assembly_and_unmapped_reads_candidates + //combine...
		      align_transcripts_to_genome + 
		      get_final_list 
		 ] + compile_all_results.using(type:".all")
}
