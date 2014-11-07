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
 ** Last Update: 2nd April 2014
 *********************************************************/

//commands="trimmomatic oases velveth velvetg R bowtie2 blat fasta_formatter samtools fastx_collapser"

codeBase = file(bpipe.Config.config.script).parentFile.absolutePath
load codeBase+"/JAFFA_stages.groovy"

//get_unmapped_as_fasta = segment { cat_reads + remove_dup + get_assembly_unmapped }
get_unmapped_as_fasta = segment { get_assembly_unmapped }

if(readLayout=="single"){ fastqInputFormat="%.gz" }
run{ run_check + fastqInputFormat * [ 
   		      make_dir_using_fastq_names +
      		      prepare_reads + 
		      run_assembly + //start the first part - assembly
		      align_transcripts_to_annotation +
		      filter_transcripts + 
		      extract_fusion_sequences + 
		      map_reads + 
		      get_spanning_reads +
		      get_unmapped_as_fasta + // start the second part - unmapped reads
		      align_reads_to_annotation + 
		      filter_transcripts + 
		      extract_fusion_sequences + 
		      make_simple_reads_table +
		      merge_assembly_and_unmapped_reads_candidates + //combine...
		      align_transcripts_to_genome + 
		      get_final_list 
		 ] + compile_all_results.using(type:".all")
}

