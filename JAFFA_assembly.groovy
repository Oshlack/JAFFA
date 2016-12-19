/***********************************************************
 ** This is the JAFFA pipeline file for fusion detection
 ** via read assembly. Run like so:
 **    bpipe run <path_to_this_file> <path_to_fastq_files>
 ** See our website for details on running options: 
 ** https://code.google.com/p/jaffa-project/.
 **
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 6th Feb 2014 
 *********************************************************/ 

//these are the commands we will check for at the start of every run.
//commands="trimmomatic oases velveth velvetg R bowtie2 blat fasta_formatter samtools"

codeBase = file(bpipe.Config.config.script).parentFile.absolutePath
load codeBase+"/JAFFA_stages.groovy"

// The actual pipeline.
if(readLayout=="single"){ fastqInputFormat="%.gz" }
run{ run_check + fastqInputFormat * [ 
		      prepare_reads +
		      run_assembly +
		      align_transcripts_to_annotation.using(tile:contigTile) + 
		      filter_transcripts + 
		      extract_fusion_sequences + 
		      map_reads + 
		      get_spanning_reads +
		      align_transcripts_to_genome + 
		      get_final_list 
		      ] + compile_all_results 
}
