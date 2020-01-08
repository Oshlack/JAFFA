/***********************************************************
 ** This is the JAFFA pipeline file for fusion detection
 ** without read assembly. Run like so:
 **    bpipe run <path_to_this_file> <path_to_fastq/fasta_files>
 ** See our website for details	on running options:
 ** https://code.google.com/p/jaffa-project/.
 **
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 2020
 *********************************************************/

codeBase = file(bpipe.Config.config.script).parentFile.absolutePath
load codeBase+"/JAFFA_stages.groovy"

common_steps = segment { filter_transcripts +
                 extract_fusion_sequences +
                 align_transcripts_to_genome }

get_fasta = {
   doc "Converting fastqs to fasta"
   output.dir=jaffa_output+branch
   produce(branch+".fasta"){
      exec "$reformat in=$input out=$output threads=$threads ;"
   }
}

minIdTrans="80"
maxIntron=100
contigTile=11
readTile=11
fastqInputFormat="%.fastq.gz"
readLayout="single"
run { run_check + fastqInputFormat * [
	    get_fasta +
	    align_reads_to_annotation +
	    common_steps +
	    make_fasta_reads_table +
	    get_final_list ] + compile_all_results 
}

