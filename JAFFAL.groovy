/***********************************************************
 ** This is the JAFFA pipeline file for fusion detection
 ** with noisy long read data. For polished long read data, 
 ** use JAFFA_direct.groovy. Run like so:
 **    bpipe run <path_to_this_file> <path_to_fastq/fasta_files>
 ** See our website for details	on running options:
 ** https://github.com/Oshlack/JAFFA/wiki.
 **
 ** Author: Nadia Davidson <nadia.davidson@petermac.org>
 ** Last Update: 2021
 *********************************************************/

codeBase = file(bpipe.Config.config.script).parentFile.absolutePath
load codeBase+"/JAFFA_stages.groovy"

get_fasta = {
   doc "Converting fastqs to fasta"
   output.dir=jaffa_output+branch
   produce(branch+".fasta"){
      exec "$reformat ignorebadquality=t in=$input out=$output threads=$threads ;"
   }
}

minimap2_transcriptome = {
   doc "Aligning candidates to transcriptome using minimap2"
   output.dir=jaffa_output+branch
   produce(branch.toString() +".paf"){
        exec """
           $minimap2 -t $threads -x map-ont -c $transFasta $input > $output1 ;
        """
   }
}

/** CODE NOT USED infer_genome_alignment = {
   doc "Bypassing genomic alignment and infering genome position from transcriptome alignments"
   output.dir=jaffa_output+branch
   produce(branch.toString() +"_genome.psl"){
      exec """
       $bypass_genomic_alignment $transTable $input.txt > $output
       """
       }
}**/

minimap2_genome = {
   doc "Aligning candidates to genome using minimap2"
   output.dir=jaffa_output+branch
   produce(branch.toString()+"_genome.paf",branch.toString()+"_genome.psl"){ 
	exec """
	   $minimap2 -t $threads -x splice -c $genomeFasta $input > $output1;
	   grep \$'\\t+\\t' $output1 | awk -F'\\t' -v OFS="\\t" '{ print \$4-\$3,0,0,0,0,0,0,0,\$5,\$1,\$2,\$3,\$4,\$6,\$7,\$8,\$9,2, 100","\$4-\$3-100",",\$3","\$3+100",",  \$8","\$9-\$4+\$3+100"," }' > $output2 ;
	   grep \$'\\t-\\t' $output1 | awk -F'\\t' -v OFS="\\t" '{ print \$4-\$3,0,0,0,0,0,0,0,\$5,\$1,\$2,\$3,\$4,\$6,\$7,\$8,\$9,2, 100","\$4-\$3-100",", \$2-\$4","\$2-\$4+100",", \$8","\$9-\$4+\$3+100"," }' >> $output2 ;
        """
   }
}

report_3_gene_fusions = {
   doc "Checking for reads that support multi-fusion transcripts"
   output.dir=jaffa_output+branch.toString()
   produce(branch.toString()+".3gene_summary",branch.toString()+".3gene_reads"){
      exec """
         $make_3_gene_fusion_table $input.summary $input.txt $output2 > $output1
      """
   }
}

reassign_dist=50

readLayout="single"
fastqInputFormat="%.gz"

common_steps = segment { 
   minimap2_transcriptome + 
   filter_transcripts +
   extract_fusion_sequences +
   minimap2_genome + 
   make_fasta_reads_table +
   get_final_list +
   report_3_gene_fusions
}


// below is the pipeline for a fasta file
if(args[0].endsWith(fastaSuffix)) {
   run { run_check + fastaInputFormat * [
	    common_steps ] + compile_all_results 
   } 
} else { //or fastq.gz will be converted to fasta.
   run { run_check + fastqInputFormat * [
	    get_fasta + common_steps ] + compile_all_results 
   }
}
