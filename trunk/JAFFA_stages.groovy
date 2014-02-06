/*********************************************************************************
 ** This file defines all the JAFFA pipeline stages (used in JAFFA.groovy and
 ** JAFFA_no_assemble.groovy). See our website for details on running, 
 ** https://code.google.com/p/jaffa-project/.
 **
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 6th Feb 2014
 ********************************************************************************/


code_base = new File(bpipe.Config.config.script).parentFile.absolutePath

/**********  Parameters that must be check by the user: *******************/

read_length=75 //Read length
threads=1 //Threads to use when running the pipeline on a single sample. ie. the total threads will be samples*threads

// Genome, Transcriptome and related data paths. You have two options:
// 1) put the full paths below. e.g. hgFasta=<path_to_genome> 
// or 2) leave as is and symlink the data files to the jaffa code directory. e.g. ln -s <path_to_genome> <path_to_jaffa_code_directory>
hgFasta=code_base+"/hg19.fa"  //genome sequence
transFasta=code_base+"/hg19_genCode.fa"  // transcript cDNA sequences
transTable=code_base+"/hg19_genCode.tab" // table of gene coordinates

// Input pattern (see bpipe documentation for how files are grouped and split )
// group on start, split on end. eg. on ReadsA_1.fastq.gz, ReadA_2.fastq.gz
// this would group the read files into pairs like we want.
fastq_filename_pattern="%_*.fastq.gz"

/************** Other configurables *************************************/

//Default output name
output_name="jaffa_results"

// trimming
scores=33
minlen=30
minQScore=10  

// for the assembly
//mem="100G"
contig_length="100"

// for aligning to known genes using blat
minId="98" //98% similar
minScore="30" // this is the minimum required flaking sequence assembled around the break-point
tile="18" //big tile size makes blat faster
maxIntron="0" //don't expect intron when mapping to the transcriptome

// filtering
gapSize="1000" //minimum distance between the two fusion candidates for the 1st filtering stage
finalGapSize="10000" //minimum distance for the final filtering FIXME
exclude="NoSupport" //fusions marked with these classifications will be thrown away. Can be a list. FIXME

//mapping and counting the coverage
MAP_COMMAND="bowtie2 -k1 --no-mixed --no-discordant --mm"
over_hang=15 //how many bases require on either side of break to count a read.

/********** Variables that shouldn't need to be changed ***********************/
//name of R scripts
R_filter_transcripts_script=code_base+"/process_transcriptome_blat_table.R"
R_get_final_list=code_base+"/make_final_table.R"
R_get_spanning_reads_script=code_base+"/get_spanning_reads.R"
R_compile_results_script=code_base+"/compile_results.R"

//A function to get the base name of a pair of reads. This is used
//for the directory name and other prefixes in the pipeline
def get_base_name(s1,s2) {
    String base = "";
    //make the directory name the part of the two
    //read end names which are in common
    for(int i=0;i<s1.length() & i<s2.length();i++){
       if(s1.charAt(i)==s2.charAt(i)){  
          base += s1.charAt(i); 
       } else { break; } 
    }  
    //remove the end of the string
    base=base.replaceAll('_R$','').replaceAll('_$','')
    //get just the filename after the path
    def stringList = base.split("/")
    base=stringList[(stringList.size()-1)]
}

/******************* Here are the pipeline stages **********************/

//lets start by checking the dependencies
run_check = {
    exec """
       echo "Using a read length of $read_length" ;
       echo "Checking for the required software..." ;
       for i in $commands ; do which $i 2>/dev/null || { echo "CAN'T FIND $i" ; 
            echo "PLEASE INSTALL IT... STOPPING NOW" ; exit 1  ; } ; done ;
       echo "Now checking for required data files..." ;
       for i in $transFasta $transTable $hgFasta ; 
            do ls $i 2>/dev/null || { echo "CAN'T FIND $i..." ; 
	    echo "PLEASE DOWNLOAD and/or FIX PATH... STOPPING NOW" ; exit 1  ; } ; done ;
       echo "All looking good"
    """
}

//Make a directory for each sample
make_dir = {
   from("*.gz"){
      def base=get_base_name(input1.toString(),input2.toString())
      output.dir=base.prefix
      produce(base+".ignore"){
         exec """ 
            if [ ! -d $output.dir ]; then mkdir $output.dir ; fi ;
      	    touch $output #this is just to get around the dir. being passed.
         """
      }
   }
}

//primarily a read trimming step
prepare_reads = {
	def base=input.split("/")[0]
	output.dir=base //set the output directory
	from("*.gz"){
	produce(base+"_trim1.fastq",base+"_trim2.fastq"){ 
	    // need to check here for whether the files are zipped - FIX
	    //trim & fix the file names so Trinity handles the paired-ends reads correctly 
           exec """
             cd $base ;
             trimmomatic PE -threads $threads -phred$scores ../$input1 ../$input2 
                         tempp1.fq tempu1.fq tempp2.fq tempu2.fq 
                         LEADING:$minQScore TRAILING:$minQScore MINLEN:$minlen ;
              
             function fix_ids { 
                    cat \$1 | awk -v app=\$2 
                               'BEGIN{ i=0 }{ 
                                  if(i==0) printf(\"%s\\/%s\\n\", \$1, app) ; 
                                  else print \$1 ; 
                                  i++ ; 
                                  if(i==4) i=0 }' 2>/dev/null 
             ; } ;
	     fix_ids tempp1.fq 1 > ${base}_trim1.fastq ;
	     fix_ids tempp2.fq 2 > ${base}_trim2.fastq ;
             rm tempu1.fq tempu2.fq tempp1.fq tempp2.fq ;
             cd ../ 
           """ 
	}
   }
}

//Cat read pairs
cat_reads = {
    output.dir=input.split("/")[0]
    exec "cat $input1.fastq $input2.fastq > $output.fastq"
}

remove_dup = { 
    output.dir=input.split("/")[0]
    exec "fastx_collapser -i $input.fastq -o $output.fastq"
}

//Remove any reads which map completely to reference transcripts
get_unmapped = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fasta"){
	exec """
	     $MAP_COMMAND --very-fast --un $output -p $threads -x $transFasta.prefix -f -U $input
                          -S ${output.dir}/temp.sam 2>&1 | tee $base/log_initial_map ;
	     rm -f ${output.dir}/temp.sam
        """
    }
}

//Run the de novo assembler (this is the most important part of the pipeline)
run_assembly = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fasta"){
       exec "my_time ./assemble.sh $base $input1 $input2 $output $threads"
    }
}

//Align the assembled transcripts against reference gene sequences using blat
align_transcripts_to_annotation = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".psl"){
       exec """
                my_time blat $transFasta $input -minIdentity=$minId -minScore=$minScore -tileSize=$tile 
                      -maxIntron=$maxIntron $output 2>&1 | tee $base/log_blat
            """
    }
}

//Parse the blat alignment table and filter for candidate fusions (uses an R script)
filter_transcripts = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".txt"){
       exec """
               time R --vanilla --args $input $output $gapSize $transTable < 
	       $R_filter_transcripts_script &> $base/log_filter 
            """
    }
}

//Extract the fasta sequences for the candidate fusions into their own fasta file
extract_fusion_sequences = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fusions.fa"){
       from("txt","fasta"){
          exec """
           cat $input1 | cut -f 1 | sed \'s/^/^>/g\' > ${output}.temp ;
           fasta_formatter -i $input2 | grep -A1 -f ${output}.temp | grep -v \"\\-\\-\" > $output ;
           rm ${output}.temp ;
          """
       }
    }
}

//Map the reads back to the candidate fusion sequences
map_reads = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".sorted.bam"){
	from("fusions.fa","1.fastq","2.fastq")
	prefix=base+"/"+base
	exec """
           bowtie2-build $input1 $input1.prefix ;
	   $MAP_COMMAND --no-unal -p $threads -x $input1.prefix -1 $input2 -2 $input3 -S ${prefix}.sam 2>&1 | tee $base/log_candidates_map;
           samtools view -S -b ${prefix}.sam > ${prefix}.bam ;
           samtools sort ${prefix}.bam ${prefix}.sorted ;
           samtools index $output ${output}.bai ;
       """
    }
}

//Calculate the number of reads which span the breakpoint of the fusions
get_spanning_reads = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".reads"){
        from("txt","bam"){
	   exec """ 
               samtools view  $input2  | cut -f 3,4,8  > $base/${base}.temp ;
               R --vanilla --args $base/$base $input1 $base/${base}.temp 
                 $output $read_length $over_hang < $R_get_spanning_reads_script 
	   """
	}
    }
}

//Used for the non-assembly pipeline
make_simple_reads_table = {
	def base=input.split("/")[0]
	output.dir=base
	produce(base+".reads"){
	   from("txt"){
           exec """
               echo  -e "transcript\tbreak_min\tbreak_max\tfusion_genes\tspanning_pairs\tspanning_reads" > $output ; 
               awk '{ print \$1"\t"\$2"\t"\$3"\t"\$4"\t"0"\t"1}' $input | sort -u  >> $output
           """
	    }
	}
}

//Align candidate fusions to the genome
align_transcripts_to_genome = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+"_genome.psl"){
	from(base+".fusions.fa"){
            exec "blat $hgFasta $input1 $output 2>&1 | tee $base/log_genome_blat"
	}
    }
}

//Do a bit more filtering and compile the final filtered list (uses an R script)
get_final_list = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".summary"){
	from("psl","reads"){
	    	exec """
	               R --vanilla --args $input1 $input2 $transTable $finalGapSize $exclude $output < $R_get_final_list 
	        """
	 }
    }
}

//Compile the results from multiple samples into an excel .csv table
//Make a fasta file with the candidates
compile_all_results = {
    produce(output_name+".fasta",output_name+".csv"){
       exec """
          R --vanilla --args $output_name $inputs < $R_compile_results_script ;
	  function get_sequence { 
	     if [ \$1 == "sample" ] ; then return ; fi ;
	     fusions_file=\$1/\$1.fusions.fa ;
	     new_id=\$1---\$2---\${13} ;
             echo ">\$new_id" >> ${output_name}.fasta ;
             grep -A1 "^>\${13}" \$fusions_file | grep -v "^>"  >> ${output_name}.fasta
	  ; } ;
	  rm -f ${output_name}.fasta ;
	  cat ${output_name}.csv | tr "," "\\t" | sed 's/\\"//g' | while read line ; do get_sequence \$line ; done ;
	  echo "Done writing ${output_name}.fasta" ;
	  echo "All Done" 
        """
       }
}



