/*********************************************************************************
 ** This file defines all the JAFFA pipeline stages (used in JAFFA.groovy, JAFFA_hybrid.groovy and
 ** JAFFA_no_assemble.groovy). See our website for details on running, 
 ** https://code.google.com/p/jaffa-project/.
 **
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 7th March 2014
 ********************************************************************************/
VERSION=0.80

codeBase = new File(bpipe.Config.config.script).parentFile.absolutePath

/**********  Parameters that must be check by the user: ***********************/
/*** Modify them below, or set them when you run bpipe (with the -p option)   */

readLength=50 //Read length

threads=1 //Threads to use when running the pipeline on a single sample. ie. the total threads will be samples*threads
          //Note that oases doesn't support threads, so this option will probably not make much difference to performance

// Genome, Transcriptome and related data paths. You have two options:
// 1) put the full paths below. e.g. hgFasta=<path_to_genome> 
// or 2) leave as is and symlink the data files to the jaffa code directory. e.g. ln -s <path_to_genome> <path_to_jaffa_code_directory>
hgFasta=codeBase+"/hg19.fa"  //genome sequence

// Input pattern (see bpipe documentation for how files are grouped and split )
// group on start, split on end. eg. on ReadsA_1.fastq.gz, ReadA_2.fastq.gz
// this would group the read files into pairs like we want.
fastqInputFormat="%_*.fastq.gz"

// Simlar to above for running JAFFA_no_assemble with a fasta file.
// You should not need to change this unless the suffix is ".fa" instead of ".fasta"
fastaSuffix="fasta"  
fastaInputFormat="%."+fastaSuffix

/***************** Other configurables *************************************/

//Default output name
outputName="jaffa_results"

// trimming
scores=33 //PHRED quality score type
minlen=30 //reads shorted that this after trimmed are thrown out.
minQScore=0 //heads and tails of reads will be trimmed if their quality score falls below this.
//note: by default the 0 above means that no trimming is done (we found this gave the
//best assembly results)

// assembly options (we founds these setting to work well on 50bp reads)
Ks="19,23,27,31,35" //kmer lengths to use for the assembly
Kmerge=27 //what kmer should Oases use to merge the assemblies.
transLength=100 //the minimum length for Oases to report an assembled contig

// for aligning to known genes using blat
minIdTrans=98 //98% similar when we blat to the human transcriptome
minScore=30 //this is the minimum matches to report an alignment - ie required flanking sequence around a breakpoint
contigTile=18 //big tile size makes blat fast
readTile=0 //This is a dummy. Gets set dynamically later. //reduce this for reads shorted than 100 bases. e.g. 15 for 75bp. 
maxIntron=0 //don't expect intron when mapping to the transcriptome

// filtering
gapSize=1000 //minimum distance between the two fusion candidates for the 1st filtering stage
finalGapSize=10000 //minimum distance for the final filtering
exclude="NoSupport" //fusions marked with these classifications will be thrown away. Can be a list. 

//mapping and counting the coverage
MapCommand="bowtie2 -k1 --no-mixed --no-discordant --mm"
overHang=15 //how many bases are require on either side of a break to count the read.

/********** Variables that shouldn't need to be changed ***********************/

//location of transcriptomic data
transFasta=codeBase+"/hg19_genCode.fa"  // transcript cDNA sequences
transTable=codeBase+"/hg19_genCode.tab" // table of gene coordinates

//name of scripts
R_filter_transcripts_script=codeBase+"/process_transcriptome_blat_table.R"
R_get_final_list=codeBase+"/make_final_table.R"
R_get_spanning_reads_script=codeBase+"/get_spanning_reads.R"
R_compile_results_script=codeBase+"/compile_results.R"
oases_assembly_script=codeBase+"/assemble.sh"

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
    produce("checks"){
    exec """
       echo "Running JAFFA version $VERSION" ;
       echo "Using a read length of $readLength" ;
       echo "Checking for the required software..." ;
       for i in $commands ; do which $i 2>/dev/null || { echo "CAN'T FIND $i" ; 
            echo "PLEASE INSTALL IT... STOPPING NOW" ; exit 1  ; } ; done ;
       echo "Now checking for required data files..." ;
       for i in $transFasta $transTable $hgFasta ; 
            do ls $i 2>/dev/null || { echo "CAN'T FIND $i..." ; 
	    echo "PLEASE DOWNLOAD and/or FIX PATH... STOPPING NOW" ; exit 1  ; } ; done ;
       echo "All looking good" ;
       echo "running JAFFA version $VERSION.. checks passed" > $output 
    """
    }
}

//Make a directory for each sample
make_dir_using_fastq_names = {
   from("*.gz"){
      def base=get_base_name(input1.toString(),input2.toString())
      output.dir=base
      produce(base+".ignore"){
         exec """ 
            if [ ! -d $output.dir ]; then mkdir $output.dir ; fi ;
      	    touch $output #this is just to get around the dir. being passed.
         """  
       }
   }
}

//Like above, but we need to do something different if there is only 1 input file.
make_dir_using_fasta_name = {
   from("*.fasta"){
      def inputPath=file(input.toString()).absolutePath
      def stringList = input.prefix.split("/")
      base=stringList[(stringList.size()-1)]
      output.dir=base
      produce(base+".fasta"){
         exec """
            if [ ! -d $output.dir ]; then mkdir $output.dir ; fi ;
            ln -s $inputPath $output
         """    
      }
   }
}

//Primarily a read trimming step - currently actually just unzips the reads.
prepare_reads = {
	def base=input.split("/")[0]
	output.dir=base //set the output directory
	from("*.gz"){
	produce(base+"_trim1.fastq",base+"_trim2.fastq"){ 
	    // need to check here for whether the files are zipped - FIX
	    //trim & fix the file names so Trinity handles the paired-ends reads correctly 
           exec """
             trimmomatic PE -threads $threads -phred$scores $input1 $input2 
                         ${base}/tempp1.fq ${base}/tempu1.fq ${base}/tempp2.fq ${base}/tempu2.fq 
                         LEADING:$minQScore TRAILING:$minQScore MINLEN:$minlen ;
              
             function fix_ids { 
                    cat \$1 | awk -v app=\$2 
                               'BEGIN{ i=0 }{ 
                                  if(i==0) printf(\"%s\\/%s\\n\", \$1, app) ; 
                                  else print \$1 ; 
                                  i++ ; 
                                  if(i==4) i=0 }' 2>/dev/null 
             ; } ;
	     fix_ids ${base}/tempp1.fq 1 > ${base}/${base}_trim1.fastq ;
	     fix_ids ${base}/tempp2.fq 2 > ${base}/${base}_trim2.fastq ;
             rm ${base}/tempu1.fq ${base}/tempu2.fq ${base}/tempp1.fq ${base}/tempp2.fq ;
           """ 
	}
   }
}

//Cat read pairs into a single file
cat_reads = {
    output.dir=input.split("/")[0]
    exec "cat $input1.fastq $input2.fastq > $output.fastq"
}

//Remove duplicated reads
remove_dup = { 
    output.dir=input.split("/")[0]
    exec "fastx_collapser -Q${scores} -i $input.fastq -o $output.fasta"
}

//Remove any reads which map completely to a reference transcript
get_unmapped = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fasta"){
	exec """
	     $MapCommand --very-fast --un $output -p $threads -x $transFasta.prefix -f -U $input
                          -S ${output.dir}/temp.sam 2>&1 | tee $base/log_initial_map_to_reference ;
	     rm -f ${output.dir}/temp.sam
        """
    }
}

//Like above: remove reads that don't map to the transcriptome, but this time use the assembled
//transcriptome as well as the reference
get_assembly_unmapped = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".unmapped.fasta"){
        exec """
             $MapCommand --very-fast --un ${base}/${base}.temp -p $threads -x $transFasta.prefix -f -U $input
                          -S ${output.dir}/temp.sam 2>&1 | tee $base/log_initial_map_to_reference ;
	     bowtie2-build ${base}/${base}.fasta ${base}/${base} ;
             $MapCommand --very-fast --un $output -p $threads -x ${base}/${base} -f -U ${base}/${base}.temp
                          -S ${output.dir}/temp.sam 2>&1 | tee $base/log_initial_map_to_assembly ;
             rm -f ${output.dir}/temp.sam ${base}/${base}.temp;
        """
    }
}

//Run the de novo assembly
run_assembly = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fasta"){
       exec "/usr/bin/time -v $oases_assembly_script $base $input1 $input2 $output $Ks $Kmerge $transLength"
    }
}

//Align the assembled transcripts against reference gene sequences using blat
align_transcripts_to_annotation = {
    def base=input.split("/")[0]
    output.dir=base
    produce(input.prefix+".psl"){ from(input.prefix+".fasta"){
       exec """
                time blat $transFasta $input -minIdentity=$minIdTrans -minScore=$minScore -tileSize=$tile 
                      -maxIntron=$maxIntron $output 2>&1 | tee $base/log_blat
            """
    } }
}

//Parse the blat alignment table and filter for candidate fusions (uses an R script)
filter_transcripts = {
    def base=input.split("/")[0]
    output.dir=base
    produce(input.prefix+".txt"){ from(input.prefix+".psl"){
       exec """
               time R --vanilla --args $input $output $gapSize $transTable < 
	       $R_filter_transcripts_script &> $base/log_filter 
            """
    }}
}

//Extract the fasta sequences for the candidate fusions into their own fasta file
extract_fusion_sequences = {
    def base=input.split("/")[0]
    output.dir=base
    produce(input.prefix+".fusions.fa"){ from(input.prefix+".psl"){
       from("txt","fasta"){
          exec """
           cat $input1 | cut -f 1 | sed \'s/^/^>/g\' > ${output}.temp ;
           fasta_formatter -i $input2 | grep -A1 -f ${output}.temp | grep -v \"\\-\\-\" > $output ;
           rm ${output}.temp ;
          """
       }
    }}
}

//Map the reads back to the candidate fusion sequences
map_reads = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".sorted.bam"){
	from("fusions.fa","1.fastq","2.fastq")
	def prefix=base+"/"+base
	exec """
           bowtie2-build $input1 $input1.prefix ;
	   $MapCommand --no-unal -p $threads -x $input1.prefix -1 $input2 -2 $input3 -S ${prefix}.sam 2>&1 | tee $base/log_candidates_map;
           samtools view -S -b ${prefix}.sam > ${prefix}.bam ;
           samtools sort ${prefix}.bam ${prefix}.sorted ;
           samtools index $output ${output}.bai ;  
       """
    }
}  // $MapCommand --no-unal -p $threads -x $input1.prefix -1 $input2 -2 $input3 -S | samtools view -bSu - | samtools sort -o - - > $output ;
 

//Calculate the number of reads which span the breakpoint of the fusions
get_spanning_reads = {
    def base=input.split("/")[0]
    output.dir=base
    produce(input.txt.prefix+".reads"){
       from("txt","bam"){
	   exec """ 
               samtools view  $input2  | cut -f 3,4,8  > $base/${base}.temp ;
               R --vanilla --args $base/$base $input1 $base/${base}.temp 
                 $output $readLength $overHang < $R_get_spanning_reads_script ;
	       rm $base/${base}.temp 
	   """
	}
    }
}

//Used for the direct and hybrid pipelines - In this case the spanning reads will be 1 for each
//read and the spanning pairs will be 0. 
make_simple_reads_table = {
	def base=input.split("/")[0]
	output.dir=base
	produce(input.txt.prefix+".reads"){
	from("txt"){
           exec """
               echo  -e "transcript\tbreak_min\tbreak_max\tfusion_genes\tspanning_pairs\tspanning_reads" > $output ; 
               awk '{ print \$1"\t"\$2"\t"\$3"\t"\$4"\t"0"\t"1}' $input | sort -u  >> $output
           """
	    }
	}
}

//This stage is only used the by hybrid mode.
//It concatenates the fusions sequence files, then the read files.
merge_assembly_and_unmapped_reads_candidates = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".all.fusions.fa",base+".all.reads"){
       from("fusions.fa",base+".fusions.fa",
            "reads",base+".reads"){
           exec """ cat $input1 $input2 > $output1 ;
                    cp $input3 $output2 ; tail -n+2 $input4 >> $output2 """
       }
    }
}


//Align candidate fusions to the genome
align_transcripts_to_genome = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+"_genome.psl"){
	from(".fusions.fa"){
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
    var type : ""
    produce(outputName+".fasta",outputName+".csv",outputName+".psl"){
       exec """
          R --vanilla --args $outputName $inputs < $R_compile_results_script ;
	  function get_sequence { 
	     if [ \$1 == "sample" ] ; then return ; fi ;
	     fusions_file=\$1/\$1${type}.fusions.fa ;
	     new_id=\$1---\$2---\${13} ;
             echo ">\$new_id" >> ${outputName}.fasta ;
	     break=\${14} ; 
	     sequence=`grep -A1 "^>\${13}" \$fusions_file | grep -v "^>"` ;
	     start=`echo \$sequence | cut -c 1-\$((\${break}-1))` ;
	     middle=`echo \$sequence | cut -c \$break-\$((\${break}+1)) | tr '[:upper:]' '[:lower:]'` ;
	     string_length=`echo \${#sequence}` ;
	     end=`echo \$sequence | cut -c \$((\$break+2))-$string_length ` ;
	     echo ${start}${middle}${end} >> ${outputName}.fasta ;
	     grep \${13} \$1/\$1_genome.psl >> ${outputName}.psl ;
	  } ;
	  rm -f ${outputName}.fasta ;	  
	  cat ${outputName}.csv | tr "," "\\t" | sed 's/\\"//g' | while read line ; do get_sequence \$line ; done ;
	  echo "Done writing ${outputName}.fasta" ;
	  echo "All Done" 
        """
       }
}
