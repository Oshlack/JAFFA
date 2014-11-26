/*********************************************************************************
 ** This file defines all the JAFFA pipeline stages (used in JAFFA_assembly.groovy, JAFFA_hybrid.groovy and
 ** JAFFA_direct.groovy). See our website for details on running, 
 ** https://code.google.com/p/jaffa-project/.
 **
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>
 ** Last Update: 24th November 2014
 ********************************************************************************/
VERSION=1.04

codeBase = file(bpipe.Config.config.script).parentFile.absolutePath

load codeBase+"/tools.groovy"

/**********  Parameters that must be check by the user: ***********************/
/*** Modify them below, or set them when you run bpipe (with the -p option)   */

readLayout="paired" //change to "single" or single-end reads

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

// Simlar to above for running JAFFA_direct with a fasta file.
// You should not need to change this unless the suffix is ".fa" instead of ".fasta"
fastaSuffix="fasta"  
fastaInputFormat="%."+fastaSuffix

/***************** Other configurables *************************************/

//Default output name
outputName="jaffa_results"

// trimming
scores=33 //PHRED quality score type
minlen=35 //reads shorted that this after trimmed are thrown out.
minQScore=0 //heads and tails of reads will be trimmed if their quality score falls below this.
//note: by default the 0 above means that no trimming is done (we found this gave the
//best assembly results)

// assembly options (we founds these setting to work well on 50bp reads)
Ks="19,36,4" //kmer lengths to use for the assembly: 19,23,27,31,15
Kmerge=27 //what kmer should Oases use to merge the assemblies.
transLength=100 //the minimum length for Oases to report an assembled contig

// for aligning to known genes using blat
minIdTrans=98 //98% similar when we blat to the human transcriptome
minScore=30 //this is the minimum matches to report an alignment - ie required flanking sequence around a breakpoint
contigTile=18 //big tile size makes blat fast
readTile=0 //This is a dummy. It gets set dynamically later. //reduce this for reads shorted than 100 bases. e.g. 15 for 75bp. 
maxIntron=0 //don't expect intron when mapping to the transcriptome

// filtering
gapSize=1000 //minimum distance between the two fusion candidates for the 1st filtering stage
finalGapSize=10000 //minimum distance for the final filtering
exclude="NoSupport,PotentialRegularTranscript" //fusions marked with these classifications will be 
					       //thrown away. Can be a comma seperated list. 

//mapping and counting the coverage
mapParams="-k1 --no-mixed --no-discordant --mm"
overHang=15 //how many bases are require on either side of a break to count the read.

/********** Variables that shouldn't need to be changed ***********************/

//location of the genome with genes masked out - used to filter the reads
maskedGenome=codeBase+"/Masked_hg19"

//location of transcriptomic data
transFasta=codeBase+"/hg19_genCode19.fa"  // transcript cDNA sequences
transTable=codeBase+"/hg19_genCode19.tab" // table of gene coordinates

//known fusions database
knownTable=codeBase+"/known_fusions.txt" //a two column table of know/recurrent fusions

//name of scripts
R_filter_transcripts_script=codeBase+"/process_transcriptome_blat_table.R"
R_get_final_list=codeBase+"/make_final_table.R"
R_get_spanning_reads_script=codeBase+"/get_spanning_reads.R"
R_get_spanning_reads_direct_script1=codeBase+"/get_spanning_reads_for_direct_1.R"
R_get_spanning_reads_direct_script2=codeBase+"/get_spanning_reads_for_direct_2.R"
R_compile_results_script=codeBase+"/compile_results.R"
oases_assembly_script=codeBase+"/assemble.sh"

//A function to get the base name of a pair of reads. This is used
//for the directory name and other prefixes in the pipeline
def get_base_name_paired(s1,s2) {
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

//A function to get the base name of single-end reads
def get_base_name_single(s1){
     def stringList = s1.split("/")
     base=stringList[(stringList.size()-1)]
}


/******************* Here are the pipeline stages **********************/

//lets start by checking the dependencies
run_check = {
    produce("checks"){
    exec """
       echo "Running JAFFA version $VERSION" ;
       echo "Checking for required data files..." ;
       for i in $transFasta $transTable $knownTable $hgFasta ${maskedGenome}.1.bt2 ; 
            do ls $i 2>/dev/null || { echo "CAN'T FIND ${i}..." ; 
	    echo "PLEASE DOWNLOAD and/or FIX PATH... STOPPING NOW" ; exit 1  ; } ; done ;
       echo "All looking good" ;
       echo "running JAFFA version $VERSION.. checks passed" > $output 
    """
    }
}

//Make a directory for each sample
make_dir_using_fastq_names = {
   from("*.gz"){
      def base=""
      if(readLayout=="single")
         base=get_base_name_single(input.prefix.prefix.toString())
      else 
         base=get_base_name_paired(input1.toString(),input2.toString())
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
      def base = get_base_name_single(input.prefix)
      output.dir=base
      produce(base+".fasta"){
         exec """
            if [ ! -d $output.dir ]; then mkdir $output.dir ; fi ;
            ln -s $inputPath $output
         """    
      }
   }
}

//Read trimming, ID fixing and filtering out reads 
//that map to chrM, introns and intergenetic regions

prepare_reads = {
	def base=input.split("/")[0]
	output.dir=base //set the output directory
	from("*.gz"){
          if(readLayout=="single"){
	     produce(base+"_filtered_reads.fastq.gz",
		     base+"_leftover_reads.fastq.gz"){
             exec """
                 $trimmomatic SE -threads $threads -phred$scores $input.gz $base/${base}_trim.fastq
                         LEADING:$minQScore TRAILING:$minQScore MINLEN:$minlen ;
              $bowtie2 $mapParams --very-fast
               	        --al-gz $output1 --un $base/temp_trans_unmap_reads.fastq
               		-p $threads -x $transFasta.prefix -U $base/${base}_trim.fastq
			-S /dev/null ;
              $bowtie2 $mapParams --very-fast --un-gz $output2 -p $threads -x $maskedGenome
               	        -U $base/temp_trans_unmap_reads.fastq -S /dev/null ;
              cat $output2 >> $output1 ;
              rm $base/temp_trans_unmap_reads.fastq ${base}/${base}_trim.fastq 
                  """
              }
	  } else {
            produce(base+"_filtered_reads.fastq.1.gz",
		    base+"_filtered_reads.fastq.2.gz",
		    base+"_leftover_reads.fastq.1.gz",
		    base+"_leftover_reads.fastq.2.gz"){
	     // need to check here for whether the files are zipped - FIX
	     //trim & fix the file names so Trinity handles the paired-ends reads correctly 
             exec """
             $trimmomatic PE -threads $threads -phred$scores $input1 $input2 
                         ${base}/tempp1.fq /dev/null ${base}/tempp2.fq /dev/null 
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
              rm ${base}/tempp1.fq ${base}/tempp2.fq ;
              $bowtie2 $mapParams --very-fast 
	         --al-conc-gz ${output1.prefix.prefix}.gz --un-conc $base/temp_trans_unmap_reads.fastq 
	       	 -p $threads -x $transFasta.prefix 
	       	 -1 ${base}/${base}_trim1.fastq -2 ${base}/${base}_trim2.fastq -S /dev/null ;
	      $bowtie2 $mapParams --very-fast --un-conc-gz ${output3.prefix.prefix}.gz -p $threads -x $maskedGenome 
               	 -1 $base/temp_trans_unmap_reads.1.fastq -2 $base/temp_trans_unmap_reads.2.fastq -S /dev/null ;
	      cat $output3 >> $output1 ;
	      cat $output4 >> $output2 ;
	      rm $base/temp_trans_unmap_reads.1.fastq $base/temp_trans_unmap_reads.2.fastq
	         ${base}/${base}_trim1.fastq ${base}/${base}_trim2.fastq ;
           """ 
	  }
       }
   }
}

//Cat read pairs into a single file
cat_reads = {
    if(readLayout=="single") return
    output.dir=input.split("/")[0]
    exec "cat $input1.fastq $input2.fastq > $output.fastq"
}


//Get read which either align discordantly or not at all
get_unmapped = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fasta",base+"_discordant_pairs.bam"){
    from("*_leftover_reads*.gz"){
    if(readLayout=="single"){
        exec """
	   $bowtie2 -k1 -p $threads --un $base/unmapped.fastq -x $transFasta.prefix -U $input |
       	   $samtools view -F 4 -S -b - | $samtools sort - $output2.prefix ;
           $samtools index $output2 ;
       """
     }else{
        exec """
	   $bowtie2 -k1 -p $threads --un $base/unmapped.fastq -x $transFasta.prefix -U $input1,$input2 |
       	   $samtools view -F 4 -S -b - | $samtools sort - $output2.prefix ;
           $samtools index $output2 ;
	"""
     }
     exec """
	$reformat in=$base/unmapped.fastq out=$base/temp.fasta threads=$threads ;
	$dedupe in=$base/temp.fasta out=$output1 threads=$threads absorbcontainment=f ; 
	rm $base/temp.fasta $base/unmapped.fastq 2> /dev/null
       """
    }
  }
}

//Like above: remove reads that don't map to the transcriptome, but this time use the assembled
//transcriptome as well as the reference
get_assembly_unmapped = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+"-unmapped.fasta",base+"_discordant_pairs.bam"){
      from("*_leftover_reads*.gz"){
        if(readLayout=="single"){
        exec """
           $bowtie2 -k1 -p $threads --un $base/unmapped_ref.fastq -x $transFasta.prefix -U $input |
           $samtools view -F 4 -S -b - | $samtools sort - $output2.prefix ;
           $samtools index $output2 ;
       """
    }else{
        exec """
           $bowtie2 -k1 -p $threads --un $base/unmapped_ref.fastq -x $transFasta.prefix -U $input1,$input2 |
           $samtools view -F 4 -S -b - | $samtools sort - $output2.prefix ;
           $samtools index $output2 ;
       """ 
    }
    exec """
           ${bowtie2}-build ${base}/${base}.fasta ${base}/${base} ;
           $bowtie2 -k1 -p $threads --un $base/unmapped_assembly.fastq -x ${base}/${base} -U ${base}/unmapped_ref.fastq
                  -S /dev/null 2>&1 | tee $base/log_initial_map_to_assembly ;

        $reformat in=$base/unmapped_assembly.fastq out=$base/temp.fasta threads=$threads ;
        $dedupe in=$base/temp.fasta out=$output1 threads=$threads absorbcontainment=f ;
	rm $base/temp.fasta $base/unmapped_assembly.fastq $base/unmapped_ref.fastq
        """
    }}
}

//Run the de novo assembly
run_assembly = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fasta"){
     from("*_filtered_reads.fastq*gz"){
       exec "time $oases_assembly_script $velveth $velvetg $oases \
                  $base $output $Ks $Kmerge $transLength $threads $inputs"
    }
    }
}

//Align the assembled transcripts against reference gene sequences using blat
align_transcripts_to_annotation = {
    def base=input.split("/")[0]
    output.dir=base
    produce(input.prefix+".psl"){
    from(".fasta"){
       exec """
          function run_blat {
                time $blat $transFasta $inputs -minIdentity=$minIdTrans -minScore=$minScore -tileSize=\$1
                      -maxIntron=$maxIntron $output 2>&1 | tee $base/log_blat ;
	  } ;
	  run_blat $contigTile;
	  `### test for the Blat tileSize bug (version 35) ###` ;
	  if [[ `cat $base/log_blat` == *"Internal error genoFind.c"* ]] ; then 
	     echo "Blat error with tileSize=$contigTile" ;
	     echo "Let's try again with tileSize=15" ;
 	     run_blat 15;
	  fi ;
          """
    } } 
}


//Align the reads to the annotation 
align_reads_to_annotation = {
    def base=input.split("/")[0]
    output.dir=base
    produce(input.prefix+".psl"){ 
    from(input.prefix+".fasta"){
       exec """
	  `#### find out the minimum read length so we can set blat's tileSize accordingly` ;
	  `#### we just use the first few thousand from the fasta file as a sample` ;
	  minReadLength=`head -n 100000 $input | awk ' NR % 2 == 0 '|  awk '{print length}' | sort -n | head -1` ;
	  if [ $readTile -eq "0" ] ; then  
               if [ \$minReadLength -le "100" ] ; then 
	           readTile=15 ; 
               else readTile=18 ; 
               fi ;
	  else readTile=$readTile ;
          fi ;
	  echo "Using tileSize of \$readTile" ;
          function run_blat {
                time $blat $transFasta $input -minIdentity=$minIdTrans -minScore=$minScore -tileSize=\$1
                      -maxIntron=$maxIntron $output 2>&1 | tee $base/log_blat ;
	  } ;
	  run_blat \$readTile;
	  `### test for the Blat tileSize bug (version 35) ###` ;
	  if [[ `cat $base/log_blat` == *"Internal error genoFind.c"* ]] ; then 
	     echo "Blat error with tileSize=\$readTile" ;
	     echo "Let's try again with tileSize=15" ;
 	     run_blat 15;
	  fi ;
          """
    } } 
}

//Parse the blat alignment table and filter for candidate fusions (uses an R script)
filter_transcripts = {
    def base=input.split("/")[0]
    output.dir=base
    produce(input.prefix+".txt"){ from(input.prefix+".psl"){
       exec """
               time $R --vanilla --args $input $output $gapSize $transTable < 
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
           cat $input1 | awk '{print \$1}' | sed \'s/^/>/g\' > ${output}.temp ;
	   $reformat in=$input2 out=stdout.fasta fastawrap=0 | awk '{print \$1}' |
	     grep -Fx -A1 -f ${output}.temp | grep -v \"^\\-\\-\" > $output ;
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
	def prefix=base+"/"+base
	def input_string=""
	from("fusions.fa","*_filtered_reads*gz"){
          if(readLayout=="single")
	     input_string="-U $input2"
	  else
	     input_string="-1 $input2 -2 $input3"
	  exec """
             ${bowtie2}-build $input1 $input1.prefix ;
	     $bowtie2 $mapParams --no-unal -p $threads -x $input1.prefix $input_string | 
             $samtools view -S -b - | $samtools sort - ${prefix}.sorted ;
             $samtools index $output
          """
          }
    }
} 
 

//Calculate the number of reads which span the breakpoint of the fusions
get_spanning_reads = {
    def base=input.split("/")[0]
    output.dir=base
    produce(input.txt.prefix+".reads"){
       from("txt","bam"){
	   exec """ 
               $samtools view $input2 | cut -f 3,4,8  > $base/${base}.temp ;
	       $samtools view $input2 | cut -f 10 | awk '{print length}' > $base/${base}.readLengths ;
               $R --vanilla --args $base/$base $input1 $base/${base}.temp 
                 $output $base/${base}.readLengths $overHang < $R_get_spanning_reads_script ;
	       rm $base/${base}.temp $base/${base}.readLengths
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
	from("txt","*_discordant_pairs.bam"){
           exec """
	       $samtools view -H $input2 | grep "@SQ" | cut -f2 | sed 's/SN://g' > $base/temp_gene_ids ;
	       R --no-save --args $input1 $transTable $base/temp_gene_ids $base/paired_contigs.temp 
	       	 	      < $R_get_spanning_reads_direct_script1 ;
	       function get_spanning_pairs {
	               gene=`echo \$1 | cut -f1 -d"?"` ;
    		       g1=`echo \$1 | cut -f2 -d"?"` ;
    		       g2=`echo \$1 | cut -f3 -d "?"` ;
    		       $samtools view $input2 \$g1 | cut -d "/" -f1 | sort -u > $base/g1 ;
    		       $samtools view $input2 \$g2 | cut -d "/" -f1 | sort -u > $base/g2 ;
    		       left=`cat $base/g1 | wc -l` ;
    		       right=`cat $base/g2 | wc -l` ;
    		       both=`cat $base/g1 $base/g2 | sort -u | wc -l` ;
    		       echo -e "\$gene\t\$(( \$left + \$right - \$both ))" ;
               } ;
	       cat $base/paired_contigs.temp | while read line ; do 
	       	   get_spanning_pairs "\$line" >> $base/spanning_pair_counts.temp ; 
	       done ;
	       R --no-save --args $input1 $base/spanning_pair_counts.temp $output < $R_get_spanning_reads_direct_script2 ;
	       rm $base/temp_gene_ids $base/spanning_pair_counts.temp $base/paired_contigs.temp $base/g1 $base/g2
           """
	    }
	}
} 

make_fasta_reads_table = {
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
            exec "$blat $hgFasta $input1 -minScore=$minScore $output 2>&1 | tee $base/log_genome_blat"
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
	               $R --vanilla --args $input1 $input2 $transTable $knownTable $finalGapSize $exclude $output < $R_get_final_list 
	        """
	 }
    }
}

//Compile the results from multiple samples into an excel .csv table
//Make a fasta file with the candidates
compile_all_results = {
    var type : ""
    produce(outputName+".fasta",outputName+".csv"){
       exec """
          $R --vanilla --args $outputName $inputs < $R_compile_results_script ;
	  function get_sequence { 
	     if [ \$1 == "sample" ] ; then return ; fi ;
	     fusions_file=\$1/\$1${type}.fusions.fa ;
	     new_id=\$1---\$2---\$3 ;
             echo ">\$new_id" >> ${outputName}.fasta ;
	     break=\$4 ; 
	     sequence=`grep -A1 "^>\$3" \$fusions_file | grep -v "^>"` ;
	     start=`echo \$sequence | cut -c 1-\$((\${break}-1))` ;
	     middle=`echo \$sequence | cut -c \$break-\$((\${break}+1)) | tr '[:upper:]' '[:lower:]'` ;
	     string_length=`echo \${#sequence}` ;
	     end=`echo \$sequence | cut -c \$((\$break+2))-$string_length ` ;
	     echo ${start}${middle}${end} >> ${outputName}.fasta ;
	     `# grep \$3 \$1/\$1_genome.psl >> ${outputName}.psl ;` ;
	  } ;
	  rm -f ${outputName}.fasta ;
	  cat ${outputName}.temp | while read line ; do get_sequence \$line ; done ;
	  rm ${outputName}.temp ;
	  echo "Done writing ${outputName}.fasta" ;
	  echo "All Done" 
        """
       }
}
