/** How to run... */

/** The following parameters must be set by the user: **/
// environment
code_base="/mnt/storage/nadiad/work_area/20121023_cancer_fusion/jaffa-project/"
trimmomatic="/usr/local/Trimmomatic-0.30/trimmomatic-0.30.jar" //"/home/nadia/Trimmomatic-0.30/trimmomatic-0.30.jar"
TRINITY="/mnt/storage/nadiad/trinityrnaseq_r20131110/Trinity.pl" //trinityrnaseq_r2013_08_14/Trinity.pl
hgFasta="/mnt/storage/shared/genomes/hg19/fasta/hg19.fa" //"/home/Shared/data/annotation/Human/genome/GRCH37/GRCH37.fa"
refSeqFasta=code_base+"hg19_genCode" //"hg19_ens" //genCode" //transcript sequences
refSeqTable=code_base+"hg19_genCode.tab" //+"hg19_ens.tab" //genCode.tab" //"refSeq_table.txt" //table of gene corrdinates
cpus="1"

// Input pattern (See bpipe documentation for how files are grouped and split )
fastq_filename_pattern="%_*.fastq.gz"
//Read length
read_length=75

/** The parameters below don't need to be set by the user, but you might want to change them **/

//Default output name
output_name="jaffa_results"

// trimming
scores=33
minlen=30
minQScore=0  //10 

// for the assembly
mem="100G"
contig_length="100"


// for aligning to known genes using blat
minId="98"
minScore="30" // this is the minimum required flaking sequence assembled around the break-point
tile="18"
maxIntron="0"

// filtering
gapSize="1000" //minimum distance between the two fusion candidates for the 1st filtering stage

//mapping
MAP_COMMAND="bowtie2 -k1 --no-mixed --no-discordant --mm"
over_hang=15

//name of R scripts
R_filter_transcripts_script=code_base+"process_transcriptome_blat_table.R"
R_get_final_list=code_base+"make_final_table.R"
R_get_spanning_reads_script=code_base+"get_spanning_reads.R"
R_compile_results_script=code_base+"compile_results.R"

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

//primarily a read trimming step
prepare_reads = {
	def base=get_base_name(input1.toString(),input2.toString())
	output.dir=base //set the output directory
	produce(base+"_trim1.fastq",base+"_trim2.fastq"){ 
	    // need to check here for whether the files are zipped.
	    //trim & fix the file names so Trinity handles the paired-ends reads correctly 
           exec """
             if [ ! -d $output.dir ]; then mkdir $output.dir ; fi ;
             cd $base ;
             java -jar $trimmomatic PE -threads $cpus -phred$scores $input1 $input2 
                         tempp1.fq tempu1.fq tempp2.fq tempu2.fq 
                         LEADING:$minQScore TRAILING:$minQScore MINLEN:$minlen ;
              
             function fix_ids { 
                    cat \$1 | awk -v app=\$2 
                               'BEGIN{ i=0 }{ 
                                  if(i==0) printf(\"%s\\/%s\\n\", \$1, app) ; 
                                  else print \$1 ; 
                                  i++ ; 
                                  if(i==4) i=0 }' 
             ; } ;
	     fix_ids tempp1.fq 1 > ${base}_trim1.fastq ;
	     fix_ids tempp2.fq 2 > ${base}_trim2.fastq ;
             rm tempu1.fq tempu2.fq tempp1.fq tempp2.fq ;
             cd ../ 
           """ 
	}
}

//Merge overlapping read pairs
/** merge_pairs = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".merged.fastq"){
	//flash MCF-7_ENCODE_20mill_trim1.fastq MCF-7_ENCODE_20mill_trim2.fastq -o ....
	cat .... > $output
    }
}

//Cat read pairs
cat_reads = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fastq"){
	exec "cat $input1 $input2 > $output"
    }
}

//Convert from fastq to fasta
fastq_to_fasta = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fasta"){
        exec "fastq_to_fasta -Q33 -i $input > $output"
    }
}**/

//Map
trans_map = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".unmapped.fasta"){
	exec """
	     $MAP_COMMAND --very-fast --un-conc ${output.prefix.prefix}.fastq -p $cpus -x $refSeqFasta -f -U $input1,$input2
                          -S temp.sam 2>&1 | tee $base/log_initial_map ;

        """
    }            // ; rm temp.sam
}

//Run the de novo assembler (this is the most important part of the pipeline)
run_assembly = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".fasta"){
       exec "my_time ./assemble.sh $base $input1 $input2 $output $cpus"
    }
}

//Align the assembled transcripts against reference gene sequences using blat
align_transcripts_to_annotation = {
    def base=input.split("/")[0]
    output.dir=base
    produce(base+".psl"){
       exec """
                my_time blat ${refSeqFasta}.fa $input -minIdentity=$minId -minScore=$minScore -tileSize=$tile 
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
               time R --slave --no-save --no-restore --no-environ --args $input $output 
                   $gapSize $refSeqTable < $R_filter_transcripts_script 2>&1 | tee $base/log_filter 
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
           bowtie2-build $output $output.prefix ;
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
	   $MAP_COMMAND --no-unal -p $cpus -x $input1.prefix -1 $input2 -2 $input3 -S ${prefix}.sam 2>&1 | tee $base/log_candidates_map;
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
               R --slave --no-save --no-restore --no-environ --args $base/$base $input1 $base/${base}.temp 
                              $output $read_length $over_hang < $R_get_spanning_reads_script 
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
	               R --slave --no-save --no-restore --no-environ --args $input1 $input2 $output $refSeqTable < $R_get_final_list
	        """
	 }
    }
}

//Remove un-needed and large files
clean_up = {
    def base=input.split("/")[0]
    prefix=base+"/"+base
    exec "rm -rf ${prefix}.sam ${prefix}.sam.bam ${prefix}.1.fastq ${prefix}.2.fastq ${prefix}.fusions.*.bt2 ${prefix}.temp"
}


//Compile the results from multiple samples into an excel .csv table
compile_all_results = {
    exec " R --slave --no-save --no-restore --no-environ --args $output_name $inputs < $R_compile_results_script"
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

