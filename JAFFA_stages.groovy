/*********************************************************************************
 ** This file defines all the JAFFA pipeline stages (used in JAFFA_assembly.groovy, JAFFA_hybrid.groovy and
 ** JAFFA_direct.groovy). See our website for details on running, 
 ** https://code.google.com/p/jaffa-project/.
 **
 ** Author: Nadia Davidson <nadia.davidson@mcri.edu.au>, Rebecca Evans <rebecca.evans@mcri.edu.au>
 ** Last Update: 8th May 2017
 ********************************************************************************/
VERSION="1.09"

codeBase = file(bpipe.Config.config.script).parentFile.absolutePath
load codeBase+"/tools.groovy"


/******** Path to reference files ********/
// leave if references are in the Jaffa folder
refBase = codeBase

// Path to reference files that are elsewhere on the file system
// refBase = "/path/to/reference/directory"

// Path to reference files 
// env = System.getenv()
// refBase = env['GENOMES']

// Should there be a folder structure for the reference files within the reference folder then
// the below variables will allow for this flexibility
fastaBase = refBase
maskedBase = refBase
transBase = refBase


/**********  Parameters that are likely to change between runs of JAFFA **************/
/*** These are usually set with the -p option in bpipe, but may also be set here    **/ 

readLayout="paired" //change to "single" or single-end reads

// Genome, Transcriptome and related data paths. 
genome="hg38"
annotation="genCode22"

// You have two options:
// 1) put the full file name (including path) below. e.g. genomeFasta=<path_to_genome>/<genome_file_name>
// or 2) leave as is and place or symlink the data files to the jaffa code directory. 
// e.g. ln -s <path_to_genome> <path_to_jaffa_code_directory>
genomeFasta=fastaBase+"/"+genome+".fa"  //genome sequence

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
jaffa_output="" // used to specify an output directory for jaffa eg. set to "jaffa/" for future use by other pipelines

// trimming
scores=33 //PHRED quality score type
minlen=35 //reads shorted that this after trimmed are thrown out.
minQScore=0 //heads and tails of reads will be trimmed if their quality score falls below this.
//note: by default the 0 above means that no trimming is done (we found this gave the
//best assembly results)

// assembly options (we founds these setting to work well on 50bp reads)
Ks="19,36,4" //kmer lengths to use for the assembly: 19,23,27,31,35
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
maskedGenome=maskedBase+"/Masked_"+genome

//location of transcriptomic data
transFasta=transBase+"/"+genome+"_"+annotation+".fa"  // transcript cDNA sequences
transTable=transBase+"/"+genome+"_"+annotation+".tab" // table of gene coordinates

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

/******************* Here are the pipeline stages **********************/

//lets start by checking the dependencies
run_check = {
    doc "check for Jaffa dependencies"
    if (jaffa_output) {
        output.dir=jaffa_output
    }
    produce("checks") {
        exec """
            echo "Running JAFFA version $VERSION" ;
            echo "Checking for required data files..." ;
            for i in $transFasta $transTable $knownTable $genomeFasta ${maskedGenome}.1.bt2 ${transFasta.prefix}.1.bt2 ; 
                 do ls $i 2>/dev/null || { echo "CAN'T FIND ${i}..." ; 
            echo "PLEASE DOWNLOAD and/or FIX PATH... STOPPING NOW" ; exit 1  ; } ; done ;
            echo "All looking good" ;
            echo "running JAFFA version $VERSION.. checks passed" > $output 
        ""","checks"
    }
}

//Read trimming, ID fixing and filtering out reads 
//that map to chrM, introns and intergenetic regions
prepare_reads = {
    doc "Prepare reads"
    output.dir=jaffa_output+branch
    if (inputs.size() == 1) {  // single reads
        produce(branch+"_filtered_reads.fastq.gz",
                branch+"_leftover_reads.fastq.gz") {
            exec """
                $trimmomatic SE -threads $threads -phred$scores $input.gz
                    ${output.dir}/${branch}_trim.fastq
                    LEADING:$minQScore TRAILING:$minQScore MINLEN:$minlen ;
                $bowtie2 $mapParams --very-fast
                    --al-gz $output1
                    --un ${output.dir}/temp_trans_unmap_reads.fastq
                    -p $threads -x $transFasta.prefix
                    -U ${output.dir}/${branch}_trim.fastq
                    -S /dev/null ; 
                $bowtie2 $mapParams --very-fast
                    --un-gz $output2 -p $threads -x $maskedGenome
                    -U ${output.dir}/temp_trans_unmap_reads.fastq -S /dev/null ;
                cat $output2 >> $output1 ;
                rm ${output.dir}/temp_trans_unmap_reads.fastq ${output.dir}/${branch}_trim.fastq
            ""","prepare_reads"
        }
    } else if (inputs.size() == 2) {  // paired reads
        produce(branch+"_filtered_reads.fastq.1.gz",
                branch+"_filtered_reads.fastq.2.gz",
                branch+"_leftover_reads.fastq.1.gz",
                branch+"_leftover_reads.fastq.2.gz") {
                // need to check here for whether the files are zipped - FIX
                //trim & fix the file names so Trinity handles the paired-ends reads correctly
            exec """
                $trimmomatic PE -threads $threads -phred$scores $input1 $input2
                    ${output.dir}/tempp1.fq /dev/null
                    ${output.dir}/tempp2.fq /dev/null
                    LEADING:$minQScore TRAILING:$minQScore MINLEN:$minlen;
                function fix_ids {
                    cat \$1 |
                    awk -v app=\$2
                        'BEGIN{ i=0 }{
                        if(i==0) print \$1 \"/\" app ;
                        else print \$1 ;
                        i++ ;
                        if(i==4) i=0 }'
                    2>/dev/null
                ; } ;
                fix_ids ${output.dir}/tempp1.fq 1 > ${output.dir}/${branch}_trim1.fastq ;
                fix_ids ${output.dir}/tempp2.fq 2 > ${output.dir}/${branch}_trim2.fastq ;
                rm ${output.dir}/tempp1.fq ${output.dir}/tempp2.fq ;

                $bowtie2 $mapParams --very-fast
                    --al-conc-gz ${output1.prefix.prefix}.gz
                    --un-conc ${output.dir}/temp_trans_unmap_reads.fastq
                    -p $threads -x $transFasta.prefix
                    -1 ${output.dir}/${branch}_trim1.fastq
                    -2 ${output.dir}/${branch}_trim2.fastq
                    -S /dev/null ;
                $bowtie2 $mapParams --very-fast
                    --un-conc-gz ${output3.prefix.prefix}.gz
                    -p $threads -x $maskedGenome
                    -1 ${output.dir}/temp_trans_unmap_reads.1.fastq
                    -2 ${output.dir}/temp_trans_unmap_reads.2.fastq
                    -S /dev/null ;
                cat $output3 >> $output1 ;
                cat $output4 >> $output2 ;
                rm  ${output.dir}/temp_trans_unmap_reads.1.fastq
                    ${output.dir}/temp_trans_unmap_reads.2.fastq
                    ${output.dir}/${branch}_trim1.fastq
                    ${output.dir}/${branch}_trim2.fastq ;
            ""","prepare_reads"
        }
    }
}

//Cat read pairs into a single file
cat_reads = {
    if (inputs.size() == 1) return
    output.dir=jaffa_output+branch
    exec "cat $input1.fastq $input2.fastq > $output.fastq"
}


//Get read which either align discordantly or not at all
get_unmapped = {
    doc "Get Unmapped"
    output.dir=jaffa_output+branch
    produce(branch+".fasta", branch+"_discordant_pairs.bam") {
        from("*_leftover_reads*.gz") {
            def input_string = ""
            if (inputs.size() == 1) {
                input_string = "$input"
            } else if (inputs.size() == 2) {
                input_string = "$input1,$input2"
            }
            exec """
               $bowtie2 -k1 -p $threads --un ${output.dir}/unmapped.fastq
                   -x $transFasta.prefix -U $input_string |
               $samtools view -s 1.0 -F 4 -S -b - |
               $samtools sort - $output2.prefix ;
               $samtools index $output2 ;
            ""","get_unmapped"
        }
        exec """
            $reformat in=${output.dir}/unmapped.fastq out=${output.dir}/temp.fasta threads=$threads ;
            $dedupe sort=d in=${output.dir}/temp.fasta out=$output1 threads=$threads absorbcontainment=f;
            rm ${output.dir}/temp.fasta ${output.dir}/unmapped.fastq 2> /dev/null
        ""","get_unmapped"
    }
}

//Like above: remove reads that don't map to the transcriptome, but this time use the assembled
//transcriptome as well as the reference
get_assembly_unmapped = {
    doc "Get assembly unmapped"
    output.dir=jaffa_output+branch
    produce(branch+"-unmapped.fasta", branch+"_discordant_pairs.bam") {
        from("*_leftover_reads*.gz") {
            def input_string = ""
            if (inputs.size() == 1) {
                input_string = "$input"
            } else if (inputs.size() == 2) {
                input_string = "$input1,$input2"
            }
            exec """
               $bowtie2 -k1 -p $threads --un ${output.dir}/unmapped_ref.fastq -x $transFasta.prefix
                   -U $input_string |
               $samtools view -s 1.0 -F 4 -S -b - |
               $samtools sort - $output2.prefix ;
               $samtools index $output2 ;
            ""","get_unmapped"
        }
        exec """
            ${bowtie2}-build ${output.dir}/${branch}.fasta ${output.dir}/${branch} ;
            $bowtie2 -k1 -p $threads --un ${output.dir}/unmapped_assembly.fastq -x ${output.dir}/${branch} -U ${output.dir}/unmapped_ref.fastq
                  -S /dev/null 2>&1 | tee ${output.dir}/log_initial_map_to_assembly ;
            $reformat in=${output.dir}/unmapped_assembly.fastq out=${output.dir}/temp.fasta threads=$threads ;
            $dedupe sort=d in=${output.dir}/temp.fasta out=$output1 threads=$threads absorbcontainment=f ;
            rm ${output.dir}/temp.fasta ${output.dir}/unmapped_assembly.fastq ${output.dir}/unmapped_ref.fastq
        ""","get_unmapped"
    }
}

//Run the de novo assembly
run_assembly = {
    doc "Align transcripts to annotation"
    output.dir=jaffa_output+branch
    produce(branch+".fasta") {
        from("*_filtered_reads.fastq*gz") {
            exec """
                time $oases_assembly_script $velveth $velvetg $oases
                ${output.dir} $output $Ks $Kmerge $transLength $threads $inputs
            ""","run_assembly"
        }
    }
}

//Align the assembled transcripts against reference gene sequences using blat
align_transcripts_to_annotation = {
    doc "Align transcripts to annotation"
    output.dir=jaffa_output+branch
    produce(input.prefix+".psl") {
        from(".fasta") {
            exec """
                function run_blat {
                    time $blat $transFasta $inputs
                        -minIdentity=$minIdTrans -minScore=$minScore -tileSize=\$1
                        -maxIntron=$maxIntron $output 2>&1 | tee ${output.dir}/log_blat ;
                } ;
                run_blat $contigTile;
                `### test for the Blat tileSize bug (version 35) ###` ;
                if [[ `cat ${output.dir}/log_blat` == *"Internal error genoFind.c"* ]] ; then 
                    echo "Blat error with tileSize=$contigTile" ;
                    echo "Let's try again with tileSize=15" ;
                    run_blat 15;
                fi ;
            ""","align_transcripts_to_annotation"
        }
    } 
}


//Align the reads to the annotation 
align_reads_to_annotation = {
    doc "Align reads to annotation"
    output.dir=jaffa_output+branch
    var CUTOFF_READ_LENGTH : 100
    var DEFAULT_TILE_SIZE : 15
    var DEFAULT_LARGE_TILE_SIZE : 18
    var SAMPLE_SIZE : 1000
    produce(input.prefix+".psl") {
        // find out the average read length so we can set blat's tileSize accordingly
        // we just use the first few thousand from the fasta file as a sample
        from(".fasta") {
            exec """
                SEQ_COUNT=`head -n $SAMPLE_SIZE $input | grep "^>" | wc -l` ;
                SUM_READ_LENGTHS=`head -n $SAMPLE_SIZE $input | grep -v ">" | tr -d "\\n" | wc --chars` ;
                AVERAGE_READ_LENGTH=`expr $SUM_READ_LENGTHS / $SEQ_COUNT` ;
                if [ $readTile -eq "0" ] ; then
                    if [ $AVERAGE_READ_LENGTH -le $CUTOFF_READ_LENGTH ] ; then
                        readTile=$DEFAULT_TILE_SIZE ;
                    else readTile=$DEFAULT_LARGE_TILE_SIZE ;
                    fi ;
                else readTile=$readTile ;
                fi ;
                echo "Using tileSize of \$readTile" ;
                function run_blat {
                    time $blat $transFasta $input -minIdentity=$minIdTrans -minScore=$minScore -tileSize=\$1
                        -maxIntron=$maxIntron $output 2>&1 | tee ${output.dir}/log_blat ;
                } ;
                run_blat \$readTile;
                `### test for the Blat tileSize bug (version 35) ###` ;
                if [[ `cat ${output.dir}/log_blat` == *"Internal error genoFind.c"* ]] ; then
                   echo "Blat error with tileSize=\$readTile" ;
                   echo "Let's try again with tileSize=15" ;
                   run_blat $DEFAULT_TILE_SIZE;
                fi ;
            ""","align_reads_to_annotation"
        }
    }
}

//Parse the blat alignment table and filter for candidate fusions (uses an R script)
filter_transcripts = {
    doc "Filter transcripts"
    output.dir=jaffa_output+branch
    produce(input.prefix+".txt") {
        from(".psl") {
            exec """
                time $R --vanilla --args $input $output $gapSize $transTable <
                    $R_filter_transcripts_script &> ${output.dir}/log_filter
            ""","filter_transcripts"
        }
    }
}

//Extract the fasta sequences for the candidate fusions into their own fasta file
extract_fusion_sequences = {
    doc "Extract fusion sequences"
    output.dir=jaffa_output+branch
    produce(input.prefix+".fusions.fa") {
        from(".txt", ".fasta") {
            exec """
                cat $input1 | awk '{print \$1}' | sed 's/^/>/g' > ${output}.temp ;
                $reformat in=$input2 out=stdout.fasta fastawrap=0 | awk '{print \$1}' |
                grep -F -A1 -f ${output}.temp | grep -v "^\\-\\-" > $output ;
                rm ${output}.temp ;
            ""","extract_fusion_sequences"
        }
    }
} 

//Map the reads back to the candidate fusion sequences
map_reads = {
    doc "Map reads back to the candidate fusion sequences"
    output.dir=jaffa_output+branch
    produce(branch+".sorted.bam") {
        from("fusions.fa","*_filtered_reads*gz") {
            def input_string=""
            if (inputs.size() == 2) {
                input_string="-U $input2"
            } else if (inputs.size() == 3) {
                input_string="-1 $input2 -2 $input3"
            }
            exec """
                ${bowtie2}-build $input1 $input1.prefix ;
                $bowtie2 $mapParams --no-unal -p $threads -x $input1.prefix $input_string | 
                $samtools view -S -b - | $samtools sort - ${output.dir}/${branch}.sorted ;
                $samtools index $output
            ""","map_reads"
        }
    }
} 
 
//Calculate the number of reads which span the breakpoint of the fusions
get_spanning_reads = {
    doc "Calculate the number of reads which span the breakpoint of the fusions"
    output.dir=jaffa_output+branch
    produce(input.txt.prefix+".reads") {
       from("txt","bam") {
           exec """ 
               $samtools view $input2 | cut -f 3,4,8  > ${output.dir}/${branch}.temp ;
               $samtools view $input2 | cut -f 10 | awk '{print length}' > ${output.dir}/${branch}.readLengths ;
               $R --vanilla --args ${output.dir}/${branch} $input1 ${output.dir}/${branch}.temp 
                   $output ${output.dir}/${branch}.readLengths $overHang < $R_get_spanning_reads_script ;
               rm ${output.dir}/${branch}.temp ${output.dir}/${branch}.readLengths
           ""","get_spanning_reads"
        }
    }
}

//Used for the direct and hybrid pipelines - In this case the spanning reads will be 1 for each
//read and the spanning pairs will be 0. 
make_simple_reads_table = {
    doc "Calculate the number of reads which span the breakpoint of the fusions"
    output.dir=jaffa_output+branch
    produce(input.txt.prefix+".reads") {
        from(".txt", "*_discordant_pairs.bam") {
            exec """
                $samtools view -H $input2 | grep "@SQ" | cut -f2 | sed 's/SN://g' > ${output.dir}/temp_gene_ids ;
                $R --no-save --args $input1 $transTable ${output.dir}/temp_gene_ids ${output.dir}/paired_contigs.temp
                    < $R_get_spanning_reads_direct_script1 ;
                function get_spanning_pairs {
                    gene=`echo \$1 | cut -f1 -d"?"` ;
                    g1=`echo \$1 | cut -f2 -d"?"` ;
                    g2=`echo \$1 | cut -f3 -d "?"` ;
                    $samtools view $input2 \$g1 | cut -d "/" -f1 | sort -u > ${output.dir}/g1 ;
                    $samtools view $input2 \$g2 | cut -d "/" -f1 | sort -u > ${output.dir}/g2 ;
                    left=`cat ${output.dir}/g1 | wc -l` ;
                    right=`cat ${output.dir}/g2 | wc -l` ;
                    both=`cat ${output.dir}/g1 ${output.dir}/g2 | sort -u | wc -l` ;
                    echo -e "\$gene\t\$(( \$left + \$right - \$both ))" ;
                } ;
                cat ${output.dir}/paired_contigs.temp | while read line ; do
                    get_spanning_pairs "\$line" >> ${output.dir}/spanning_pair_counts.temp ;
                done ;
                $R --no-save --args $input1 ${output.dir}/spanning_pair_counts.temp $output < $R_get_spanning_reads_direct_script2 ;
                rm ${output.dir}/temp_gene_ids ${output.dir}/spanning_pair_counts.temp ${output.dir}/paired_contigs.temp ${output.dir}/g1 ${output.dir}/g2
            ""","make_simple_reads_table"
        }
    }
} 

make_fasta_reads_table = {
    doc "Make fasta reads table"
    output.dir=jaffa_output+branch
    produce(input.txt.prefix+".reads") {
        from("txt") {
            exec """
                echo  -e "transcript\tbreak_min\tbreak_max\tfusion_genes\tspanning_pairs\tspanning_reads" > $output ; 
                awk '{ print \$1"\t"\$2"\t"\$3"\t"\$4"\t"0"\t"1}' $input | sort -u  >> $output
            ""","make_fasta_reads_table"
        }
    }
}


//This stage is only used the by hybrid mode.
//It concatenates the fusions sequence files, then the read files.
merge_assembly_and_unmapped_reads_candidates = {
    doc "Concatenate fusion sequence files and reads files (hybrid only)"
    output.dir=jaffa_output+branch
    produce(branch+".all.fusions.fa", branch+".all.reads") {
        from("fusions.fa", branch+".fusions.fa", "reads", branch+".reads") {
            exec """
                cat $input1 $input2 > $output1 ;
                cp $input3 $output2 ; tail -n+2 $input4 >> $output2
            ""","merge_assembly_and_unmapped_reads_candidates"
        }
    }
}


//Align candidate fusions to the genome
align_transcripts_to_genome = {
    doc "Align candidate fusions to the genome"
    output.dir=jaffa_output+branch
    produce(branch+"_genome.psl") {
        from(".fusions.fa") {
            exec """
                set -o pipefail; $blat $genomeFasta $input1 -minScore=$minScore $output 2>&1 | tee ${output.dir}/log_genome_blat
            ""","align_transcripts_to_genome"
        }
    }
}

//Do a bit more filtering and compile the final filtered list (uses an R script)
get_final_list = {
    doc "Get final list"
    output.dir=jaffa_output+branch
    produce(branch+".summary") {
        from(".psl", ".reads") {
            exec """
                $R --vanilla --args $input1 $input2 $transTable $knownTable $finalGapSize $exclude $output < $R_get_final_list
            ""","get_final_list"
        }
    }
}

//Compile the results from multiple samples into an excel .csv table
//Make a fasta file with the candidates
compile_all_results = {
    doc "Compile all results"
    var type : "" 
    if (jaffa_output) {
        output.dir=jaffa_output
    }
    produce(outputName+".fasta",outputName+".csv") {
        // change to the jaffa output directory
        exec """
            cd ${output.dir} ;
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
        ""","compile_all_results"
    }
}
