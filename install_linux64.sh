#!/bin/bash

## This script will install the tools required for the JAFFA pipeline.
## It will fetched each tool from the web and placed into the tools/ subdirectory.
## Paths to all installed tools can be found in the file tools.groovy at the
## end of execution of this script. These paths can be changed if a different
## version of software is required. Note that R must be installed manually
##
## Last Modified: Sep. 2021 by Nadia Davidson

mkdir -p tools/bin 
cd tools 

#a list of which programs need to be installed
commands="bpipe velveth velvetg oases trimmomatic samtools bowtie2 blat dedupe reformat extract_seq_from_fasta make_simple_read_table blastn minimap2 process_transcriptome_align_table make_3_gene_fusion_table"

#installation methods

function minimap2_install {
   wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2
   wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2
   tar -xvf minimap2-2.17_x64-linux.tar.bz2 ; rm minimap2-2.17_x64-linux.tar.bz2
#   make -C minimap2-2.17
   cp minimap2-2.17_x64-linux/minimap2 bin/
}


function blastn_install {
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz
    tar -xvf ncbi-blast-2.9.0+-x64-linux.tar.gz
    rm ncbi-blast-2.9.0+-x64-linux.tar.gz
    cp ncbi-blast-2.9.0+/bin/blastn bin/blastn
}

function make_3_gene_fusion_table_install {
    g++ -std=c++11 -O3 -o bin/make_3_gene_fusion_table ../src/make_3_gene_fusion_table.c++
}

function extract_seq_from_fasta_install {
    g++ -std=c++11 -O3 -o bin/extract_seq_from_fasta ../src/extract_seq_from_fasta.c++
}

function make_simple_read_table_install {
    g++ -std=c++11 -O3 -o bin/make_simple_read_table ../src/make_simple_read_table.c++
}

function process_transcriptome_align_table_install {
    g++ -std=c++11 -O3 -o bin/process_transcriptome_align_table ../src/process_transcriptome_align_table.c++
}

function make_count_table_install {
    g++ -O3 -o bin/make_count_table ../src/make_count_table.c++
}

#function bypass_genomic_alignment_install {
#    g++ -std=c++11 -O3 -o bin/bypass_genomic_alignment ../src/bypass_genomic_alignment.c++
#}

function bpipe_install {
   wget -O bpipe-0.9.9.2.tar.gz https://github.com/ssadedin/bpipe/releases/download/0.9.9.2/bpipe-0.9.9.2.tar.gz
   tar -zxvf bpipe-0.9.9.2.tar.gz ; rm bpipe-0.9.9.2.tar.gz
   ln -s $PWD/bpipe-0.9.9.2/bin/* $PWD/bin/
}

function velveth_install {
    wget -O velvet-1.2.10.tar.gz https://github.com/dzerbino/velvet/archive/v1.2.10.tar.gz
    tar -zxvf velvet-1.2.10.tar.gz ; rm velvet-1.2.10.tar.gz
    make -C velvet-1.2.10/ MAXKMERLENGTH=37 LONGSEQUENCES=1 #OPENMP=1
    ln -s $PWD/velvet-1.2.10/velvetg $PWD/bin/
    ln -s $PWD/velvet-1.2.10/velveth $PWD/bin/
}

function oases_install {
    wget -O oases_0.2.09.tgz https://github.com/dzerbino/oases/archive/0.2.09.tar.gz
    tar -zxvf oases_0.2.09.tgz ; rm oases_0.2.09.tgz
    make -C oases-0.2.09/ MAXKMERLENGTH=37 LONGSEQUENCES=1 'VELVET_DIR=../velvet-1.2.10'
    ln -s $PWD/oases-0.2.09/oases $PWD/bin/
}

function trimmomatic_install {
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
    unzip Trimmomatic-0.32.zip ; rm Trimmomatic-0.32.zip
    echo "java -jar $PWD/Trimmomatic-0.32/trimmomatic-0.32.jar \$*"  > Trimmomatic-0.32/trimmomatic.sh
    chmod +x Trimmomatic-0.32/trimmomatic.sh
    ln -s $PWD/Trimmomatic-0.32/trimmomatic.sh $PWD/bin/trimmomatic
}

function samtools_install {
   wget --no-check-certificate http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2
   tar -jxvf samtools-1.1.tar.bz2
   rm samtools-1.1.tar.bz2
   make prefix=$PWD install -C samtools-1.1/
}

function bowtie2_install {
    wget --no-check-certificate http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip
    unzip bowtie2-2.4.4-linux-x86_64.zip ; rm bowtie2-2.4.4-linux-x86_64.zip
    ln -s $PWD/bowtie2-2.4.4-linux-x86_64/bowtie2 $PWD/bin
    ln -s $PWD/bowtie2-2.4.4-linux-x86_64/bowtie2-build $PWD/bin
}

function blat_install {
   wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
   mv blat $PWD/bin
   chmod +x $PWD/bin/blat
}

function fasta_formatter_install {
    wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    tar -jxvf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    rm fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
}

function dedupe_install {
    wget --no-check-certificate https://sourceforge.net/projects/bbmap/files/BBMap_36.59.tar.gz
    tar -zxvf BBMap_36.59.tar.gz
    rm BBMap_36.59.tar.gz
    for script in `ls $PWD/bbmap/*.sh` ; do
	s=`basename $script`
	s_pre=`echo $s | sed 's/.sh//g'`
	echo "$PWD/bbmap/$s \$@" > $PWD/bin/$s_pre
	chmod +x $PWD/bin/$s_pre
    done
}

#Check if the version of gcc is >= 4.9
gcc_version=`gcc -dumpversion`
gcc_check=`echo -e "$gcc_version\n4.9" | sort -n | tail -n1`
if [[ $gcc_chek = "4.9" ]] 
then 
   echo "Your version of gcc is $gcc_version."
   echo "gcc must be >= 4.9 to install JAFFA. Exiting..."
   exit 1
fi

echo "gcc check passed"

echo "// Path to tools used by the JAFFA pipeline" > ../tools.groovy

for c in $commands ; do 
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then 
	echo "$c not found, fetching it"
	${c}_install
	c_path=`which $PWD/bin/$c 2>/dev/null`
    fi
    echo "$c=\"$c_path\"" >> ../tools.groovy
done

#finally check that R is install
R_path=`which R 2>/dev/null`
if [ -z $R_path ] ; then
    echo "R not found!"
    echo "Please go to http://www.r-project.org/ and follow the installation instructions."
    echo "Note that the IRanges R package must be installed."
fi
echo "R=\"$R_path\"" >> ../tools.groovy

#loop through commands to check they are all installed
echo "Checking that all required tools were installed:"
Final_message="All commands installed successfully!"
for c in $commands ; do
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then
	echo -n "WARNING: $c could not be found!!!! " 
	echo "You will need to download and install $c manually, then add its path to tools.groovy"
	Final_message="WARNING: One or more command did not install successfully. See warning messages above. \
                       You will need to correct this before running JAFFA."
    else 
        echo "$c looks like it has been installed"
    fi
done
echo "**********************************************************"
echo $Final_message



