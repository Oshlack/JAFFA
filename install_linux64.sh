#!/bin/bash

## This script will install the tools required for the JAFFA pipeline.
## The script will first check whether the required tool is install and if
## not it will be fetched from the web and placed into the tools/ subdirectory
## paths to all installed tools can be found in the file tools.groovy at the
## end of execution of this script. These paths can be changed if a different
## version of software is required.
##
## Last Modified: 16th September by Nadia Davidson

mkdir -p tools/bin 
cd tools 

#export PATH=$PATH:$PWD/bin
PATH=$PWD/bin:$PWD/../my_usrbin/:/bin:/usr/local/sbin:/usr/sbin:/sbin
echo $PATH

#a list of which programs need to be installed
commands="bpipe velveth velvetg oases trimmomatic samtools bowtie2 blat fasta_formatter fastx_collapser R"

#installation method
function bpipe_install {
   wget http://download.bpipe.org/versions/bpipe-0.9.8.6_rc2.tar.gz
   gunzip bpipe-0.9.8.6_rc2.tar.gz ; tar -xvf bpipe-0.9.8.6_rc2.tar ; rm bpipe-0.9.8.6_rc2.tar
   ln -s $PWD/bpipe-0.9.8.6_rc2/bin/* $PWD/bin/
}

function velveth_install {
    wget http://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
    gunzip velvet_1.2.10.tgz ; tar -xvf velvet_1.2.10.tar ; rm velvet_1.2.10.tar
    make -C velvet_1.2.10/ MAXKMERLENGTH=37 OPENMP=1
    ln -s $PWD/velvet_1.2.10/velvetg $PWD/bin/
    ln -s $PWD/velvet_1.2.10/velveth $PWD/bin/
}

function oases_install {
    wget http://www.ebi.ac.uk/~zerbino/oases/oases_0.2.08.tgz
    gunzip oases_0.2.08.tgz ; tar -xvf oases_0.2.08.tar ; rm oases_0.2.08.tar
    make -C oases_0.2.8/ MAXKMERLENGTH=35 'VELVET_DIR=../velvet_1.2.10'
    ln -s $PWD/oases_0.2.8/oases $PWD/bin/
}

function trimmomatic_install {
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
    unzip Trimmomatic-0.32.zip ; rm Trimmomatic-0.32.zip
    echo "java -jar $PWD/Trimmomatic-0.32/trimmomatic-0.32.jar \$*"  > Trimmomatic-0.32/trimmomatic.sh
    chmod +x Trimmomatic-0.32/trimmomatic.sh
    ln -s $PWD/Trimmomatic-0.32/trimmomatic.sh $PWD/bin/trimmomatic
}

function samtools_install {
    wget http://sourceforge.net/projects/samtools/files/samtools/1.0/samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
    bunzip2 samtools-bcftools-htslib-1.0_x64-linux.tar.bz2
    tar -xvf samtools-bcftools-htslib-1.0_x64-linux.tar
    rm samtools-bcftools-htslib-1.0_x64-linux.tar
    ln -s $PWD/samtools-bcftools-htslib-1.0_x64-linux/bin/* $PWD/bin/
}

function bowtie2_install {
    wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip
    unzip bowtie2-2.2.3-linux-x86_64.zip ; rm bowtie2-2.2.3-linux-x86_64.zip
    ln -s $PWD/bowtie2-2.2.3/bowtie2 $PWD/bin
    ln -s $PWD/bowtie2-2.2.3/bowtie2-build $PWD/bin
}

function blat_install {
   wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
   mv blat $PWD/bin
}

function fasta_formatter_install {
    wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    bunzip2 fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 
    tar -xvf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar
    rm fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar
}

function R_install {
    echo  "Please go to http://www.r-project.org/ and follow the installation instructions"
}

echo "// Path to tools used by the JAFFA pipeline" > ../tools.groovy

for c in $commands ; do 
    echo "checking if $c is installed..." 
    c_path=`which $c 2>/dev/null`
    if [ -z $c_path ] ; then 
	echo "$c not found, fetching it"
	${c}_install
	c_path=`which $c 2>/dev/null`
    fi
    if [ -z $c_path ] ; then 
	echo ""
	echo "WARNING!!! Could not install $c. You will need to install it manually!" ; 
	echo "Then edit tools.groovy with the tool's path"
	echo ""
    fi
    echo "$c=\"$c_path\"" >> ../tools.groovy
done

echo "All done. Please check that the file, tools.groovy, lists all paths correctly."




