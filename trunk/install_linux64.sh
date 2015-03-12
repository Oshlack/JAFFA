#!/bin/bash

## This script will install the tools required for the JAFFA pipeline.
## It will fetched each tool from the web and placed into the tools/ subdirectory.
## Paths to all installed tools can be found in the file tools.groovy at the
## end of execution of this script. These paths can be changed if a different
## version of software is required. Note that R must be installed manually
##
## Last Modified: 22nd October by Nadia Davidson

mkdir -p tools/bin 
cd tools 

#a list of which programs need to be installed
commands="bpipe velveth velvetg oases trimmomatic samtools bowtie2 blat dedupe reformat"

#installation method
function bpipe_install {
   wget http://download.bpipe.org/versions/bpipe-0.9.8.6_rc2.tar.gz
   tar -zxvf bpipe-0.9.8.6_rc2.tar.gz ; rm bpipe-0.9.8.6_rc2.tar.gz
   ln -s $PWD/bpipe-0.9.8.6_rc2/bin/* $PWD/bin/
}

function velveth_install {
    wget http://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
    tar -zxvf velvet_1.2.10.tgz ; rm velvet_1.2.10.tgz
    make -C velvet_1.2.10/ MAXKMERLENGTH=37 LONGSEQUENCES=1 #OPENMP=1
    ln -s $PWD/velvet_1.2.10/velvetg $PWD/bin/
    ln -s $PWD/velvet_1.2.10/velveth $PWD/bin/
}

function oases_install {
    wget http://www.ebi.ac.uk/~zerbino/oases/oases_0.2.08.tgz
    tar -zxvf oases_0.2.08.tgz ; rm oases_0.2.08.tgz
    make -C oases_0.2.8/ MAXKMERLENGTH=37 LONGSEQUENCES=1 'VELVET_DIR=../velvet_1.2.10'
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
   wget http://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2
   tar -jxvf samtools-1.1.tar.bz2
   rm samtools-1.1.tar.bz2
   make prefix=$PWD install -C samtools-1.1/
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
   chmod +x $PWD/bin/blat
}

function fasta_formatter_install {
    wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    tar -jxvf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    rm fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
}

function dedupe_install {
    wget http://sourceforge.net/projects/bbmap/files/BBMap_33.41_java7.tar.gz
    tar -zxvf BBMap_33.41_java7.tar.gz
    rm BBMap_33.41_java7.tar.gz
    for script in `ls $PWD/bbmap/*.sh` ; do
	s=`basename $script`
	s_pre=`echo $s | sed 's/.sh//g'`
	echo "$PWD/bbmap/$s \$@" > $PWD/bin/$s_pre
	chmod +x $PWD/bin/$s_pre
    done
}


echo "// Path to tools used by the JAFFA pipeline" > ../tools.groovy

for c in $commands ; do 
    c_path=`which bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then 
	echo "$c not found, fetching it"
	${c}_install
	c_path=`which bin/$c 2>/dev/null`
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
    c_path=`which bin/$c 2>/dev/null`
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



