FROM ubuntu:18.04

MAINTAINER Rebecca Louise Evans (rebecca.louise.evans@gmail.com)

LABEL Description="This image is used to run JAFFA" Version="2.1"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y \
    && apt-get install -y \
        bowtie2 \
        bzip2 \
        g++ \
        git \
        gzip \
        libncurses5-dev \
        libpng-dev \
        libtool \
        lmod \
        make \
        openjdk-8-jdk \
        python \
        r-base \
        r-base-dev \
        time \
        trimmomatic \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

VOLUME /data/batch

WORKDIR /opt

# Set Standard settings
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64
ENV CLASSPATH .
ENV CP ${CLASSPATH}
ENV BASH_ENV /usr/share/lmod/lmod/init/bash
ENV PATH /usr/bin:/bin:/usr/local/bin:/opt/bin:/opt/bbmap
#ENV MODULEPATH

# setup lmod
RUN ln -s /usr/share/lmod/lmod/init/profile /etc/profile.d/modules.sh
RUN ln -s /usr/share/lmod/lmod/init/cshrc /etc/profile.d/modules.csh

# install minimap2
RUN git clone https://github.com/lh3/minimap2
RUN make -C minimap2
RUN cp minimap2/minimap2 /usr/local/bin
RUN rm -rf minimap2

# install blastn
ENV BLAST_BASE_URL https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+
#RUN wget ${BLAST_BASE_URL}/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz -O - | tar xz
RUN wget ${BLAST_BASE_URL}/LATEST/ncbi-blast-$(wget -O - ${BLAST_BASE_URL}/VERSION)+-x64-linux.tar.gz -O - | tar -xz
RUN cp ncbi-blast-*+/bin/blastn /usr/local/bin/blastn
RUN rm -rf ncbi-blast-*+

# install samtools/1.1 (due to backwards incompatibility)
RUN wget --max-redirect 5 https://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2 -O - | tar -xj
RUN make prefix=/usr/local install -C samtools-1.1
RUN rm -rf samtools-1.1

# install bbmap
RUN wget --max-redirect 5 https://sourceforge.net/projects/bbmap/files/latest/download?source=files -O - | tar -xz
RUN make -C bbmap/jni -f makefile.linux
RUN find bbmap -name \*.sh -type f -exec ln -s '{}' /usr/local/bin/ \;

# install oases and velvet
RUN git clone --recursive https://github.com/dzerbino/oases.git
RUN make -C oases/velvet/ MAXKMERLENGTH=37 LONGSEQUENCES=1
RUN make -C oases/ MAXKMERLENGTH=37 LONGSEQUENCES=1 'VELVET_DIR=velvet'
RUN cp oases/velvet/velvetg /usr/bin/
RUN cp oases/velvet/velveth /usr/bin/
RUN cp oases/oases /usr/bin
RUN rm -rf oases

# install fastx_toolkit, creates ${WORKDIR}/bin directory
#RUN wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 -O - | tar -xj
RUN git clone https://github.com/agordon/libgtextutils.git
RUN cd libgtextutils && ./reconf && ./configure && make && make install
RUN rm -rf libgtextutils
RUN git clone https://github.com/agordon/fastx_toolkit.git
# --disable-wall (due to -Werror)
#RUN cd fastx_toolkit && ./reconf && ./configure && make && make install
RUN cd fastx_toolkit && ./reconf && ./configure --disable-wall && make && make install

# install blat
RUN wget http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip
RUN unzip blatSrc35.zip
RUN rm blatSrc35.zip
ENV MACHTYPE=x86_64
RUN mkdir -p ${HOME}/bin/${MACHTYPE}
RUN make -C blatSrc
RUN mv ${HOME}/bin/${MACHTYPE}/* /usr/local/bin
RUN rmdir ${HOME}/bin/${MACHTYPE}
RUN rmdir ${HOME}/bin
RUN rm -rf blatSrc

# install bpipe
RUN git clone https://github.com/ssadedin/bpipe.git
# gradle properties holds http.proxyHost and http.proxyPort
# COPY gradle.properties bpipe/gradle.properties
# RUN bpipe/gradlew -p bpipe dist
RUN cd bpipe; ./gradlew dist
RUN mv bpipe/build/stage/bpipe* .
RUN rm -rf bpipe
RUN mv bpipe* bpipe
RUN chmod 755 /opt/bpipe/bin/*
RUN find /opt/bpipe/bin -type f -exec ln -s '{}' /usr/local/bin/ \;

# install R dependencies (required by JAFFA)
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("IRanges")'

# install jaffa
RUN git clone https://github.com/Oshlack/JAFFA.git -b master
RUN mkdir -p JAFFA/tools/bin
RUN g++ -std=c++11 -O3 -o JAFFA/tools/bin/process_transcriptome_align_table JAFFA/src/process_transcriptome_align_table.c++
RUN g++ -O3 -o JAFFA/tools/bin/extract_seq_from_fasta JAFFA/src/extract_seq_from_fasta.c++
RUN g++ -std=c++11 -O3 -o JAFFA/tools/bin/make_simple_read_table JAFFA/src/make_simple_read_table.c++
#RUN g++ -std=c++11 -O3 -o JAFFA/bin/bypass_genomic_alignment JAFFA/src/bypass_genomic_alignment.c++
ENV PATH ${PATH}:/opt/JAFFA/tools/bin

# set the tools
COPY tools.groovy JAFFA/tools.groovy
RUN chmod 644 JAFFA/tools.groovy

COPY convert_jaffa_to_bedpe.py /usr/local/bin/convert_jaffa_to_bedpe.py
RUN chmod 755 /usr/local/bin/convert_jaffa_to_bedpe.py

WORKDIR /data/batch

CMD ["bpipe", "run", "-p", "fastqInputFormat='%_*.fastq.gz'", "-p", "refBase=/data/reference", "-p", "genome=hg38", "-p", "annotation=genCode22", "/opt/JAFFA/JAFFA_direct.groovy", "/data/example/BT474-demo_1.fastq.gz", "/data/example/BT474-demo_2.fastq.gz"]
