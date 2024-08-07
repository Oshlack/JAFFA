# syntax=docker/dockerfile:1

#######################################
#              BASE IMAGE             #
#######################################
#
# This image contains basic packages which are needed in both
# the build image and the runtime image. In this case,
# it sets the locale information to en_US.UTF8 for R,
# installs a modern version of R, and also the Java runtime.


# bullseye is stable and contains openjdk-11-jre in stable by default
FROM debian:bullseye-slim AS base

# get latest version of R: https://cran.r-project.org/bin/linux/debian/
RUN apt-get update
RUN apt-get install -y gnupg2
RUN apt-get install -y perl

# set locale
RUN apt-get install -y locales
RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

RUN gpg --keyserver keyserver.ubuntu.com \
    --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' | \
    tee /etc/apt/trusted.gpg.d/cran_debian_key.asc

RUN echo 'deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/' >> /etc/apt/sources.list

RUN apt-get update

# fix version of Java for bpipe 0.9.9.2
RUN apt-get install -y r-base-core r-bioc-iranges openjdk-11-jre --no-install-recommends

# # install bioconductor
# RUN R -e 'install.packages("BiocManager")'

# # install needed packages
# RUN R -e 'BiocManager::install("IRanges")'

#######################################
#            BUILD PACKAGE            #
#######################################
#
# This is the build image, which will download and install
# all the dependencies using the bundled ./install_linux64.sh
# script.


FROM base AS build
# install useful tools for building non-R software
# this is placed after the R installs so that they can
# be cached
RUN apt-get install -y build-essential wget make cmake unzip zip python3 zlib1g-dev libncurses5-dev libncursesw5-dev


WORKDIR /JAFFA

COPY ./install_linux64.sh .
COPY ./src ./src

RUN ./install_linux64.sh

COPY . .

# set a reference file path to /ref for easy binding
RUN sed -i 's/refBase = codeBase/refBase = "\/ref"/' JAFFA_stages.groovy


#######################################
#              RUN IMAGE              #
#######################################
#
# This is the runtime image, which is what will be executed.
# Ensure that all necessary files are copied into this image!


FROM base

# add authorship information
LABEL org.opencontainers.image.title "JAFFA" \
      org.opencontainers.image.description "High sensitivity transcriptome-focused fusion gene detection." \
      org.opencontainers.image.authors "Davidson, N.M., Majewski, I.J. & Oshlack, A." \
      org.opencontainers.image.source "https://github.com/Oshlack/JAFFA" \
      org.opencontainers.image.documentation "https://github.com/Oshlack/JAFFA/wiki/HowToSetUpJAFFA#docker"

# copy packages from the build repository
COPY --from=build /JAFFA /JAFFA

# clean up any packages lying around, if they are installed accidentially for some reason
RUN apt remove -y gcc
RUN apt autoremove

ENV PATH=${PATH}:/JAFFA/tools/bin:/JAFFA

# this flag disables the "WARNING: An illegal reflective access operation has occurred" message.
# this can currently be safely disabled, as the Java runtime version is fixed at openjdk-11-jre,
# where this is safely permitted. in the future, if Java or Groovy is updated, this variable
# should be removed.
ENV JDK_JAVA_OPTIONS="$JDK_JAVA_OPTIONS --illegal-access=permit"

ENTRYPOINT ["/JAFFA/tools/bin/bpipe", "run"]
