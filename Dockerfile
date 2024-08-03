# syntax=docker/dockerfile:1

# must use bookworm for directly compatible R version >=4.4.0
FROM debian:bullseye-slim AS base

# get latest version of R: https://cran.r-project.org/bin/linux/debian/
RUN apt-get update && apt-get install -y gnupg2

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



# we only need this to run the image
FROM base

COPY --from=build /JAFFA /JAFFA
# RUN strip --strip-debug /usr/local/lib/R/site-library/*/libs/*.so

RUN apt remove -y gcc
RUN apt autoremove

ENV PATH=${PATH}:/JAFFA/tools/bin:/JAFFA

# set locale
ENV export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

ENTRYPOINT ["/JAFFA/tools/bin/bpipe", "run"]
