# syntax=docker/dockerfile:1

# ================================================
# BASE BUILD - everything will use this
# ================================================

# must use bookworm for directly compatible R version >=4.4.0
FROM debian:bullseye AS base

# get latest version of R: https://cran.r-project.org/bin/linux/debian/
RUN apt-get update
RUN apt-get install -y gnupg2
RUN apt-get install -y perl python3

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

# install procps to get ps, needed in Nextflow
RUN apt-get install -y procps

RUN apt-get install -y --no-install-recommends \
    r-base-core=4.5.1-1~bullseyecran.0

# fix version of Java for bpipe 0.9.9.2
RUN apt-get install -y --no-install-recommends openjdk-11-jre

# install libncurses.so.6 for minimap2
RUN apt-get install -y libncurses6 --no-install-recommends

# ================================================
# R BUILD: BUILD AND INSTALL ALL R PACKAGES
# ================================================
FROM base AS r-build
RUN apt-get install -y \
    r-base-core=4.5.1-1~bullseyecran.0 \
    r-base-dev=4.5.1-1~bullseyecran.0

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && R -e "BiocManager::install('IRanges')"


# ================================================
# SOFTWARE BUILD: INSTALL ALL NON-R SOFTWARE
# ================================================
FROM base AS software-build
# install useful tools for building non-R software
# this is placed after the R installs so that they can
# be cached
RUN apt-get install -y build-essential wget make cmake unzip zip python3 zlib1g-dev libncurses-dev

WORKDIR /JAFFA

COPY ./install_linux64.sh .
COPY ./src ./src

RUN ./install_linux64.sh

COPY . .

# set a reference file path to /ref for easy binding
RUN sed -i 's/refBase = codeBase/refBase = "\/ref"/' JAFFA_stages.groovy


# ================================================
# RUNNER
# ================================================
FROM base

LABEL org.opencontainers.image.title "JAFFA" \
    org.opencontainers.image.description "High sensitivity transcriptome-focused fusion gene detection." \
    org.opencontainers.image.authors "Davidson, N.M., Majewski, I.J. & Oshlack, A." \
    org.opencontainers.image.source "https://github.com/Oshlack/JAFFA" \
    org.opencontainers.image.documentation "https://github.com/Oshlack/JAFFA/wiki/HowToSetUpJAFFA#docker"


# copy over software
COPY --from=software-build /JAFFA /JAFFA

# copy over R packages
COPY --from=r-build /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=r-build /usr/lib/R/site-library /usr/lib/R/site-library
COPY --from=r-build /usr/lib/R/library /usr/lib/R/library

# RUN strip --strip-debug /usr/local/lib/R/site-library/*/libs/*.so

RUN apt remove -y gcc
RUN apt autoremove -y

ENV PATH=${PATH}:/JAFFA/tools/bin:/JAFFA
# ENV _JAVA_OPTIONS="--add-opens java.base/java.lang=ALL-UNNAMED"

ENTRYPOINT ["/JAFFA/tools/bin/bpipe", "run"]
