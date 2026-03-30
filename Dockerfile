# syntax=docker/dockerfile:1

# ================================================
# BASE BUILD - everything will use this
# ================================================

FROM debian:bookworm AS base

# get latest version of R: https://cran.r-project.org/bin/linux/debian/
RUN apt-get update
RUN apt-get install -y gnupg2
RUN apt-get install -y perl python3

# set locale
RUN apt-get install -y locales
RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# install procps to get ps, needed in Nextflow
RUN apt-get install -y procps

# install libncurses.so.6 for minimap2
RUN apt-get install -y libncurses6 --no-install-recommends


# ================================================
# SOFTWARE BUILD: INSTALL ALL SOFTWARE
# ================================================
FROM base AS software-build
RUN apt-get install -y build-essential wget make cmake unzip zip python3 zlib1g-dev libncurses-dev libgomp1

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

# fix version of Java for bpipe 0.9.9.2
# while bpipe requires openjdk 11, JAFFA testing is in Java 17 so
# openjdk-17 is preferred
RUN apt-get install -y --no-install-recommends openjdk-17-jre

# dependencies for individual tools
RUN apt-get install -y libgomp1
    
LABEL org.opencontainers.image.title "JAFFA" \
    org.opencontainers.image.description "High sensitivity transcriptome-focused fusion gene detection." \
    org.opencontainers.image.authors "Davidson, N.M., Majewski, I.J. & Oshlack, A." \
    org.opencontainers.image.source "https://github.com/Oshlack/JAFFA" \
    org.opencontainers.image.documentation "https://github.com/Oshlack/JAFFA/wiki/HowToSetUpJAFFA#docker"


# copy over software
COPY --from=software-build /JAFFA /JAFFA

RUN apt remove -y gcc
RUN apt autoremove -y

ENV PATH=${PATH}:/JAFFA/tools/bin:/JAFFA

ENTRYPOINT ["/JAFFA/tools/bin/bpipe", "run"]
