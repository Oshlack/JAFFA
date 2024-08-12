# syntax=docker/dockerfile:1

# must use bookworm for directly compatible R version >=4.4.0
FROM debian:bullseye-slim AS base

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

# fix version of Java for bpipe 0.9.9.2
RUN apt-get install -y r-base-core r-bioc-iranges openjdk-11-jre --no-install-recommends

# install libncurses.so.6 for minimap2
RUN apt-get install -y libncurses6 --no-install-recommends


FROM base AS build
# install useful tools for building non-R software
# this is placed after the R installs so that they can
# be cached
RUN apt-get install -y build-essential wget make cmake unzip zip python3 zlib1g-dev libncurses-dev


WORKDIR /JAFFA

COPY ./install_linux64.sh .
COPY ./src ./src

RUN ./install_linux64.sh

# apply a patch to the issue at https://github.com/ssadedin/bpipe/pull/293,
# which seems to disproportionately affect Docker images
RUN patch ./tools/bpipe-0.9.9.2/bin/bpipe ./src/pid-error.patch

COPY . .

# set a reference file path to /ref for easy binding
RUN sed -i 's/refBase = codeBase/refBase = "\/ref"/' JAFFA_stages.groovy

# we only need this to run the image
FROM base

LABEL org.opencontainers.image.title "JAFFA" \
      org.opencontainers.image.description "High sensitivity transcriptome-focused fusion gene detection." \
      org.opencontainers.image.authors "Davidson, N.M., Majewski, I.J. & Oshlack, A." \
      org.opencontainers.image.source "https://github.com/Oshlack/JAFFA" \
      org.opencontainers.image.documentation "https://github.com/Oshlack/JAFFA/wiki/HowToSetUpJAFFA#docker"


COPY --from=build /JAFFA /JAFFA
# RUN strip --strip-debug /usr/local/lib/R/site-library/*/libs/*.so

RUN apt remove -y gcc
RUN apt autoremove

ENV PATH=${PATH}:/JAFFA/tools/bin:/JAFFA

ENTRYPOINT ["/JAFFA/tools/bin/bpipe", "run"]
