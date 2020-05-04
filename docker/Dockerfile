FROM r-base:3.6.2



### make required directories
RUN mkdir /home/runDir && \
    mkdir /home/runDir/bin && \
    mkdir /usr/local/packages

### move to this directory to install programs
WORKDIR /usr/local/packages



##########
# Install programs

### get java, cpanminus, curl and mafft
RUN apt-get update -y && \
    apt-get install -y default-jre \
    mafft \
    curl \
    cpanminus \
    libdbi-perl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

### get blast
ARG BLAST_version=2.10.0

RUN wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_version}/ncbi-blast-${BLAST_version}+-x64-linux.tar.gz && \
    tar -zxf ncbi-blast-${BLAST_version}+-x64-linux.tar.gz && \
    rm ncbi-blast-${BLAST_version}+-x64-linux.tar.gz

### get fasttree
ARG FT_version=2.1.11

RUN wget -q http://microbesonline.org/fasttree/FastTree-${FT_version}.c && \
    gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree-${FT_version}.c -lm && \
    rm FastTree-${FT_version}.c && \
    chmod 777 FastTree

### get dendroscope
COPY dendroInstallInput.txt dendroInstallInput.txt
ARG DEND_version=3_7_2

RUN wget -q https://software-ab.informatik.uni-tuebingen.de/download/dendroscope/Dendroscope_unix_${DEND_version}.sh && \
     chmod 777 Dendroscope_unix_${DEND_version}.sh && \
     ./Dendroscope_unix_${DEND_version}.sh < dendroInstallInput.txt && \
     rm Dendroscope_unix_${DEND_version}.sh dendroInstallInput.txt

### get imageMagick
ARG IM_version=7.0.9-27

RUN wget -q https://imagemagick.org/download/releases/ImageMagick-${IM_version}.tar.gz && \
     tar -zxf ImageMagick-${IM_version}.tar.gz && \
     cd ImageMagick-${IM_version} && \
     ./configure && \
     make && \
     make install && \
     ldconfig /usr/local/lib && \
     cd ../ && \
     rm ImageMagick-${IM_version}.tar.gz


#####################
### Install R packages
RUN R -e 'install.packages("tidyverse", dependencies = TRUE, quiet = TRUE, Ncpus = 16)'
RUN R -e 'install.packages(c("knitrBootstrap", "gridExtra", "ggwordcloud"), dependencies = TRUE, quiet = TRUE, Ncpus = 16)'


### get MVesuviusC stuff from github
ARG PE_version=v0.1.1
ARG GT_version=V2.2.2
ARG LPBA_version=v0.1.1

RUN mkdir primerEvaluation && mkdir getTaxa && mkdir localPrimerBlastAlternative

RUN wget -q https://github.com/MVesuviusC/primerEvaluation/archive/${PE_version}.zip && \
    unzip -jd primerEvaluation/ ${PE_version}.zip && \
    rm ${PE_version}.zip && \
    wget -q https://github.com/MVesuviusC/getTaxa/archive/${GT_version}.zip && \
    unzip -jd getTaxa/ ${GT_version}.zip && \
    rm ${GT_version}.zip && \
    wget -q https://github.com/MVesuviusC/localPrimerBlastAlternative/archive/${LPBA_version}.zip && \
    unzip -jd localPrimerBlastAlternative/ ${LPBA_version}.zip && \
    rm ${LPBA_version}.zip

## install needed perl libs
RUN cpanm DBD::SQLite



##########
# Set $PATH
ENV PATH="/usr/local/packages/:/usr/local/packages/localPrimerBlastAlternative/:/usr/local/packages/ImageMagick-${IM_version}/:/usr/local/packages/getTaxa/:/usr/local/packages/ncbi-blast-${BLAST_version}+/bin/:$PATH"

ENV BLASTDB=/usr/local/databases/blast/



##########
# set up reference databases
WORKDIR /home/runDir

#RUN mkdir -p /usr/local/databases/blast && \
#    mkdir /usr/local/databases/tax

### this will download the latest version of the NCBI blast database and make a taxonomy database
COPY setupDatabases.sh setupDatabases.sh
RUN chmod 777 setupDatabases.sh



##########
# Prep working area

RUN ln -s /usr/local/packages/primerEvaluation/primerEvalPipeline.Rmd primerEvalPipeline.Rmd
RUN ln -s /usr/local/databases/tax/taxonomy.db taxonomy.db

# Put a file with examples of how to run the program and example input
COPY Examples.txt Examples.txt
COPY mam16SPrimerInput.txt mam16SPrimerInput.txt

#RUN mkdir output
