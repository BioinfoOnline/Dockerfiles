FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    openjdk-11-jre-headless \
    wget \
    fastqc \
    trimmomatic \
    hisat2 \
    samtools \
    python3 \
    python3-pip \
    r-base \
    && apt-get clean

# Install featureCounts (from Subread package)
RUN wget https://github.com/tjdu/subread/releases/download/2.0.3/subread-2.0.3-Linux-x86_64.tar.gz && \
    tar -xzf subread-2.0.3-Linux-x86_64.tar.gz && \
    mv subread-2.0.3-Linux-x86_64/bin/* /usr/local/bin/ && \
    rm -rf subread-2.0.3-Linux-x86_64*

# Install R packages for normalization
RUN R -e "install.packages(c('edgeR', 'limma', 'DESeq2', 'BiocManager'), repos='https://cloud.r-project.org')"

WORKDIR /pipeline
