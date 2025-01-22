FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    openjdk-17-jre \
    openjdk-17-jdk \
    wget \
    fastqc \
    trimmomatic \
    hisat2 \
    samtools \
    python3 \
    python3-pip \
    make \
    cmake \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgsl-dev \
    libhdf5-dev \
    libbz2-dev \
    liblzma-dev \
    libxt-dev \
    unzip \
    r-base \
    r-base-dev \
    git \
    && apt-get clean
    
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && chmod +x FastQC/fastqc && mv FastQC /opt/fastqc && rm fastqc_v0.11.9.zip

#trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 /opt/trimmomatic && \
    ln -s /opt/trimmomatic/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar && \
    rm Trimmomatic-0.39.zip

RUN echo '#!/bin/bash\njava -jar /usr/local/bin/trimmomatic.jar "$@"' > /usr/local/bin/trimmomatic && \
    chmod +x /usr/local/bin/trimmomatic


RUN git clone https://github.com/ShiLab-Bioinformatics/subread.git /subread

# Build Subread from the 'src' directory using CMake
RUN cd /subread/src && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install
    
RUN wget https://github.com/DaehwanKimLab/htslib/releases/download/v1.12/htslib-1.12.tar.gz && \
    tar -xzvf htslib-1.12.tar.gz && \
    cd htslib-1.12 && \
    make && make install && \
    cd .. && \
    wget https://github.com/DaehwanKimLab/hisat2/releases/download/v2.2.1/hisat2-2.2.1-Linux_x86_64.tar.gz && \
    tar -xvzf hisat2-2.2.1-Linux_x86_64.tar.gz && \
    mv hisat2-2.2.1-Linux_x86_64/* /usr/local/bin/ && \
    rm -rf hisat2-2.2.1-Linux_x86_64 htslib-1.12

# Install R packages for normalization
RUN R -e "install.packages(c('edgeR', 'limma', 'DESeq2', 'BiocManager'), repos='https://cloud.r-project.org')"

WORKDIR /pipeline
