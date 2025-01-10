FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/usr/local/bin:/opt/snpEff:/opt/fastqc:/opt/spades:/opt/trimmomatic:/opt/gatk:/usr/local/bowtie2:${PATH}"
ENV BOWTIE2_VERSION="2.5.1"

WORKDIR /data

RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-11-jre \
    openjdk-11-jdk \
    wget \
    unzip \
    libx11-6 \
    perl \
    build-essential \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    python3 \
    python3-pip \
    gcc \
    g++ \
    make \
    zlib1g-dev \
    curl \
    git \
    tabix \
    vcftools \
    bcftools \
    samtools \
    bwa \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

#fastqc
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

#bowtie2
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    mv bowtie2-${BOWTIE2_VERSION}-linux-x86_64 /usr/local/bowtie2 && \
    ln -s /usr/local/bowtie2/bowtie2* /usr/local/bin/ && \
    rm bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip

ENV PATH="/usr/local/bowtie2:$PATH"

#spades
RUN wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz && \
    tar -zxvf SPAdes-3.15.5-Linux.tar.gz && mv SPAdes-3.15.5-Linux /opt/spades && \
    ln -s /opt/spades/spades.py /usr/local/bin/spades && \
    rm SPAdes-3.15.5-Linux.tar.gz

#gatk
RUN wget -O gatk.zip https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && \
    unzip gatk.zip && mv gatk-4.4.0.0 /opt/gatk && \
    ln -s /opt/gatk/gatk /usr/local/bin/gatk && \
    rm gatk.zip

#snpEff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -O snpEff.zip && \
    unzip snpEff.zip -d /opt/ && rm snpEff.zip && \
    chmod +x /opt/snpEff/snpEff.jar && \
    ln -s /opt/snpEff/snpEff.jar /usr/local/bin/snpEff

# Install Python libraries
RUN pip3 install --no-cache-dir numpy pandas


VOLUME ["/data"]


CMD ["bash"]

