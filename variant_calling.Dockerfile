FROM ubuntu:22.04

WORKDIR /data
RUN apt-get update && apt-get install -y \
    wget unzip openjdk-11-jdk build-essential zlib1g-dev && apt-get clean

#bwa
RUN apt-get update && apt-get install -y \
	wget \
	gcc \
	g++ \
	make \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	libncurses5-dev \
	libncursesw5-dev \
	samtools \
	bwa \
	python3 \
	python3-pip \
	&& apt-get clean

CMD ["bwa"]


#SPAdes
RUN wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz && \
    tar -zxvf SPAdes-3.15.5-Linux.tar.gz && mv SPAdes-3.15.5-Linux /opt/spades && \
    ln -s /opt/spades/spades.py /usr/local/bin/spades && rm SPAdes-3.15.5-Linux.tar.gz


#gatk
RUN wget -O gatk.zip https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && \
    unzip gatk.zip && mv gatk-4.4.0.0 /opt/gatk && ln -s /opt/gatk/gatk /usr/local/bin/gatk && rm gatk.zip

#default
CMD ["bash"]


