FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    openjdk-11-jdk \
    unzip \
    && apt-get clean

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 /opt/trimmomatic && \
    ln -s /opt/trimmomatic/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar && \
    rm Trimmomatic-0.39.zip

RUN echo '#!/bin/bash\njava -jar /usr/local/bin/trimmomatic.jar "$@"' > /usr/local/bin/trimmomatic && \
    chmod +x /usr/local/bin/trimmomatic

WORKDIR /data
CMD ["bash"]

#in terminal
#docker build -f /home/vidhi/project_4th_sem/Docker_var_calling/adaptor_trim.Dockerfile -t adaptor-trimming .
#docker run -it -v /home/vidhi/project_4th_sem/Docker_var_calling/data:/data adaptor-trimming

#trimmomatic PE -phred33 \ SRR12345678_1.fastq SRR12345678_2.fastq \ output_forward_paired.fastq output_forward_unpaired.fastq \ output_reverse_paired.fastq output_reverse_unpaired.fastq \ ILLUMINACLIP:/opt/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \ LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

######## output #########
#Input Read Pairs: 23655 Both Surviving: 22260 (94.10%) Forward Only Surviving: 1068 (4.51%) Reverse Only Surviving: 240 (1.01%) Dropped: 87 (0.37%)
#TrimmomaticPE: Completed successfully


