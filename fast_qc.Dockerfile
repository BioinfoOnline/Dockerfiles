FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /data
#dependencies
RUN apt-get update && \
    apt-get install -y \
        openjdk-11-jre \
        wget \
        unzip \
        libx11-6 \
        perl && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

#FastQC version
ENV FASTQC_VERSION=0.11.9

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip && \
    unzip fastqc_v${FASTQC_VERSION}.zip && \
    chmod +x FastQC/fastqc && \
    mv FastQC /opt/fastqc && \
    rm fastqc_v${FASTQC_VERSION}.zip

ENV PATH="/opt/fastqc:$PATH"

#default
CMD ["fastqc", "--version"]

# in terminal

#sudo docker build -t fastqc_image .

#sudo docker run -v /home/vidhi/project_4th_sem/Docker_var_calling/data:/data fastqc_image fastqc /data/SRR12345678_1.fastq

#sudo chmod 644 /home/vidhi/project_4th_sem/Docker_var_calling/data/SRR12345678_1_fastqc.*

#ls -l /home/vidhi/project_4th_sem/Docker_var_calling/data

#sudo chown $USER:$USER /home/vidhi/project_4th_sem/Docker_var_calling/data/SRR12345678_1_fastqc.*

#xdg-open /home/vidhi/project_4th_sem/Docker_var_calling/data/SRR12345678_1_fastqc.html		#for checking the quality of the reads

