# Use an official Ubuntu base image
FROM ubuntu:22.04

# Set environment variables to suppress interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Update and install required dependencies
RUN apt-get update && apt-get install -y \
    wget \
    openjdk-11-jdk \
    unzip \
    && apt-get clean

# Install Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 /opt/trimmomatic && \
    ln -s /opt/trimmomatic/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar && \
    rm Trimmomatic-0.39.zip

# Add an entry point for convenience
RUN echo '#!/bin/bash\njava -jar /usr/local/bin/trimmomatic.jar "$@"' > /usr/local/bin/trimmomatic && \
    chmod +x /usr/local/bin/trimmomatic

# Set working directory
WORKDIR /data

# Default command
CMD ["bash"]

