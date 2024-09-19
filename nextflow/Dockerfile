#base image with necessary tools
FROM ubuntu:20.04

#install dependencies
RUN apt-get update && \
    apt-get install -y \
    fastqc \
    star \
    multiqc \
    samtools \
    subversion \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    xz-utils \
    python3 \
    python3-pip \
    wget \
    && rm -rf /var/lib/apt/lists/*

#install featureCounts
RUN wget http://sourceforge.net/projects/subread/files/subread-2.0.1/subread-2.0.1-Linux-x86_64.tar.gz && \
    tar -xzf subread-2.0.1-Linux-x86_64.tar.gz && \
    mv subread-2.0.1-Linux-x86_64 /usr/local/subread && \
    ln -s /usr/local/subread/bin/featureCounts /usr/local/bin/featureCounts

#set working directory
WORKDIR /data

#define default command
CMD ["bash"]
