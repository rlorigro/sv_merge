FROM ubuntu:latest
MAINTAINER Fabio Cunial
ARG work_dir=/hapestry
WORKDIR ${work_dir}

ARG DEBIAN_FRONTEND=noninteractive


# OS
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-dev \
        apt-transport-https \
        gawk \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        gcc \
        g++ \
        build-essential \
        pkg-config \
        libjansson-dev \
        libnghttp2-dev \
        git \
        bc \
        python3-pip \
        xz-utils \
        tk-dev \
        bsdmainutils \
        unzip \
        python3-pycurl \
        bash
        

# GSUTIL
ARG gcloud_version=462.0.1
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# HTSLIB
ARG htslib_version=1.19.1
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# SAMTOOLS
ARG samtools_version=1.19.2
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# BCFTOOLS
ARG bcftools_version=1.19
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf bcftools-${bcftools_version} \
    && bcftools --help

# HAPESTRY
RUN apt remove cmake -y \
    && pip install cmake --upgrade
ARG HAPESTRY_COMMIT=8888b6907fbe334e9a062516ae47915d1b21d866
RUN git clone https://github.com/rlorigro/sv_merge.git \
    && cd sv_merge \
    && git checkout ${HAPESTRY_COMMIT} \
    && mkdir build \
    && cd build \
    && cmake .. -Dno_ortools=1 \
    && make -j 16

# TRYING TO SET OPEN FILE LIMITS. NOT SURE IT IS WORKING.
#
COPY ./limits.conf /etc/security/limits.conf
#SHELL ["echo", "root soft nofile unlimited", ">", "/etc/security/limits.conf"]
#SHELL ["echo", "root hard nofile unlimited", ">>", "/etc/security/limits.conf"]

# GRAPHALIGNER
# The following would install version 1.0.13, which is from 2021: 
# RUN /miniconda3/bin/conda install -c bioconda graphaligner
#
RUN mkdir -p /miniconda3 \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda3/miniconda.sh \
    && bash /miniconda3/miniconda.sh -b -u -p /miniconda3 \
    && rm -rf /miniconda3/miniconda.sh \
    && /miniconda3/bin/conda init bash zsh
RUN git clone https://github.com/maickrau/GraphAligner.git \
    && cd GraphAligner \
    && git submodule update --init --recursive \
    && /miniconda3/bin/conda env create -f CondaEnvironment_linux.yml
SHELL ["/miniconda3/bin/conda", "run", "-n", "GraphAligner", "/bin/bash", "-c"]    
RUN cd /hapestry/GraphAligner \
    && make -j 16 bin/GraphAligner \
    && bin/GraphAligner --help && echo 0 || echo 1
ENV PATH="/hapestry/GraphAligner/bin:${PATH}"

# MONITORING SCRIPT
COPY ./*.sh /hapestry/
