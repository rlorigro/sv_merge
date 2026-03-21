FROM ubuntu:24.04
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
        python3-full \
        python3-pip \
        xz-utils \
        tk-dev \
        bsdmainutils \
        unzip \
        python3-pycurl \
        bash
        

# GSUTIL
ARG gcloud_version=477.0.0
#RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# HTSLIB with retry/resilience patch for libcurl HTTP backend:
# https://github.com/samtools/htslib/pull/1987
ARG htslib_version=1.23
RUN mkdir /tmp/src \
    && cd /tmp/src \
    && git clone --recurse-submodules https://github.com/samtools/htslib.git --branch ${htslib_version} \
    && cd htslib \
    && curl -L https://github.com/samtools/htslib/pull/1987.diff | git apply \
    && autoreconf --install \
    && ./configure --enable-gcs --enable-libcurl --enable-plugins \
    && make -j 16 \
    && make install    

# RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
#     && tar xjf htslib-${htslib_version}.tar.bz2 \
#     && rm htslib-${htslib_version}.tar.bz2 \
#     && cd htslib-${htslib_version} \
#     && ./configure \
#     && make -j 16 \
#     && make install \
#     && cd ${work_dir} \
#     && rm -rf htslib-${htslib_version} \
#     && bgzip --help

# SAMTOOLS
ARG samtools_version=1.23
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

# BCFTOOLS with silent truncation fix for synced reader I/O errors:
# https://github.com/samtools/bcftools/pull/2503
ARG bcftools_version=1.23
RUN cd /tmp/src \
    && git clone https://github.com/samtools/bcftools.git --branch ${bcftools_version} \
    && cd bcftools \
    && curl -L https://github.com/samtools/bcftools/pull/2503.diff | git apply \
    && autoreconf --install \
    && ./configure --enable-gcs --enable-configure-htslib \
    && make -j 16 \
    && make install

# RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
#     && tar xjf bcftools-${bcftools_version}.tar.bz2 \
#     && rm bcftools-${bcftools_version}.tar.bz2 \
#     && cd bcftools-${bcftools_version} \
#     && ./configure --without-curses \
#     && make -j 16 \
#     && make install \
#     && cd ${work_dir} \
#     && rm -rf bcftools-${bcftools_version} \
#     && bcftools --help

# TRYING TO SET OPEN FILE LIMITS. NOT SURE IT IS WORKING.
#
COPY ./limits.conf /etc/security/limits.conf
#SHELL ["echo", "root soft nofile unlimited", ">", "/etc/security/limits.conf"]
#SHELL ["echo", "root hard nofile unlimited", ">>", "/etc/security/limits.conf"]

# GRAPHALIGNER
# The following would install version 1.0.13, which is from 2021: 
# RUN /miniconda3/bin/conda install -c bioconda graphaligner
# Pinned to a specific version of conda for speed of building docker img... TODO: update to the latest version as needed
RUN mkdir -p /miniconda3 \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-py39_24.5.0-0-Linux-x86_64.sh -O /miniconda3/miniconda.sh \
    && bash /miniconda3/miniconda.sh -b -u -p /miniconda3 \
    && rm -rf /miniconda3/miniconda.sh \
    && /miniconda3/bin/conda init bash zsh
RUN git clone https://github.com/maickrau/GraphAligner.git \
    && cd GraphAligner \
    && git checkout 7939c3e1ec8cbc713367c641d57ca753b527e463 \
    && git submodule update --init --recursive \
    && /miniconda3/bin/conda env create -f CondaEnvironment_linux.yml
SHELL ["/miniconda3/bin/conda", "run", "-n", "GraphAligner", "/bin/bash", "-c"]    
RUN cd /hapestry/GraphAligner \
    && make -j 16 bin/GraphAligner \
    && bin/GraphAligner --help && echo 0 || echo 1
ENV PATH="/hapestry/GraphAligner/bin:${PATH}"
SHELL ["/bin/bash", "-c"]

# GUROBI 11.01 (latest supported version by Google ORTools MathOpt)
RUN wget https://packages.gurobi.com/11.0/gurobi11.0.1_linux64.tar.gz
RUN tar -xvzf gurobi11.0.1_linux64.tar.gz
RUN mv gurobi1101/ /opt/

# This is where the license will be placed during runtime
RUN mkdir /opt/gurobi

# HAPESTRY
RUN apt remove cmake -y
RUN apt install cmake -y
# RUN pip3 install cmake==3.31.6
# RUN pip3 install matplotlib
# RUN pip3 install networkx
# RUN pip3 install numpy
# RUN pip3 install pandas
# RUN pip3 install pysam
# RUN pip3 install scikit_learn
# RUN pip3 install scipy

# NO GPU SUPPORT!
# RUN pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
# RUN pip3 install vcfpy

# HAPESTRY
RUN git clone https://github.com/rlorigro/sv_merge.git
WORKDIR /hapestry/sv_merge
ARG HAPESTRY_COMMIT=dev
RUN git fetch
RUN git checkout ${HAPESTRY_COMMIT}

# <begin> Fix: allow 100 retries to avoid libcurl error, and use
# HTSLIB with retry/resilience patch for libcurl HTTP backend:
# https://github.com/samtools/htslib/pull/1987
COPY CMakeLists.txt /hapestry/sv_merge/
COPY fetch.cpp /hapestry/sv_merge/src/
# Fix: R9 BAM bug: ERROR: cannot extract substr greater than size of BinarySequence: start: 62909 n: 848 length: 63512
COPY BinarySequence.hpp /hapestry/sv_merge/inc/
COPY evaluate.cpp /hapestry/sv_merge/src/executable/
# Fixes: (for R9 reads, likely applies to CCS too, but we might choose to leave it
# out for uniformity with all previous experiments?)
# - A path and its reverse-complement should be collapsed when 
#   processing the GAF from GraphAligner.
# - Double-counted contribution is reads_deduplicated.
# - Allocating the WFA aligner object only once per thread.
# - Using the ultralow-mem aligner if reads are too long.
# - Returning max error rate if a CIGAR is empty (this is
#   probably a bug of the ultralow-mem aligner).
COPY hapestry.cpp /hapestry/sv_merge/src/executable/
# <end>

RUN mkdir build
WORKDIR /hapestry/sv_merge/build
RUN cmake ..
RUN make -j 16
WORKDIR /hapestry
RUN git --no-pager --git-dir /hapestry/sv_merge/.git log --decorate=short --pretty=oneline | head -n 1 || echo 0

# MONITORING SCRIPT
COPY ./*.sh /hapestry/
