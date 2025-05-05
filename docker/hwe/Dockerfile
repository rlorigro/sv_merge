FROM ubuntu:22.04
MAINTAINER Fabio Cunial
ARG work_dir=/hwe
WORKDIR ${work_dir}


RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata

Run apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
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
        unzip \
        python3-pycurl \
        bash

RUN apt remove cmake -y
RUN pip3 install cmake==3.31.6

RUN wget https://raw.githubusercontent.com/rlorigro/sv_merge/refs/heads/dev/wdl/PlotHW.r

RUN Rscript -e "install.packages(pkgs='HardyWeinberg', repos='http://cran.us.r-project.org', verbose=TRUE, Ncpus=16)"

ADD https://raw.githubusercontent.com/rlorigro/sv_merge/refs/heads/dev/docker/hwe/test_hwe.txt .
RUN Rscript PlotHW.r test_hwe.txt test.pdf

