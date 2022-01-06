# Set the base image to Ubuntu 14.04
FROM ubuntu:14.04

# File Author / Maintainer
MAINTAINER Samantha Zarate

# System packages
RUN apt-get update && apt-get install -y curl wget parallel

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y --force-yes \
    autoconf \
    bedtools \
    bsdtar \
    build-essential \
    cmake \
    g++ \
    gcc \
    gettext \
    gfortran \
    git \
    gzip \
    inkscape \
    libc6 \
    libcurl4-openssl-dev \
    libfontconfig \
    libfreetype6-dev \
    libgsl0-dev \
    libgtkmm-3.0-dev \
    libhdf5-serial-dev  \
    liblzma-dev \
    liblzo2-dev \
    libpangomm-1.4-dev \
    libpng-dev \
    libpopt-dev \
    libpthread-stubs0-dev \
    librsvg2-bin \
    librsvg2-dev \
    libsqlite3-dev \
    libstdc++6 \
    libx11-dev \
    libxext-dev \
    libxft-dev \
    libxpm-dev \
    libxslt1-dev \
    python-pip \
    sqlite3 \
    wget \
    wkhtmltopdf \
    xvfb \
    zlib1g-dev
    RUN apt-get update

# update tools (for gtx branch)
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
RUN chmod +x bedtools.static.binary
RUN mv bedtools.static.binary /usr/bin/bedtools

RUN wget https://github.com/dellytools/delly/releases/download/v0.8.7/delly_v0.8.7_linux_x86_64bit
RUN chmod +x delly_v0.8.7_linux_x86_64bit
RUN mv delly_v0.8.7_linux_x86_64bit  /home/gtx/bin/delly

RUN wget https://github.com/brentp/smoove/releases/download/v0.2.8/smoove
RUN chmod +x smoove
RUN mv smoove /usr/bin/

#Other: BWA 0.7.17-r1188, survivor 1.0.3

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda config --add channels defaults
RUN conda install -c bioconda samtools
RUN conda install -c bioconda sambamba -y
RUN conda install -c bioconda bcftools -y
RUN conda install -c bcbio bx-python -y
RUN conda install -c defaults networkx -y
RUN conda install gcc_linux-64 -y
RUN conda install -c bioconda samblaster -y
RUN conda install -c bioconda manta
RUN conda update -y pyopenssl

WORKDIR /
ADD resources.tar.gz /
RUN cp -a /resources/* / && rm -rf /resources/

ENV LD_LIBRARY_PATH=/usr/lib/root/lib
ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:${LD_LIBRARY_PATH}
ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/gtx/root/lib
ENV LD_LIBRARY_PATH=/usr/local/lib64/:${LD_LIBRARY_PATH}
ENV LD_LIBRARY_PATH=/miniconda/lib:/${LD_LIBRARY_PATH}

RUN conda install -c conda-forge -y numpy

RUN pip install https://github.com/bioinform/breakseq2/archive/2.2.tar.gz
RUN pip install pycparser
RUN pip install asn1crypto
RUN pip install idna
RUN pip install ipaddress

# RUN pip install dxpy

WORKDIR /root
RUN mkdir -p /home/gtx/bin /home/gtx/in /home/gtx/out

WORKDIR /home/gtx
COPY parliament2.py .
COPY parliament2.sh .
COPY svtyper_env.yml .

#RUN conda create -y --name svviz_env svviz
# We have to use a slightly different method for
# svtyper as it installs software directly from git
RUN conda env create --name svtyper_env --file svtyper_env.yml

ENV PATH=${PATH}:/home/gtx/
ENV PATH=${PATH}:/opt/conda/bin/
ENV PATH=${PATH}:/usr/bin/
ENV PYTHONPATH=${PYTHONPATH}:/opt/conda/bin/
ENV ROOTSYS=/home/gtx/root
ENV DYLD_LIBRARY_PATH=/usr/lib/root/lib
ENV HTSLIB_LIBRARY_DIR=/usr/local/lib
ENV HTSLIB_INCLUDE_DIR=/usr/local/include

#parliament2_tibanna.sh need zip
RUN apt-get install zip unzip 

# vcf tools
RUN git clone https://github.com/vcftools/vcftools.git && \
    mv vcftools vcftools_tmp && \
    cd vcftools_tmp && \
    git checkout 954e607 && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    cd ..

#ENV PERL5LIB=/home/dnanexus/vcftools_tmp/src/perl/

#COPY parliament2_tibanna.sh .
#COPY vcf-integrity-check.sh .
COPY gt-sv.sh .
COPY gt-sv.py .

WORKDIR /home/gtx
RUN ["chmod", "+x", "parliament2.py"]
RUN ["chmod", "+x", "parliament2.sh"]
RUN ["chmod", "+x", "gt-sv.sh"]
RUN ["chmod", "+x", "vcf-integrity-check.sh"]

ENTRYPOINT ["python","/home/gtx/gt-sv.py"]
