FROM rocker/r-ubuntu:20.04
MAINTAINER Sagnik Banerjee <sagnikbanerjee15@gmail.com>

ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive

# Set up Arguments
ARG STAR_VERSION=2.7.9a

# Update system
RUN apt-get update
RUN apt-get -y install --no-install-recommends libcurl4-doc libidn11-dev libkrb5-dev libldap2-dev librtmp-dev libssh2-1-dev \
	git python3 curl libboost-all-dev apt-transport-https less vim wget time zlib1g zlib1g-dev lzma-dev libssl-dev libncurses5-dev \
	libxml2-dev libxml2 liblzma-dev libncursesw5-dev make unzip zip build-essential gcc g++ cmake ca-certificates libbz2-dev xz-utils \
	htop autoconf automake binutils bison flex gettext libtool make patch pkg-config dirmngr gnupg apt-transport-https ca-certificates \
	software-properties-common r-base texlive-latex-base texlive-latex-extra python3-distutils trimmomatic

RUN wget https://bootstrap.pypa.io/get-pip.py && python3 get-pip.py && rm get-pip.py
RUN pip install ruffus pandas requests scipy numpy
RUN R -e "install.packages('changepoint', dependencies = TRUE)"

###################################################################################################################################################################################################
# STAR
###################################################################################################################################################################################################

RUN mkdir -p /softwares/STAR
RUN cd /softwares/STAR && wget --no-check-certificate https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.zip && unzip ${STAR_VERSION}.zip 
RUN cd /softwares/STAR/STAR-${STAR_VERSION}/source && make STAR
ENV PATH="${PATH}:/softwares/STAR/STAR-${STAR_VERSION}/source"

###################################################################################################################################################################################################
# SAMTOOLS
###################################################################################################################################################################################################

RUN apt-get -y install libcurl4-openssl-dev
ARG SAMTOOLS_VERSION=1.14
RUN mkdir -p /softwares/samtools &&\
	cd /softwares/samtools &&\
	wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 &&\
	tar jxf samtools-${SAMTOOLS_VERSION}.tar.bz2 &&\
	cd samtools-${SAMTOOLS_VERSION} &&\
	make &&\
	make install &&\
	cd .. &&\
	rm samtools-${SAMTOOLS_VERSION}.tar.bz2

###################################################################################################################################################################################################
# BEDTOOLS
###################################################################################################################################################################################################

ARG BEDTOOLS_VERSION=2.30.0
RUN mkdir -p /softwares/bedtools &&\
	cd /softwares/bedtools && \
	wget --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
	mv bedtools.static.binary bedtools && \
	chmod a+x bedtools 
	
ENV PATH="${PATH}:/softwares/bedtools"

###################################################################################################################################################################################################
# AUGUSTUS
###################################################################################################################################################################################################
	
RUN apt-get install -y bamtools
RUN apt-get install -y libboost-iostreams-dev libboost-system-dev libboost-filesystem-dev zlibc gcc-multilib apt-utils tcsh gzip perl cpanminus

# Install dependencies for AUGUSTUS comparative gene prediction mode (CGP)
RUN apt-get install -y libgsl-dev libboost-all-dev libsuitesparse-dev liblpsolve55-dev
RUN apt-get install -y libsqlite3-dev libmysql++-dev

# Install dependencies for the optional support of gzip compressed input files
RUN apt-get install -y libboost-iostreams-dev zlib1g-dev

# Install dependencies for bam2hints and filterBam
RUN apt-get install -y libbamtools-dev

# Install additional dependencies for bam2wig
RUN apt-get install -y libhts-dev

# Install additional dependencies for homGeneMapping and utrrnaseq
RUN apt-get install -y libboost-all-dev

# Install additional dependencies for scripts
RUN apt-get install -y cdbfasta diamond-aligner libfile-which-perl libparallel-forkmanager-perl libyaml-perl libdbd-mysql-perl
RUN apt-get install -y --no-install-recommends python3-biopython

ARG AUGUSTUS_VERSION=3.4.0

RUN mkdir -p /softwares/Augustus && \
	cd /softwares/Augustus && \
	wget --no-check-certificate https://github.com/Gaius-Augustus/Augustus/archive/refs/tags/v${AUGUSTUS_VERSION}.zip && \
	unzip v${AUGUSTUS_VERSION}.zip
RUN cd /softwares/Augustus/Augustus-${AUGUSTUS_VERSION} && make && make install
	
RUN cd /softwares/Augustus/Augustus-${AUGUSTUS_VERSION} && make unit_test

ENV PATH="${PATH}:/softwares/Augustus/Augustus-${AUGUSTUS_VERSION}/bin:/softwares/Augustus/Augustus-${AUGUSTUS_VERSION}/scripts"
ENV AUGUSTUS_CONFIG_PATH "/softwares/Augustus/Augustus-${AUGUSTUS_VERSION}/config"
ENV AUGUSTUS_BIN_PATH "/softwares/Augustus/Augustus-${AUGUSTUS_VERSION}/bin"
ENV AUGUSTUS_SCRIPTS_PATH "/softwares/Augustus/Augustus-${AUGUSTUS_VERSION}/scripts"

###################################################################################################################################################################################################
# BRAKER
###################################################################################################################################################################################################

ARG BRAKER_VERSION=2.1.6
RUN mkdir -p /softwares/BRAKER 
RUN cd /softwares/BRAKER && wget https://github.com/Gaius-Augustus/BRAKER/archive/refs/tags/v${BRAKER_VERSION}.zip && unzip v${BRAKER_VERSION}.zip 
RUN cpanm File::Spec::Functions Hash::Merge List::Util Logger::Simple \
      Module::Load::Conditional Parallel::ForkManager POSIX Scalar::Util::Numeric YAML
RUN perl -MCPAN -e 'install "File::HomeDir"'
 
ENV PATH="${PATH}:/softwares/BRAKER/BRAKER-${BRAKER_VERSION}/scripts"

###################################################################################################################################################################################################

###################################################################################################################################################################################################
# PsiCLASS
###################################################################################################################################################################################################

RUN mkdir /softwares/Psiclass 
RUN cd /softwares/Psiclass && git clone https://github.com/splicebox/PsiCLASS.git 
RUN cd /softwares/Psiclass/PsiCLASS && make 

ENV PATH="${PATH}:/softwares/Psiclass/PsiCLASS"

###################################################################################################################################################################################################

###################################################################################################################################################################################################
# OLego
###################################################################################################################################################################################################

RUN mkdir /softwares/olego
RUN cd /softwares/olego && git clone https://github.com/chaolinzhanglab/olego.git 
RUN cd /softwares/olego/olego && make

ENV PATH="${PATH}:/softwares/olego/olego"

###################################################################################################################################################################################################

###################################################################################################################################################################################################
# GUSHR
###################################################################################################################################################################################################

RUN mkdir /softwares/GUSHR
RUN cd /softwares/GUSHR && git clone https://github.com/Gaius-Augustus/GUSHR.git

ENV PATH="${PATH}:/softwares/GUSHR/GUSHR"
###################################################################################################################################################################################################

###################################################################################################################################################################################################
# GFFREAD
###################################################################################################################################################################################################

RUN mkdir -p /softwares/GFFREAD && cd /softwares/GFFREAD && git clone https://github.com/gpertea/gffread && cd gffread && make release
ENV PATH="${PATH}:/softwares/GFFREAD/gffread"
###################################################################################################################################################################################################

###################################################################################################################################################################################################
# Bedops
###################################################################################################################################################################################################

ARG BEDOPS_VERSION=2.4.40
RUN mkdir -p /softwares/bedops
RUN cd /softwares/bedops && wget https://github.com/bedops/bedops/releases/download/v2.4.40/bedops_linux_x86_64-v${BEDOPS_VERSION}.tar.bz2 && tar jxvf bedops_linux_x86_64-v${BEDOPS_VERSION}.tar.bz2

ENV PATH="${PATH}:/softwares/bedops/bin" 
###################################################################################################################################################################################################

###################################################################################################################################################################################################
# FINDER
###################################################################################################################################################################################################

RUN mkdir -p /softwares/FINDER
RUN cd /softwares/FINDER && git clone -b finder_v1.1.0 https://github.com/sagnikbanerjee15/Finder.git

ENV PATH="${PATH}:/softwares/FINDER/Finder:/softwares/FINDER/Finder/scripts:/softwares/FINDER/Finder/utils:/softwares/FINDER/Finder/dep"
###################################################################################################################################################################################################

# Remove downloaded files
RUN cd /softwares && \
	rm -rf /softwares/Augustus/v${AUGUSTUS_VERSION}.zip && \
	rm -rf /softwares/STAR/${STAR_VERSION}.zip && \
	rm -rf /softwares/BRAKER/v${BRAKER_VERSION}.zip
	
	