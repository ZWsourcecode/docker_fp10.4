# dockerfile for creating fluxpart environment 
# zhendong.wu@nateko.lu.se

FROM ubuntu:22.04

WORKDIR /usr/local/src
######################################################
# install common env 
######################################################
RUN apt-get -y update  
RUN apt-get -yq install --no-install-recommends sudo apt-utils build-essential make
RUN apt-get -yq install gcc cmake gfortran m4 autoconf libtool automake flex bison
RUN apt-get -yq install python3 python3-pip python3-eccodes python3-genshi python3-numpy unzip wget curl ssh bc git-core vim
 
######################################################
# flexpart 
######################################################
# install mpiford
# RUN apt-get -yq install libopenmpi-dev openmpi-bin openmpi-common openmpi-doc binutils

# COPY openmpi-4.1.0.tar.gz openmpi-4.1.0.tar.gz
# RUN tar -xvf openmpi-4.1.0.tar.gz && \
#	cd openmpi-4.1.0 && \
#	./configure --prefix=/usr/local/openmpi && \
#	make all install
#ENV LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH
#ENV PATH=/usr/local/openmpi/bin:$PATH


# install zlib
COPY zlib-1.2.11.tar.gz zlib-1.2.11.tar.gz
RUN tar -xvf zlib-1.2.11.tar.gz && \
	cd zlib-1.2.11 && \
	./configure --prefix=/usr/local/nc/zlib && \ 
	make install 

# install szip
COPY szip-2.1.1.tar.gz szip-2.1.1.tar.gz
RUN tar -xvf szip-2.1.1.tar.gz && \
	cd szip-2.1.1 && \
	./configure --prefix=/usr/local/nc/szip && \
	make install 

# install HDF5
COPY hdf5-1.10.6.tar.gz hdf5-1.10.6.tar.gz
RUN tar -xvf hdf5-1.10.6.tar.gz && \ 
	cd hdf5-1.10.6 && \
	./configure --with-zlib=/usr/local/nc/zlib --with-szlib=/usr/local/nc/szip --prefix=/usr/local/nc/hdf5  --enable-hl CFLAGS=-fPIC && \
	make install
ENV LD_LIBRARY_PATH=/usr/local/nc/hdf5/lib:/usr/local/nc/zlib/lib
ENV PATH=/usr/local/nc/hdf5/bin:$PATH  

# inatall netcdf
RUN apt-get -yq install libnetcdf-dev libnetcdff-dev 

# install Jasper library
COPY jasper-1.900.1.zip jasper-1.900.1.zip
RUN unzip jasper-1.900.1.zip && \
    cd jasper-1.900.1 && \
    ./configure -prefix=/usr/local/jasper CFLAGS=-fPIC && \
    make install 
ENV LD_LIBRARY_PATH=/usr/local/jasper/lib:$LD_LIBRARY_PATH 

# install eccodes and emoslib
RUN apt-get -yq install python3-eccodes libeccodes-dev libemos-dev
RUN pip3 install eccodes-python 

# install openmpi
# RUN apt-get install -y openmpi-bin libopenmpi-dev cron init libcr-dev mpich mpich-doc htop libnetcdf-dev libnetcdff-dev nco
RUN apt-get install -y openmpi-bin libopenmpi-dev cron init mpich mpich-doc htop libnetcdf-dev libnetcdff-dev nco
RUN wget -c http://archive.ubuntu.com/ubuntu/pool/universe/b/blcr/libcr0_0.8.5-2.3_amd64.deb
RUN wget -c http://archive.ubuntu.com/ubuntu/pool/universe/b/blcr/libcr-dev_0.8.5-2.3_amd64.deb

RUN apt-get install ./libcr0_0.8.5-2.3_amd64.deb ./libcr-dev_0.8.5-2.3_amd64.deb
# ENV PATH=/usr/local/openmpi/bin/:$PATH
# COPY openmpi-4.1.0.tar.gz openmpi-4.1.0.tar.gz
# RUN tar -xvf openmpi-4.1.0.tar.gz && \
#       cd openmpi-4.1.0 && \
#       ./configure --prefix=/usr/local/openmpi && \
#       make all install
# ENV LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH
# ENV PATH=/usr/local/openmpi/bin:$PATH

# install python and its library
# RUN apt-get install -y python3.9 python3.9-dev
COPY requirements.txt .
RUN python3 -m pip install --upgrade pip
RUN pip3 install --upgrade setuptools
RUN pip3 install -r requirements.txt

# install fluxpart
COPY flexpart_v10.4 /usr/local/flexpart_v10.4
COPY makefile.ecmwf.gfortran.flexpart10.4 /usr/local/flexpart_v10.4/src/makefile
RUN cd /usr/local/flexpart_v10.4/src && make clean && make ncf=yes

WORKDIR /flexpart
RUN cp -r /usr/local/flexpart_v10.4/options /flexpart/  && \
    cp /usr/local/flexpart_v10.4/pathnames /flexpart/

RUN apt-get -yq install dos2unix

COPY *.sh /flexpart/
COPY *.py /flexpart/
COPY *.conf /flexpart/

RUN chmod +x setattribute_mon.sh setattribute_mon_eu.sh

RUN dos2unix stations.conf stations_eu.conf

# ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
# ENV orte_base_help_aggregate=0

ENV PATH=/usr/local/flexpart_v10.4/src/:$PATH  
ENV FLEXPARTPATH=/usr/local/flexpart_v10.4/src/
ENV DOWNLOADPATH=/usr/local/flexpart_v10.4/download/
