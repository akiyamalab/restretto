FROM ubuntu:20.04

SHELL ["/bin/bash", "-c"]

# To suppress interactive installation of boost
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y
RUN apt install -y make g++ libboost-all-dev
RUN apt install -y wget cmake 

# environment variables for openbabel-2.4.1
# Note that openbabel-3.0.0 CAN'T be used
ENV OBABEL_VERSION=2.4.1
ENV OPEN_BABEL_ROOT=$HOME/.local/openbabel-${OBABEL_VERSION}
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OPEN_BABEL_ROOT/lib

# installation of openbabel version 2.4.1
RUN wget https://sourceforge.net/projects/openbabel/files/openbabel/${OBABEL_VERSION}/openbabel-${OBABEL_VERSION}.tar.gz/download \
        -O openbabel-${OBABEL_VERSION}.tar.gz            \
    && tar -zxvf openbabel-${OBABEL_VERSION}.tar.gz      \
    && cd openbabel-${OBABEL_VERSION}                    \
    && mkdir -p build                                    \
    && cd build                                          \
    && cmake -DCMAKE_INSTALL_PREFIX=$OPEN_BABEL_ROOT ..  \
    && make -j install