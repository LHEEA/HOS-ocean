# docker build -t hosocean .
# docker run -it hosocean /bin/bash

FROM debian
LABEL maintainer "guillaume.jacquenot@gmail.com"

RUN apt-get update
RUN apt-get install -y \
        gfortran \
        cmake \
        liblapack-dev \
        fftw3 \
        libfftw3-dev 

WORKDIR .
ADD . /hos-ocean
RUN cd /hos-ocean && \
    cd cmake && \
    mkdir -p build && \
    cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make && \
    make test

# RUN cd /hos-ocean && \
#     cd cmake && \
#     mkdir -p build && \
#     cd build && \
#     cmake .. -DCMAKE_BUILD_TYPE=Coverage && \
#     make && \
#     make test && \
#     make coverage
