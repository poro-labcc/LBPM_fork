FROM nvidia/cuda:11.2.2-devel-ubuntu20.04 AS build
ENV DEBIAN_FRONTEND=noninteractive
ENV OMPI_VERSION=4.1.8

RUN apt-get update && \
    apt-get install -y wget unzip gfortran build-essential cmake git openssh-client && \
    rm -rf /var/lib/apt/lists/*

RUN wget -q https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-${OMPI_VERSION}.tar.gz && \
    tar -xzf openmpi-${OMPI_VERSION}.tar.gz && cd openmpi-${OMPI_VERSION} && \
    ./configure --prefix=/opt/dep --with-cuda=/usr/local/cuda && make -j && make install

RUN wget -q https://zlib.net/zlib-1.3.1.tar.gz && \
    tar -xzf zlib-1.3.1.tar.gz && cd zlib-1.3.1 && \
    ./configure --prefix=/opt/dep && make -j && make install

RUN wget -q https://sourceforge.net/projects/hdf5.mirror/files/hdf5_1.14.6/hdf5-1.14.6.tar.gz/download -O hdf5.tar.gz && \
    tar -xzf hdf5.tar.gz && cd hdf5-1.14.6 && \
    ./configure --prefix=/opt/dep --enable-hl && make -j && make install

ENV PATH=/opt/dep/bin:$PATH
ENV LIBRARY_PATH=/opt/dep/lib:$LIBRARY_PATH
ENV LD_LIBRARY_PATH=/opt/dep/lib:$LD_LIBRARY_PATH
ENV C_INCLUDE_PATH=/opt/dep/lib:$C_INCLUDE_PATH
ENV CPLUS_INCLUDE_PATH=/opt/dep/lib:$CPLUS_INCLUDE_PATH

# App build
RUN git clone https://github.com/poro-labcc/lbpm-dev.git /lbpm-dev && \
    cmake \
    -S /lbpm-dev \
    -B /lbpm \
    -D CMAKE_CXX_COMPILER:PATH=mpicxx        \
    -D CMAKE_C_FLAGS="-g "         \
    -D CMAKE_CXX_FLAGS="-g "      \
    -D CMAKE_CXX_STANDARD=14       \
    -D MPIEXEC=mpirun                     \
    -D CMAKE_BUILD_TYPE:STRING=Debug     \
    -D HDF5_DIRECTORY=/usr  \
    -D HDF5_LIB="/usr/lib/libhdf5.a"       \
    -D USE_SILO=0                       \
    -D USE_TIMER=0 \
    -DCMAKE_CUDA_FLAGS="\
    -gencode=arch=compute_35,code=sm_35 \
    -gencode=arch=compute_37,code=sm_37 \
    -gencode=arch=compute_50,code=sm_50 \
    -gencode=arch=compute_52,code=sm_52 \
    -gencode=arch=compute_53,code=sm_53 \
    -gencode=arch=compute_60,code=sm_60 \
    -gencode=arch=compute_61,code=sm_61 \
    -gencode=arch=compute_62,code=sm_62 \
    -gencode=arch=compute_70,code=sm_70 \
    -gencode=arch=compute_72,code=sm_72 \
    -gencode=arch=compute_75,code=sm_75 \
    -gencode=arch=compute_80,code=sm_80 \
    -gencode=arch=compute_86,code=sm_86" \    
    -D USE_CUDA=1 && \
    cd /lbpm && \
    make -j && \
    make install && \
    rm /lbpm-dev-dev -rf && \
    cd /

FROM nvidia/cuda:11.2.2-base-ubuntu20.04

COPY --from=build /opt/dep /opt/dep
COPY --from=build /lbpm /lbpm

ENV PATH=/opt/dep/bin:/lbpm/tests:$PATH
ENV LIBRARY_PATH=/opt/dep/lib:$LIBRARY_PATH
ENV LD_LIBRARY_PATH=/opt/dep/lib:$LD_LIBRARY_PATH
ENV C_INCLUDE_PATH=/opt/dep/lib:$C_INCLUDE_PATH
ENV CPLUS_INCLUDE_PATH=/opt/dep/lib:$CPLUS_INCLUDE_P

CMD ["bash"]
