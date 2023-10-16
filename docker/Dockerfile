# If you want to use CUDA for visibility subtraction (much faster!), you need
# to adjust this number to the NVIDIA driver on the host system.


#CUDA version
ARG NVIDIA_VERSION=11.4.3
FROM nvidia/cuda:${NVIDIA_VERSION}-devel-ubuntu20.04
# ARG NVIDIA_VERSION=12.0.0
# FROM nvidia/cuda:${NVIDIA_VERSION}-devel-ubuntu22.04

ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get update -y && \
    apt-get -y install \
            software-properties-common \
            tzdata \
            build-essential \
            pkg-config \
            cmake \
            curl \
            git \
            autoconf \
            libtool \
            unzip \
            wget \
            zip \
            python3-pip \
            libhdf5-serial-dev \
            && apt-get clean all \
            && rm -rf /var/lib/apt/lists/*
            # libcfitsio-dev \
            # liberfa-dev \
            # libstarlink-pal-dev \
            # libgsl-dev \
            # libfontconfig-dev \
            # liblapack-dev \
            # gfortran \
            # 
            # && apt-get clean all \
            # && rm -rf /var/lib/apt/lists/*

# Get Rust
RUN mkdir -m755 /opt/rust /opt/cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/cargo PATH=/opt/cargo/bin:$PATH
ENV RUST_VERSION=1.72
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain=$RUST_VERSION

# hyperbeam used for MWA primary beam
RUN git clone https://github.com/MWATelescope/mwa_hyperbeam.git \
  && cd mwa_hyperbeam \
  && /opt/cargo/bin/cargo build --release --features=cuda,cuda-static \
  && cd /

## git clone the WODEN library and complile the C/CUDA code
RUN git clone https://github.com/JLBLine/WODEN.git \
  && cd WODEN \
  && mkdir build && cd build \
  && cmake .. -DHBEAM_INC=/mwa_hyperbeam/include/ \
           -DHBEAM_LIB=/mwa_hyperbeam/target/release/libmwa_hyperbeam.so \
  && make
  
##Now do the pip install
RUN cd /WODEN && pip install -r requirements.txt && pip install .

# unity - used for testing
RUN cd / && git clone https://github.com/ThrowTheSwitch/Unity.git \
  && cd Unity \
  && mkdir build && cd build \
  && cmake .. \
  && make -j 4 \
  && make install

##install the tests. no sense in running the tests of course because
##that is user GPU dependent. But we can install them.
RUN cd /WODEN/build && cmake .. -DTARGET_GROUP=test -DUNITY_ROOT=/Unity \
    && make -j4

##==============================================================================

RUN cd / && mkdir MWA_beam_files && cd MWA_beam_files \
    && wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5 \
    && wget http://ws.mwatelescope.org/static/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5


ENV MWA_FEE_HDF5=/MWA_beam_files/mwa_full_embedded_element_pattern.h5
ENV MWA_FEE_HDF5_INTERP=/MWA_beam_files/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5