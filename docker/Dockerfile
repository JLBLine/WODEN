#Base image including CUDA to write everything on to
FROM nvidia/cuda:11.4.3-devel-ubuntu20.04

##install the dependencies
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

# Get Rust
RUN mkdir -m755 /opt/rust /opt/cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/cargo PATH=/opt/cargo/bin:$PATH
ENV RUST_VERSION=1.72
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain=$RUST_VERSION

# hyperbeam used for MWA primary beam
RUN wget https://github.com/MWATelescope/mwa_hyperbeam/archive/refs/tags/v0.6.0.tar.gz \
  && tar -xvf v0.6.0.tar.gz \
  && cd mwa_hyperbeam-0.6.0 \
  && /opt/cargo/bin/cargo build --release --features=cuda,cuda-static \
  && cd /

## clone the release verion of WODEN and complile the C/CUDA code
RUN git clone -b v2.0.0 https://github.com/JLBLine/WODEN.git \
  && cd WODEN \
  && mkdir build && cd build \
  && cmake .. -DHBEAM_INC=/mwa_hyperbeam-0.6.0/include/ \
           -DHBEAM_LIB=/mwa_hyperbeam-0.6.0/target/release/libmwa_hyperbeam.so \
  && make
  
##Now do the pip install
RUN cd /WODEN && pip install -r requirements.txt && pip install .

# unity - used for testing
RUN cd / && git clone https://github.com/ThrowTheSwitch/Unity.git

##install the tests. no sense in running the tests of course because
##that is user GPU dependent. But we can install them.
RUN cd /WODEN/build && cmake .. -DTARGET_GROUP=test -DUNITY_ROOT=/Unity \
    && make -j4

##copy the MWA primary beam files
RUN cd / && mkdir MWA_beam_files && cd MWA_beam_files \
    && wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5 \
    && wget http://ws.mwatelescope.org/static/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5

##make some environment variables
ENV MWA_FEE_HDF5=/MWA_beam_files/mwa_full_embedded_element_pattern.h5
ENV MWA_FEE_HDF5_INTERP=/MWA_beam_files/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5