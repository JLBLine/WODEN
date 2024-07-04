#!/bin/bash
# salloc --nodes=1 --partition=gpu --account=${PAWSEY_PROJECT}-gpu -t 00:30:00 --gres=gpu:1

##what to call WODEN directory
WODEN_NAME="woden-hip"

##Load some dependencies
module load rust
module load rocm/5.7.3
module load python/3.11.6

cd $MYSOFTWARE

##Do HYPERBEAM------------------------------------------------------------------

export HBEAM_GIT=$MYSOFTWARE/hyperbeam # where to put hyperbeam
git clone --branch v0.9.2 https://github.com/MWATelescope/mwa_hyperbeam.git $HBEAM_GIT && cd $HBEAM_GIT

export HYPERBEAM_HIP_ARCH=gfx90a
cargo build --locked --release --features=hip,hdf5-static


##Do WODEN----------------------------------------------------------------------

export WODEN_GIT=${MYSOFTWARE}/${WODEN_NAME} ##where to put WODEN

##git clone and switch to appropriate branch
git clone https://github.com/JLBLine/WODEN.git $WODEN_GIT && cd $WODEN_GIT
cd $MYSOFTWARE/${WODEN_NAME}
git checkout compile_both_CUDA_and_HIP

##cd in and cmake to compile the HIP version
mkdir build && cd build
cmake .. -DUSE_HIP=ON -DHIP_ARCH=gfx90a \
    -DHBEAM_INC=${HBEAM_GIT}/include \
    -DHBEAM_LIB=${HBEAM_GIT}/target/release/libmwa_hyperbeam.so
make

##get a python environment and install in the WODEN git repo
cd $WODEN_GIT
python -m venv pyenv-${WODEN_NAME}
source ${MYSOFTWARE}/${WODEN_NAME}/pyenv-${WODEN_NAME}/bin/activate

##install dependencies
pip install -U pip
pip install -r requirements.txt
##install actual WODEN
pip install .

##get the beam files
cd $MYSOFTWARE
mkdir MWA_beam_files && cd MWA_beam_files \
    && wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5 \
    && wget http://ws.mwatelescope.org/static/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5

export MWA_FEE_HDF5=${MYSOFTWARE}/MWA_beam_files/mwa_full_embedded_element_pattern.h5
export MWA_FEE_HDF5_INTERP=${MYSOFTWARE}/MWA_beam_files/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5


##THIS IS NEEDED FOR CTEST------------------------------------------------------
##NOTE you need to have MWA_FEE_HDF5 and MWA_FEE_HDF5_INTERP defined for all
##tests to run. FEE beam tests will just skip, not error, if missing

##this is used to get helpful testing and coverage framework
export UNITY_GIT=$MYSOFTWARE/Unity
cd $MYSOFTWARE
git clone https://github.com/ThrowTheSwitch/Unity.git

##switches to making tests. If you're on a cluster, sometimes there are
##multiple python executables. Let's point cmake directly at the one we want
cd $WODEN_GIT/build
cmake .. -DTARGET_GROUP=test -DUNITY_ROOT=$UNITY_GIT \
  -DPYTHON_EXECUTABLE=${MYSOFTWARE}/${WODEN_NAME}/pyenv-${WODEN_NAME}/bin/python
make

##you can now run tests via
#ctest