module purge
module load cuda/12.0.0
module load gcc/12.3.0
module load rust/1.70.0
module load cmake/3.26.3
module load python/3.11.3

MYSOFTWARE=/fred/oz048/jline/software_nt
WODEN_NAME="woden-2.1.0"
HYPER_NAME="hyperbeam-0.9.2"

cd $MYSOFTWARE

##Do HYPERBEAM------------------------------------------------------------------

export HYPERDRIVE_CUDA_COMPUTE=80
export HBEAM_GIT=$MYSOFTWARE/$HYPER_NAME # where to put hyperbeam
git clone --branch v0.9.2 https://github.com/MWATelescope/mwa_hyperbeam.git $HBEAM_GIT && cd $HBEAM_GIT
cargo build --locked --release --features=cuda,hdf5-static


##Do WODEN----------------------------------------------------------------------

export WODEN_GIT=${MYSOFTWARE}/${WODEN_NAME} ##where to put WODEN

##git clone and switch to appropriate branch
git clone https://github.com/JLBLine/WODEN.git $WODEN_GIT && cd $WODEN_GIT
git checkout v2.1.0

# ##cd in and cmake to compile the HIP version
mkdir build && cd build
cmake .. -DUSE_CUDA=ON \
    -DHBEAM_INC=${HBEAM_GIT}/include \
    -DHBEAM_LIB=${HBEAM_GIT}/target/release/libmwa_hyperbeam.so
make

cd $WODEN_GIT
python -m venv pyenv-${WODEN_NAME}
source ${MYSOFTWARE}/${WODEN_NAME}/pyenv-${WODEN_NAME}/bin/activate

# ##install dependencies
pip install -U pip
pip install -r requirements.txt
# ##install actual WODEN

pip install .


# ##this is used to get helpful testing and coverage framework
export UNITY_GIT=$MYSOFTWARE/Unity
cd $MYSOFTWARE
git clone https://github.com/ThrowTheSwitch/Unity.git


##DO this if you want to run some tests, easy interactive jobbie-----------------
# sinteractive --time 1:00:00 --nodes 1 --cpus-per-task 2 --mem=50gb --gres=gpu:1 --account=oz048

# module load cuda/12.0.0
# module load gcc/12.3.0
# module load cmake/3.26.3
# module load python/3.11.3
# module load mwa_beams/latest

# MYSOFTWARE=/fred/oz048/jline/software_nt
# WODEN_NAME="woden-2.1.0"

# export PATH=$PATH:"${MYSOFTWARE}/${WODEN_NAME}/pyenv-${WODEN_NAME}/bin"
# export PYTHONPATH=$PYTHONPATH:"${MYSOFTWARE}/${WODEN_NAME}/pyenv-${WODEN_NAME}/lib/python3.11/site-packages/"

# export MWA_FEE_HDF5=/fred/oz048/achokshi/software/mwa_beams/mwa_full_embedded_element_pattern.h5
# export MWA_FEE_HDF5_INTERP=/fred/oz048/achokshi/software/mwa_beams/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5

# # ##switches to making tests. If you're on a cluster, sometimes there are
# # ##multiple python executables. Let's point cmake directly at the one we want
# cd $WODEN_GIT/build
# cmake .. -DTARGET_GROUP=test -DUNITY_ROOT=$UNITY_GIT \
#   -DPYTHON_EXECUTABLE=${MYSOFTWARE}/${WODEN_NAME}/pyenv-${WODEN_NAME}/bin/python
# make

# # ##you can now run tests via
# ctest