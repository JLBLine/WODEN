# WODEN installation
WODEN is built around the `atomicAdd` function so you will need an NVIDIA GPU to run it. Currently, WODEN has only been tested and run on linux, specifically Ubuntu 18.04 and 16.04, and the OzStar super cluster of Swinburne University. If you're mad keen to run on Windows or Mac, please contact Jack at jack.line@curtin.edu.au and we can give it a go.

## Dependencies
CMake (https://cmake.org) version >= 3.10 \
NVIDIA CUDA (https://developer.nvidia.com/cuda-downloads) \
CFITSIO (https://heasarc.gsfc.nasa.gov/fitsio/) \
json-c (https://github.com/json-c/json-c)

The wrapper `run_woden.py` relies on the python modules:
astropy
numpy
struct
subprocess
jdcal

Note in future work I plan removing the bespoke uvfits behaviour and linking with pyuvdata (https://pyuvdata.readthedocs.io/en/latest/).

## Installation
Installation is run via CMake. In the WODEN repo, create a `build` dir:
```sh
cd WODEN && mkdir build && cd build
```
The basic way to run CMake from within the build dir is
```sh
cmake ..
make
```
CMake will try and find all the necessary libraries. If it's successful you can then run `make` to compile the code.

### Machine specifics
Depending on your machine and installation, CMake may or may not find certain libraries. CMake should be pretty good at finding CUDA, but I've included the variables `CFITSIO_ROOT` and `JSONC_ROOT` variables which you can set as options to the CMake command to find them (all CMake options have a -D in front):
```
cmake .. \
    -DJSONC_ROOT=/path/to/json-c/installation \
    -DCFITSIO_ROOT=/path/to/cfitsio/installation
```
All NVIDIA GPUs have a specfic compute capability, which relates to their internal architecture. You can tell the compiler which architecture to compile for, which in theory should make compilation quicker, and ensure the code runs correctly on your GPU. You can find out the compute value here (https://developer.nvidia.com/cuda-gpus), and pass it to CMake via:
```
cmake .. -DCUDA_ARCH=6.0
```
If you need to pass extra flags to your CUDA compiler, you can do so by adding something like the following
```sh
-DCMAKE_CUDA_FLAGS="-someflag"
```

### Example
Using the above, an example CMake command used on OzStar was
```sh
cmake .. -DCUDA_ARCH=6.0 \
    -DJSONC_ROOT=/fred/oz048/jline/software/json-c/install/ \
    -DCFITSIO_ROOT=/fred/oz048/MWA/CODE/rts_libs/cfitsio/

make
```

## Post compilation
The python wrapper `run_woden.py` uses an environment variable to locate the `woden` executable. Once you have compiled WODEN, create the following variable (e.g. for `bash`):
```
export WODEN_DIR=/path/to/your/WODEN/build
```
which you can put in your `~/.bashrc` and then forget about it.

The `json-c` library is dynamic, so you need to ensure that WODEN can see it when running. It may already be in your `$LD_LIBRARY_PATH`, but if it isn't, you can add it to your `~/.bashrc` if you are a `bash` user via
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/dir/containing/json-c/lib/
```
