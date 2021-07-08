# WODEN installation
WODEN is built on CUDA so you will need an NVIDIA GPU to run it. Currently, WODEN has only been tested and run on linux, specifically Ubuntu 16.04 up to 20.04,the OzStar super cluster of Swinburne University, and Garrwarla cluster of Pawsey. If you're mad keen to run on Windows or Mac, please contact Jack at jack.line@curtin.edu.au and we can give it a go.

## Dependencies
 - CMake https://cmake.org version >= 3.10 \
 - NVIDIA CUDA https://developer.nvidia.com/cuda-downloads \
 - json-c https://github.com/json-c/json-c \
 - ERFA (Essential Routines for Fundamental Astronomy) https://github.com/liberfa/erfa \
 - HDF5 https://www.hdfgroup.org/downloads/hdf5/ - if on Ubuntu, do something like `sudo apt install libhdf5-serial-dev` \
 - PAL (Positional Astronomy Library) https://github.com/Starlink/pal/releases - Easiest thing is to download a release version from this page, rather than trying to install from git. Then do
```sh
./configure --prefix=/usr/local --without-starlink
    make
    make install
```
which if something like Ubuntu will put the `include` files in `/usr/local/include/star` and `lib` in `/usr/local/lib`. Make sure `/usr/local/lib` is in your `LD_LIBRARY_PATH` so `WODEN` can find the PAL library.

The wrapper `run_woden.py` relies on the python modules:
astropy
numpy
struct
subprocess

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
All NVIDIA GPUs have a specific compute capability, which relates to their internal architecture. You can tell the compiler which architecture to compile for, which in theory should make compilation quicker, and ensure the code runs correctly on your GPU. You can find out the compute value here (https://developer.nvidia.com/cuda-gpus), and pass it to CMake via:
```
cmake .. -DCUDA_ARCH=6.0
```
>WARNING - for newer CUDA versions, some compute capabilities are deprecated, so the compiler leaves them out by default. For example, using CUDA version 11.2, compute capabilities 3.5 to 5.0 are ignored. If you card has a compute capability of 5.0, you _must_ include the flag `-DCUDA_ARCH=5.0`, otherwise the `nvcc` compiler will not create an executable capable of running on your device.

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
The python wrapper `run_woden.py` uses an environment variable to locate the `woden` executable. Once you have compiled WODEN, run the following script (e.g. for `bash`):
```
source /path/to/your/WODEN/build/init_WODEN.sh
```
This will create the variable `$WODEN_DIR`, and add it to your `$PATH`. You can put in your `~/.bashrc` and then forget about it.

The `json-c` library is dynamic, so you need to ensure that WODEN can see it when running. It may already be in your `$LD_LIBRARY_PATH`, but if it isn't, you can add it to your `~/.bashrc` if you are a `bash` user via
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/dir/containing/json-c/lib/
```

# Installing tests
You can run the unit/itegration tests I use for developing `WODEN` to check your system / the code runs as expected. The tests are compiled using the same `cmake` script as the main code. The tests use the following C library:

 - Unity https://github.com/ThrowTheSwitch/Unity/releases

You don't need to install Unity anywhere, `WODEN` uses the `C` library directly. You just need to tell `WODEN` where Unity lives in your system (in the example below I have downloaded the Unity release version 2.5.2). You also need to add `TARGET_GROUP=test` to your command to tell CMake to build tests instead of the main code:

```sh
cd WODEN/build
cmake .. -DTARGET_GROUP=test -DUNITY_ROOT=/usr/local/Unity-2.5.2
make
```

Once that compiles, you can run the tests by running
```
ctest
```
NOTE - To test MWA Fully Embedded Element beam code, you must have the environment variable
```
MWA_FEE_HDF5=/path/to/mwa_full_embedded_element_pattern.h5
```
declared. If you don't have that, the MWA FEE beam code tests will just be skipped. If you want more detail of what the tests are doing, run `ctest --verbose`.



> WARNING - once you have done this, to go back to compiling the main `woden`
executable, you need to run
```
cmake .. -DTARGET_GROUP=production
```
otherwise you'll just keep building the tests.
