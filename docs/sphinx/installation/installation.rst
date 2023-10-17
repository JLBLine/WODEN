*************
Installation
*************

WODEN is built on CUDA so you will need an NVIDIA GPU to run it. Currently, WODEN has only been tested and run on linux, specifically Ubuntu 16.04 up to 20.04, the OzStar super cluster of Swinburne University, and Garrawarla cluster of Pawsey. You have two options for installation:

 - More work, but tailored to your system: :ref:`install manual`
 - Less work, but less flexibility/performance: :ref:`install docker`

Both options are described below, jump to whatever suits you.

.. _install manual:

Manual Installation
######################

Dependencies
-----------------

``WODEN`` has a number of dependencies so it doesn't reinvent the wheel. A brief list of them here is followed by detailed instructions on how I installed them in the following subsection.

- **CMake** - https://cmake.org version >= 3.10
- **NVIDIA CUDA** - https://developer.nvidia.com/cuda-downloads
- **HDF5** - https://www.hdfgroup.org/downloads/hdf5/ (needed for ``mwa_hyperbeam``)
- **rust** - https://www.rust-lang.org/tools/install (needed for ``mwa_hyperbeam``)
- **mwa_hyperbeam** - https://github.com/MWATelescope/mwa_hyperbeam
- **Python >= 3.8** (as well as a number of Python modules, see below)

How to install dependencies
****************************

These instructions are for Ubuntu 20.04, but can be used as a guide for other
linux-like systems.

+ **CMake** - https://cmake.org version >= 3.10::

   $ sudo snap install cmake

+ **NVIDIA CUDA** - https://developer.nvidia.com/cuda-downloads. I typically download the runfile option, which you run as::

  $ sudo sh cuda_11.2.2_460.32.03_linux.run

  but I do NOT install the drivers at this point, as I'll already have drivers. Up to you and how your system works. Also, don't ignore the step of adding something like ``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-11.2/lib64`` to your ``~/.bashrc``, or your system won't find ``CUDA``.
+ **HDF5** - https://www.hdfgroup.org/downloads/hdf5/ - just do::

  $ sudo apt install libhdf5-serial-dev
+ **mwa_hyperbeam** - https://github.com/MWATelescope/mwa_hyperbeam - ``mwa_hyperbeam`` is the go-to package for calculating the MWA Fully Embedded Element (FEE) primary beam model. At the time of writing (23/03/2022), we'll have to install and compile from source to get the CUDA code that we want to link to. We should be able to install release versions in the future. For now, you'll first need to install ``rust``, the language the library is written in. I followed the installation guide at https://www.rust-lang.org/tools/install, which for me on Ubuntu just means running::

  $ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

  Once that's installed, I ran the following commands (you can choose where to install it, I'm just putting where I happended to do it this time round)::

  $ cd /home/jline/software
  $ git clone https://github.com/MWATelescope/mwa_hyperbeam.git
  $ cd mwa_hyperbeam
  $ cargo build --release --features=cuda,cuda-static

  That's it! I'll show you how to link to it later when we install ``WODEN``. If you don't want to have to tell ``CMake`` where to look for the libraries, you'll need to link/copy ``libmwa_hyperbeam.so`` somewhere your compiler can see, as well as ``mwa_hyperbeam.h``.
+ **python >= 3.8** - The requirements can be found in ``WODEN/requirements.txt``, which you can install via something like::

  $ pip3 install -r requirements_testing.txt

For completeness, those packages are::

  sphinx_argparse
  breathe
  astropy
  numpy
  pyerfa
  palpy
  importlib_resources
  sphinx-math-dollar
  matplotlib
  pyuvdata
  python-casacore

Phew! That's it for now.

Compiling ``WODEN`` ``C/CUDA`` code
**************************************

In an ideal world, if the installation of your dependencies went perfectly and
you have a newer NVIDIA GPU, you should be able to simply run::

  $ git clone https://github.com/JLBLine/WODEN.git
  $ cd WODEN
  $ mkdir build && cd build
  $ cmake ..
  $ make -j 4

et voila, your code is compiled. Keep reading to see how to install ``WODEN`` so you can run it from anywhere.

.. warning:: Even if the code compiled, if your GPU has a compute capability < 5.1, newer versions of ``nvcc`` won't compile code that will work. You'll get error messages like "No kernel image available". Check out how to fix that in 'Machine specifics' below.

Machine specifics
**********************
It's almost a guarantee ``cmake`` won't be able to find ``mwa_hyperbeam``, so you'll have to point it to where things are installed. You can use two keywords in the following way to achieve that::

  $ cmake .. -DHBEAM_INC=/home/jline/software/mwa_hyperbeam/include \
             -DHBEAM_LIB=/home/jline/software/mwa_hyperbeam/target/release/libmwa_hyperbeam.so

Obviously you'll need to point to where you have installed things. If *you* have a library with my name in the path I'd be concerned, so edit it as appropriate.

All NVIDIA GPUs have a specific compute capability, which relates to their internal architecture. You can tell the compiler which architecture to compile for, which in theory should make compilation quicker, and ensure the code runs correctly on your GPU. You can find out the compute value here (https://developer.nvidia.com/cuda-gpus), and pass it to CMake via::

  $ cmake .. -DCUDA_ARCH=6.0

(for a compute capability of 6.0, for example).

.. warning:: For newer ``CUDA`` versions, some compute capabilities are deprecated, so the compiler leaves them out by default. For example, using ``CUDA`` version 11.2, compute capabilities 3.5 to 5.0 are ignored. If you card has a compute capability of 5.0, you **must** include the flag ``-DCUDA_ARCH=5.0``, otherwise the `nvcc` compiler will not create an executable capable of running on your device.

If you need to pass extra flags to your CUDA compiler, you can do so by adding something like the following (noting that all CMake flags start with ``-D``)::

  -DCMAKE_CUDA_FLAGS="-Dsomeflag"


Installing ``WODEN``
*****************************

We've compiled the C/CUDA libraries; now to install the ``WODEN`` Python package and executables. You can do this by running::

  $ cd WODEN
  $ pip3 install .

That's it. You should be able to run ``run_woden.py --help`` on the command line.

Post compilation (optional)
*****************************

If you want to use the MWA FEE primary beam model, you must have the stored spherical harmonic coefficients hdf5 file ``mwa_full_embedded_element_pattern.h5``. You can then define this environment variable in your ``~/.bash_rc``::

  export MWA_FEE_HDF5=/path/to/your/location/mwa_full_embedded_element_pattern.h5

so ``run_woden.py`` can find it. There is a command line option ``--hdf5_beam_path`` in ``run_woden.py`` which you can use instead of this environment variable if you want.

If you don't have the spherical harmonic file you can obtain it via the command::

  $ wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5

To use the interpolated MWA FEE beam model, do similarly::

  $ wget http://ws.mwatelescope.org/static/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5
  $ export MWA_FEE_HDF5_INTERP=/path/to/your/location/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5


.. _install docker:

Use the ``Docker`` image
##########################
Fair warning, this is a new option, and hasn't been heavily tested. I have successfully run it on my desktop and the Garrawarla supercluster of Pawsey, but that's it. This docker image is built upon the ``nvidia/cuda:11.4.3-devel-ubuntu20.04`` image, and so your local NVIDIA cards / drives **have to work with CUDA 11.4.3; docker uses the local NVIDIA drivers**. You can pull the image from Docker Hub via::

  $ docker pull jlbline/woden-2.0

Then in theory, you can just run WODEN commands by doing something like this::

  $ docker run -it --gpus all woden-2.0 \
    --env XDG_CONFIG_HOME=/somewhere/astropy_storage \
    --env XDG_CACHE_HOME=/somewhere/astropy_storage \
    run_woden.py --help

where the ``--gpus all`` means the docker instance can see your GPUs. The environment variables point to somewhere to keep your ``astropy`` outputs, which is useful if you're running somewhere you're not admin (like on a clsuter). There must be a better way to do this but I'm a ``docker`` noob.

Using singularity
******************
If your system has ``singularity`` and not docker, you can convert the docker image to a singularity image via::

  $ singularity build woden-2.0.sif docker://jlbline/woden-2.0

with an example of running the help looking something like::

  $ singularity exec --nv --home=/astro/mwaeor/jline \
    woden-2.0.sif run_woden.py --help

Similarly to the ``docker`` image, ``--nv`` means use the GPUs, and ``--home`` sets a specific location to treat as home if you're not on a local machine.