.. _`Windows Subsystem for Linux 2 (WSL 2)`: https://docs.microsoft.com/en-us/windows/wsl/
.. _everybeam insallation page: https://everybeam.readthedocs.io/en/latest/build-instructions.html
.. _casacore installation: https://github.com/casacore/casacore#installation
.. _`table of NVIDIA compute capabilities`: https://developer.nvidia.com/cuda-gpus
.. _`CUDAARCHS explanation`: https://cmake.org/cmake/help/latest/envvar/CUDAARCHS.html

###############
Installation
###############

``WODEN`` is built for speed and intended to work with a GPU. Currently, you need either an NVIDIA GPU to use ``CUDA`` functionality, or something else that can use ``HIP`` (likely an AMD GPU). ``CUDA`` is tried and tested in ``WODEN``, whereas ``HIP`` is new in version 2.2 and not well tested. Furthermore, ``WODEN`` has only been tested to run on linux, specifically Ubuntu 16.04 up to 24.04. This does however include the `Windows Subsystem for Linux 2 (WSL 2)`_., so you can technically run on Windows kinda.

You have two options for installation:

- More work, but tailored to your system: :ref:`install manual`
- Less work, but less flexibility/performance: :ref:`install docker`

Both options are described below, jump to whatever suits you.

``WODEN`` has been tested to run on the following Australian super computers:

- Garrawarla (Pawsey) CUDA (only tested up to ``WODEN`` v2.3)
- OzStar (Swinburne University) CUDA (tested ``WODEN`` v2.6, EveryBeam does not work from singularity)
- Ngarrgu Tindebeek (Swinburne University) CUDA (tested ``WODEN`` v2.6, including EveryBeam)
- Setonix (Pawsey) HIP (only tested up to ``WODEN`` v2.3)

.. _install manual:

Manual Installation
**********************

If you don't want to read, an example installation command (assuming you've installed all dependencies) is::

   HYPER_VERSION=0.10.1
   git clone https://github.com/JLBLine/WODEN.git &&
   cd WODEN &&
   mkdir build && cd build &&
   cmake .. -DHBEAM_INC=/mwa_hyperbeam-${HYPER_VERSION}/include/ \
            -DHBEAM_LIB=/mwa_hyperbeam-${HYPER_VERSION}/target/release/libmwa_hyperbeam.so \
            -DEBEAM_INSTALL=/everybeam_install/ \
            -DEBEAM_ROOT=/EveryBeam/ &&
   make -j8 &&
   cd .. &&
   pip install -r requirements.txt && pip install .

Obviously you'll need to change the paths to where you've installed things. The following sections will guide you through the installation process.

Dependencies
==========================

``WODEN`` has a number of dependencies so it doesn't reinvent the wheel. A brief list of them here is followed by detailed instructions on how I installed them in the following subsection.

- **CMake** - https://cmake.org version >= 3.21
- Either **NVIDIA CUDA** - https://developer.nvidia.com/cuda-downloads
- or **AMD ROCm** - https://rocm.docs.amd.com/projects/install-on-linux/en/latest/
- **rust** - https://www.rust-lang.org/tools/install (needed for ``mwa_hyperbeam``)
- **mwa_hyperbeam** - https://github.com/MWATelescope/mwa_hyperbeam
- **Python >= 3.8** (as well as a number of Python modules, see below)
- *Optional* **EveryBeam** - https://everybeam.readthedocs.io/en/latest/build-instructions.html

How to install dependencies
------------------------------

These instructions are for Ubuntu 24.04, but can be used as a guide for other
linux-like systems.

+ **CMake** - https://cmake.org version >= 3.21::

   $ sudo snap install cmake

+ **NVIDIA CUDA** - https://developer.nvidia.com/cuda-downloads. Best used if you have an NVIDIA GPU. I typically download the runfile option, which you run as::

  $ sudo sh cuda_11.2.2_460.32.03_linux.run ##your version will likely be different

  but I do NOT install the drivers at this point, as I'll already have drivers. Up to you and how your system works. Also, don't ignore the step of adding something like ``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-11.2/lib64`` to your ``~/.bashrc``, or your system won't find ``CUDA``.
+ **AMD ROCm** - https://rocm.docs.amd.com/projects/install-on-linux/en/latest/::

  I don't have an AMD GPU, so I've never done this. Fingers crossed the linked instructions work for you!
+ **mwa_hyperbeam** - https://github.com/MWATelescope/mwa_hyperbeam - ``mwa_hyperbeam`` is the go-to package for calculating the MWA Fully Embedded Element (FEE) primary beam model. At the time of writing (23/03/2022), we'll have to install and compile from source to get the CUDA code that we want to link to. We should be able to install release versions in the future. For now, you'll first need to install ``rust``, the language the library is written in. I followed the installation guide at https://www.rust-lang.org/tools/install, which for me on Ubuntu just means running::

  $ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

  Once that's installed, I run the following commands for a CUDA installation (you can choose where to install it, I'm just putting where I happened to do it this time round)::

  $ cd /home/jline/software
  $ git clone https://github.com/MWATelescope/mwa_hyperbeam.git
  $ cd mwa_hyperbeam
  $ export HYPERDRIVE_CUDA_COMPUTE=60 ##your compute capability
  $ cargo build --locked --release --features=cuda,hdf5-static

  .. note:: ``export HYPERDRIVE_CUDA_COMPUTE=60`` is not essential as the compiler should be smart enough, but you *might* get a speed boost but setting the correct architecture. This of course depends on your GPU; see 'Machine specifics' below on how to work out your architecture.

  If you have an AMD GPU, replace the last two lines with something like::

  $ export HYPERBEAM_HIP_ARCH=gfx90a
  $ cargo build --locked --release --features=hip,hdf5-static

  where again the value of ``HYPERBEAM_HIP_ARCH`` depends on what kind of GPU you have.

  That's it! I'll show you how to link to it later when we install ``WODEN``. If you don't want to have to tell ``CMake`` where to look for the libraries, you'll need to link/copy ``libmwa_hyperbeam.so`` somewhere your compiler can see, as well as ``mwa_hyperbeam.h``.
+ **python >= 3.8** - 3.8 should work, but I'd suggest going with 3.11 or 3.12. The requirements can be found in ``WODEN/requirements.txt``, which you can install via something like::

  $ pip3 install -r requirements_testing.txt

+ **everybeam** - You only need to do this if you want to use EveryBeam primary beams in your simulations. ``WODEN`` will check for the ``everybeam`` Python package at runtime, and run fine if it's missing (unless you ask it to run with EveryBeam, then it will complain). The build instructions for EveryBeam live on the `everybeam insallation page`_. You'll have to install a bunch of dependencies, including ``casacore``. This is sad, but not too hard on Ubuntu. You can follow them with *one important difference* - you need to currently install my branch, as it's the only one that has MWA Python bindings. You can do this by running::

  $ git clone -b mwa_python_wrapper --recursive -j4 https://git.astron.nl/RD/EveryBeam.git

  If you follow the EveryBeam instructions, you'll install a system-wide version. If this isn't suitable for your system, i.e. you're not running in a container, this was may approach. When I ran ``cmake``, I did::

  $ cmake .. -DCMAKE_INSTALL_PREFIX=/some/path/software/installed/ 

+ **casacore** - If you installed ``EveryBeam`` via the instructions above, you should have ``casacore`` installed. If not, you'll need to install it. You can do this via::

  $ sudo apt install libcasacore-dev

  or if you want to install from source, follow the instructions at `casacore installation`_.  If you installed it in a non-standard location, you'll need to point ``CMake`` to it when you compile ``WODEN`` (see below).


Compiling ``WODEN`` ``C/CUDA`` code
=========================================
.. note:: The installation files used to make docker images live in ``WODEN/docker``. You can use them as a template to install on Ubuntu. If you have another OS, you should be able to use them as a guide.

In an ideal world, if the installation of your dependencies went perfectly,
you have a newer NVIDIA GPU, and you don't care about EveryBeam, you should be able to simply run::

  $ git clone https://github.com/JLBLine/WODEN.git
  $ cd WODEN
  $ mkdir build && cd build
  $ cmake ..
  $ make -j 4

et voila, your code is compiled. Keep reading to see how to install ``WODEN`` so you can run it from anywhere.

.. warning:: Even if the code compiled, if your GPU has a compute capability < 5.1, newer versions of ``nvcc`` won't compile code that will work. You'll get error messages like "No kernel image available". Check out how to fix that in 'NVIDIA specifics' below.

HyperBeam
------------------------------
It's almost a guarantee ``cmake`` won't be able to find ``mwa_hyperbeam``, so you'll have to point it to where things are installed. You can use two keywords in the following way to achieve that::

  $ cmake .. -DHBEAM_INC=/home/jline/software/mwa_hyperbeam/include \
             -DHBEAM_LIB=/home/jline/software/mwa_hyperbeam/target/release/libmwa_hyperbeam.so

Obviously you'll need to point to where you have installed things. If *you* have a library with my name in the path I'd be concerned, so edit it as appropriate.

EveryBeam
------------------------------
``EveryBeam`` relies on ``casacore``. As ``WODEN`` compiles some EveryBeam code directly,
you have to link ``WODEN`` against it. If you've installed via Ubuntu, ``CMake`` will be able to find it. If you've installed it in a non-standard location, you'll need to point ``CMake`` to it. You can do this via::

  $ cmake .. -DCASACORE_ROOT_DIR=/path/to/casacore

Next, you not only have to point ``CMake`` to where you installed the compiled ``EveryBeam``, but also to the source code as well. Parts of ``EveryBeam`` rely on ``aocommon`` headers, which doesn't get copied during ``EveryBeam``  installation. So you need to point ``EBEAM_INSTALL`` to where you installed the compiled
``EveryBeam`` outputs, and ``EBEAM_ROOT`` to where you installed the source code. You can do this via::

  $ cmake .. -DEBEAM_INSTALL=/some/path/software/installed/ \
             -DEBEAM_ROOT=/path/to/EveryBeam


NVIDIA/CUDA specifics
------------------------------
All NVIDIA GPUs have a specific compute capability, which relates to their internal architecture. CMake seems to be pretty good at figuring things out itself, so the default it to use the CMake defaults. You can tell the compiler however which architecture to compile for, which in theory should make compilation quicker, and ensure the code runs correctly on your GPU. You can find out the compute value via this `table of NVIDIA compute capabilities`_ , and pass it to CMake via setting the ``CUDAARCHS`` environment variable (`CUDAARCHS explanation`_) BEFORE you run the call to ``cmake``::

  $ export CUDAARCHS=60

(for a compute capability of 6.0, for example).

You can pass multiple compute capabilities by separating them with a comma if you want to suuport multiple architectures. This will certainly make the compilation take longer, and *possibly* slow things down on certain cards (GPUs can be fickle). If you want to do that, you can pass something like the following::

  $ export CUDAARCHS="60;61;70"

.. warning:: For newer ``CUDA`` versions, some compute capabilities are deprecated, so the compiler leaves them out by default. For example, using ``CUDA`` version 11.2, compute capabilities 3.5 to 5.0 are ignored. If you card has a compute capability of 5.0, you **must** set ``export CUDAARCHS=50``, otherwise the ``nvcc`` compiler will not create an executable capable of running on your device.

If you need to pass extra flags to your CUDA compiler, you can do so by adding something like the following (noting that all CMake flags start with ``-D``)::

  -DCMAKE_CUDA_FLAGS="-Dsomeflag"


AMD/HIP specifics
------------------------------

If you have an AMD GPU, you can compile the ``HIP`` code instead of the ``CUDA`` code. This is a new feature in ``WODEN`` and not as well tested. You can compile the ``HIP`` code by setting the ``USE_HIP`` flag to ``ON`` when you run ``cmake`` (you'll still need to link )::

  $ cmake .. -DUSE_HIP=ON \
      -DHBEAM_INC=/home/jline/software/mwa_hyperbeam/include \
      -DHBEAM_LIB=/home/jline/software/mwa_hyperbeam/target/release/libmwa_hyperbeam.so

Similarly to ``CUDA``, you can set a ``HIP`` architecture. To find out which one you need, try::
  
  $ offload-arch

which spat out ``gfx90a`` for me. You pass that onto ``cmake`` via the ``HIP_ARCH`` flag::

  $ cmake .. -DUSE_HIP=ON -DHIP_ARCH=gfx90a \
      -DHBEAM_INC=/home/jline/software/mwa_hyperbeam/include \
      -DHBEAM_LIB=/home/jline/software/mwa_hyperbeam/target/release/libmwa_hyperbeam.so

Fair warning, I *had* to include the ``HIP_ARCH`` flag. The code would compile fine but not work at runtime, so a bit nasty.

Installing ``wodenpy``
==========================

OK, we've compiled the C/GPU libraries; now to install the ``WODEN`` Python package and executables. You can do this by running::

  $ cd WODEN
  $ pip3 install .

That's it. You should be able to run ``run_woden.py --help`` on the command line.

Post compilation (optional)
==============================

If you want to use the MWA FEE primary beam model, you must have the stored spherical harmonic coefficients hdf5 file ``mwa_full_embedded_element_pattern.h5``. You can then define this environment variable in your ``~/.bash_rc``::

  export MWA_FEE_HDF5=/path/to/your/location/mwa_full_embedded_element_pattern.h5

so ``run_woden.py`` can find it. There is a command line option ``--hdf5_beam_path`` in ``run_woden.py`` which you can use instead of this environment variable if you want.

If you don't have the spherical harmonic file you can obtain it via the command::

  $ wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5

To use the interpolated MWA FEE beam model, do similarly::

  $ wget http://ws.mwatelescope.org/static/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5
  $ export MWA_FEE_HDF5_INTERP=/path/to/your/location/MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5


.. _install docker:

Use a ``Docker`` image
**************************

.. note:: All the images listed here were created with the script ``WODEN/docker/make_docker_image.sh``. If a particular image doesn't work for you, you can edit the source to hopefully get it working.

For CUDA
--------------

Fair warning, this is a new option, and hasn't been heavily tested. I have successfully run it on a number of clusters (via singularity). Which version you pull depends on your GPU. If you have an NVIDIA GPU, you need to work out what your compute capability is, and pull the appropriate image. Say you have an NVIDIA V100 card, you have a compute capacity of 7.0, so you'd pull the image like this::

  $ docker pull jlbline/woden-2.6:cuda-70

I have made images for computes ``60,61,70,75,80,86``. If you need another compute, either run the Docker script to make a new docker image, or just compile from source as instructed above. 

Alternatively, if you really want an image which can in theory run on many architectures, I've built an image with all ``60,61,70,75,80,86`` computes, called

$ docker pull jlbline/woden-2.6:cuda-multi

This should run on any NVIDIA GPU with a compute capability in that list. However, depending on your compute environment, it may not run as fast as a specific compute image.

In theory, you can just run ``WODEN`` commands by doing something like this::

  $ docker run -it --gpus all woden-2.5:cuda-70 \
    --env XDG_CONFIG_HOME=/somewhere/astropy_storage \
    --env XDG_CACHE_HOME=/somewhere/astropy_storage \
    run_woden.py --help

where the ``--gpus all`` means the docker instance can see your GPUs. The environment variables point to somewhere to keep your ``astropy`` outputs, which is useful if you're running somewhere you're not admin (like on a cluster). There must be a better way to do this but I'm a ``docker`` noob.

For HIP
--------------

The only HIP image I've made is for the Setonix cluster, and is based on a Pawsey specific image https://quay.io/repository/pawsey/rocm-mpich-base?tab=tags&tag=latest. You can pull it like this::

  $ docker pull jlbline/woden-2.3:setonix

It is *highly* unlikely it won't work anywhere else.

Using singularity
**************************

For CUDA
--------------

If your system has ``singularity`` and not docker, you can convert the docker image to a singularity image via::

  $ singularity build woden-2.5-70.sif docker://jlbline/woden-2.3:cuda-70

with an example of running the help looking something like::

  $ singularity exec --nv --home=/astro/mwaeor/jline \
    woden-2.5-70.sif run_woden.py --help

Similarly to the ``docker`` image, ``--nv`` means use the NVIDIA GPUs, and ``--home`` sets a specific location to treat as home if you're not on a local machine.

For HIP
--------------

Again, the only HIP image I've made is for the Setonix cluster where you can do::
  
    $ singularity build woden-2.3:setonix.sif docker://jlbline/woden-2.3:setonix

and run it like::

  $ singularity exec --home=/scratch/mwaeor/jline \
  ${MYSOFTWARE}/woden-2.3-setonix.sif run_woden.py --help

.. warning:: EVERYTHING on the internet will tell you to use the ``--rocm`` flag. This WILL NOT WORK with the Setonix based image, because of shenanigans. So leave it be.

Updating ``astropy`` IERS
**************************

If you're working on a cluster, often you won't have access to the internet. If you try and a recent date, you'll likely need to update the IERS data. This normally happens automatically, but obviously without internet you're cooked. If you're using a singularity image, I've included a little helper script, which you can launch on the head node (which will have internet) and then ``astropy`` will store a new IERS file in a location it thinks is sensible (which should be accessible to the compute nodes). You can run it like this::

  $ apptainer shell --home=/path/to/cluster/home /path/to/sing_file/woden-2.5-80.sif
  $ Apptainer> python3 /WODEN/docker/fetch_iers_data.py

Note here I'm using ``apptainer``, but you can swap ``singularity`` in for ``apptainer``. It just depends on what system you're using. Also note that you need to change ``/path/to/cluster/home`` and ``/path/to/sing_file`` to the appropriate locations, but leave ``/WODEN/docker/fetch_iers_data.py`` as is, as it's in the image.