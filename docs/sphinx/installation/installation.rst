*************
Installation
*************

WODEN is built on CUDA so you will need an NVIDIA GPU to run it. Currently, WODEN has only been tested and run on linux, specifically Ubuntu 16.04 up to 20.04, the OzStar super cluster of Swinburne University, and Garrawarla cluster of Pawsey. If you're mad keen to run on Windows or Mac, please contact Jack at jack.l.b.line@gmail.com and we can give it a go.

Dependencies
##############

``WODEN`` has a number of dependencies so it doesn't reinvent the wheel. A brief list of them here is followed by detailed instructions on how I installed them in the following subsection. Note that the explicit installation instructions I have included for ``json-c``, ``erfa``, and ``pal`` are the only way I have reliably managed to install these packages - the package installation manager sometimes does whacky things for them.

- **CMake** - https://cmake.org version >= 3.10
- **NVIDIA CUDA** - https://developer.nvidia.com/cuda-downloads
- **json-c** - https://github.com/json-c/json-c
- **ERFA** - https://github.com/liberfa/erfa/releases
- **HDF5** - https://www.hdfgroup.org/downloads/hdf5/
- **PAL** - https://github.com/Starlink/pal/releases
- **python >= 3.6**

How to install dependencies
****************************

These instructions are for Ubuntu 20.04, but can be used as a guide for other
linux-like systems.

+ **CMake** - https://cmake.org version >= 3.10::

   $ sudo snap install cmake

+ **NVIDIA CUDA** - https://developer.nvidia.com/cuda-downloads. I typically download the runfile option, which you run as::

  $ sudo sh cuda_11.2.2_460.32.03_linux.run

  but I do NOT install the drivers at this point, as I'll already have drivers. Up to you and how your system works. Also, don't ignore the step of adding something like ``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-11.2/lib64`` to your ``~/.bashrc``, or your system won't find ``CUDA``.
+ **json-c** - https://github.com/json-c/json-c. This is a typical ``cmake`` installation::

  $ git clone https://github.com/json-c/json-c.git
  $ cd json-c
  $ mkdir build && cd build
  $ cmake ..
  $ make -j 4
  $ sudo make install

  When you run ``cmake ..`` you should find out what dependencies you are missing and can install them as needed.
+ **ERFA** - https://github.com/liberfa/erfa/releases. I think it's best to install a release version of ``ERFA``. Comes with a ``configure`` file, while the ``git`` repo doesn't. An installation route would look like::

  $ wget https://github.com/liberfa/erfa/releases/download/v2.0.0/erfa-2.0.0.tar.gz
  $ tar -xvf erfa-2.0.0.tar.gz
  $ cd erfa-2.0.0
  $ ./configure
  $ make -j 4
  $ sudo make install
+ **HDF5** - https://www.hdfgroup.org/downloads/hdf5/ - just do::

  $ sudo apt install libhdf5-serial-dev
+ **PAL** - https://github.com/Starlink/pal/releases - ``PAL`` is a little mental with it's default installation paths. I *HIGHLY* recommend downloading a release version, and then using the ``--without-starlink`` option::

  $ wget https://github.com/Starlink/pal/releases/download/v0.9.8/pal-0.9.8.tar.gz
  $ tar -xvf pal-0.9.8.tar.gz
  $ cd pal-0.9.8
  $ ./configure --prefix=/usr/local --without-starlink
  $ make
  $ sudo make install

  Doing it this way installs things in normal locations, making life easier during linking.
+ **python >= 3.6** - the best way to run ``WODEN`` is through the script ``run_woden.py``, which has a number of package dependencies. One of these is ``pyerfa``, which uses f-strings during installation, so you have to use a python version >= 3.6. Sorry. The requirements can be found in ``WODEN/docs/sphinx/sphinx/requirements_testing.txt``, which you can install via something like::

  $ pip3 install -r requirements_testing.txt

For completeness, those packages are::

  sphinx_argparse
  breathe
  astropy
  numpy
  pyerfa
  palpy
  matplotlib

The ``sphinx_argparse, breathe`` packages are used for the documentation. Further packages of ``palpy, matplotlib`` are only used in the ``test_installation/absolute_accuracy`` test, so if you're aiming for a minimal installation, you only need ``numpy, astropy, and pyerfa``.

Phew! That's it for now.

Compiling ``WODEN``
######################

In an ideal world, if the installation of your dependencies went perfectly and
you have a newer NVIDIA GPU, you should be able to simply run::

  $ git clone https://github.com/JLBLine/WODEN.git
  $ cd WODEN
  $ mkdir build && cd build
  $ cmake ..
  $ make -j 4

et voila, your code is compiled. If this worked, and you're happy to install ``WODEN`` into the system default location, just run::

  $ sudo make install

(usually the default is something like ``/usr/local`` hence you need admin privileges). If complilation fails or you're not used to ``cmake``, check out the 'Machine specifics' for help. If you don't want to install or don't have admin rights, head to the 'Post Compilation' section below to finish off your installation.

.. warning:: Even if the code compiled, if your GPU has a compute capability < 5.1, newer versions of ``nvcc`` won't compile code that will work. You'll get error messages like "No kernel image available". Check out how to fix that in 'Machine specifics' below.

Machine specifics
######################
``cmake`` is pretty good at trying to find all the necessary libraries, but every machine is unique, so often you'll need to point ``cmake`` in the correct direction. To that end, I've include 4 keywords: ``JSONC_ROOT``, ``ERFA_ROOT``, ``HDF5_ROOT``, ``PAL_ROOT`` that you can pass to ``cmake``. When passing an option to ``cmake``, you add ``-D`` to the front. For example, on ``OzStar``, I used the command::

  $ cmake ..  -DJSONC_ROOT=/fred/oz048/jline/software/json-c/install/

which tells ``cmake`` to look for ``libjson-c.so`` in paths like ``${JSONC_ROOT}/lib`` or ``${JSONC_ROOT}/lib64``, and ``json.h`` in paths like ``${JSONC_ROOT}/include`` and ``${JSONC_ROOT}/include/json-c``. Read the errors out of ``cmake`` to see which libraries it can't find and add whatever you need to your ``cmake`` command to point to the correct libraries.

.. note:: If you install a dependency in an unusual place on you machine, you have to make sure ``woden`` can find it at run time. So if you compiled with the ``json-c`` library in the ``cmake`` example above, you'd need to call ``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fred/oz048/jline/software/json-c/install/lib64`` before you call ``woden`` (or put that line in your ``~/.bashrc`` or equivalent).

All NVIDIA GPUs have a specific compute capability, which relates to their internal architecture. You can tell the compiler which architecture to compile for, which in theory should make compilation quicker, and ensure the code runs correctly on your GPU. You can find out the compute value here (https://developer.nvidia.com/cuda-gpus), and pass it to CMake via::

  $ cmake .. -DCUDA_ARCH=6.0

(for a compute capability of 6.0, for example).

.. warning:: For newer ``CUDA`` versions, some compute capabilities are deprecated, so the compiler leaves them out by default. For example, using ``CUDA`` version 11.2, compute capabilities 3.5 to 5.0 are ignored. If you card has a compute capability of 5.0, you **must** include the flag ``-DCUDA_ARCH=5.0``, otherwise the `nvcc` compiler will not create an executable capable of running on your device.

If you need to pass extra flags to your CUDA compiler, you can do so by adding something like the following::

  -DCMAKE_CUDA_FLAGS="-Dsomeflag"


Post compilation (required if you don't run ``make install``)
###############################################################

If you don't run ``make install``, ``run_woden.py`` won't be able to find the ``woden`` executable. Default installation locations often need admin privileges. If you can't install to them (or just want to keep ``WODEN`` contained inside a single directory), you can instead just add::

  source /path/to/your/location/WODEN/build/init_WODEN.sh

to your ``~/.bash_rc`` (where you replace ``/path/to/your/location`` to wherever you installed ``WODEN``). This will create the variable ``$WODEN_DIR``, and add it to your ``$PATH``. Furthermore, ``init_WODEN.sh`` is generated by the script ``src/update_init_WODEN.py``, which looks through ``CMakeCache.txt`` for the locations of ``ERFA``, ``HDF5``, ``JSONC``, ``PAL``. It then appends lines to ``init_WODEN.sh`` to add these locations to ``LD_LIBRARY_PATH``, so ``woden`` can find these libraries at run time. For example, on my machine, ``init_WODEN.sh`` ends up looking like::

  ##This line finds the current directory at sets the env variable WODEN_DIR
  export WODEN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  ##This adds the line to PATH
  export PATH=$WODEN_DIR:$PATH
  ##Add library paths to LD_LIBRARY_PATH so the can be found at runtime
  export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial/:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH

.. note:: Every time you run ``make``, ``init_WODEN.sh`` is regenerated, so any edits you make will be overwritten. I suggest any other customisation of you ``LD_LIBRARY_PATH`` happens in your ``~/.bashrc`` or equivalent.

Post compilation (optional)
###############################

If you want to use the MWA FEE primary beam model, you must have the stored spherical harmonic coefficients hdf5 file ``mwa_full_embedded_element_pattern.h5``. You can then define this environment variable in your ``~/.bash_rc``::

  export MWA_FEE_HDF5=/path/to/your/location/mwa_full_embedded_element_pattern.h5

so again ``run_woden.py`` can find it. There is a command line option ``--hdf5_beam_path`` in ``run_woden.py`` which you can use instead of this environment variable if you want.

If you don't have the spherical harmonic file you can obtain it via the command::

  $ wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5
