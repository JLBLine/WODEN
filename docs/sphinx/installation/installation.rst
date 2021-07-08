Installation
=============

WODEN is built on CUDA so you will need an NVIDIA GPU to run it. Currently, WODEN has only been tested and run on linux, specifically Ubuntu 16.04 up to 20.04, the OzStar super cluster of Swinburne University, and Garrawarla cluster of Pawsey. If you're mad keen to run on Windows or Mac, please contact Jack at jack.l.b.line@gmail.com and we can give it a go.

Dependencies
-------------

``WODEN`` has a number of dependencies so it doesn't reinvent the wheel. I've linked detailed instructions on how I installed them on Ubuntu 20.04 here, and listed them below.

.. toctree::
  :maxdepth: 2

  dependencies

- **CMake** - https://cmake.org version >= 3.10
- **NVIDIA CUDA** - https://developer.nvidia.com/cuda-downloads
- **json-c** - https://github.com/json-c/json-c
- **ERFA** - https://github.com/liberfa/erfa/releases
- **HDF5** - https://www.hdfgroup.org/downloads/hdf5/
- **PAL** - https://github.com/Starlink/pal/releases
- **python >= 3.6**

Compiling ``WODEN``
---------------------

In an ideal world, if the installation of your dependencies went perfectly and
you have a newer NVIDIA GPU, you should be able to simply run::

  $ git clone https://github.com/JLBLine/WODEN.git
  $ cd WODEN
  $ mkdir build && cd build
  $ cmake ..
  $ make -j 4

et voila, your code is compiled. If this worked, head to the 'Post Compilation' section below to finish off your installation. If it didn't, and you're not used to ``cmake``, check out the 'Machine specifics' for help.

.. warning:: Even if the code compiled, if your GPU has a compute capability < 5.1, newer versions of ``nvcc`` won't compile code that will work. You'll get error messages like "No kernel image available". Check out how to fix that in 'Machine specifics' below.

Machine specifics
------------------
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


Post compilation (required)
----------------------------

Once compiled, just add::

  source /path/to/your/location/WODEN/build/init_WODEN.sh

to your ``~/.bash_rc`` (where you replace ``/path/to/your/location`` to wherever you installed ``WODEN``). This will create the variable ``$WODEN_DIR``, and add it to your ``$PATH``. This allows ``run_woden.py`` to find the ``woden`` executable.

Post compilation (optional)
----------------------------

If you want to use the MWA FEE primary beam model, you must have the stored spherical harmonic coefficients hdf5 file ``mwa_full_embedded_element_pattern.h5``. You can then define this environment variable in your ``~/.bash_rc``::

  export MWA_FEE_HDF5=/path/to/your/location/mwa_full_embedded_element_pattern.h5

so again ``run_woden.py`` can find it. There is a command line option ``--hdf5_beam_path`` in ``run_woden.py`` which you can use instead of this environment variable if you want.

If you don't have the spherical harmonic file you can obtain it via the command::

  $ wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5
