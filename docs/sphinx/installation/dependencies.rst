How to install dependencies
============================

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
+ **PAL** - https://github.com/Starlink/pal/releases - ``PAL`` is a little mental with it's default installation paths. I *HIGHLY* recommmend downloading a release version, and then using the ``--without-starlink`` option::

  $ wget https://github.com/Starlink/pal/releases/download/v0.9.8/pal-0.9.8.tar.gz
  $ tar -xvf pal-0.9.8.tar.gz
  $ cd pal-0.9.8
  $ ./configure --prefix=/usr/local --without-starlink
  $ make
  $ sudo make install

  Doing it this way installs things in normal locations, making life easier during linking.
+ **python >= 3.6** - the best way to run ``WODEN`` is through the script ``run_woden.py``, which has a number of package dependencies. One of these is ``pyerfa``, which uses f-strings during installation, so you have to use a python version >= 3.6. Sorry. The requirements can be found in ``WODEN/docs/sphinx/sphinx/requirements.txt``, which you can install via something like::

  $ pip3 install -r requirements.txt

For completeness, those packages are::

  sphinx_argparse
  breathe
  astropy
  numpy
  pyerfa

The first two packages are used for the documentation.

Phew! That's it for now.
