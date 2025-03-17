# WODEN

[![status](https://joss.theoj.org/papers/bbc90ec4cd925ade93ed0781e571d247/status.svg)](https://joss.theoj.org/papers/bbc90ec4cd925ade93ed0781e571d247)

[![Documentation Status](https://readthedocs.org/projects/woden/badge/?version=latest)](https://woden.readthedocs.io/en/latest/?badge=latest) [![codecov](https://codecov.io/gh/JLBLine/WODEN/branch/master/graph/badge.svg?token=Q3JFCI5GOC)](https://codecov.io/gh/JLBLine/WODEN) _*note code coverage only applies to `python3`, `C`, `C++` code, as `CUDA` code is not currently supported by free code coverage software. However, every GPU function now has a CPU equivalent which is tested, so this number is indicative of the entire package_

`WODEN` is Python / C / C++ / GPU code designed to be able to simulate low-frequency radio interferometric data. It is written to be simplistic and *fast* to allow all-sky simulations. Although `WODEN` was primarily written to simulate Murchinson Widefield Array ([MWA, Tingay et al. 2013](https://doi.org/10.1017/pasa.2012.007)) visibilities, it is capable of bespoke array layouts, and has a number of primary beam options. Currently, the ``MWA`` primary beam is supported via [hyperbeam](https://github.com/MWATelescope/mwa_hyperbeam); ``LOFAR`` and ``OSKAR``-style primary beams are supported via the [EveryBeam](https://everybeam.readthedocs.io/en/latest/index.html) package; ``EDA2``, an analytic ``MWA`` beam, and a simple Gaussian beam are supported natively to ``WODEN``.

`WODEN` outputs `uvfits` files. The `WODEN` documentation lives [here on readthedocs](https://woden.readthedocs.io/en/latest/). If your internet has broken and you have already installed `WODEN`, you can build a local copy by navigating into `WODEN/docs/sphinx` and running `make html` (you'll need to have `doxygen` installed).

> **Note**
> Jul 2024: Jack is back on a temporary contract to do some `WODEN` development. Keep your eyes peeled for new features!

<!-- > Although ``WODEN`` is still very much a great tool for MWA interferomteric simulations, I (Jack Line) am no longer working in astronomy. I'll drop in to advise and/or fix bugs from time to time, but I can't commit to developing new features I'm afraid. If you want new features (or even better, want to write new features), please reach out to the [Epoch of Reionisation group at Curtin University](https://astronomy.curtin.edu.au/research/epoch-of-reionisation/). They should know if anyone is actively working on/using the software. If you end up taking over this project, feel free to delete this message! -->

> New in version 2.5-alpha:
 - CPU-only mode via `--cpu_mode`! You should use this when running with an `EveryBeam` primary beam, as around 70% of all computation time is spent in the primary beam calculation. On a cluster, the GPU just sits there doing nothing but costing you and environment.
 - `EveryBeam` is now called directly from `C++` rather than the `Python` wrapper. It's a good deal faster that way.

> New in version 2.4-alpha:
 - You can now use EveryBeam primary beams in `WODEN`. This allows you to simulate `LOFAR` and `SKA` like instruments.
 - As a consequence of the above, the sky model reading part of `WODEN` is now multi-threaded. If you want exact behaviour to match, you should run `WODEN` with `--num_threads=1`. 
 - See the [EveryBeam testing and developing](https://woden.readthedocs.io/en/everybeam/latest/everybeam_testing.html) of the documentation for more information.
 - You don't have to install `EveryBeam` to run `WODEN` if you're not going to use it, so you can ignore this feature if you're only interested in MWA simulations. `WODEN` will just check if the `EveryBeam` Python library is installed and use it if it is.

> New in version 2.3-alpha:
 - Full polarisation is back! You can now simulate QUV visibilities, and they can be independent to Stokes I. In fact, you can set Stokes I to zero in the sky model, and still enter QUV information.
 - Stokes V can either be a power-law, curved power-law, be a polarisation fraction of Stokes I, or be a list of flux densities.
 - Linear polarisation (Q and U) can either be a power-law, curved power-law, be a polarisation fraction of Stokes I, or be a list of flux densities. The first three are always used in conjunction with a rotation measure. List-type fluxes can be given separately for Q and U, which makes Q and U independent of each other. Alternatively, a single list of linear polarisation fluxes can be given, which is used in conjunction with a rotation measure to determine Q and U.

> New in version 2.2-alpha:
 - Thanks to Marcin Sokolowski, and code from the [PaCER Blink project](https://github.com/PaCER-BLINK-Project), there is nominal support for AMD GPUs via ROCm. I've managed to compile and run on the Pawsey Setonix cluster, but this is the only AMD GPU I have access to. If you have an AMD GPU, please test and let me know how you go!
 - There are some installation examples in ``WODEN/templates/`` for a CUDA and an AMD installation on superclusters. Hopefully you can adapt these to your needs if the Docker images don't work for you.
 - CUDA arch flag on compilation is changed to setting an environment variable `CUDAARCHS` to use `CMakes` `CMAKE_CUDA_ARCHITECTURES` feature.

> New in version 2.1:
 - Dipole flagging for the MWA FEE beam is now enabled (`--use_MWA_dipflags`)
 - Dipole amplitudes for the MWA FEE beam is now enabled (`--use_MWA_dipamps`). Needs a modified metafits file
 - Things are still all Stokes I; plan is to implement QUV in version 2.3

> New in version 2.0:
 - Large swathes of the code have been transferred from `C` in `python`
 - The `C/CUDA` code that remains is called directly from `python` now, meaning the `.json` and `.dat` files are no longer needed
 - A new FITS sky model format is supported, meaning you can run simulations with greater than 25 million components via lazy loading
 - I've done away with using Stokes QUV to speed things up. The plan is to implement some kind of rotation measure model to include Q/U, and reinstate Stokes V in a future release

Please be aware, before ``WODEN`` version 1.4.0, in the output `uvfits` files, the first polarisation (usually called XX) was derived from North-South dipoles, as is the labelling convention according to the IAU. However, most `uvfits` users I've met, as well as the data out of the MWA telescope, define XX as East-West. So although the internal labelling and mathematics within the C/GPU code is to IAU spec, by default, ``run_woden.py`` now writes out XX as East-West and YY as North-South. From version 1.4.0, a header value of ``IAUORDER = F`` will appear, with ``F`` meaning IAU ordering is False, so the polarisations go EW-EW, NS-NS, EW-NS, NS-EW. If ``IAUORDER = T``, the order is NS-NS, EW-EW, NS-EW, EW-NS. If there is no ``IAUORDER`` at all, assume ``IAUORDER = T``.

If you have feature requests or want to contribute to `WODEN`, have a read of the
[guide to contributing](CONTRIBUTION_GUIDE.md) to get started. I welcome the feedback and/or help!

Jack Line \
July 2024


## 1. Installation
Read the comprehensive [installation guide on readthedocs](https://woden.readthedocs.io/en/latest/installation/installation.html#dependencies). 

The quickest way is to use a docker image (for a CUDA arch of 7.5 for example):

```
docker pull docker://jlbline/woden-2.3:cuda-75
```

and then run things through Docker or Singularity (more detail on that on [ReadTheDocs](https://woden.readthedocs.io/en/latest/installation/installation.html#dependencies)). These versions come bundled with all the MWA FEE primary beam files which is nice.

Docker versions tested to run on various Australian superclusters are listed below:

- `jlbline/woden-2.3:cuda-60` - tested on Swinburne OzStar
- `jlbline/woden-2.3:cuda-80` - tested on Swinburne Ngarrgu Tindebeek
- `jlbline/woden-2.3:setonix` - tested on Pawsey Setonix (this is *very* experimental)

To get the most control, install yourself. Again, go to [ReadTheDocs](https://woden.readthedocs.io/en/latest/installation/installation.html#dependencies) for details, but in short, you will need the dependencies:

- CMake - https://cmake.org version >= 3.21
- Either NVIDIA CUDA - https://developer.nvidia.com/cuda-downloads
- or AMD ROCm - https://rocm.docs.amd.com/projects/install-on-linux/en/latest/
- rust - https://www.rust-lang.org/tools/install
- mwa_hyperbeam - https://github.com/MWATelescope/mwa_hyperbeam
- python >= 3.8


Once you have those, compilation is done via `CMake`. Ideally, this will be enough:
```bash
$ cd WODEN
$ mkdir build && cd build
$ cmake ..
$ make -j 4
```
You then install the `python3` code via something like
```bash
$ cd WODEN
$ pip install -r requirements.txt
$ pip install .
```
with a couple of post-compilation environment variables needed to use the MWA FEE primary beam model. See the full installation guide for more details, including how to use the `EveryBeam` functionality.

## 2. Testing
There are two routes to test WODEN:
- via the unit/integration tests run by `ctest`, which will test your dependencies (see [cmake testing on readthedocs](https://woden.readthedocs.io/en/latest/testing/cmake_testing.html))
- via simple example scripts, which will test your installation (see [testing via scripts on readthedocs](https://woden.readthedocs.io/en/latest/testing/script_testing.html))

## 3. Usage

`WODEN` is run via the script `run_woden.py`, which will launch the `woden` executable. All [ `run_woden.py` arguments are explained here on readthedocs](https://woden.readthedocs.io/en/latest/API_reference/python_code/run_woden.html) (or you can type `run_woden.py --help`).

As `WODEN` is written primarily for the MWA, the simplest way to feed in observational properties is via a `metafits` file. You can [obtain metafits files here](https://asvo.mwatelescope.org/). A _very_ minimalistic example command using a metafits file looks like

```bash
run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.fits \
    --metafits_filename=1202815152_metafits_ppds.fits \
    --primary_beam=MWA_FEE
```
where the array layout, observational settings, frequency and time resolution are all read in from the metafits file. All you have to specify is a phase centre (`ra0`, `dec0`), a sky model (`cat_filename`), and a primary beam (in this case the MWA Fully Embedded Element beam - the delays are also read from the metafits). For full control, you can also specify all observational settings explicitly:

```bash
run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.fits \
    --primary_beam=MWA_FEE \
    --MWA_FEE_delays=[6,4,2,0,8,6,4,2,10,8,6,4,12,10,8,6] \
    --lowest_channel_freq=169.6e+6 \
    --freq_res=10e+3 \
    --num_time_steps=240 \
    --time_res=0.5 \
    --date=2018-02-28T08:47:06 \
    --array_layout=MWA_phase2_extended.txt
```

`WODEN` can read in an array layout as specified in local east, north, and height, and can simulate for any location on Earth (defaults to the MWA but you can specify Earth location with `--longitude` / `--latitude`).

## 4. LoBES/hyperdrive/WODEN sky model
This sky model uses point sources, elliptical Gaussians, and shapelet models. All sources can either have power-law, curved power-law, or list of flux densities for Stokes I spectral behaviour. It can optionally also have a polarisation fraction, power-law, or curved power-law for linear (Q and U) and circular (V) polarisation.

 A full [breakdown of the sky model lives here on readthedocs](https://woden.readthedocs.io/en/latest/operating_principles/skymodel.html), along with the possible formats to feed into `WODEN`.
