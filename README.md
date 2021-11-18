# WODEN

[![status](https://joss.theoj.org/papers/bbc90ec4cd925ade93ed0781e571d247/status.svg)](https://joss.theoj.org/papers/bbc90ec4cd925ade93ed0781e571d247)
[![codecov](https://codecov.io/gh/JLBLine/WODEN/branch/joss_review/graph/badge.svg?token=Q3JFCI5GOC)](https://codecov.io/gh/JLBLine/WODEN) _*note code coverage only applies to `python3` and `C` code, `CUDA` code is not currently supported by code coverage software_

> The `WODEN` documentation lives [here on readthedocs](https://woden.readthedocs.io/en/latest/). If your internet has broken and you have already installed `WODEN`, you can build a local copy by navigating into `WODEN/docs/sphinx` and running `make html` (you'll need to have `doxygen` installed).

`WODEN` is C / CUDA code designed to be able to simulate low-frequency radio interferometric data. It is written to be simplistic and *fast* to allow all-sky simulations. Although `WODEN` was primarily written to simulate Murchinson Widefield Array ([MWA, Tingay et al. 2013](https://doi.org/10.1017/pasa.2012.007)) visibilities, it is becoming less instrument-specific as time goes on. `WODEN` outputs `uvfits` files.

The unique part of `WODEN` is that it can simulate shapelet model sources (along with point and Gaussian) that are compatible with the `RTS` ([Mitchell et al. 2008](https://ieeexplore.ieee.org/document/4703504?arnumber=4703504 "IEEExplorer")). These models are generated with SHApelet Modelling For Interferometers ([SHAMFI](https://github.com/JLBLine/SHAMFI)), specified with the `--woden_srclist` SHAMFI option. It also includes a script to convert a multi-scale CLEAN component list out of [WSClean](https://sourceforge.net/projects/wsclean/) into a `WODEN`-style srclist (when running `WSClean` use the `-save-source-list` option). `WODEN` can also produce visibilities that can be fed directly into the `RTS` to allow testing of calibration and modelling methodologies.

If you have feature requests or want to contribute to `WODEN`, have a read of the
[guide to contributing](CONTRIBUTION_GUIDE.md) to get started. I welcome the feedback and/or help!

Jack Line \
November 2021


## 1. Installation
Read the comprehensive [installation guide on readthedocs](https://woden.readthedocs.io/en/latest/installation/installation.html#dependencies). In short, you will need the dependencies:

- CMake - https://cmake.org version >= 3.10
- NVIDIA CUDA - https://developer.nvidia.com/cuda-downloads
- json-c - https://github.com/json-c/json-c
- ERFA - https://github.com/liberfa/erfa/releases
- HDF5 - https://www.hdfgroup.org/downloads/hdf5/
- PAL - https://github.com/Starlink/pal/releases
- python >= 3.6

Once you have those, installation is done via `CMake`. Ideally, this will be enough:
```bash
$ cd WODEN
$ mkdir build && cd build
$ cmake ..
$ make -j 4
$ sudo make install #(this is optional)
```
with a couple of post-compilation environment variables needed (if you don't want to `make install`). Checkout the [installation guide on readthedocs](https://woden.readthedocs.io/en/latest/installation/installation.html#dependencies) for full details.

## 2. Testing
There are two routes to test WODEN:
- via the unit/integration tests run by `ctest`, which will test your dependencies (see [cmake testing on readthedocs](https://woden.readthedocs.io/en/latest/testing/cmake_testing.html))
- via simple example scripts, which will test your installation (see [testing via scripts on readthedocs](https://woden.readthedocs.io/en/latest/testing/script_testing.html))

## 3. Usage

### 3.1 Python wrapper

The recommended way to run `WODEN` is via the wrapper script `run_woden.py`, which will launch the `woden` executable. All [ `run_woden.py` arguments are explained here on readthedocs](https://woden.readthedocs.io/en/latest/API_reference/python_code/run_woden.html) (or you can type `run_woden.py --help`).

As `WODEN` is written primarily for the MWA, the simplest way to feed in observational properties is via a `metafits` file. You can [obtain metafits files here](https://asvo.mwatelescope.org/). A minimalistic example command using a metafits file looks like

```bash
run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.txt \
    --metafits_filename=1202815152_metafits_ppds.fits \
    --primary_beam=MWA_FEE
```
where the array layout, observational settings, frequency and time resolution are all read in from the metafits file. All you have to specify is a phase centre (`ra0`, `dec0`), a sky model (`cat_filename`), and a primary beam (in this case the MWA Fully Embedded Element beam - the delays are also read from the metafits). For full control, you can also specify all observational settings explicitly:

```bash
run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.txt \
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

### 3.2 Direct command line

Alternatively, one can run WODEN directly via the command line, using a `.json` file such as
```sh
woden input_file.json
```
where an `input_file.json` reads as
```json
{
  "ra0": 50.6700000000,
  "dec0": -37.2000000000,
  "num_freqs": 128,
  "num_time_steps": 240,
  "cat_filename": "srclist_msclean_fornaxA_phase1+2.txt",
  "time_res": 0.50000,
  "frequency_resolution": 10000.000,
  "chunking_size": 0,
  "jd_date": 2458165.9714583335444331,
  "LST": 72.79744734,
  "array_layout": "WODEN_array_layout.txt",
  "lowest_channel_freq": 1.6959500000e+08,
  "latitude": -26.70331944,
  "coarse_band_width": 1.2800000000e+06,
  "band_nums": [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
}
```
with the input values matching the options listed in the table above. This will only create binary output files however, and not `uvfits` files.

### 3.3 Using multi-scale CLEAN components from WSCLEAN
To generate a `WODEN`-style source catalogue from a `WSClean` multi-scale CLEAN, simply use the script `convert_WSClean_list_to_WODEN.py`, which will be linked to the build directory upon compilation. An exmaple command is
```sh
convert_WSClean_list_to_WODEN.py \
 --file=msclean_output_from_WSClean-sources.txt \
 --outname=srclist-woden_wsclean-model.txt
```

## 4. WODEN sky model
The WODEN sky model uses point sources, elliptical Gaussians, and shapelet models. All sources are currently given a power-law spectral behaviour. A full [breakdown of the sky model lives here on readthedocs](https://woden.readthedocs.io/en/latest/operating_principles/skymodel.html). As an idea, here is a simple sky model of thee point sources:

```
SOURCE multi_point P 3 G 0 S 0 0
COMPONENT POINT 4.0 -27.0
LINEAR 1.8e+08 10.0 0 0 0 -0.4
ENDCOMPONENT
COMPONENT POINT 3.0 -37.0
LINEAR 1.3e+08 1.0.0 0 0 0 -0.786
ENDCOMPONENT
COMPONENT POINT 5.0 -47.0
LINEAR 3.9e+08 0.04 0 0 0 .02
ENDCOMPONENT
ENDSOURCE
```
