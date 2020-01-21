# WODEN

Jack Line, 15/01/2020 \
Any issues, questions, wisdom, feel free to email me at jack.line@curtin.edu.au

WODEN is C / CUDA code designed to compare different visibility generating methods. WODEN was written to simulate Murchinson Widefield Array ([MWA, Tingay et al. 2013](https://doi.org/10.1017/pasa.2012.007)) visibilities. WODEN creates a simulated phase tracked set of visibilities for the MWA (currently only the layout, the simulations do _not_ include the primary beam, other instrumental effects, or noise). It generates them in 1.28MHz chunks (bands, like the real MWA).

WODEN reads parameters in from a `.json` file, relying heavily on an MWA metafits file for observational settings, as well as the array layout. You can obtain metafits files from the [MWA ASVO](https://asvo.mwatelescope.org/dashboard). WODEN dumps all data to binary files. This is all currently controlled using `run_woden.py`, which also converts binary outputs into uvfits files, via a uvfits template.

Code borrows the gaussian and shapelet methodology from the RTS [Mitchell et al. 2008](https://ieeexplore.ieee.org/document/4703504?arnumber=4703504 "IEEExplorer"), as well as directly uses the metafits file reader.

WODEN is able to run using shapelet model source catalogues generated with SHApelet Modelling For Interferometers ([SHAMFI](https://github.com/JLBLine/SHAMFI)), specified with the `--woden_srclist` SHAMFI option. It also includes a script to convert a multi-scale CLEAN component list out of [WSClean](https://sourceforge.net/projects/wsclean/) into a WODEN-style srclist (when running WSClean use the `-save-source-list` option).

## 1. Installation
See INSTALL.md for details. Only supported on linux so far, and with only limited installation testing. Requires an NVIDIA GPU.

## 2. Testing
Limited testing of simulating a single point source, gaussian source, and shapelet source are included in the `tests` directory. Using bash, simply run
```sh
source run_all_tests.sh
```
to generate example uvfits files and test your installation. As an added step, if you have `casa` ([download casa](https://casa.nrao.edu/casa_obtaining.shtml)) and `wsclean` ([download WSClean](https://sourceforge.net/p/wsclean/wiki/Installation/)) installed, I have included a script to convert the output uvfits files into measurement sets, and CLEAN them. The script `image_all_tests.sh` assumes `casa` and `wsclean` are available in the command line, so you may need to add `casa` to your `$PATH`. Alternatively, I use an `alias` inside my `~/.bashrc`, e.g. `alias casa='/usr/local/casa-pipeline-release-5.6.2-2.el7/bin/casa --nologger'`. You can then run
```sh
source image_all_tests.sh
```
to image the test simulations and check they make sense given the input source catalogues (see Section 4 for more about source catalogues)

## 3. Usage
Run using the python wrapper script `run_woden.py`, which generates the necessary `.json` input file. The arguments for the script are detailed below (also see `run_woden.py --help`)

Argument  | Type | Description
--|--|--
`--ra0`  | Mandatory  | RA of the desired phase centre (deg)
`--dec0`  | Mandatory  | Dec of the desired phase centre (deg)
`--cat_filename`  | Mandatory  | A WODEN style source catalogue to simulate
`--metafits_filename`  | Mandatory  | Metafits file to base the simulation on. WODEN takes the array layout, LST, observing frequency, and time and frequency resolution from the metafits
`--band_nums`  | Optional | Defaults to running all 24 course bands. Alternatively, enter required numbers delineated by commas, e.g. --band_nums=1,7,9
`--output_uvfits_prepend`  | Optional  | Prepend name for uvfits - will append `band%02d.uvfits %band_num` at the end. Can include the path to another directory.
`--num_freq_channels`  | Optional  | The number of frequency channels within a band to simulate - defaults to 1.28MHz / --freq_res
`--num_time_steps`  | Optional  | The number of time steps to simulate - defaults to how many are in the metafits
`--freq_res`  | Optional | Fine channel frequency resolution (Hz) - will default to what is in the metafits
`--time_res`  | Optional  | Time resolution (s) - will default to what is in the metafits
`--no_tidy`  | Optional  | Defaults to deleting output binary files from WODEN and `.json` files. Add this flag to not delete these files
`--nvprof`  | Optional  | Add to switch on the nvidia profiler when running WODEN and save output profile results to a `.txt` file, in the same location as the output uvfits. This may or may not work depending on how your CUDA installation works
`--template_uvfits`  | Optional | The template uvfits file to base outputs on - defaults to `/build/template_MWA_128T.uvfits`

Alternatively, one can run WODEN directly via the command line, using a `.json` file such as
```sh
woden input_file.json
```
where an `input_file.json` reads as
```json
{
  "ra0": 50.6700000000,
  "dec0": -37.2000000000,
  "num_freqs": 16,
  "num_time_steps": 14,
  "cat_filename": "srclist_point_source.txt",
  "metafits_filename": "1102865128_metafits_ppds.fits",
  "time_res": 8.00000,
  "frequency_resolution": 80000.000,
  "band_nums": [1,2,3]
}
```
with the input values matching the options listed in the table about. This will only create binary files however and not uvfits files. In the future I plan to use [pyuvdata](https://pyuvdata.readthedocs.io/en/latest/) to read directly from the binary files to then generate output products.

### 3.1 Using multi-scale CLEAN components from WSCLEAN
To generate a WODEN-style source catalogue from a WSClean multi-scale CLEAN, simply use the script `convert_WSClean_list_to_WODEN.py`, which will be linked to the build directory upon compilation. An exmaple command is
```sh
convert_WSClean_list_to_WODEN.py \
 --file=msclean_output_from_WSClean-sources.txt \
 --outname=srclist-woden_wsclean-model.txt
```

## 4. WODEN source catalogue
The WODEN source catalogue is a modified version of the RTS srclist. In the `v0.1` release of WODEN, you create one single SOURCE which can include as many COMPONENTS (which can be point, gaussian or shapelet) as you want. In future releases I plan to implement multiple sources and have `--num_sources` as an option in `run_woden.py`.

### 4.1 Basic point source
An example of a single SOURCE with a single point source COMPONENT is below.
```
SOURCE source_name P 1 G 0 S 0 0
COMPONENT POINT 4.0 -27.0
FREQ 1.8e+08 10.0 0 0 0
ENDCOMPONENT
ENDSOURCE
```
An explanation of each line and value:
```
SOURCE source_name P 1 G 0 S 0 0
##Initialises the SOURCE, giving it the name `source_name`, and specifying the
##number and type of components (P = point, G = gaussian, S = shapelet).
##For shapelet, the two numbers are total number of coefficients and total
##number of components
```
```
COMPONENT POINT 4.0 -27.0
##Initialises a component, specifying the type (either POINT, GAUSSIAN, SHAPELET)
##and the RA and DEC (hours, deg). So this line means a point source at
##RA,DEC = 4h, -27deg.
```
```
FREQ 1.8e+08 10.0 0 0 0
##Specifies a reference Stokes flux density as `FREQ Hz I Q U V`. In release v0.1
##all components are given an SI of -0.8, and there is no support for Q,U,V, so
##always set the to zero
```
```
ENDCOMPONENT
##Ends the component
```
```
ENDSOURCE
##Ends the source
```
To add multiple point sources, simply repeat the COMPONENT / ENDCOMPONENT sections with new details, i.e.
```
SOURCE multi_point P 3 G 0 S 0 0
COMPONENT POINT 4.0 -27.0
FREQ 1.8e+08 10.0 0 0 0
ENDCOMPONENT
COMPONENT POINT 3.0 -37.0
FREQ 1.3e+08 1.0.0 0 0 0
ENDCOMPONENT
COMPONENT POINT 5.0 -47.0
FREQ 3.9e+08 0.04 0 0 0
ENDCOMPONENT
ENDSOURCE
```
noting that at the very top line, I have updated `P 3` to reflect there are now three point sources. These numbers are used to quickly allocate memory so things will go wrong if they are wrong.

### 4.2 Gaussian sources
An example srclist containing a single gaussian:
```
SOURCE gaussian_source P 0 G 1 S 0 0
COMPONENT GAUSSIAN 3.378 -37.2
FREQ 1.8e+08 10.0 0 0 0
GPARAMS 45.0000000000 6.0 3.0
ENDCOMPONENT
ENDSOURCE
```
where all lines have the same meaning as the point source, and the meaning of the extra line:
```
GPARAMS 45.0000000000 6.0 3.0
##Specify the gaussian parameters as `GPARAMS pa(deg) major_axis(arcmin) minor_axis(arcmin)`
##The major and minor axes are specified as FWHM
##Note this line needs to sit in between the lines starting with
##`COMPONENT GAUSSIAN` and `ENDCOMPONENT`
```

### 4.3 Shapelet sources
To generate shapelet models compatible with WODEN, simply use [SHAMFI](https://github.com/JLBLine/SHAMFI) to fit an image with the `--woden_srclist` option. This will ensure all normalisations are correct. An example srclist (made by hand so the normalisations _won't_ be correct) is:
```
SOURCE shapelet_source P 0 G 0 S 1 3
COMPONENT SHAPELET 3.378 -37.2
FREQ 1.8e+08 10.0 0 0 0
SPARAMS 45.0000000000 6.0 3.0
SCOEFF 0 0 0.92342
SCOEFF 1 10 0.0002354
SCOEFF 4 5 0.004567
ENDCOMPONENT
ENDSOURCE
```
which generates a single shapelet component, including 3 coffecients, hence `S 1 3` in the first line. The `SPARAMS` line again details `SPARAMS pa(deg) major_axis(arcmin) minor_axis(arcmin)` similarly to the gaussian source. The extra lines detail:
```
SCOEFF 0 0 0.92342
##Details the order of the shapelet basis function (see Line et al. 2020
##for details, link to come once published) and fitted coefficient as
##`SCOEFF p1 p2 coeff`
```
You can add as many `SCOEFF` lines as necessary, with a maximum order < 100. If you use SHAMFI, the coefficients will be scaled such that the I flux density of the full source will be 10Jy at 180MHz.

### 4.4 Putting it all together
An example srclist with all component types would look something like this:
```
SOURCE multi_sources P 3 G 1 S 2 7
COMPONENT SHAPELET 3.378 -37.2
FREQ 1.8e+08 10.0 0 0 0
SPARAMS 45.0000000000 6.0 3.0
SCOEFF 0 0 0.92342
SCOEFF 1 10 0.0002354
SCOEFF 4 5 0.004567
ENDCOMPONENT
COMPONENT SHAPELET 3.12 -32.2
FREQ 1.8e+08 3.1 0 0 0
SPARAMS 56.0000000000 9.0 3.0
SCOEFF 0 0 0.02345
SCOEFF 3 0 -0.234234
SCOEFF 21 34 0.82342
SCOEFF 31 5 -0.00876234
ENDCOMPONENT
COMPONENT GAUSSIAN 3.378 -37.2
FREQ 1.8e+08 10.0 0 0 0
GPARAMS 45.0000000000 6.0 3.0
ENDCOMPONENT
COMPONENT POINT 4.0 -27.0
FREQ 1.8e+08 10.0 0 0 0
ENDCOMPONENT
COMPONENT POINT 3.0 -37.0
FREQ 1.3e+08 1.0.0 0 0 0
ENDCOMPONENT
COMPONENT POINT 5.0 -47.0
FREQ 3.9e+08 0.04 0 0 0
ENDCOMPONENT
ENDSOURCE
```
