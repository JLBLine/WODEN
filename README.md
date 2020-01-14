# WODEN

Jack Line, 08/01/2020

CUDA code designed to compare different visibility generating methods.

Reads params in from `.json` file, relying heavily on an MWA metafits file for observational settings. Dumps all data to a binary file output_visi_band%02d.dat. Currently controlled using `run_woden.py`, which also converts binary outputs into uvfits files, via a uvfits template.

Code borrows the gaussian and shapelet methodology from the RTS [Mitchell et al. 2008](https://ieeexplore.ieee.org/document/4703504?arnumber=4703504 "IEEExplorer"), as well as the metafits file reader.
