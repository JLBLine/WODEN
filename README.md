# WODEN

Jack Line, 17/09/2019

CUDA code designed to compare different visibility generating methods.

Reads params in from .json file, relying heavily on metafits files for observational settings. Dumps all data to a binary file output_visi_band%02d.dat. Needs a new-style RTS srclist. Will add documentation once options are added / if this is going to be published.

Currently controlled using `run_woden.py`, which also converts binary outputs into uvfits files, via a uvfits template.

Code borrows heavily / sometimes outright steals functions from the RTS.
