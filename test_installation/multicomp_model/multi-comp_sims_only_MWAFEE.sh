#!/bin/sh

mkdir -p data

run_woden.py \
   --ra0=60.0 --dec0=-40.0 \
   --num_freq_channels=5 --num_time_steps=5 \
   --freq_res=80e+3 --time_res=8.0 \
   --cat_filename=../skymodels/srclist_multi-comp_grid.txt \
   --metafits_filename=../metafits/1102865128_metafits_ppds.fits \
   --band_nums=1,2 --output_uvfits_prepend=./data/multi-comp_grid_MWA_FEE \
   --primary_beam=MWA_FEE
