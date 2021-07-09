#!/bin/sh

mkdir -p data

for name in "point" "gauss" "shapelet"
do

  run_woden.py \
     --ra0=60.0 --dec0=-40.0 \
     --num_freq_channels=5 --num_time_steps=5 \
     --freq_res=80e+3 --time_res=8.0 \
     --cat_filename=../skymodels/srclist_${name}_source_grid.txt \
     --metafits_filename=../metafits/1102865128_metafits_ppds.fits \
     --band_nums=1,2 --output_uvfits_prepend=./data/grid_${name} \
     --primary_beam=None

done
