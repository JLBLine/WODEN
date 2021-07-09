#!/bin/sh

mkdir -p data

for name in "point" "gauss" "shapelet"
do

  run_woden.py \
     --ra0=50.67 --dec0=-37.2 \
     --num_freq_channels=5 --num_time_steps=5 \
     --freq_res=80e+3 --time_res=8.0 \
     --cat_filename=../skymodels/srclist_${name}_source.txt \
     --metafits_filename=../metafits/1102865128_metafits_ppds.fits \
     --band_nums=1,2 --output_uvfits_prepend=./data/single_${name} \
     --primary_beam=None

done
