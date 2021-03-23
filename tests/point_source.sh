#!/bin/sh

mkdir -p data

run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --num_freq_channels=16 --num_time_steps=14 \
    --freq_res=80e+3 --time_res=8.0 \
    --cat_filename=srclist_point_source.txt \
    --metafits_filename=1102865128_metafits_ppds.fits \
    --band_nums=1,2,3 --output_uvfits_prepend=./data/point_zenith \
    --primary_beam=none
