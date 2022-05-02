#!/bin/sh

mkdir -p data

for beam in "None" "MWA_analy" "Gaussian"
do

  run_woden.py \
     --ra0=60.0 --dec0=-40.0 \
     --num_freq_channels=10 --num_time_steps=10 \
     --freq_res=80e+3 --time_res=8.0 \
     --cat_filename=../skymodels/srclist_multi-comp_grid.txt \
     --metafits_filename=../metafits/1102865128_metafits_ppds.fits \
     --band_nums=1,2,3 --output_uvfits_prepend=./data/multi-comp_grid_${beam} \
     --primary_beam=${beam} \
     --precision=float \

done

run_woden.py \
   --ra0=60.0 --dec0=-40.0 \
   --num_freq_channels=5 --num_time_steps=5 \
   --freq_res=80e+3 --time_res=8.0 \
   --cat_filename=../skymodels/srclist_multi-comp_grid.txt \
   --metafits_filename=../metafits/1102865128_metafits_ppds.fits \
   --band_nums=1,2,3 --output_uvfits_prepend=./data/multi-comp_grid_EDA2 \
   --primary_beam=EDA2 \
   --array_layout=../array_layouts/EDA2_layout_255.txt \
   --precision=float \
   --chunking_size=1e9
