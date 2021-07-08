#!/bin/sh

mkdir -p data

for beam in "None" "Gaussian"
do

  run_woden.py \
     --ra0=50.67 --dec0=-37.2 \
     --num_freq_channels=5 --num_time_steps=5 \
     --freq_res=80e+3 --time_res=8.0 \
     --cat_filename=../skymodels/srclist_point_source.txt \
     --metafits_filename=../metafits/1102865128_metafits_ppds.fits \
     --band_nums=1,2 --output_uvfits_prepend=./data/point_${beam} \
     --primary_beam=${beam}

  run_woden.py \
     --ra0=60.0 --dec0=-40.0 \
     --num_freq_channels=5 --num_time_steps=5 \
     --freq_res=80e+3 --time_res=8.0 \
     --cat_filename=../skymodels/srclist_point_source_grid.txt \
     --metafits_filename=../metafits/1102865128_metafits_ppds.fits \
     --band_nums=1,2 --output_uvfits_prepend=./data/point_grid_${beam} \
     --primary_beam=${beam}

done

beam="EDA2"

run_woden.py \
   --ra0=50.67 --dec0=-37.2 \
   --num_freq_channels=5 --num_time_steps=5 \
   --freq_res=80e+3 --time_res=8.0 \
   --cat_filename=../skymodels/srclist_point_source.txt \
   --metafits_filename=../metafits/1102865128_metafits_ppds.fits \
   --band_nums=1,2 --output_uvfits_prepend=./data/point_${beam} \
   --primary_beam=${beam} \
   --array_layout=../array_layouts/EDA2_layout_255.txt

run_woden.py \
   --ra0=60.0 --dec0=-40.0 \
   --num_freq_channels=5 --num_time_steps=5 \
   --freq_res=80e+3 --time_res=8.0 \
   --cat_filename=../skymodels/srclist_point_source_grid.txt \
   --metafits_filename=../metafits/1102865128_metafits_ppds.fits \
   --band_nums=1,2 --output_uvfits_prepend=./data/point_grid_${beam} \
   --primary_beam=${beam} \
   --array_layout=../array_layouts/EDA2_layout_255.txt
