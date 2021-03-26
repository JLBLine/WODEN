#!/bin/sh

mkdir -p data

#run_woden.py \
#    --ra0=50.67 --dec0=-37.2 \
#    --num_freq_channels=16 --num_time_steps=14 \
#    --freq_res=80e+3 --time_res=8.0 \
#    --cat_filename=srclist_point_source.txt \
#    --metafits_filename=1102865128_metafits_ppds.fits \
#    --band_nums=1,2,3 --output_uvfits_prepend=./data/point_zenith \
#    --primary_beam=Gaussian \
#    --latitude=20.0


#run_woden.py \
#    --ra0=50.67 --dec0=20 \
#    --num_freq_channels=16 --num_time_steps=14 \
#    --freq_res=80e+3 --time_res=8.0 \
#    --cat_filename=srclist_point_source_lat20.txt \
#    --metafits_filename=1102865128_metafits_ppds.fits \
#    --band_nums=1,2,3 --output_uvfits_prepend=./data/point_zenith \
#    --primary_beam=Gaussian \
#    --gauss_ra_point=50.67 \
#    --gauss_dec_point=20.0 \
#    --latitude=20.0 \
#    --no_tidy

#run_woden.py \
#    --ra0=50.67 --dec0=20 \
#    --num_freq_channels=16 --num_time_steps=14 \
#    --freq_res=80e+3 --time_res=8.0 \
#    --cat_filename=srclist_point_source_lat20.txt \
#    --metafits_filename=1102865128_metafits_ppds.fits \
#    --band_nums=1,2,3 --output_uvfits_prepend=./data/point_zenith \
#    --primary_beam=MWA_FEE \
#    --latitude=20.0 \
#    --no_tidy
    
#run_woden.py \
#    --ra0=50.67 --dec0=-37.2 \
#    --num_freq_channels=16 --num_time_steps=14 \
#    --freq_res=80e+3 --time_res=8.0 \
#    --cat_filename=srclist_point_source.txt \
#    --metafits_filename=1102865128_metafits_ppds.fits \
#    --band_nums=1,2,3 --output_uvfits_prepend=./data/point_zenith \
#    --primary_beam=MWA_FEE \
#    --no_tidy

run_woden.py \
    --ra0=56.0 --dec0=-39.0 \
    --num_freq_channels=16 --num_time_steps=14 \
    --freq_res=80e+3 --time_res=8.0 \
    --cat_filename=srclist_grid.txt \
    --metafits_filename=1102865128_metafits_ppds.fits \
    --band_nums=1,2,3 --output_uvfits_prepend=./data/grid_pointing \
    --primary_beam=MWA_FEE \
    --no_tidy \
    --chunking_size=9000
