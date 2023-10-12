#!/bin/sh

mkdir -p images
##See if we can find where the scripts live - annoyingly, have to give an
##absolute path to casa to run uv2ms.py, so need to know where it lives

for beam in "MWA_FEE" "MWA_FEE_interp"
do

  ##Convert uvfits to measurement sets. --nologger stops a logging popup
  ##that can slow things down. Remove it if you want a graphical logger
  woden_uv2ms.py \
      --uvfits_prepend=./data/multi-comp_grid_${beam}_band \
      --band_nums=1,2,3

  ##CLEAN them up
  wsclean -name ./images/multi-comp_grid_${beam} -size 1250 1250 -niter 30000 \
      -auto-threshold 0.5 -auto-mask 3 \
      -pol I -multiscale -weight briggs 0 -scale 0.02 -j 8 -mgain 0.85 \
      -no-update-model-required \
      ./data/multi-comp_grid_${beam}_band*.ms


done

rm casa*.log
rm importuvfits.last
