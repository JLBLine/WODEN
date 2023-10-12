#!/bin/sh

mkdir -p images

for name in "point" "gauss" "shapelet"
do
  ##Convert uvfits to measurement sets. --nologger stops a logging popup
  ##that can slow things down. Remove it if you want a graphical logger
  woden_uv2ms.py \
      --uvfits_prepend=./data/single_${name}_band \
      --band_nums=1,2

  ##CLEAN them up
  wsclean -name ./images/single_${name} -size 150 150 -niter 10000 \
      -auto-threshold 0.5 -auto-mask 3 \
      -pol I -multiscale -weight briggs 0 -scale 0.007 -j 8 -mgain 0.85 \
      -no-update-model-required \
      ./data/single_${name}_band*.ms
done

rm casa*.log
rm importuvfits.last
