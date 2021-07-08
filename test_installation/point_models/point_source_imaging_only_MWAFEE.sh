#!/bin/sh

mkdir -p images

for name in "point" "point_grid"
  do
  ##Convert uvfits to measurement sets. --nologger stops a logging popup
  ##that can slow things down. Remove it if you want a graphical logger
  ${CASA_DIR}/casa --nologger -c  ${WODEN_DIR}/uv2ms.py \
      --uvfits_prepend=./data/${name}_MWA_FEE_band \
      --band_nums=1,2

  ##CLEAN them up
  wsclean -name ./images/${name}_MWA_FEE -size 750 750 -niter 10000 \
      -auto-threshold 0.5 -auto-mask 3 \
      -pol I -multiscale -weight briggs 0 -scale 0.007 -j 8 -mgain 0.85 \
      -no-update-model-required \
      ./data/${name}_MWA_FEE_band*.ms
done

rm casa*.log
rm importuvfits.last
