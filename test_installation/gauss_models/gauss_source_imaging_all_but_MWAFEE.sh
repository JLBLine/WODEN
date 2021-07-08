#!/bin/sh

mkdir -p images

for name in "gauss" "gauss_grid"
  do
  for beam in "None" "Gaussian" "EDA2"
  do
    ##Convert uvfits to measurement sets. --nologger stops a logging popup
    ##that can slow things down. Remove it if you want a graphical logger
    ${CASA_DIR}/casa --nologger -c  ${WODEN_DIR}/uv2ms.py \
        --uvfits_prepend=./data/${name}_${beam}_band \
        --band_nums=1,2

    ##CLEAN them up
    wsclean -name ./images/${name}_${beam} -size 750 750 -niter 10000 \
        -auto-threshold 0.5 -auto-mask 3 \
        -pol I -multiscale -weight briggs 0 -scale 0.007 -j 8 -mgain 0.85 \
        -no-update-model-required \
        ./data/${name}_${beam}_band*.ms
  done
done

rm casa*.log
rm importuvfits.last
