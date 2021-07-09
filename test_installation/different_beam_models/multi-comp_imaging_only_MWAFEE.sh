#!/bin/sh

mkdir -p images

##Convert uvfits to measurement sets. --nologger stops a logging popup
##that can slow things down. Remove it if you want a graphical logger
${CASA_DIR}/casa --nologger -c  ${WODEN_DIR}/uv2ms.py \
    --uvfits_prepend=./data/multi-comp_grid_MWA_FEE_band \
    --band_nums=1,2

##CLEAN them up
wsclean -name ./images/multi-comp_grid_MWA_FEE -size 1250 1250 -niter 30000 \
    -auto-threshold 0.5 -auto-mask 3 \
    -pol I -multiscale -weight briggs 0 -scale 0.02 -j 8 -mgain 0.85 \
    -no-update-model-required \
    ./data/multi-comp_grid_MWA_FEE_band*.ms


rm casa*.log
rm importuvfits.last
