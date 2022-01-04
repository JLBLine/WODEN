#!/bin/sh

mkdir -p images
##See if we can find where the scripts live - annoyingly, have to give an
##absolute path to casa to run uv2ms.py, so need to know where it lives

##This checks if the ENV VARIABLE WODEN_DIR does not exists
if [ -z ${WODEN_DIR+x} ]; then

echo "WODEN_DIR is unset. Searching for uv2ms.py instead"
##This searches for uv2ms.py apparently
uv2ms=$(command -v uv2ms.py)
echo "Will attempt to use ${uv2ms}"

##If it exists, use it
else
uv2ms="${WODEN_DIR}/uv2ms.py"
echo "WODEN_DIR is is set to '${WODEN_DIR}'.";
echo "Will use ${uv2ms}"

fi

##Convert uvfits to measurement sets. --nologger stops a logging popup
##that can slow things down. Remove it if you want a graphical logger
${CASA_DIR}/casa --nologger -c ${uv2ms} \
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
