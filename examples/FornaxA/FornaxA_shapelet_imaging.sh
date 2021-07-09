mkdir -p images
#rm images/*

##Convert uvfits to measurement sets. --nologger stops a logging popup
##that can slow things down. Remove it if you want a graphical logger
${CASA_DIR}/casa --nologger -c  ${WODEN_DIR}/uv2ms.py \
  --uvfits_prepend=./data/FornaxA_shapelets_band \
  --band_nums=1,2,3,4,5

wsclean -name ./images/FornaxA_shapelets -size 350 350 -niter 80000 \
  -auto-threshold 0.5 -auto-mask 3 \
  -pol I -multiscale -weight briggs 0 -scale 0.004 -j 12 -mgain 0.85 \
  -no-update-model-required \
  ./data/FornaxA_shapelets_band*.ms

rm casa*.log
