mkdir -p images
#rm images/*

##Convert uvfits to measurement sets. --nologger stops a logging popup
##that can slow things down. Remove it if you want a graphical logger
${CASA_DIR}/casa --nologger -c  ${WODEN_DIR}/uv2ms.py \
  --uvfits_prepend=./data/EDA2_haslam_band \
  --band_nums=1

time wsclean -name ./images/EDA2_haslam -size 1000 1000 -niter 80000 \
  -auto-threshold 0.5 -auto-mask 3 \
  -pol I -multiscale -weight briggs 0 -scale 0.12 -j 12 -mgain 0.85 \
  -no-update-model-required \
  ./data/EDA2_haslam_band*.ms

rm casa*.log
