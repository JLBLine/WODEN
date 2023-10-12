mkdir -p images

woden_uv2ms.py \
  --uvfits_prepend=./data/EDA2_haslam_band \
  --band_nums=1

time wsclean -name ./images/EDA2_haslam -size 1000 1000 -niter 80000 \
  -auto-threshold 0.5 -auto-mask 3 \
  -pol I -multiscale -weight briggs 0 -scale 0.12 -j 12 -mgain 0.85 \
  -no-update-model-required \
  ./data/EDA2_haslam_band*.ms