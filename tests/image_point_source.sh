mkdir -p images

casa -c uv2ms.py ./data/point_zenith

wsclean -name ./images/point_zenith -size 512 512 -niter 80000 -auto-threshold 0.5 -auto-mask 3 \
  -pol I -multiscale -weight briggs 0 -scale 0.004 -j 12 -mgain 0.85 \
  -no-update-model-required \
  ./data/point_zenith_band*.ms

rm casa*.log