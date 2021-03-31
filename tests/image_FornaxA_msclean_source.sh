mkdir -p images

casa -c uv2ms.py ./data/FornaxA_msclean_zenith

wsclean -name ./images/FornaxA_msclean_zenith -size 350 350 -niter 80000 -auto-threshold 0.5 -auto-mask 3 \
  -pol I -multiscale -weight briggs 0 -scale 0.004 -j 12 -mgain 0.85 \
  -no-update-model-required \
  ./data/FornaxA_msclean_zenith_band*.ms

rm casa*.log
