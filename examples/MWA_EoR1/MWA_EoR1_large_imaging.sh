mkdir -p images

for precision in "float" "double"
do

  ##Convert uvfits to measurement sets. --nologger stops a logging popup
  ##that can slow things down. Remove it if you want a graphical logger
  woden_uv2ms.py \
   --uvfits_prepend=./data/MWA_EoR1_${precision}_band \
   --band_nums=1,2,3,4,5

  time wsclean -name ./images/MWA_EoR1_larger_${precision} -size 6000 6000 -niter 80000 \
   -auto-threshold 0.5 -auto-mask 3 \
   -pol I -multiscale -weight briggs 0 -scale 0.015 -j 12 -mgain 0.85 \
   -no-update-model-required \
   -abs-mem 16 \
   ./data/MWA_EoR1_${precision}_band*.ms

done