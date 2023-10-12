mkdir -p images

fi

for precision in "float" "double"

do

  woden_uv2ms.py \
   --uvfits_prepend=./data/MWA_EoR1_${precision}_band \
   --band_nums=1,2,3,4,5

  time wsclean -name ./images/MWA_EoR1_smaller_${precision} -size 3000 3000 -niter 80000 \
    -auto-threshold 0.5 -auto-mask 3 \
    -pol I -multiscale -weight briggs 0 -scale 0.01 -j 12 -mgain 0.85 \
    -no-update-model-required \
    ./data/MWA_EoR1_${precision}_band*.ms

done

#rm casa*.log
