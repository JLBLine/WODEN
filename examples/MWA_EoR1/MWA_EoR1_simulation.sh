mkdir -p data

for precision in "float" "double"
do

  time run_woden.py \
    --ra0=60.0 --dec0=-27.0 \
    --num_freq_channels=16 --num_time_steps=14 \
    --freq_res=80e+3 --time_res=8.0 \
    --cat_filename=srclist_pumav3_EoR0LoBES_EoR1pietro_CenA-GP_2023-11-07.fits \
    --metafits_filename=../metafits/1136380296_metafits_ppds.fits \
    --band_nums=1,2,3,4,5 \
    --output_uvfits_prepend=./data/MWA_EoR1_${precision} \
    --primary_beam=MWA_FEE \
    --sky_crop_components \
    --precision=${precision}

done
