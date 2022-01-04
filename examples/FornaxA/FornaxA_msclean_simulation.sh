mkdir -p data

for precision in "float" "double"
do
  time run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --num_freq_channels=16 --num_time_steps=14 \
    --freq_res=80e+3 --time_res=8.0 \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.txt \
    --metafits_filename=../metafits/1202815152_metafits_ppds.fits \
    --band_nums=1,2,3,4,5 \
    --output_uvfits_prepend=./data/FornaxA_msclean_${precision} \
    --primary_beam=MWA_FEE \
    --precision=${precision}
done
