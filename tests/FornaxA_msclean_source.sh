mkdir -p data

run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --num_freq_channels=64 --num_time_steps=5 \
    --freq_res=20e+3 --time_res=0.5 \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.txt \
    --metafits_filename=1202815152_metafits_ppds.fits \
    --band_nums=1,2,3 --output_uvfits_prepend=./data/FornaxA_msclean_zenith \
    --primary_beam=MWA_FEE \
    --chunking_size=2000
