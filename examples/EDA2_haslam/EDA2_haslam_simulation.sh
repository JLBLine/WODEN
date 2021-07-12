mkdir -p data

time run_woden.py \
    --ra0=74.79589467 --dec0=-27.0 \
    --num_freq_channels=10 --num_time_steps=10 \
    --freq_res=10e+3 --time_res=10.0 \
    --lowest_channel_freq=99e+6 \
    --cat_filename=pygsm_woden-list_100MHz_n256.txt \
    --array_layout=../../test_installation/array_layouts/EDA2_layout_255.txt \
    --date=2020-02-01T12:27:45.900 \
    --output_uvfits_prepend=./data/EDA2_haslam \
    --primary_beam=EDA2 \
    --sky_crop_components \
    --band_nums=1,2,3,4,5
