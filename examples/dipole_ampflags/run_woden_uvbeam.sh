SRCLIST=srclist_singlepoint_power.yaml
META=../metafits/1126115208_metafits.fits

# time run_woden.py \
#     --ra0=0.0 --dec0=-27.0 \
#     --num_freq_channels=32 --num_time_steps=14 \
#     --freq_res=40e+3 --time_res=8.0 \
#     --cat_filename=srclist_singlepoint_power.yaml \
#     --metafits_filename=$META \
#     --band_nums=1 \
#     --output_uvfits_prepend=woden_uvbeam_dipflags \
#     --primary_beam=uvbeam_MWA \
#     --lowest_channel_freq=181.775e+6 \
#     --use_MWA_dipflags

# time run_woden.py \
#     --ra0=0.0 --dec0=-27.0 \
#     --num_freq_channels=32 --num_time_steps=14 \
#     --freq_res=40e+3 --time_res=8.0 \
#     --cat_filename=srclist_singlepoint_power.yaml \
#     --metafits_filename=$META \
#     --band_nums=1 \
#     --output_uvfits_prepend=woden_uvbeam_noflags \
#     --primary_beam=uvbeam_MWA \
#     --lowest_channel_freq=181.775e+6

META=../metafits/1088285600_DipAmps.metafits

time run_woden.py \
    --ra0=0.0 --dec0=-27.0 \
    --num_freq_channels=32 --num_time_steps=14 \
    --freq_res=40e+3 --time_res=8.0 \
    --cat_filename=srclist_singlepoint_power.yaml \
    --metafits_filename=$META \
    --band_nums=1 \
    --output_uvfits_prepend=woden_uvbeam_dipamps \
    --primary_beam=uvbeam_MWA \
    --lowest_channel_freq=181.775e+6 \
    --use_MWA_dipamps


time run_woden.py \
    --ra0=0.0 --dec0=-27.0 \
    --num_freq_channels=32 --num_time_steps=14 \
    --freq_res=40e+3 --time_res=8.0 \
    --cat_filename=srclist_singlepoint_power.yaml \
    --metafits_filename=$META \
    --band_nums=1 \
    --output_uvfits_prepend=woden_uvbeam_default \
    --primary_beam=uvbeam_MWA \
    --lowest_channel_freq=181.775e+6