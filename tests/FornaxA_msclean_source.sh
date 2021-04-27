mkdir -p data

time run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --num_freq_channels=16 --num_time_steps=14 \
    --freq_res=80e+3 --time_res=8.0 \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.txt \
    --metafits_filename=1202815152_metafits_ppds.fits \
    --band_nums=1 --output_uvfits_prepend=./data/FornaxA_msclean_zenith \
    --primary_beam=None \
    --hdf5_beam_path=/home/jline/software/useful/mwa_full_embedded_element_pattern.h5 \
    --chunking_size=2000
    
    
#       --hdf5_beam_path=/home/jline/software/useful/mwa_full_embedded_element_pattern.h5 \
