##This runs hyperdrive and stores the results in text files
python compare_to_hyperdrive_multi.py

##I run this command after to turn some of those text files into the
##header file test_multifreq_RTS_FEE_beam.h. These values are then
##used in testing in test_run_and_map_multifreq_calc_CUDA_FEE_beam.c
python make_values_header_multi.py


##Next, you have to run the ctest command in build otherwise there are no
##WODEN outputs to link. If you have your build dir where I suggest, this
##command should symlink things to here:
ln -s ../../build/cmake_testing/FEE_primary_beam_cuda/MWA_FEE_multifreq_gains_*.txt .

##Now make the comparison plots
python plot_multi_freq_results.py
