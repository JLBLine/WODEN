
##Need to run ctest first, then you should link in the results like
ln -s ../../build/cmake_testing/FEE_primary_beam_cuda/*float*.txt .
ln -s ../../build/cmake_testing/FEE_primary_beam_cuda/*double*.txt .

##Then you can make hyperdrive outputs like this
python compare_to_hyperdrive.py

##Then can run the plotting scripts to compare the jones matrices
##from hyperdrive to
mkdir -p beam_plots

for tag in "zenith" #"offzen1" "offzen2"
do
  for freq in "100" #"150" "200"
  do

    echo "Doing ${tag} for FLOAT precision"
    python plot_beam_results_diff.py hyperbeam_${tag}_${freq}.txt ${tag}_${freq}_float.txt
    python plot_beam_results_diff.py hyperbeam_${tag}_${freq}_rot.txt ${tag}_${freq}_rot_float.txt

    echo "Doing ${tag} for DOUBLE precision"
    python plot_beam_results_diff.py hyperbeam_${tag}_${freq}.txt ${tag}_${freq}_double.txt
    python plot_beam_results_diff.py hyperbeam_${tag}_${freq}_rot.txt ${tag}_${freq}_rot_double.txt
  done
done
