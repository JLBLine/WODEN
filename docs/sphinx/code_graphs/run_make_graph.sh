python make_graph.py

for image in all array_layout observational skymodel use_libwoden uvfits wodenpy_setup
do
    dot -Tsvg wodenpy_${image}.dot > wodenpy_${image}.svg
done

##then go to docs/doxygen and run:
# doxygen Doxyfile_make_graphs
##and then look in the html directory for appropriate graphs