for flux in curve list power
  do
  for comp in singlepoint singlegauss singleshape threecomponents threesources
    do
        python convert_skymodel_yaml_to_FITS.py \
            --input_yaml ../yaml_models/srclist_${comp}_${flux}.yaml \
            --output_fits srclist_${comp}_${flux}.fits
  done
done

python convert_skymodel_yaml_to_FITS.py \
            --input_yaml ../yaml_models/srclist_mulitple_source-components.yaml \
            --output_fits srclist_mulitple_source-components.fits

