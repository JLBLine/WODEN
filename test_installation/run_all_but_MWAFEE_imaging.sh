#!/bin/sh

come_home=$PWD

cd ${come_home}/single_component_models
./single_component_images.sh

cd ${come_home}/grid_component_models
./grid_component_images.sh

cd ${come_home}/different_beam_models
./multi-comp_imaging_all_but_MWAFEE.sh

cd ${come_home}
