#!/bin/sh

come_home=$PWD

cd ${come_home}/single_component_models
./single_component_sims.sh

cd ${come_home}/grid_component_models
./grid_component_sims.sh

cd ${come_home}/different_beam_models
./multi-comp_sims_all_but_MWAFEE.sh

cd ${come_home}/absolute_accuracy
./run_the_absolute_accuracy_test.sh

cd ${come_home}
