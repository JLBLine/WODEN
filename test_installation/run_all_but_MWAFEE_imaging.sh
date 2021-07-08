#!/bin/sh

come_home=$PWD

cd ${come_home}/point_models
./point_source_imaging_all_but_MWAFEE.sh

cd ${come_home}/gauss_models
./gauss_source_imaging_all_but_MWAFEE.sh

cd ${come_home}/shapelet_models
./shapelet_source_imaging_all_but_MWAFEE.sh

cd ${come_home}/multicomp_model
./multi-comp_imaging_all_but_MWAFEE.sh

cd ${come_home}
