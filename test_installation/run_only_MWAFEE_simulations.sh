#!/bin/sh

come_home=$PWD

cd ${come_home}/point_models
./point_source_sims_only_MWAFEE.sh

cd ${come_home}/gauss_models
./gauss_source_sims_only_MWAFEE.sh

cd ${come_home}/shapelet_models
./shapelet_source_sims_only_MWAFEE.sh

cd ${come_home}/multicomp_model
./multi-comp_sims_only_MWAFEE.sh

cd ${come_home}
