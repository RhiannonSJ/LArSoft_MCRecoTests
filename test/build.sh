#!/bin/bash

export FHICL_DIR=/sbnd/app/users/rsjones/LArSoft_reco_validation/LArSoft-v06_45_00/srcs/recoperformance/recoperformance
cd $MRB_BUILDDIR
mrbsetenv
mrb i -j4
mrbslp
cd -
export FHICL_FILE_PATH=${FHICL_DIR}:${FHICL_FILE_PATH}
