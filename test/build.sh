#!/bin/bash

export FHICL_DIR=/sbnd/app/users/rsjones/LArSoft_v06_49_03/LArSoft-v06_49_03/srcs/recoperformance/recoperformance/fcl
cd $MRB_BUILDDIR
mrbsetenv
mrb i -j4
mrbslp
cd -
export FHICL_FILE_PATH=${FHICL_DIR}:${FHICL_FILE_PATH}
