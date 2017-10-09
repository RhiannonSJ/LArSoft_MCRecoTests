#!/bin/bash

WORKING_DIR=/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/working_dir
RECO_COMP_DIR=/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance

TEST_ALL=/pnfs/sbnd/scratch/users/rsjones/anatree_files_mcc_0_75.txt

cd $WORKING_DIR

lar -c run_vtxParameters.fcl -S $TEST_ALL

cd $RECO_COMP_DIR
