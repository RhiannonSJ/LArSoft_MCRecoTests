#!/bin/bash

WORKING_DIR=/sbnd/app/users/rsjones/LArSoft_v06_49_03/LArSoft-v06_49_03/srcs/recoperformance/recoperformance/working_dir

rm -rf /sbnd/app/users/rsjones/LArSoft_v06_49_03/LArSoft-v06_49_03/srcs/recoperformance/recoperformance/distances.txt
cd $WORKING_DIR

lar -c run_mcRecoComp.fcl -S /pnfs/sbnd/scratch/users/rsjones/reco_files_mcc_0_75.txt

cd -
