#!/bin/bash

# test exclude

WORKING_DIR=/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/working_dir
RECO_COMP_DIR=/sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance

TEST_1=/pnfs/sbnd/persistent/sbndpro/mcc0.75/v06_48_00_MCC/prodgenie_nu_singleinteraction_cryostat_gsimple-configb-v1/reco/19889636_0/prodgenie_sbnd_GenieGen-20170905T141930_32611fa6-040c-4904-bae8-93e5abf1ae2f_G4-20170905T143213_DetSim-20170905T151617_Reco-20170905T163304.root
TEST_20=/pnfs/sbnd/scratch/users/rsjones/cc0pi_files_mcc_0_75_20.txt
TEST_ALL=/pnfs/sbnd/scratch/users/rsjones/reco_files_mcc_0_75.txt

rm -rf /sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/distances.txt

rm -rf /sbnd/app/users/rsjones/LArSoft_v06_52_00/LArSoft-v06_52_00/srcs/recoperformance/recoperformance/plots/nt_distances.root

cd $WORKING_DIR

lar -c run_mcRecoComp.fcl -S $TEST_ALL

cd $RECO_COMP_DIR
