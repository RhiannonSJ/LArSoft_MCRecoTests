#!/bin/bash

INPUT_FILE=/sbnd/data/users/rsjones/hit_filtered_cc0pi_sample/19366014_12/prodgenie_sbnd_GenieGen-20170808T233908_5dd110c9-c622-4eb5-b5ed-8b800d0cb0d5_topologyFilteredEvents_TopologyFilter-20170809T000505_G4-20170809T012154__b7bc1705-c469-4b40-940c-eb5e548df128.root

N_EVENTS=10

WORKING_DIR=/sbnd/app/users/rsjones/LArSoft_reco_validation/LArSoft-v06_45_00/srcs/recoperformance/recoperformance/working_dir

cd $WORKING_DIR

lar -c run_mcRecoComp.fcl -s ${INPUT_FILE} -n ${N_EVENTS}

cd -
