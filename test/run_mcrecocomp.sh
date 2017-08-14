#!/bin/bash

INPUT_FILE=/sbnd/data/users/rsjones/cczeropi_outputs/cczeropi_events_rerun_reco_1.root
N_EVENTS=10

lar -c run_mcRecoComp.fcl -s ${INPUT_FILE} -n ${N_EVENTS}
