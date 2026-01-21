#!/bin/bash
# Executable script for CERN LXPLUS HTCondor
# Job: sim_0.0deg_particleID_0_288GeV_setup_1_Run999

# Setup environment using LCG views
cd /afs/cern.ch/user/s/sritter/prototype_simulation/DD4SHiP
currDir=`pwd`
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
cd $currDir
unset currDir

source install/bin/thisdd4hep.sh

cd DDDetectors/shipdetectors

ddsim --compactFile=./Caloprototype.xml \
      --runType=batch \
      -G \
      -N=10000 \
      --steeringFile steering.py \
      --outputFile=/eos/user/s/sritter/ship_sim_data/sim_0.0deg_particleID_0_288GeV_setup_1_Run999_v1.root \
      --gun.position "0.0 0.0 -110.0*cm" \
      --gun.direction "0.0 0.0 1.0" \
      --gun.energy "288*GeV" \
      --part.userParticleHandler="" \
      --gun.particle "e-"

echo "Job completed successfully"