#!/bin/bash
# Executable script for CERN LXPLUS HTCondor
# Job: sps_3GeV_0.0deg_particleID_0_setup_1_116

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
      -N=10 \
      --steeringFile steering.py \
      --outputFile=testProtoECAL.root \
      --gun.position "0.0 0.0 -110.0*cm" \
      --gun.direction "0.0 0.0 1.0" \
      --gun.energy "10*GeV" \
      --part.userParticleHandler="" \
      --gun.particle "mu-"

echo "Job completed successfully"