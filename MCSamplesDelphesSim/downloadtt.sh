#!/bin/bash


cd ../../..

##################Downloading 
# Base path for the remote directory
scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/archive_pythia8308_16_15_25May2024/run/aa_fixed_energy_tt_hepmcout.hepmc analysis/filesPostDelphes/ttHepMC
########################Downloading

##################Delphes 
# Specify the Delphes card
# CARD_NAME="cards/delphes_card_ILD_DSiDi.tcl"
CARD_NAME="cards/delphes_card_SiD_2024_XCC.tcl"

# ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/ttHepMC/GammaGammattAllILDDSiDi.root analysis/filesPostDelphes/ttHepMC/aa_fixed_energy_tt_hepmcout.hepmc
./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/ttHepMC/GammaGammattAllSiD2024XCC.root analysis/filesPostDelphes/ttHepMC/aa_fixed_energy_tt_hepmcout.hepmc
######################Delphes

