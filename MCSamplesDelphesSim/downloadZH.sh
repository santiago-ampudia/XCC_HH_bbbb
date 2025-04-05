#!/bin/bash


cd ../../..

###################Downloading 
# Base path for the remote directory
remote_base="ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/"

# Local directory where files will be downloaded
local_base="analysis/filesPostDelphes/ZHHepMC"

# List of file patterns to download
file_patterns=(
    "AA_HZ_{i}_ps_had_SM_ac_CKM_v06395.++.gg/AA_HZ_380.hepmc"
)

for i in {0..9}; do
    for file_pattern in "${file_patterns[@]}"; do
        file=$(echo $file_pattern | sed "s/{i}/$i/g")
        remote_file="${remote_base}/${file}"
        local_file="${local_base}/AA_HZ_380_${i}.hepmc"
        scp $remote_file $local_file
    done
done
########################Downloading

##################Delphes 
# Specify the Delphes card
# CARD_NAME="cards/delphes_card_ILD_DSiDi.tcl"
CARD_NAME="cards/delphes_card_SiD_2024_XCC.tcl"

for i in {0..9}; do
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/ZHHepMC/AA_HZ_380_${i}.root analysis/filesPostDelphes/ZHHepMC/AA_HZ_380_${i}.hepmc
done
######################Delphes

#######root file addition
# Define the output file name
# output_file="analysis/filesPostDelphes/ZHHepMC/GammaGammaZHESpreadAllILDDSiDi.root"
output_file="analysis/filesPostDelphes/ZHHepMC/GammaGammaZHESpreadAllSiD2024XCC.root"

# Initialize an empty string to hold all the root files
root_files=""

for i in {0..9}; do
    root_files+="analysis/filesPostDelphes/ZHHepMC/AA_HZ_380_${i}.root "
done


hadd $output_file $root_files
########addition

####### Remove all files except for GammaGammaZHESpreadAllILDDSiDi80BTag.root and AA_HZ_380_0.hepmc
#find analysis/filesPostDelphes/ZHHepMC -type f ! -name "GammaGammaZHESpreadAllILDDSiDi80BTag.root" ! -name "AA_HZ_380_0.hepmc" -exec rm {} +
#######removal
