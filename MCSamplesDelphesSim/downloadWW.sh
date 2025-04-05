#!/bin/bash


cd ../../..

###################Downloading 
# Base path for the remote directory
remote_base="ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/archive_pythia8308_16_15_25May2024/run"

# Local directory where files will be downloaded
local_base="analysis/filesPostDelphes/WWHepMC"

# List of file patterns to download
file_patterns=(
    "e1E1_WW0_v05741_iEcm{i}_scale676.hepmc"
)

for i in {0..23}; do
    for file_pattern in "${file_patterns[@]}"; do
        file=$(echo $file_pattern | sed "s/{i}/$i/g")
        remote_file="${remote_base}/${file}"
        scp $remote_file $local_base
    done
done

# List of file patterns to download
file_patterns=(
    "e1E1_WW1_v05741_iEcm{i}_scale676.hepmc"
)

for i in {0..23}; do
    for file_pattern in "${file_patterns[@]}"; do
        file=$(echo $file_pattern | sed "s/{i}/$i/g")
        remote_file="${remote_base}/${file}"
        scp $remote_file $local_base
    done
done
########################Downloading

##################Delphes 
# Specify the Delphes card
# CARD_NAME="cards/delphes_card_ILD_DSiDi.tcl"
CARD_NAME="cards/delphes_card_SiD_2024_XCC.tcl"

for i in {0..23}; do
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/WWHepMC/e1E1_WW0_v05741_iEcm${i}_scale676.root analysis/filesPostDelphes/WWHepMC/e1E1_WW0_v05741_iEcm${i}_scale676.hepmc
done
for i in {0..23}; do
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/WWHepMC/e1E1_WW1_v05741_iEcm${i}_scale676.root analysis/filesPostDelphes/WWHepMC/e1E1_WW1_v05741_iEcm${i}_scale676.hepmc
done
######################Delphes

#######root file addition
# Define the output file name
# output_file="analysis/filesPostDelphes/WWHepMC/GammaGammaWWESpreadAllILDDSiDi.root"
output_file="analysis/filesPostDelphes/WWHepMC/GammaGammaWWESpreadAllSiD2024XCC.root"

# Initialize an empty string to hold all the root files
root_files=""

for i in {0..23}; do
    root_files+="analysis/filesPostDelphes/WWHepMC/e1E1_WW0_v05741_iEcm${i}_scale676.root "
done
for i in {0..23}; do
    root_files+="analysis/filesPostDelphes/WWHepMC/e1E1_WW1_v05741_iEcm${i}_scale676.root "
done


hadd $output_file $root_files
