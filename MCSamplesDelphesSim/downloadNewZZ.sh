#!/bin/bash


cd ../../..

###################Downloading 
scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_ZZ_00_ps_had_SM_ac_CKM_v06395.++.gg/AA_ZZ_380.hepmc analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_0.hepmc
scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_ZZ_00_ps_had_SM_ac_CKM_v06395.+-.gg/AA_ZZ_380.hepmc analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_1.hepmc
scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_ZZ_00_ps_had_SM_ac_CKM_v06395.-+.gg/AA_ZZ_380.hepmc analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_2.hepmc
scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_ZZ_00_ps_had_SM_ac_CKM_v06395.--.gg/AA_ZZ_380.hepmc analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_3.hepmc

scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_ZZ_01_ps_had_SM_ac_CKM_v06395.++.gg/AA_ZZ_380.hepmc analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_4.hepmc
scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_ZZ_01_ps_had_SM_ac_CKM_v06395.+-.gg/AA_ZZ_380.hepmc analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_5.hepmc
scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_ZZ_01_ps_had_SM_ac_CKM_v06395.-+.gg/AA_ZZ_380.hepmc analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_6.hepmc

# Base path for the remote directory
remote_base="ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee"

# Local directory where files will be downloaded
local_base="analysis/filesPostDelphes/ZZHepMC/newZZ"

# File pattern template
file_pattern="AA_ZZ_0{i}_ps_had_SM_ac_CKM_v06395.++.gg/AA_ZZ_380.hepmc"

# Loop over the desired range
for i in {2..9}; do
    # Replace {i} in the file pattern with the loop variable
    file=$(echo "$file_pattern" | sed "s/{i}/$i/g")
    new_index=$((i + 5))

    # Construct the remote and local file paths
    remote_file="${remote_base}/${file}"
    local_file="${local_base}/AA_ZZ_380_${new_index}.hepmc"

    # Download the file
    echo "Downloading $remote_file to $local_file..."
    scp "$remote_file" "$local_file"
done

# File pattern template
file_pattern="AA_ZZ_{i}_ps_had_SM_ac_CKM_v06395.++.gg/AA_ZZ_380.hepmc"

# Loop over the desired range
for i in {1..19}; do
    # Replace {i} in the file pattern with the loop variable
    file=$(echo "$file_pattern" | sed "s/{i}/$i/g")
    new_index=$((i + 5))

    # Construct the remote and local file paths
    remote_file="${remote_base}/${file}"
    local_file="${local_base}/AA_ZZ_380_${new_index}.hepmc"

    # Download the file
    echo "Downloading $remote_file to $local_file..."
    scp "$remote_file" "$local_file"
done

########################Downloading

##################Delphes 
# Specify the Delphes card
# CARD_NAME="cards/delphes_card_ILD_DSiDi.tcl"
CARD_NAME="cards/delphes_card_SiD_2024_XCC.tcl"

for i in {0..24}; do
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_${i}.root analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_${i}.hepmc
done
######################Delphes

#######root file addition
# Define the output file name
# output_file="analysis/filesPostDelphes/ZZHepMC/newZZ/GammaGammaZZESpreadAllILDDSiDi.root"
output_file="analysis/filesPostDelphes/ZZHepMC/newZZ/GammaGammaZZESpreadAllSiD2024XCC.root"

# Initialize an empty string to hold all the root files
root_files=""

for i in {0..24}; do
    root_files+="analysis/filesPostDelphes/ZZHepMC/newZZ/AA_ZZ_380_${i}.root "
done

hadd $output_file $root_files
