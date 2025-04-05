#!/bin/bash


cd ../../..

###################Downloading 
#scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_HH_00_ps_had_SM_ac_CKM_v06395.++.gg/AA_HH_380_nofilter.hepmc analysis/filesPostDelphes/HHHepMC/newHH/AA_HH_380_nofilter_0.hepmc
#scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_HH_00_ps_had_SM_ac_CKM_v06395.+-.gg/AA_HH_380_nofilter.hepmc analysis/filesPostDelphes/HHHepMC/newHH/AA_HH_380_nofilter_1.hepmc
#scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_HH_00_ps_had_SM_ac_CKM_v06395.-+.gg/AA_HH_380_nofilter.hepmc analysis/filesPostDelphes/HHHepMC/newHH/AA_HH_380_nofilter_2.hepmc
#scp ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/AA_HH_00_ps_had_SM_ac_CKM_v06395.--.gg/AA_HH_380_nofilter.hepmc analysis/filesPostDelphes/HHHepMC/newHH/AA_HH_380_nofilter_3.hepmc

# Base path for the remote directory
#remote_base="ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee"

# Local directory where files will be downloaded
#local_base="analysis/filesPostDelphes/HHHepMC/newHH"

# List of file patterns to download
#file_patterns=(
#    "AA_HH_0{i}_ps_had_SM_ac_CKM_v06395.++.gg/AA_HH_380_nofilter.hepmc"
#)

# Loop over the desired range
#for i in {1..9}; do
#    file=$(echo $file_pattern | sed "s/{i}/$i/g")
#    new_index=$((i + 3))
#    remote_file="${remote_base}/${file}"
#    local_file="${local_base}/AA_HH_380_nofilter_${new_index}.hepmc"
#    echo "Downloading $remote_file to $local_file..."
#    scp "$remote_file" "$local_file"
#done
########################Downloading

##################Delphes 
# Specify the Delphes card
# CARD_NAME="cards/delphes_card_ILD_DSiDi.tcl"
CARD_NAME="cards/delphes_card_SiD_2024_XCC.tcl"

for i in {0..12}; do
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/HHHepMC/newHH/AA_HH_380_nofilter_${i}.root analysis/filesPostDelphes/HHHepMC/newHH/AA_HH_380_nofilter_${i}.hepmc
done
######################Delphes

#######root file addition
# Define the output file name
# output_file="analysis/filesPostDelphes/HHHepMC/newHH/GammaGammaHHESpreadAllILDDSiDi.root"
output_file="analysis/filesPostDelphes/HHHepMC/newHH/GammaGammaHHESpreadAllSiD2024XCC.root"

# Initialize an empty string to hold all the root files
root_files=""

for i in {0..12}; do
    root_files+="analysis/filesPostDelphes/HHHepMC/newHH/AA_HH_380_nofilter_${i}.root "
done


hadd $output_file $root_files
