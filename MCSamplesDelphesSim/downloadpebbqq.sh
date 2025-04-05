#!/bin/bash


cd ../../..


###################Downloading 
# Base path for the remote directory
remote_base="ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/whizard-3.1.5/run/eA_eA_eee/"

# Local directory where files will be downloaded
local_base="analysis/filesPostDelphes/pebbqqHepMC"

# List of file patterns to download

file_patterns=(
    ep_bbqq_00_ps_had_allhad250_v06888.+-.ep/ep_bbqq_ps_had_allhad250_nofilter_+-.hepmc
    pe_bbqq_00_ps_had_allhad250_v06888.+-.pe/pe_bbqq_ps_had_allhad250_nofilter_+-.hepmc
    ep_bbqq_00_ps_had_allhad250_v06888.-+.ep/ep_bbqq_ps_had_allhad250_nofilter_-+.hepmc
    pe_bbqq_00_ps_had_allhad250_v06888.-+.pe/pe_bbqq_ps_had_allhad250_nofilter_-+.hepmc
    ep_bbqq_01_ps_had_allhad250_v06888.+-.ep/ep_bbqq_ps_had_allhad250_nofilter_+-.hepmc
    pe_bbqq_01_ps_had_allhad250_v06888.+-.pe/pe_bbqq_ps_had_allhad250_nofilter_+-.hepmc
    ep_bbqq_01_ps_had_allhad250_v06888.-+.ep/ep_bbqq_ps_had_allhad250_nofilter_-+.hepmc
    pe_bbqq_01_ps_had_allhad250_v06888.-+.pe/pe_bbqq_ps_had_allhad250_nofilter_-+.hepmc
    ep_bbqq_02_ps_had_allhad250_v06888.+-.ep/ep_bbqq_ps_had_allhad250_nofilter_+-.hepmc
    pe_bbqq_02_ps_had_allhad250_v06888.+-.pe/pe_bbqq_ps_had_allhad250_nofilter_+-.hepmc
    ep_bbqq_02_ps_had_allhad250_v06888.-+.ep/ep_bbqq_ps_had_allhad250_nofilter_-+.hepmc
    pe_bbqq_02_ps_had_allhad250_v06888.-+.pe/pe_bbqq_ps_had_allhad250_nofilter_-+.hepmc
)

i=0
for file_pattern in "${file_patterns[@]}"; do
    remote_file="${remote_base}/${file_pattern}"
    local_file="${local_base}/pe_bbqq_${i}.hepmc"
    echo "Downloading file $i to $local_file"
    scp $remote_file $local_file
    ((i++))
done
########################Downloading



##################Delphes 
# Specify the Delphes card
# CARD_NAME="cards/delphes_card_ILD_DSiDi.tcl"
CARD_NAME="cards/delphes_card_SiD_2024_XCC.tcl"

for i in {0..11}; do
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/pebbqqHepMC/pe_bbqq_${i}.root analysis/filesPostDelphes/pebbqqHepMC/pe_bbqq_${i}.hepmc
done
######################Delphes

#######root file addition
# Define the output file name
# output_file="analysis/filesPostDelphes/pebbqqHepMC/pebbqqESpreadAllILDDSiDi.root"
output_file="analysis/filesPostDelphes/pebbqqHepMC/pebbqqESpreadAllSiD2024XCC.root"

# Initialize an empty string to hold all the root files
root_files=""

for i in {0..11}; do
    root_files+="analysis/filesPostDelphes/pebbqqHepMC/pe_bbqq_${i}.root "
done

hadd $output_file $root_files
########addition

