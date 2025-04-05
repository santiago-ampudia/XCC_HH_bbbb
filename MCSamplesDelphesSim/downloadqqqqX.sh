cd ../../..

###################Downloading 
# Base path for the remote directory
remote_base="ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/archive_pythia8308_16_15_25May2024/run"

# Local directory where files will be downloaded
local_base="analysis/filesPostDelphes/qqqqXHepMC"

# List of file patterns to download
file_patterns=(
    "Ae_qqqqe_ps_had_{i}_allhad250_v06395.++.ge/Ae_qqqqe_ps_had_{i}_allhad250_++.hepmc"
    "Ae_qqqqe_ps_had_{i}_allhad250_v06395.+-.ge/Ae_qqqqe_ps_had_{i}_allhad250_+-.hepmc"
    "Ae_qqqqe_ps_had_{i}_allhad250_v06395.-+.ge/Ae_qqqqe_ps_had_{i}_allhad250_-+.hepmc"
    "Ae_qqqqe_ps_had_{i}_allhad250_v06395.--.ge/Ae_qqqqe_ps_had_{i}_allhad250_--.hepmc"
    "Ae_qqqqv_ps_had_{i}_allhad250_v06395.+-.ge/Ae_qqqqv_ps_had_{i}_allhad250_+-.hepmc"
    "Ae_qqqqv_ps_had_{i}_allhad250_v06395.--.ge/Ae_qqqqv_ps_had_{i}_allhad250_--.hepmc"
    "eA_eqqqq_ps_had_{i}_allhad250_v06395.++.eg/eA_eqqqq_ps_had_{i}_allhad250_++.hepmc"
    "eA_eqqqq_ps_had_{i}_allhad250_v06395.+-.eg/eA_eqqqq_ps_had_{i}_allhad250_+-.hepmc"
    "eA_eqqqq_ps_had_{i}_allhad250_v06395.-+.eg/eA_eqqqq_ps_had_{i}_allhad250_-+.hepmc"
    "eA_eqqqq_ps_had_{i}_allhad250_v06395.--.eg/eA_eqqqq_ps_had_{i}_allhad250_--.hepmc"
    "eA_vqqqq_ps_had_{i}_allhad250_v06395.-+.eg/eA_vqqqq_ps_had_{i}_allhad250_-+.hepmc"
    "eA_vqqqq_ps_had_{i}_allhad250_v06395.--.eg/eA_vqqqq_ps_had_{i}_allhad250_--.hepmc"
    "ee_eqqqq_ps_had_{i}_allhad250_v06395.++.ee/ee_eqqqq_ps_had_{i}_allhad250_++.hepmc"
    "ee_eqqqq_ps_had_{i}_allhad250_v06395.+-.ee/ee_eqqqq_ps_had_{i}_allhad250_+-.hepmc"
    "ee_eqqqq_ps_had_{i}_allhad250_v06395.-+.ee/ee_eqqqq_ps_had_{i}_allhad250_-+.hepmc"
    "ee_qqqqe_ps_had_{i}_allhad250_v06395.++.ee/ee_qqqqe_ps_had_{i}_allhad250_++.hepmc"
    "ee_qqqqe_ps_had_{i}_allhad250_v06395.+-.ee/ee_qqqqe_ps_had_{i}_allhad250_+-.hepmc"
    "ee_qqqqe_ps_had_{i}_allhad250_v06395.-+.ee/ee_qqqqe_ps_had_{i}_allhad250_-+.hepmc"
    "ee_qqqqv_ps_had_{i}_allhad250_v06395.+-.ee/ee_qqqqv_ps_had_{i}_allhad250_+-.hepmc"
    "ee_vqqqq_ps_had_{i}_allhad250_v06395.-+.ee/ee_vqqqq_ps_had_{i}_allhad250_-+.hepmc"
)

#for i in {1..25}; do
#    for file_pattern in "${file_patterns[@]}"; do
#        # Replace {i} with the current number in the loop
#        file=$(echo $file_pattern | sed "s/{i}/$i/g")
#        remote_file="${remote_base}/${file}"
#        scp $remote_file $local_base
#    done
#done
########################Downloading

#######################Delphes
# Specify the Delphes card
# CARD_NAME="cards/delphes_card_ILD_DSiDi.tcl"
CARD_NAME="cards/delphes_card_SiD_2024_XCC.tcl"

for i in {1..25}; do
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_++.root analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_++.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqv_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqv_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqv_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqv_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_++.root analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_++.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/eA_vqqqq_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqqqXHepMC/eA_vqqqq_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/eA_vqqqq_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqqqXHepMC/eA_vqqqq_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/ee_eqqqq_ps_had_${i}_allhad250_++.root analysis/filesPostDelphes/qqqqXHepMC/ee_eqqqq_ps_had_${i}_allhad250_++.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/ee_eqqqq_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqqqXHepMC/ee_eqqqq_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/ee_eqqqq_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqqqXHepMC/ee_eqqqq_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqe_ps_had_${i}_allhad250_++.root analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqe_ps_had_${i}_allhad250_++.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqe_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqe_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqe_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqe_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqv_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqv_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqqqXHepMC/ee_vqqqq_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqqqXHepMC/ee_vqqqq_ps_had_${i}_allhad250_-+.hepmc
done
##########################Delphes

#######root file addition
# Define the output file name
# output_file="analysis/filesPostDelphes/qqqqXHepMC/eGammaqqqqXAllILDDSiDi.root"
output_file="analysis/filesPostDelphes/qqqqXHepMC/eGammaqqqqXAllSiD2024XCC.root"

# Initialize an empty string to hold all the root files
root_files=""
for i in {1..25}; do
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqe_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqv_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/Ae_qqqqv_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/eA_eqqqq_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/eA_vqqqq_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/eA_vqqqq_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/ee_eqqqq_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/ee_eqqqq_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/ee_eqqqq_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqe_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqe_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqe_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/ee_qqqqv_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqqqXHepMC/ee_vqqqq_ps_had_${i}_allhad250_-+.root "
done

hadd $output_file $root_files

###########root file addition