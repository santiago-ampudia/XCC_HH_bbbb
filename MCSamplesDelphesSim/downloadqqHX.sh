cd ../../..

# ###################Downloading 
# # Base path for the remote directory
# remote_base="ampudia@s3dfdtn.slac.stanford.edu:/fs/ddn/sdf/group/atlas/d/timb/ubuntu/archive_pythia8308_16_15_25May2024/run"
# 
# # Local directory where files will be downloaded
# local_base="analysis/filesPostDelphes/qqHXHepMC"
# 
# # List of file patterns to download
# file_patterns=(
#     "Ae_qqHe_ps_had_{i}_allhad250_v06395.++.ge/Ae_qqHe_ps_had_{i}_allhad250_++.hepmc"
#     "Ae_qqHe_ps_had_{i}_allhad250_v06395.+-.ge/Ae_qqHe_ps_had_{i}_allhad250_+-.hepmc"
#     "Ae_qqHe_ps_had_{i}_allhad250_v06395.-+.ge/Ae_qqHe_ps_had_{i}_allhad250_-+.hepmc"
#     "Ae_qqHe_ps_had_{i}_allhad250_v06395.--.ge/Ae_qqHe_ps_had_{i}_allhad250_--.hepmc"
#     "Ae_qqHv_ps_had_{i}_allhad250_v06395.+-.ge/Ae_qqHv_ps_had_{i}_allhad250_+-.hepmc"
#     "Ae_qqHv_ps_had_{i}_allhad250_v06395.--.ge/Ae_qqHv_ps_had_{i}_allhad250_--.hepmc"
#     "eA_eqqH_ps_had_{i}_allhad250_v06395.++.eg/eA_eqqH_ps_had_{i}_allhad250_++.hepmc"
#     "eA_eqqH_ps_had_{i}_allhad250_v06395.+-.eg/eA_eqqH_ps_had_{i}_allhad250_+-.hepmc"
#     "eA_eqqH_ps_had_{i}_allhad250_v06395.-+.eg/eA_eqqH_ps_had_{i}_allhad250_-+.hepmc"
#     "eA_eqqH_ps_had_{i}_allhad250_v06395.--.eg/eA_eqqH_ps_had_{i}_allhad250_--.hepmc"
#     "eA_vqqH_ps_had_{i}_allhad250_v06395.-+.eg/eA_vqqH_ps_had_{i}_allhad250_-+.hepmc"
#     "eA_vqqH_ps_had_{i}_allhad250_v06395.--.eg/eA_vqqH_ps_had_{i}_allhad250_--.hepmc"
#     "ee_eqqH_ps_had_{i}_allhad250_v06395.++.ee/ee_eqqH_ps_had_{i}_allhad250_++.hepmc"
#     "ee_eqqH_ps_had_{i}_allhad250_v06395.+-.ee/ee_eqqH_ps_had_{i}_allhad250_+-.hepmc"
#     "ee_eqqH_ps_had_{i}_allhad250_v06395.-+.ee/ee_eqqH_ps_had_{i}_allhad250_-+.hepmc"
#     "ee_eqqH_ps_had_{i}_allhad250_v06395.--.ee/ee_eqqH_ps_had_{i}_allhad250_--.hepmc"
#     "ee_qqHe_ps_had_{i}_allhad250_v06395.++.ee/ee_qqHe_ps_had_{i}_allhad250_++.hepmc"
#     "ee_qqHe_ps_had_{i}_allhad250_v06395.+-.ee/ee_qqHe_ps_had_{i}_allhad250_+-.hepmc"
#     "ee_qqHe_ps_had_{i}_allhad250_v06395.-+.ee/ee_qqHe_ps_had_{i}_allhad250_-+.hepmc"
#     "ee_qqHe_ps_had_{i}_allhad250_v06395.--.ee/ee_qqHe_ps_had_{i}_allhad250_--.hepmc"
#     "ee_qqHv_ps_had_{i}_allhad250_v06395.+-.ee/ee_qqHv_ps_had_{i}_allhad250_+-.hepmc"
#     "ee_qqHv_ps_had_{i}_allhad250_v06395.--.ee/ee_qqHv_ps_had_{i}_allhad250_--.hepmc"
#     "ee_vqqH_ps_had_{i}_allhad250_v06395.-+.ee/ee_vqqH_ps_had_{i}_allhad250_-+.hepmc"
#     "ee_vqqH_ps_had_{i}_allhad250_v06395.--.ee/ee_vqqH_ps_had_{i}_allhad250_--.hepmc"
# )
# 
# for i in {1..5}; do
#     for file_pattern in "${file_patterns[@]}"; do
#         # Replace {i} with the current number in the loop
#         file=$(echo $file_pattern | sed "s/{i}/$i/g")
#         remote_file="${remote_base}/${file}"
#         scp $remote_file $local_base
#     done
# done
# ########################Downloading

#######################Delphes
# Specify the Delphes card
# CARD_NAME="cards/delphes_card_ILD_DSiDi.tcl"
CARD_NAME="cards/delphes_card_SiD_2024_XCC.tcl"

for i in {1..5}; do
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_++.root analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_++.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/Ae_qqHv_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqHXHepMC/Ae_qqHv_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/Ae_qqHv_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqHXHepMC/Ae_qqHv_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_++.root analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_++.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/eA_vqqH_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqHXHepMC/eA_vqqH_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/eA_vqqH_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqHXHepMC/eA_vqqH_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_++.root analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_++.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_++.root analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_++.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_qqHv_ps_had_${i}_allhad250_+-.root analysis/filesPostDelphes/qqHXHepMC/ee_qqHv_ps_had_${i}_allhad250_+-.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_qqHv_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqHXHepMC/ee_qqHv_ps_had_${i}_allhad250_--.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_vqqH_ps_had_${i}_allhad250_-+.root analysis/filesPostDelphes/qqHXHepMC/ee_vqqH_ps_had_${i}_allhad250_-+.hepmc
    ./DelphesHepMC3 $CARD_NAME analysis/filesPostDelphes/qqHXHepMC/ee_vqqH_ps_had_${i}_allhad250_--.root analysis/filesPostDelphes/qqHXHepMC/ee_vqqH_ps_had_${i}_allhad250_--.hepmc
done
##########################Delphes

#######root file addition
# Define the output file name
# output_file="analysis/filesPostDelphes/qqHXHepMC/eGammaqqHXAllILDDSiDi.root"
output_file="analysis/filesPostDelphes/qqHXHepMC/eGammaqqHXAllSiD2024XCC1.root"

# Initialize an empty string to hold all the root files
root_files=""
for i in {1..3}; do
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHv_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHv_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_vqqH_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_vqqH_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHv_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHv_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_vqqH_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_vqqH_ps_had_${i}_allhad250_--.root "
done

hadd $output_file $root_files

output_file="analysis/filesPostDelphes/qqHXHepMC/eGammaqqHXAllSiD2024XCC2.root"

# Initialize an empty string to hold all the root files
root_files=""
for i in {4..5}; do
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHe_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHv_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/Ae_qqHv_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_eqqH_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_vqqH_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/eA_vqqH_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_eqqH_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_++.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHe_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHv_ps_had_${i}_allhad250_+-.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_qqHv_ps_had_${i}_allhad250_--.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_vqqH_ps_had_${i}_allhad250_-+.root "
    root_files+="analysis/filesPostDelphes/qqHXHepMC/ee_vqqH_ps_had_${i}_allhad250_--.root "
done

hadd $output_file $root_files

###########root file addition