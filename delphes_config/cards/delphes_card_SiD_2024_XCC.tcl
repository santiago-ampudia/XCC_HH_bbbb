############################################################
# DSiD: Delphes card with SiD performance parameters
# Based on papers from: Chris Potter, Dimitris Ntounis
# Reference: ILC Technical Design Report Volume 4: Detectors, https://arxiv.org/pdf/2110.09965
# Adapted by Santiago Ampudia March 2025 for XCC study
############################################################


set B 5.0
set R 2.493
set HL 3.018

set ExecutionPath {

  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  TrackMergerPre
  TrackSmearing
  TrackMerger

  ECal
  HCal

  EFlowTrackMerger
  EFlowMerger
  EFlowFilter

  PhotonEfficiency
  PhotonIsolation

  MuonFilter
  TowerMerger

  ElectronFilter
  ElectronEfficiency
  ElectronIsolation

  ChargedHadronFilter

  MuonEfficiency
  MuonIsolation

  MissingET

  NeutrinoFilter
  GenJetFinderAntiKt
  GenJetFinder10
  GenMissingET

  FastJetFinderAntiKt
  FastJetFinder0
  FastJetFinder5
  FastJetFinder10
  FastJetFinder15
  FastJetFinder20
  FastJetFinder25
  FastJetFinder30

  JetEnergyScaleAntiKt
  JetEnergyScale0
  JetEnergyScale5
  JetEnergyScale10
  JetEnergyScale15
  JetEnergyScale20
  JetEnergyScale25
  JetEnergyScale30

  JetFlavorAssociationAntiKt
  JetFlavorAssociation0
  JetFlavorAssociation5
  JetFlavorAssociation10
  JetFlavorAssociation15
  JetFlavorAssociation20
  JetFlavorAssociation25
  JetFlavorAssociation30

  BTaggingAntiKt
  BTagging0
  BTagging5
  BTagging10
  BTagging15
  BTagging20
  BTagging25
  BTagging30

  TauTaggingAntiKt
  TauTagging0
  TauTagging5
  TauTagging10
  TauTagging15
  TauTagging20
  TauTagging25
  TauTagging30

  ScalarHT

  UniqueObjectFinderAntiKt
  UniqueObjectFinder0
  UniqueObjectFinder5
  UniqueObjectFinder10
  UniqueObjectFinder15
  UniqueObjectFinder20
  UniqueObjectFinder25
  UniqueObjectFinder30

  TreeWriter
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius $R
  # half-length of the magnetic field coverage, in m
  set HalfLength $HL
  # CP outer radii of the SiD HCAL 
  # CP reference Table II-1.1
  # magnetic field
  set Bz  $B
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  
  set OutputArray chargedHadrons

  set UseMomentumVector true

  source delphes_card_SiD_2024_XCC_params/SiD_ChargedHadronTrackingEfficiency.tcl
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons

  set OutputArray electrons

  set UseMomentumVector true

  # CP reference Figure 11-3.5 left (only muon efficiencies are available)
  source delphes_card_SiD_2024_XCC_params/SiD_ChargedHadronTrackingEfficiency.tcl
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons

  set OutputArray muons

  set UseMomentumVector true

  # CP reference Figure 11-3.5 left (only muon efficiencies are available)
  source delphes_card_SiD_2024_XCC_params/SiD_ChargedHadronTrackingEfficiency.tcl
}

##############
# Track merger
##############

module Merger TrackMergerPre {
  # add InputArray InputArray
  # add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  # add InputArray ElectronMomentumSmearing/electrons
  # add InputArray MuonMomentumSmearing/muons
  add InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  add InputArray ElectronTrackingEfficiency/electrons
  add InputArray MuonTrackingEfficiency/muons  

  set OutputArray tracks
}

########################################
# Smearing for charged tracks
########################################

module TrackCovariance TrackSmearing {
  set InputArray TrackMergerPre/tracks
  
  set OutputArray tracks

  ## magnetic field
  set Bz $B

  source delphes_card_SiD_2024_XCC_params/SiD_TrackCovariance.tcl

}

##############
# Track merger
##############

module Merger TrackMerger {
  # add InputArray InputArray
  # add InputArray TimeOfFlight/tracks
  add InputArray TrackSmearing/tracks

  set OutputArray tracks
}

#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray ecalTowers
  set PhotonOutputArray photons
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons


  set IsEcal true
  set EnergyMin 0.5
  set EnergySignificanceMin 1.0
  set SmearTowerCenter true

  source delphes_card_SiD_2024_XCC_params/SiD_ECal_Binning.tcl

  source delphes_card_SiD_2024_XCC_params/SiD_ECal_EnergyFractions.tcl

  source delphes_card_SiD_2024_XCC_params/SiD_ECal_Resolution.tcl

}

#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ECal/eflowTracks
  # set TrackInputArray TrackMerger/tracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false
  set EnergyMin 1.0
  set EnergySignificanceMin 1.0
  set SmearTowerCenter true


  source delphes_card_SiD_2024_XCC_params/SiD_HCal_Binning.tcl

  source delphes_card_SiD_2024_XCC_params/SiD_HCal_EnergyFractions.tcl

  source delphes_card_SiD_2024_XCC_params/SiD_HCal_Resolution.tcl
                 
}

############################
# Jim: Energy flow track merger
############################

module Merger EFlowTrackMerger {
  # add InputArray InputArray
  # add InputArray ECal/eflowTracks #Jim test: comment this out
  add InputArray HCal/eflowTracks

  set OutputArray eflowTracks
}

####################
# Energy flow merger
####################

module Merger EFlowMerger {
  # add InputArray InputArray
  # add InputArray HCal/eflowTracks
  add InputArray EFlowTrackMerger/eflowTracks
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons

  set OutputArray eflow
}

######################
# EFlowFilter
######################

module PdgCodeFilter EFlowFilter {
  set InputArray EFlowMerger/eflow

  set OutputArray eflow
  
  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}

###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray ECal/eflowPhotons

  set OutputArray photons

  source delphes_card_SiD_2024_XCC_params/SiD_PhotonEfficiency.tcl

}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  # set IsolationInputArray EFlowMerger/eflow
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray photons

  set DeltaRMax 0.5
  set PTMin 0.5
  set PTRatioMax 0.12 
}

#################
# Muon filter
#################

module PdgCodeFilter MuonFilter {
  set InputArray EFlowTrackMerger/eflowTracks

  set OutputArray muons

  set Invert true
  add PdgCode {13}
  add PdgCode {-13}
}




###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger TowerMerger {
# add InputArray InputArray
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  
  set OutputArray towers
}

#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray EFlowTrackMerger/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}

#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray ElectronFilter/electrons

  set OutputArray electrons

  source delphes_card_SiD_2024_XCC_params/SiD_ElectronEfficiency.tcl
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons
  
  set DeltaRMax 0.5
  set PTMin 0.5
  set PTRatioMax 0.12
}

######################
# ChargedHadronFilter
######################

module PdgCodeFilter ChargedHadronFilter {
  # set InputArray HCal/eflowTracks
  set InputArray EFlowTrackMerger/eflowTracks

  set OutputArray chargedHadrons
  
  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray MuonFilter/muons

  set OutputArray muons

  source delphes_card_SiD_2024_XCC_params/SiD_MuonEfficiency.tcl

}

################
# Muon isolation
################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray muons
  
  set DeltaRMax 0.5
  set PTMin 0.5
  set PTRatioMax 0.25 
}

###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow

  set MomentumOutputArray momentum
}


#################
# Neutrino Filter
#################

module PdgCodeFilter NeutrinoFilter {
  set InputArray Delphes/stableParticles

  set OutputArray filteredParticles

  set PTMin 0.0
  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}
}



#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinderAntiKt {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets
  
  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  # set excl_ymerge34 400.0
  set ExclusiveClustering false
  set rtd_min 0.0
}

module FastJetFinder GenJetFinder10 {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets
  
  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 11
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  #set excl_ymerge34 400.0
  set ExclusiveClustering true
  set rtd_min 10.0
}

#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {
  # add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}

module FastJetFinder FastJetFinderAntiKt {
#  set InputArray TowerMerger/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  #set excl_ymerge34 400.0
  set ExclusiveClustering false
  set rtd_min 0.0
}

############
# Jet finder
############

module FastJetFinder FastJetFinder0 {
#  set InputArray TowerMerger/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 11
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  #set excl_ymerge34 400.0
  set ExclusiveClustering true
  set rtd_min 0.0
}

module FastJetFinder FastJetFinder5 {
#  set InputArray TowerMerger/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 11
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  #set excl_ymerge34 400.0
  set ExclusiveClustering true
  set rtd_min 5.0
}

module FastJetFinder FastJetFinder10 {
#  set InputArray TowerMerger/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 11
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  #set excl_ymerge34 400.0
  set ExclusiveClustering true
  set rtd_min 10.0
}

module FastJetFinder FastJetFinder15 {
#  set InputArray TowerMerger/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 11
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  #set excl_ymerge34 400.0
  set ExclusiveClustering true
  set rtd_min 15.0
}

module FastJetFinder FastJetFinder20 {
#  set InputArray TowerMerger/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 11
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  #set excl_ymerge34 400.0
  set ExclusiveClustering true
  set rtd_min 20.0
}

module FastJetFinder FastJetFinder25 {
#  set InputArray TowerMerger/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 11
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  #set excl_ymerge34 400.0
  set ExclusiveClustering true
  set rtd_min 25.0
}

module FastJetFinder FastJetFinder30 {
#  set InputArray TowerMerger/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 11
  set ParameterR 0.5
  set NJets 4
  set JetPTMin 10.0
  #set excl_ymerge34 400.0
  set ExclusiveClustering true
  set rtd_min 30.0
}

##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScaleAntiKt {
  set InputArray FastJetFinderAntiKt/jets

  set OutputArray jets

  # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScale0 {
  set InputArray FastJetFinder0/jets

  set OutputArray jets

  # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScale5 {
  set InputArray FastJetFinder5/jets

  set OutputArray jets

  # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScale10 {
  set InputArray FastJetFinder10/jets

  set OutputArray jets

  # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScale15 {
  set InputArray FastJetFinder15/jets

  set OutputArray jets

  # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScale20 {
  set InputArray FastJetFinder20/jets

  set OutputArray jets

  # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScale25 {
  set InputArray FastJetFinder25/jets

  set OutputArray jets

  # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScale30 {
  set InputArray FastJetFinder30/jets

  set OutputArray jets

  # scale formula for jets
  set ScaleFormula {1.00}
}

########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociationAntiKt {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScaleAntiKt/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5
}

module JetFlavorAssociation JetFlavorAssociation0 {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale0/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5
}

module JetFlavorAssociation JetFlavorAssociation5 {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale5/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5
}

module JetFlavorAssociation JetFlavorAssociation10 {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale10/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5
}

module JetFlavorAssociation JetFlavorAssociation15 {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale15/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5
}

module JetFlavorAssociation JetFlavorAssociation20 {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale20/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5
}

module JetFlavorAssociation JetFlavorAssociation25 {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale25/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5
}

module JetFlavorAssociation JetFlavorAssociation30 {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale30/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5
}

###########
# b-tagging
###########

module BTagging BTaggingAntiKt {
  set JetInputArray JetEnergyScaleAntiKt/jets

  set BitNumber 0

  # based on arXiv:2501.16584 -- working point for b-tagging effi. 0.85
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.00045+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.007+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.85+0.0} 
}

module BTagging BTagging0 {
  set JetInputArray JetEnergyScale0/jets

  set BitNumber 0

  # based on arXiv:2501.16584 -- working point for b-tagging effi. 0.85
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.00045+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.007+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.85+0.0} 
}

module BTagging BTagging5 {
  set JetInputArray JetEnergyScale5/jets

  set BitNumber 0

  # based on arXiv:2501.16584 -- working point for b-tagging effi. 0.85
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.00045+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.007+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.85+0.0} 
}

module BTagging BTagging10 {
  set JetInputArray JetEnergyScale10/jets

  set BitNumber 0

  # based on arXiv:2501.16584 -- working point for b-tagging effi. 0.85
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.00045+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.007+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.85+0.0} 
}

module BTagging BTagging15 {
  set JetInputArray JetEnergyScale15/jets

  set BitNumber 0

  # based on arXiv:2501.16584 -- working point for b-tagging effi. 0.85
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.00045+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.007+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.85+0.0} 
}

module BTagging BTagging20 {
  set JetInputArray JetEnergyScale20/jets

  set BitNumber 0

  # based on arXiv:2501.16584 -- working point for b-tagging effi. 0.85
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.00045+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.007+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.85+0.0} 
}

module BTagging BTagging25 {
  set JetInputArray JetEnergyScale25/jets

  set BitNumber 0

  # based on arXiv:2501.16584 -- working point for b-tagging effi. 0.85
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.00045+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.007+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.85+0.0} 
}

module BTagging BTagging30 {
  set JetInputArray JetEnergyScale30/jets

  set BitNumber 0

  # based on arXiv:2501.16584 -- working point for b-tagging effi. 0.85
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.00045+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.007+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.85+0.0} 
}

#############
# tau-tagging
#############

module TauTagging TauTaggingAntiKt {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScaleAntiKt/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTagging0 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale0/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTagging5 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale5/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTagging10 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale10/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTagging15 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale15/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTagging20 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale20/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTagging25 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale25/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTagging30 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale30/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

##################
# Scalar HT merger
##################

module Merger ScalarHT {
  # add InputArray InputArray
  add InputArray EFlowMerger/eflow

  set EnergyOutputArray energy
}

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinderAntiKt {
  # earlier arrays take precedence over later ones
  # add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScaleAntiKt/jets jets
}

module UniqueObjectFinder UniqueObjectFinder0 {
  # earlier arrays take precedence over later ones
  # add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale0/jets jets
}

module UniqueObjectFinder UniqueObjectFinder5 {
  # earlier arrays take precedence over later ones
  # add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale5/jets jets
}

module UniqueObjectFinder UniqueObjectFinder10 {
  # earlier arrays take precedence over later ones
  # add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale10/jets jets
}

module UniqueObjectFinder UniqueObjectFinder15 {
  # earlier arrays take precedence over later ones
  # add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale15/jets jets
}

module UniqueObjectFinder UniqueObjectFinder20 {
  # earlier arrays take precedence over later ones
  # add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale20/jets jets
}

module UniqueObjectFinder UniqueObjectFinder25 {
  # earlier arrays take precedence over later ones
  # add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale25/jets jets
}

module UniqueObjectFinder UniqueObjectFinder30 {
  # earlier arrays take precedence over later ones
  # add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale30/jets jets
}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
  # add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle

  add Branch EFlowTrackMerger/eflowTracks EFlowTrack Track
  add Branch TrackSmearing/tracks Track Track
  add Branch ECal/eflowPhotons EFlowPhoton Tower
  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch EFlowMerger/eflow ParticleFlowCandidate ParticleFlowCandidate
  add Branch TowerMerger/towers Tower Tower

  add Branch UniqueObjectFinder10/photons Photon Photon
  add Branch UniqueObjectFinder10/electrons Electron Electron
  add Branch UniqueObjectFinder10/muons Muon Muon
  add Branch UniqueObjectFinder10/jets JetAntiKt Jet
  add Branch UniqueObjectFinder0/jets Jet0 Jet
  add Branch UniqueObjectFinder5/jets Jet5 Jet
  add Branch UniqueObjectFinder10/jets Jet10 Jet
  add Branch UniqueObjectFinder15/jets Jet15 Jet
  add Branch UniqueObjectFinder20/jets Jet20 Jet
  add Branch UniqueObjectFinder25/jets Jet25 Jet
  add Branch UniqueObjectFinder30/jets Jet30 Jet

  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
  add Branch GenMissingET/momentum GenMissingET MissingET

  add Branch GenJetFinderAntiKt/jets GenJetAntiKt Jet
  add Branch GenJetFinder10/jets GenJet10 Jet

  # add Info InfoName InfoValue
  add Info Bz $B
}


