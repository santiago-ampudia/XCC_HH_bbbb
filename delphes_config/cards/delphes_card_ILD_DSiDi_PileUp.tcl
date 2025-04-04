# ILD card modified to have the same parameters as the DSiDi card

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  PileUpMerger
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  MuonMomentumSmearing

  TrackMerger
  
  ECal
  HCal

  Calorimeter
  ElectronFilter
  TrackPileUpSubtractor
  
  EFlowMerger
  EFlowFilter
  
  Rho
  
  PhotonEfficiency
  PhotonIsolation

  ElectronEfficiency
  ElectronIsolation

  ChargedHadronFilter

  MuonEfficiency
  MuonIsolation

  NeutrinoFilter
  
  GenJetFinderAntiKt
  GenJetFinder0
  GenJetFinder5
  GenJetFinder10
  GenJetFinder15
  GenJetFinder20
  GenJetFinder25
  GenJetFinder30
  
  FastJetFinderAntiKt
  FastJetFinder0
  FastJetFinder5
  FastJetFinder10
  FastJetFinder15
  FastJetFinder20
  FastJetFinder25
  FastJetFinder30
  
  JetPileUpSubtractorAntiKt
  JetPileUpSubtractor0
  JetPileUpSubtractor5
  JetPileUpSubtractor10
  JetPileUpSubtractor15
  JetPileUpSubtractor20
  JetPileUpSubtractor25
  JetPileUpSubtractor30

  MissingET
  GenMissingET
  GenPileUpMissingET


  JetEnergyScaleAntiKt
  JetEnergyScale0
  JetEnergyScale5
  JetEnergyScale10
  JetEnergyScale15
  JetEnergyScale20
  JetEnergyScale25
  JetEnergyScale30
  
  JetEnergyScalePileUpSubtractionAntiKt
  JetEnergyScalePileUpSubtraction0
  JetEnergyScalePileUpSubtraction5
  JetEnergyScalePileUpSubtraction10
  JetEnergyScalePileUpSubtraction15
  JetEnergyScalePileUpSubtraction20
  JetEnergyScalePileUpSubtraction25
  JetEnergyScalePileUpSubtraction30

  JetFlavorAssociationAntiKt
  JetFlavorAssociation0
  JetFlavorAssociation5
  JetFlavorAssociation10
  JetFlavorAssociation15
  JetFlavorAssociation20
  JetFlavorAssociation25
  JetFlavorAssociation30
  
  JetFlavorAssociationPileUpSubtractionAntiKt
  JetFlavorAssociationPileUpSubtraction0
  JetFlavorAssociationPileUpSubtraction5
  JetFlavorAssociationPileUpSubtraction10
  JetFlavorAssociationPileUpSubtraction15
  JetFlavorAssociationPileUpSubtraction20
  JetFlavorAssociationPileUpSubtraction25
  JetFlavorAssociationPileUpSubtraction30

  BTaggingAntiKt
  BTagging0
  BTagging5
  BTagging10
  BTagging15
  BTagging20
  BTagging25
  BTagging30
  
  BTaggingPileUpSubtractionAntiKt
  BTaggingPileUpSubtraction0
  BTaggingPileUpSubtraction5
  BTaggingPileUpSubtraction10
  BTaggingPileUpSubtraction15
  BTaggingPileUpSubtraction20
  BTaggingPileUpSubtraction25
  BTaggingPileUpSubtraction30
  
  TauTaggingAntiKt
  TauTagging0
  TauTagging5
  TauTagging10
  TauTagging15
  TauTagging20
  TauTagging25
  TauTagging30
  
  TauTaggingPileUpSubtractionAntiKt
  TauTaggingPileUpSubtraction0
  TauTaggingPileUpSubtraction5
  TauTaggingPileUpSubtraction10
  TauTaggingPileUpSubtraction15
  TauTaggingPileUpSubtraction20
  TauTaggingPileUpSubtraction25
  TauTaggingPileUpSubtraction30
    
  ScalarHT

  UniqueObjectFinderAntiKt
  UniqueObjectFinder0
  UniqueObjectFinder5
  UniqueObjectFinder10
  UniqueObjectFinder15
  UniqueObjectFinder20
  UniqueObjectFinder25
  UniqueObjectFinder30
  
  UniqueObjectFinderPileUpSubtractionAntiKt
  UniqueObjectFinderPileUpSubtraction0
  UniqueObjectFinderPileUpSubtraction5
  UniqueObjectFinderPileUpSubtraction10
  UniqueObjectFinderPileUpSubtraction15
  UniqueObjectFinderPileUpSubtraction20
  UniqueObjectFinderPileUpSubtraction25
  UniqueObjectFinderPileUpSubtraction30

  TreeWriter
}

###############
# PileUp Merger
###############

module PileUpMerger PileUpMerger {
  set InputArray Delphes/stableParticles

  set ParticleOutputArray stableParticles
  set VertexOutputArray vertices

  # pre-generated minbias input file
  set PileUpFile MinBias.pileup

  # average expected pile up
  set MeanPileUp 6.81

   # maximum spread in the beam direction in m
  set ZVertexSpread 3.3E-5

  # maximum spread in time in s
  set TVertexSpread 0.11E-12

  # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s)
  set VertexDistributionFormula {exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))}


}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray PileUpMerger/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  #set Radius 1.8 ###Original
  set Radius 2.493  
  #DSiDi 
  # half-length of the magnetic field coverage, in m
  #set HalfLength 2.4 ###Original
  set HalfLength 3.018
  #DSiDi
  # magnetic field
  #set Bz 3.5  ###Original
  set Bz 5.0 
  #DSiDi
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # add EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for charged hadrons
  #set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
  #                                         (abs(eta) <= 2.4)               * (pt > 0.1)    * (0.99) +
  #                                         (abs(eta) >  2.4)                               * (0.00)}
  
  set EfficiencyFormula { (pt<=0.1)*0.0+
                          (abs(eta)<=1.32)*(pt>0.1&&pt<=0.6)*0.90+
                          (abs(eta)<=1.32)*(pt>0.6&&pt<=2.0)*0.98+
                          (abs(eta)<=1.32)*(pt>2.0&&pt<=4.0)*0.99+
                          (abs(eta)<=1.32)*(pt>4.0&&pt<=10000.)*0.99+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>0.1&&pt<=0.6)*0.95+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>0.6&&pt<=2.0)*0.99+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>2.0&&pt<=4.0)*0.98+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>4.0&&pt<=10000.)*(0.99-0.00021*(pt-4.))+
                          (abs(eta)>2.44)*0.0 } 
                          #DSiDi 
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for electrons
  #set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
  #                                         (abs(eta) <= 2.4)               * (pt > 0.1)    * (0.99) +
  #                                         (abs(eta) >  2.4)                               * (0.00)}
  
  set EfficiencyFormula { (pt<=0.1)*0.0+
                          (abs(eta)<=1.32)*(pt>0.1&&pt<=0.6)*0.90+
                          (abs(eta)<=1.32)*(pt>0.6&&pt<=2.0)*0.98+
                          (abs(eta)<=1.32)*(pt>2.0&&pt<=4.0)*0.99+
                          (abs(eta)<=1.32)*(pt>4.0&&pt<=10000.)*0.99+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>0.1&&pt<=0.6)*0.95+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>0.6&&pt<=2.0)*0.99+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>2.0&&pt<=4.0)*0.98+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>4.0&&pt<=10000.)*(0.99-0.00021*(pt-4.))+
                          (abs(eta)>2.44)*0.0 } 
                          #DSiDi
                          
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for muons
  #set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
  #                                         (abs(eta) <= 2.4)               * (pt > 0.1)    * (0.99) +
  #                                        (abs(eta) >  2.4)                               * (0.00)}
  
  set EfficiencyFormula { (pt<=0.1)*0.0+
                          (abs(eta)<=1.32)*(pt>0.1&&pt<=0.6)*0.90+
                          (abs(eta)<=1.32)*(pt>0.6&&pt<=2.0)*0.98+
                          (abs(eta)<=1.32)*(pt>2.0&&pt<=4.0)*0.99+
                          (abs(eta)<=1.32)*(pt>4.0&&pt<=10000.)*0.99+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>0.1&&pt<=0.6)*0.95+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>0.6&&pt<=2.0)*0.99+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>2.0&&pt<=4.0)*0.98+
                          (abs(eta)<=2.44&&abs(eta)>1.32)*(pt>4.0&&pt<=10000.)*(0.99-0.00021*(pt-4.))+
                          (abs(eta)>2.44)*0.0 } 
                          #DSiDi
                          
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for charged hadrons
  #set ResolutionFormula {    (abs(eta) <= 1.0)                   * sqrt(0.001^2 + pt^2*1.e-5^2) +
   #                          (abs(eta) > 1.0 && abs(eta) <= 2.4) * sqrt(0.01^2 + pt^2*1.e-4^2)}
  set ResolutionFormula {(abs(eta)<=1.32)*sqrt(0.0000146^2*pt^2+0.00217^2)+
                         (abs(eta)>1.32)*sqrt(0.0000237^2*pt^2+0.00423^2) }  
                         #DSiDi   


}

###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

   # resolution formula for charged hadrons
  #set ResolutionFormula {    (abs(eta) <= 1.0)                   * sqrt(0.001^2 + pt^2*1.e-5^2) +
   #                          (abs(eta) > 1.0 && abs(eta) <= 2.4) * sqrt(0.01^2 + pt^2*1.e-4^2)}
   
   set ResolutionFormula {(abs(eta)<=1.32)*sqrt(0.0000146^2*pt^2+0.00217^2)+
                         (abs(eta)>1.32)*sqrt(0.0000237^2*pt^2+0.00423^2) } 
                         #DSiDi
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

   # resolution formula for charged hadrons
  #set ResolutionFormula {    (abs(eta) <= 1.0)                   * sqrt(0.001^2 + pt^2*1.e-5^2) +
   #                          (abs(eta) > 1.0 && abs(eta) <= 2.4) * sqrt(0.01^2 + pt^2*1.e-4^2)}
                             
   set ResolutionFormula {(abs(eta)<=1.32)*sqrt(0.0000146^2*pt^2+0.00217^2)+
                       (abs(eta)>1.32)*sqrt(0.0000237^2*pt^2+0.00423^2) } 
                       #DSiDi                       

}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronMomentumSmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}

#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true 
 
  set EnergyMin 0.5
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 0.5 degree towers (5x5 mm^2)
  set PhiBins {}
  for {set i -360} {$i <= 360} {incr i} {
    add PhiBins [expr {$i * $pi/360.0}]
  }

  # 0.01 unit in eta up to eta = 3.0
  #for {set i -300} {$i <= 300} {incr i} {
  #  set eta [expr {$i * 0.01}]
  #  add EtaPhiBins $eta $PhiBins
  #}
  
  # 0.01 unit in eta up to eta = 2.5
  for {set i -500} {$i <= 500} {incr i} {
    set eta [expr {$i * 0.005}]
    add EtaPhiBins $eta $PhiBins
  } 
  #DSiDi

  # default energy fractions {abs(PDG code)} {fraction of energy deposited in ECAL}

  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0}
  add EnergyFraction {22} {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3}
  add EnergyFraction {3122} {0.3}

  # set ECalResolutionFormula {resolution formula as a function of eta and energy}

  #set ResolutionFormula { (abs(eta) <= 3.0)                   * sqrt(energy^2*0.01^2 + energy*0.15^2) }
  
  set ResolutionFormula {sqrt(energy^2*0.01^2 + energy*0.17^2)} 
  #DSiDi

}

#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ECal/eflowTracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false 
 
  set EnergyMin 1.0
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower


  # 6 degree towers
  set PhiBins {}
  for {set i -60} {$i <= 60} {incr i} {
    add PhiBins [expr {$i * $pi/60.0}]
  }

  # 0.5 unit in eta up to eta = 3
  for {set i -60} {$i <= 60} {incr i} {
    set eta [expr {$i * 0.05}]
    add EtaPhiBins $eta $PhiBins
  }


  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {0.0}
  add EnergyFraction {22} {0.0}
  add EnergyFraction {111} {0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.7}
  add EnergyFraction {3122} {0.7}

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}

  #set ResolutionFormula {                  (abs(eta) <= 3.0) * sqrt(energy^2*0.015^2 + energy*0.50^2)}
  
  set ResolutionFormula {sqrt(energy^2*0.094^2 + energy*0.559^2)} 
  #DSiDi

}

#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}

##########################
# Track pile-up subtractor
##########################

module TrackPileUpSubtractor TrackPileUpSubtractor {
# add InputArray InputArray OutputArray
  add InputArray HCal/eflowTracks eflowTracks
  add InputArray ElectronFilter/electrons electrons
  add InputArray MuonMomentumSmearing/muons muons

  set VertexInputArray PileUpMerger/vertices
  # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution
  # Z vertex resolution in m
  # 0.1 micron ideal for DSiDi --- 1 micron more realistic
  set ZVertexResolution {0.0000001}
}

######################
# ChargedHadronFilter
######################

module PdgCodeFilter ChargedHadronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray chargedHadrons
  
  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}



###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger Calorimeter {
# add InputArray InputArray
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  set OutputArray towers
}


####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray TrackPileUpSubtractor/eflowTracks
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
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}


##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set EnergyOutputArray energy
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
  #set excl_ymerge34 400.0
  set ExclusiveClustering false
  set rtd_min 0.0
}

module FastJetFinder GenJetFinder0 {
  set InputArray NeutrinoFilter/filteredParticles

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

module FastJetFinder GenJetFinder5 {
  set InputArray NeutrinoFilter/filteredParticles

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

module FastJetFinder GenJetFinder15 {
  set InputArray NeutrinoFilter/filteredParticles

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

module FastJetFinder GenJetFinder20 {
  set InputArray NeutrinoFilter/filteredParticles

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

module FastJetFinder GenJetFinder25 {
  set InputArray NeutrinoFilter/filteredParticles

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

module FastJetFinder GenJetFinder30 {
  set InputArray NeutrinoFilter/filteredParticles

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

#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {
# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}

###################
# Gen PileUp Missing ET
###################

module Merger GenPileUpMissingET {
# add InputArray InputArray
#  add InputArray RunPUPPI/PuppiParticles
  add InputArray ParticlePropagator/stableParticles
  set MomentumOutputArray momentum
}



############
# Jet finder
############

module FastJetFinder FastJetFinderAntiKt {
#  set InputArray Calorimeter/towers
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

module FastJetFinder FastJetFinder0 {
#  set InputArray Calorimeter/towers
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
#  set InputArray Calorimeter/towers
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
#  set InputArray Calorimeter/towers
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
#  set InputArray Calorimeter/towers
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
#  set InputArray Calorimeter/towers
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
#  set InputArray Calorimeter/towers
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
#  set InputArray Calorimeter/towers
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

module EnergyScale JetEnergyScalePileUpSubtractionAntiKt {
  set InputArray JetPileUpSubtractorAntiKt/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScalePileUpSubtraction0 {
  set InputArray JetPileUpSubtractor0/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScalePileUpSubtraction5 {
  set InputArray JetPileUpSubtractor5/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScalePileUpSubtraction10 {
  set InputArray JetPileUpSubtractor10/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScalePileUpSubtraction15 {
  set InputArray JetPileUpSubtractor15/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScalePileUpSubtraction20 {
  set InputArray JetPileUpSubtractor20/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScalePileUpSubtraction25 {
  set InputArray JetPileUpSubtractor25/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}

module EnergyScale JetEnergyScalePileUpSubtraction30 {
  set InputArray JetPileUpSubtractor30/jets
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

module JetFlavorAssociation JetFlavorAssociationPileUpSubtractionAntiKt {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScalePileUpSubtractionAntiKt/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}

module JetFlavorAssociation JetFlavorAssociationPileUpSubtraction0 {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScalePileUpSubtraction0/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}

module JetFlavorAssociation JetFlavorAssociationPileUpSubtraction5 {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScalePileUpSubtraction5/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}

module JetFlavorAssociation JetFlavorAssociationPileUpSubtraction10 {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScalePileUpSubtraction10/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}

module JetFlavorAssociation JetFlavorAssociationPileUpSubtraction15 {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScalePileUpSubtraction15/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}

module JetFlavorAssociation JetFlavorAssociationPileUpSubtraction20 {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScalePileUpSubtraction20/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}

module JetFlavorAssociation JetFlavorAssociationPileUpSubtraction25 {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScalePileUpSubtraction25/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}

module JetFlavorAssociation JetFlavorAssociationPileUpSubtraction30 {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScalePileUpSubtraction30/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}


#############
# Rho pile-up
#############

module FastJetGridMedianEstimator Rho {

  set InputArray EFlowMerger/eflow
  set RhoOutputArray rho

  # add GridRange rapmin rapmax drap dphi
  # rapmin - the minimum rapidity extent of the grid
  # rapmax - the maximum rapidity extent of the grid
  # drap - the grid spacing in rapidity
  # dphi - the grid spacing in azimuth

  add GridRange -5.0 -2.5 1.0 1.0
  add GridRange -2.5 2.5 1.0 1.0
  add GridRange 2.5 5.0 1.0 1.0

}

###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray ECal/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons
  #set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) +
  #                                         (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) +
  #                       (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.95) +
  #                       (abs(eta) > 2.5)                                   * (0.00)}
                         
  set EfficiencyFormula { (pt <= 10.0) * (0.00) +
      (abs(eta)<=1.01)*(pt > 10.0)*0.99+
      (abs(eta)>1.01&&abs(eta)<=1.32)*(pt > 10.0)*0.95+
      (abs(eta)>1.32&&abs(eta)<=2.44)*(pt > 10.0)*0.99+
      (abs(eta)>2.44)*0.0} 
      #DSiDi
}


##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowFilter/eflow
  set RhoInputArray Rho/rho

  set OutputArray photons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}

#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray TrackPileUpSubtractor/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
  #set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) +
  #                                         (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) +
  #                       (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.95) +
  #                       (abs(eta) > 2.5)                                   * (0.00)}
                         
  set EfficiencyFormula { (pt <= 10.0) * (0.00) +
      (abs(eta)<=1.01)*(pt > 10.0)*0.98+
      (abs(eta)>1.01&&abs(eta)<=1.32)*(pt > 10.0)*0.75+
      (abs(eta)>1.32&&abs(eta)<=2.44)*(pt > 10.0)*0.95+
      (abs(eta)>2.44)*0.0} 
      #DSiDi
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowFilter/eflow
  set RhoInputArray Rho/rho

  set OutputArray electrons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray TrackPileUpSubtractor/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency as a function of eta and pt}

  # efficiency formula for muons
  #set EfficiencyFormula {                                      (pt <= 10.0)               * (0.00) +
  #                                         (abs(eta) <= 1.5) * (pt > 10.0 && pt <= 1.0e3) * (0.95) +
  #                                         (abs(eta) <= 1.5) * (pt > 1.0e3)               * (0.95 * exp(0.5 - pt*5.0e-4)) +
  #                       (abs(eta) > 1.5 && abs(eta) <= 2.4) * (pt > 10.0 && pt <= 1.0e3) * (0.95) +
  #                       (abs(eta) > 1.5 && abs(eta) <= 2.4) * (pt > 1.0e3)               * (0.95 * exp(0.5 - pt*5.0e-4)) +
  #                       (abs(eta) > 2.4)                                                 * (0.00)}
                         
  set EfficiencyFormula { (pt <= 10.0) * (0.00) +
    (abs(eta)<=2.44)*(pt > 10.0)*0.98+0.0}
     #DSiDi
}

################
# Muon isolation
################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set IsolationInputArray EFlowFilter/eflow
  set RhoInputArray Rho/rho

  set OutputArray muons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.25
}


###########
# b-tagging
###########

module BTagging BTaggingAntiKt {
  set JetInputArray JetEnergyScaleAntiKt/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTagging0 {
  set JetInputArray JetEnergyScale0/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTagging5 {
  set JetInputArray JetEnergyScale5/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTagging10 {
  set JetInputArray JetEnergyScale10/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTagging15 {
  set JetInputArray JetEnergyScale15/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTagging20 {
  set JetInputArray JetEnergyScale20/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTagging25 {
  set JetInputArray JetEnergyScale25/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTagging30 {
  set JetInputArray JetEnergyScale30/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTaggingPileUpSubtractionAntiKt {
  set JetInputArray JetEnergyScalePileUpSubtractionAntiKt/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTaggingPileUpSubtraction0 {
  set JetInputArray JetEnergyScalePileUpSubtraction0/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTaggingPileUpSubtraction5 {
  set JetInputArray JetEnergyScalePileUpSubtraction5/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTaggingPileUpSubtraction10 {
  set JetInputArray JetEnergyScalePileUpSubtraction10/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTaggingPileUpSubtraction15 {
  set JetInputArray JetEnergyScalePileUpSubtraction15/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTaggingPileUpSubtraction20 {
  set JetInputArray JetEnergyScalePileUpSubtraction20/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTaggingPileUpSubtraction25 {
  set JetInputArray JetEnergyScalePileUpSubtraction25/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
}

module BTagging BTaggingPileUpSubtraction30 {
  set JetInputArray JetEnergyScalePileUpSubtraction30/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # based on arXiv:1211.4462
  
    # default efficiency formula (misidentification rate)
  #add EfficiencyFormula {0} {0.01+0.000038*pt}

  # efficiency formula for c-jets (misidentification rate)
  #add EfficiencyFormula {4} {0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))}

  # efficiency formula for b-jets
  #add EfficiencyFormula {5} {0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))}
  
  add EfficiencyFormula {0} {(abs(eta)<2.17)*0.003+0.0}
  add EfficiencyFormula {4} {(abs(eta)<2.17)*0.02+0.0}
  add EfficiencyFormula {5} {(abs(eta)<2.17)*0.80+0.0} 
  #DSiDi
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

module TauTagging TauTaggingPileUpSubtractionAntiKt {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScalePileUpSubtractionAntiKt/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTaggingPileUpSubtraction0 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScalePileUpSubtraction0/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTaggingPileUpSubtraction5 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScalePileUpSubtraction5/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTaggingPileUpSubtraction10 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScalePileUpSubtraction10/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTaggingPileUpSubtraction15 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScalePileUpSubtraction15/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTaggingPileUpSubtraction20 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScalePileUpSubtraction20/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTaggingPileUpSubtraction25 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScalePileUpSubtraction25/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

module TauTagging TauTaggingPileUpSubtraction30 {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScalePileUpSubtraction30/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
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

module UniqueObjectFinder UniqueObjectFinderPileUpSubtractionAntiKt {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScalePileUpSubtractionAntiKt/jets jets
}

module UniqueObjectFinder UniqueObjectFinderPileUpSubtraction0 {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScalePileUpSubtraction0/jets jets
}

module UniqueObjectFinder UniqueObjectFinderPileUpSubtraction5 {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScalePileUpSubtraction5/jets jets
}

module UniqueObjectFinder UniqueObjectFinderPileUpSubtraction10 {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScalePileUpSubtraction10/jets jets
}

module UniqueObjectFinder UniqueObjectFinderPileUpSubtraction15 {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScalePileUpSubtraction15/jets jets
}

module UniqueObjectFinder UniqueObjectFinderPileUpSubtraction20 {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScalePileUpSubtraction20/jets jets
}

module UniqueObjectFinder UniqueObjectFinderPileUpSubtraction25 {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScalePileUpSubtraction25/jets jets
}

module UniqueObjectFinder UniqueObjectFinderPileUpSubtraction30 {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScalePileUpSubtraction30/jets jets
}


##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle
  
  add Branch GenJetFinderAntiKt/jets GenJetAntiKt Jet
  add Branch GenJetFinder0/jets GenJet0 Jet
  add Branch GenJetFinder5/jets GenJet5 Jet
  add Branch GenJetFinder10/jets GenJet10 Jet
  add Branch GenJetFinder15/jets GenJet15 Jet
  add Branch GenJetFinder20/jets GenJet20 Jet
  add Branch GenJetFinder25/jets GenJet25 Jet
  add Branch GenJetFinder30/jets GenJet30 Jet
  add Branch GenMissingET/momentum GenMissingET MissingET
  add Branch GenPileUpMissingET/momentum GenPileUpMissingET MissingET

  add Branch TrackMerger/tracks Track Track
  add Branch Calorimeter/towers Tower Tower

  add Branch HCal/eflowTracks EFlowTrack Track
  add Branch ECal/eflowPhotons EFlowPhoton Tower
  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower
  
  add Branch EFlowMerger/eflow EFlowCandidate ParticleFlowCandidate
  
  add Branch UniqueObjectFinderAntiKt/photons Photon Photon
  add Branch UniqueObjectFinderAntiKt/electrons Electron Electron
  add Branch UniqueObjectFinderAntiKt/muons Muon Muon
  add Branch UniqueObjectFinderAntiKt/jets JetAntiKt Jet
  add Branch UniqueObjectFinder0/jets Jet0 Jet
  add Branch UniqueObjectFinder5/jets Jet5 Jet
  add Branch UniqueObjectFinder10/jets Jet10 Jet
  add Branch UniqueObjectFinder15/jets Jet15 Jet
  add Branch UniqueObjectFinder20/jets Jet20 Jet
  add Branch UniqueObjectFinder25/jets Jet25 Jet
  add Branch UniqueObjectFinder30/jets Jet30 Jet
  
  add Branch UniqueObjectFinderPileUpSubtractionAntiKt/jets JetPileUpSubtractionAntiKt Jet
  add Branch UniqueObjectFinderPileUpSubtraction0/jets JetPileUpSubtraction0 Jet
  add Branch UniqueObjectFinderPileUpSubtraction5/jets JetPileUpSubtraction5 Jet
  add Branch UniqueObjectFinderPileUpSubtraction10/jets JetPileUpSubtraction10 Jet
  add Branch UniqueObjectFinderPileUpSubtraction15/jets JetPileUpSubtraction15 Jet
  add Branch UniqueObjectFinderPileUpSubtraction20/jets JetPileUpSubtraction20 Jet
  add Branch UniqueObjectFinderPileUpSubtraction25/jets JetPileUpSubtraction25 Jet
  add Branch UniqueObjectFinderPileUpSubtraction30/jets JetPileUpSubtraction30 Jet
  
  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
  add Branch Rho/rho Rho Rho
  add Branch PileUpMerger/vertices Vertex Vertex
}
