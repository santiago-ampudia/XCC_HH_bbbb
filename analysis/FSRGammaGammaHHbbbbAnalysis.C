/////Preselection code for the analysis

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TMatrixDSymEigen.h"
#include "EventShape/Class/interface/EventShape.h"
#include "vector"
#include <vector>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <cstdlib>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <random>
#include <set>
#include <queue>
#endif

#include "epConstrainHH.h"
#include "LorentzVectorWithErrors.h"

using namespace std;

void functionLinking()
{	
	cout<<5+9-7<<endl;
	cout<<"functionLinking working"<<endl;
}

//Function that receives the true particles branch, a particle, and a number 1 or 2 saying if it should return the PID of the mother 1 or the mother 2. If it has no mother for the number requested, it reutrns -999.
int findMother(TClonesArray *branchParticle, GenParticle *particle, int m) 
{	
  	int mIndex;
  	
  	if(m==1) 
  	{
  		if(particle->M1 == -1) return -999;
  		else mIndex=particle->M1;
  	}
  	if(m==2) 
  	{
  		if(particle->M2 == -1) return -999;
  		else mIndex=particle->M2;
  	}
  	
  	return abs(static_cast<GenParticle*>(branchParticle->At(mIndex))->PID);
}

//Function that receives the true particles branch, a particle, and a number 1 or 2 saying if it should return the index of the mother 1 or the mother 2. If it has no mother for the number requested, it reutrns -999.
int findMotherIndex(TClonesArray *branchParticle, GenParticle *particle, int m) 
{	
  	int mIndex;
  	
  	if(m==1) 
  	{
  		if(particle->M1 == -1) return -999;
  		else mIndex=particle->M1;
  	}
  	if(m==2) 
  	{
  		if(particle->M2 == -1) return -999;
  		else mIndex=particle->M2;
  	}
  	
  	return mIndex;
}

//Function that receives the true particles branch, a particle, and a number 1 or 2 saying if it should return the PID of the daughter 1 or the daughter 2. If it has no daughter for the number requested, it reutrns -999.
int findDaughter(TClonesArray *branchParticle, GenParticle *particle, int d) 
{	
  	int dIndex;
  	
  	if(d==1) 
  	{
  		if(particle->D1 == -1) return -999;
  		else dIndex=particle->D1;
  	}
  	if(d==2) 
  	{
  		if(particle->D2 == -1) return -999;
  		else dIndex=particle->D2;
  	}
  	
  	return abs(static_cast<GenParticle*>(branchParticle->At(dIndex))->PID);
}

//Function that receives the true particles branch, a particle, and a number 1 or 2 saying if it should return the index of the daughter 1 or the daughter 2. If it has no daughter for the number requested, it reutrns -999.
int findDaughterIndex(TClonesArray *branchParticle, GenParticle *particle, int d) 
{	
  	int dIndex;
  	
  	if(d==1) 
  	{
  		if(particle->D1 == -1) return -999;
  		else dIndex=particle->D1;
  	}
  	if(d==2) 
  	{
  		if(particle->D2 == -1) return -999;
  		else dIndex=particle->D2;
  	}
  	
  	return dIndex;
}

//Function that recursively goes over the daughters of a particle until it finds a quark, a lepton or a b. It also finds if the top decays into a W. Used only for ttbar.
void findFinalDecaysttbar(TClonesArray *branchParticle, GenParticle *particle, vector<int>& finalDecays, int& contWEvent, int& contLeptonEvent, int& contHadronEvent, int index) 
{	
  	int pid = abs(particle->PID);
  	finalDecays.push_back(pid);
  	int m1 = abs(findMother(branchParticle, particle, 1));
  	//cout<<endl<<endl;
  	//cout<<"Index: "<<index<<endl<<"Pid func: "<<pid<<endl<<"m1 func: "<<m1<<endl;
  	
  	
  	if((pid == 1 || pid == 2 || pid == 3 || pid == 4 || pid == 5) && m1 == 24) //lookinf for u, d, s, c quarks from W
  	{
  		contHadronEvent++;
  		//cout<<"Entra hadron"<<endl;
  		return;
  	}
  	if((pid == 11 || pid == 12 || pid == 13 || pid == 14 || pid == 15 || pid == 16 || pid == 17 || pid == 18) && m1 == 24)  //looking for leptons from W
  	{
  		contLeptonEvent++;
  		//cout<<"Entra lepton"<<endl;
  		return;
  	}
  	if(pid == 5 && m1 == 6)  //looking for bs from tops
  	{
  		//cout<<"Entra b"<<endl;
  		return;
  	}
  	if(pid == 24 && m1 == 6)  //looking for Ws from tops
  	{
  		contWEvent++;
  		//cout<<"Entra w"<<endl;
  	}
  	
  	
  	int d1Index = findDaughterIndex(branchParticle, particle, 1);
  	if(d1Index != -999)
  	{
  		GenParticle *d1 = (GenParticle*) branchParticle->At(d1Index);
  		//cout<<"Entra d1"<<endl;
  		findFinalDecaysttbar(branchParticle, d1, finalDecays, contWEvent, contLeptonEvent, contHadronEvent, d1Index);
  	}
  	int d2Index = findDaughterIndex(branchParticle, particle, 2);
  	if(d2Index != -999 && d1Index != d2Index)
  	{
  		GenParticle *d2 = (GenParticle*) branchParticle->At(d2Index);
  		//cout<<"Entra d2"<<endl;
  		findFinalDecaysttbar(branchParticle, d2, finalDecays, contWEvent, contLeptonEvent, contHadronEvent, d2Index);
  	}
  		
  	return;
}

///Function that checks if a particle is coming from a photon.
void isGammaDaughter(TClonesArray *branchParticle, GenParticle *particle, bool& isGammaDaughterN)
{
	int pid = abs(particle->PID);
  	int m1 = abs(findMother(branchParticle, particle, 1));
  	
  	int m1Index = findMotherIndex(branchParticle, particle, 1);
  	//cout<<endl<<"isGammaDaughter func: "<<endl<<"pid: "<<pid<<"    m1: "<<m1;
  	if(pid != m1 && m1 != 22)
  	{
  		isGammaDaughterN = false;
  		return;
  	}
  	else if(m1 == 22) return;
  	else if(pid == m1)
  	{
  		if(m1Index != -999)
	  	{
	  		GenParticle *m1 = (GenParticle*) branchParticle->At(m1Index);
	  		//cout<<"Entra m1"<<endl;
	  		isGammaDaughter(branchParticle, m1, isGammaDaughterN);
	  	}
  	}
  	
  	return;
}

//Function that recursively goes over the daughters of a particle until it finds a quark, a lepton or a b. It also finds if the top decays into a W. Used only for ZZ.
void findFinalDecaysZZ(TClonesArray *branchParticle, GenParticle *particle, vector<int>& finalDecays, int& contNeutrinoEvent, int& contLeptonEvent, int& contHadronEvent, int index) 
{	
  	int pid = abs(particle->PID);
  	finalDecays.push_back(pid);
  	int m1 = abs(findMother(branchParticle, particle, 1));
  	//cout<<endl<<endl;
  	//cout<<"Index: "<<index<<endl<<"Pid func: "<<pid<<endl<<"m1 func: "<<m1<<endl;
  	
  	
  	if((pid == 1 || pid == 2 || pid == 3 || pid == 4 || pid == 5) && m1 == 23) //lookinf for u, d, s, c quarks from Z
  	{
  		contHadronEvent++;
  		//cout<<"Entra hadron"<<endl;
  		return;
  	}
  	if((pid == 11 || pid == 13 || pid == 15 || pid == 17) && m1 == 23)  //looking for leptons from Z
  	{
  		contLeptonEvent++;
  		//cout<<"Entra lepton"<<endl;
  		return;
  	}
  	if((pid == 12 || pid == 14 || pid == 16 || pid == 18) && m1 == 23)  //looking for neutrino from Z
  	{
  		contNeutrinoEvent++;
  		//cout<<"Entra neutrino"<<endl;
  		return;
  	}
  	
  	
  	int d1Index = findDaughterIndex(branchParticle, particle, 1);
  	if(d1Index != -999)
  	{
  		GenParticle *d1 = (GenParticle*) branchParticle->At(d1Index);
  		//cout<<"Entra d1"<<endl;
  		findFinalDecaysZZ(branchParticle, d1, finalDecays, contNeutrinoEvent, contLeptonEvent, contHadronEvent, d1Index);
  	}
  	int d2Index = findDaughterIndex(branchParticle, particle, 2);
  	if(d2Index != -999 && d1Index != d2Index)
  	{
  		GenParticle *d2 = (GenParticle*) branchParticle->At(d2Index);
  		//cout<<"Entra d2"<<endl;
  		findFinalDecaysZZ(branchParticle, d2, finalDecays, contNeutrinoEvent, contLeptonEvent, contHadronEvent, d2Index);
  	}
  		
  	return;
}

//Function that recursively goes over the daughters of a particle until it finds a quark, a lepton or a b. It also finds if the top decays into a W. Used only for WW.
void findFinalDecaysWW(TClonesArray *branchParticle, GenParticle *particle, vector<int>& finalDecays, int& contNeutrinoEvent, int& contLeptonEvent, int& contHadronEvent, int index) 
{	
  	int pid = abs(particle->PID);
  	finalDecays.push_back(pid);
  	int m1 = abs(findMother(branchParticle, particle, 1));
  	//cout<<endl<<endl;
  	//cout<<"Index: "<<index<<endl<<"Pid func: "<<pid<<endl<<"m1 func: "<<m1<<endl;
  	
  	
  	if((pid == 1 || pid == 2 || pid == 3 || pid == 4 || pid == 5) && m1 == 24) //lookinf for u, d, s, c quarks from Z
  	{
  		contHadronEvent++;
  		//cout<<"Entra hadron"<<endl;
  		return;
  	}
  	if((pid == 11 || pid == 13 || pid == 15 || pid == 17) && m1 == 24)  //looking for leptons from Z
  	{
  		contLeptonEvent++;
  		//cout<<"Entra lepton"<<endl;
  		return;
  	}
  	if((pid == 12 || pid == 14 || pid == 16 || pid == 18) && m1 == 24)  //looking for neutrino from Z
  	{
  		contNeutrinoEvent++;
  		//cout<<"Entra neutrino"<<endl;
  		return;
  	}
  	
  	
  	int d1Index = findDaughterIndex(branchParticle, particle, 1);
  	if(d1Index != -999)
  	{
  		GenParticle *d1 = (GenParticle*) branchParticle->At(d1Index);
  		//cout<<"Entra d1"<<endl;
  		findFinalDecaysWW(branchParticle, d1, finalDecays, contNeutrinoEvent, contLeptonEvent, contHadronEvent, d1Index);
  	}
  	int d2Index = findDaughterIndex(branchParticle, particle, 2);
  	if(d2Index != -999 && d1Index != d2Index)
  	{
  		GenParticle *d2 = (GenParticle*) branchParticle->At(d2Index);
  		//cout<<"Entra d2"<<endl;
  		findFinalDecaysWW(branchParticle, d2, finalDecays, contNeutrinoEvent, contLeptonEvent, contHadronEvent, d2Index);
  	}
  		
  	return;
}



//Function that recursively goes over the mothers of a particle until it finds a H, a q, or a t. Then it returns the pid of the first decay of those. 
void findInitialMother(TClonesArray *branchParticle, GenParticle *particle, int index, int& initialMotherPid) 
{	
  	int pid = abs(particle->PID);
  	int m1 = abs(findMother(branchParticle, particle, 1));
  	
  	int m1Index = findMotherIndex(branchParticle, particle, 1);
  	GenParticle *m1Part = (GenParticle*) branchParticle->At(m1Index);
  	bool isGammaDaughter1=true;
  	isGammaDaughter(branchParticle, m1Part, isGammaDaughter1);
  	
  	//cout<<endl<<endl;
  	//cout<<"Index: "<<index<<"    Pid func: "<<pid<<endl<<"m1 index: "<<m1Index<<"    m1 pid: "<<m1<<endl<<"isGammaDaughter1: "<<isGammaDaughter1;
  	
  	if((m1 == 25 || m1 == 36 || m1 == 1 || m1 == 2 || m1 == 3 || m1 == 4 || m1 == 5 || m1 == 6 || m1 == 7 || m1 == 8) && isGammaDaughter1 == true)
  	{
  		initialMotherPid = pid;
  		return;
  	}
  	
  	
  	if(initialMotherPid == -999)
  	{
	  	if(m1Index != -999)
	  	{
	  		//cout<<"Entra m1"<<endl;
	  		findInitialMother(branchParticle, m1Part, m1Index, initialMotherPid);
	  	}
  	}
  		
  	return;
}

//Function that finds b-tagged and non-b-tagged jets and associates them to TLorentzVectors. It also finds their index in the jet branch.
void jetBTagging(TClonesArray *branchJet, int nJets, int& contBJets, int& contNBJets, TLorentzVector& jetB1, TLorentzVector& jetB2, TLorentzVector& jetB3, TLorentzVector& jetB4, int& indexb1, int& indexb2, int& indexb3, int& indexb4, int& contEtaCut, int& contCosThetaCut, float& pt1, float& pt2, float& pt3, float& pt4, float& ptSum, float& jetB1M, float& jetB2M, float& jetB3M, float& jetB4M, float& constSizeB1, float& constSizeB2, float& constSizeB3, float& constSizeB4, float& jetB1NCharged, float& jetB2NCharged, float& jetB3NCharged, float& jetB4NCharged, float& jetB1NNeutrals, float& jetB2NNeutrals, float& jetB3NNeutrals, float& jetB4NNeutrals)
{
	//cout<<"nJets: "<<nJets<<endl;
	for(int i=0;i<nJets;i++)
	  {
		Jet *jet = (Jet*) branchJet->At(i);
		TRefArray* constituents = &(jet->Constituents);
		if(jet->BTag==1)
		{
			contBJets++;
			if(i==0)
			{
				jetB1=jet->P4();
				indexb1=i;
				jetB1M=jet->Mass;
				constSizeB1 = constituents->GetSize();
				jetB1NCharged = jet->NCharged;
				jetB1NNeutrals = jet->NNeutrals;
			}	
			else if(i==1)
			{
				jetB2=jet->P4();
				indexb2=i;
				jetB2M=jet->Mass;
				constSizeB2 = constituents->GetSize();
				jetB2NCharged = jet->NCharged;
				jetB2NNeutrals = jet->NNeutrals;
			}
			else if(i==2)
			{
				jetB3=jet->P4();
				indexb3=i;
				jetB3M=jet->Mass;
				constSizeB3 = constituents->GetSize();
				jetB3NCharged = jet->NCharged;
				jetB3NNeutrals = jet->NNeutrals;
			}
			else if(i==3)
			{
				jetB4=jet->P4();
				indexb4=i;
				jetB4M=jet->Mass;
				constSizeB4 = constituents->GetSize();
				jetB4NCharged = jet->NCharged;
				jetB4NNeutrals = jet->NNeutrals;
			}
		}
		else
		{
			contNBJets++;
			if(i==0)
			{
				jetB1=jet->P4();
				indexb1=i;
				jetB1M=jet->Mass;
				constSizeB1 = constituents->GetSize();
				jetB1NCharged = jet->NCharged;
				jetB1NNeutrals = jet->NNeutrals;
			}	
			else if(i==1)
			{
				jetB2=jet->P4();
				indexb2=i;
				jetB2M=jet->Mass;
				constSizeB2 = constituents->GetSize();
				jetB2NCharged = jet->NCharged;
				jetB2NNeutrals = jet->NNeutrals;
			}
			else if(i==2)
			{
				jetB3=jet->P4();
				indexb3=i;
				jetB3M=jet->Mass;
				constSizeB3 = constituents->GetSize();
				jetB3NCharged = jet->NCharged;
				jetB3NNeutrals = jet->NNeutrals;
			}
			else if(i==3)
			{
				jetB4=jet->P4();
				indexb4=i;
				jetB4M=jet->Mass;
				constSizeB4 = constituents->GetSize();
				jetB4NCharged = jet->NCharged;
				jetB4NNeutrals = jet->NNeutrals;
			}
		}
		if(contBJets>4 || contNBJets>1) 
		{
			return;
		}
		/*
		//////////jet pt
		ptSum+=jet->PT;
		if(i==0)
		{ 
			pt1 = jet->PT;
		}
		else if(i==1)
		{
			pt2 = jet->PT;
		}
		else if(i==2)
		{ 
			pt3 = jet->PT;
		}
		else if(i==3)
		{
			pt4 = jet->PT;
		}
		/////////jet pt
			
		////////angles
		float eta, eEta, theta, cosTheta, jetEta;
		eta = jet->Eta;
		eEta = pow(TMath::E(), -eta);
		theta = 2*TMath::ATan(eEta);
		cosTheta = TMath::Cos(theta);
      		jetEta=abs(jet->Eta);
		if(jetEta<1.5) contEtaCut++;
        	if(cosTheta>-0.91 && cosTheta<0.91) contCosThetaCut++;
        	////////angles
        	*/
	  }
	  return;
}

/////function that receives eta and returns cosTheta
float findCosTheta(float eta)
{
	float cosTheta = pow(TMath::E(), -eta);
	cosTheta = 2*TMath::ATan(cosTheta);
	cosTheta = TMath::Cos(cosTheta);
	
	return cosTheta;
}

//////function that finds the mass of the jet-pair comb. with the least inv. mass
float findMinJetM(TLorentzVector jetB1, TLorentzVector jetB2, TLorentzVector jetB3, TLorentzVector jetB4)
{
	TLorentzVector jetPairBB1=jetB1+jetB2;
	TLorentzVector jetPairBB2=jetB3+jetB4;
	TLorentzVector jetPairBB3=jetB1+jetB3;
	TLorentzVector jetPairBB4=jetB2+jetB4;
	TLorentzVector jetPairBB5=jetB1+jetB4;
	TLorentzVector jetPairBB6=jetB2+jetB3;
	float minMass = TMath::Min(static_cast<float>(jetPairBB1.M()), static_cast<float>(jetPairBB2.M()));
	minMass = TMath::Min(minMass, static_cast<float>(jetPairBB3.M()));
	minMass = TMath::Min(minMass, static_cast<float>(jetPairBB4.M()));
	minMass = TMath::Min(minMass, static_cast<float>(jetPairBB5.M()));
	minMass = TMath::Min(minMass, static_cast<float>(jetPairBB6.M()));
	
	return minMass;
}

float findEventShape(TLorentzVector j1, TLorentzVector j2, TLorentzVector j3, TLorentzVector j4, string shape)
{
	vector<TVector3> vectorCollection = {
		TVector3(j1.Px(), j1.Py(), j1.Pz()),
		TVector3(j2.Px(), j2.Py(), j2.Pz()),
		TVector3(j3.Px(), j3.Py(), j3.Pz()),
		TVector3(j4.Px(), j4.Py(), j4.Pz())
	};
	TMatrixDSym momentumTensor(3);
	for(std::vector<TVector3>::const_iterator p=vectorCollection.begin(); p!=vectorCollection.end(); p++){
		for(int k=0; k<3; k++){
			for(int m=0; m<=k; m++){
				momentumTensor[k][m] += (*p)[k]*(*p)[m];
			}
		}
	}
	momentumTensor *= 1/(momentumTensor[0][0] + momentumTensor[1][1] + momentumTensor[2][2]);
	TMatrixDSymEigen eigenSystem(momentumTensor);
	TVectorD eigenvalues = eigenSystem.GetEigenValues();
	float eigenvalue1_ = eigenvalues[0];
	float eigenvalue2_ = eigenvalues[1];
	float eigenvalue3_ = eigenvalues[2];
	if(shape == "sphericity") return 1.5*(eigenvalue2_ + eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
	else if(shape == "aplanarity") return 1.5*(eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
	else return 0;
}

//////Function that finds the jet pair combinations that minimizes chiSquared. 
void findJetPairs(TLorentzVector& jetPairB1, TLorentzVector& jetPairB2, TLorentzVector jetB1, TLorentzVector jetB2, TLorentzVector jetB3, TLorentzVector jetB4, double& jetPairB1Index1, double& jetPairB1Index2, double& jetPairB2Index1, double& jetPairB2Index2)
{
	double errorJetPairMass=1; ///////Should eventually measure this (? -- CHECK)
	jetPairB1=jetB1+jetB2;
	jetPairB2=jetB3+jetB4;
	double jetChiS12=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	jetPairB1=jetB1+jetB3;
	jetPairB2=jetB2+jetB4;
	double jetChiS13=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	jetPairB1=jetB1+jetB4;
	jetPairB2=jetB2+jetB3;
	double jetChiS14=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	double minJetChiS=TMath::Min(TMath::Min(jetChiS12, jetChiS13), jetChiS14);
	if(minJetChiS==jetChiS12)
	{
		jetPairB1=jetB1+jetB2;
		jetPairB2=jetB3+jetB4;
		jetPairB1Index1=0;
		jetPairB1Index2=1;
		jetPairB2Index1=2;
		jetPairB2Index2=3;
	}
	else if(minJetChiS==jetChiS13)
	{
		jetPairB1=jetB1+jetB3;
		jetPairB2=jetB2+jetB4;
		jetPairB1Index1=0;
		jetPairB1Index2=2;
		jetPairB2Index1=1;
		jetPairB2Index2=3;
	}
	else if(minJetChiS==jetChiS14)
	{
		jetPairB1=jetB1+jetB4;
		jetPairB2=jetB2+jetB3;
		jetPairB1Index1=0;
		jetPairB1Index2=3;
		jetPairB2Index1=1;
		jetPairB2Index2=2;
	}
	
	return;
}
///////

//////// Function to calculate chi-squared
double calculateChiSquared(TClonesArray *branchJet, const vector<int>& group1, const vector<int>& group2) {
    TLorentzVector jetPairB1, jetPairB2;
    for (const auto& index : group1)
    { 
    	Jet *jet = (Jet*) branchJet->At(index);
    	jetPairB1 += jet->P4();
    }
    for (const auto& index : group2)
    { 
    	Jet *jet = (Jet*) branchJet->At(index);
    	jetPairB2 += jet->P4();
    }
    double errorJetPairMass=1;
    double chiSquared=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
    return chiSquared;
}
///////

///////Function to try all jet pair combinations
void findBestCombination(TClonesArray *branchJet, int n, vector<int>& bestGroup1, vector<int>& bestGroup2, vector<double>& chiSquares, vector<vector<int>>& group1s, vector<vector<int>>& group2s) 
{
    	double minChiSquared = numeric_limits<double>::max();

    	// Generate all possible combinations of jets into two groups
    	for (int i = 1; i < (1 << n) - 1; ++i) {
        	vector<int> group1, group2;
        	group1.push_back(0);
        	for (int j = 1; j < n; ++j) {
            		if (i & (1 << (j-1))) group1.push_back(j);
            		else group2.push_back(j);
        	}

        double chiSquared = calculateChiSquared(branchJet, group1, group2);
        if (chiSquared < minChiSquared) {
            	minChiSquared = chiSquared;
            	bestGroup1 = group1;
            	bestGroup2 = group2;
        	}
        chiSquares.push_back(chiSquared);
        group1s.push_back(group1);
        group2s.push_back(group2);
    	}
}
///////
 
////Function that, given three vectors, sorts the first one by increasing order and rearranges the entries of the second, third ones to follow the first one. For example, if I have vector1 3,1,2, vector2 9,7,8, and vector3 4,2,3, I need vector1 to be 1,2,3 and vector2 to be 7,8,9 and vector3 to be 2,3,4. However, this modified version is for one vector of double and two vectors of vectors of int.
void sortVectorTriple(std::vector<double>& vector1, std::vector<vector<int>>& vector2, std::vector<vector<int>>& vector3)
{
    // Check that all vectors are the same size
    if (vector1.size() != vector2.size() || vector1.size() != vector3.size()) {
        throw std::invalid_argument("All vectors must have the same size.");
    }
    // Step 1: Pair all vectors into a vector of tuples
    std::vector<std::tuple<double, vector<int>, vector<int>>> pairedVectors(vector1.size());
    for (size_t i = 0; i < vector1.size(); ++i) {
        pairedVectors[i] = std::make_tuple(vector1[i], vector2[i], vector3[i]);
    }

    // Step 2: Sort the vector of tuples based on the first element of the tuple
    std::sort(pairedVectors.begin(), pairedVectors.end(), 
        [](const std::tuple<double, vector<int>, vector<int>>& a, const std::tuple<double, vector<int>, vector<int>>& b) {
            return std::get<0>(a) < std::get<0>(b);
        });

    // Step 3: Unpack the sorted tuples back into the four separate vectors
    for (size_t i = 0; i < pairedVectors.size(); ++i) {
        vector1[i] = std::get<0>(pairedVectors[i]);
        vector2[i] = std::get<1>(pairedVectors[i]);
        vector3[i] = std::get<2>(pairedVectors[i]);
    }

    return;
}
////////////

/*gamma gamma -> xx analysis. For topology: 
1) HH
2) qqbar
3) ttbar
4) ZZ*/ 
void analysis(const char *inputFile, int topology, string size, float weight, string jetAlgoText, string jetAlgo, string genJetAlgo, TTree& TreeTrain, TTree& TreeTest)
{
	if(topology == 1) cout<<"HHAnalysis working!"<<endl;
	else if(topology == 2) cout<<"qqAnalysis working!"<<endl;
	else if(topology == 3) cout<<"ttAnalysis working!"<<endl;
	else if(topology == 4) cout<<"ZZAnalysis working!"<<endl;
	else if(topology == 5) cout<<"WWAnalysis working!"<<endl;
	//else if(topology == -999) return;
	
	
	functionLinking();
	
	gSystem->Load("libDelphes");

  	TChain chain("Delphes");
  	chain.Add(inputFile);
  	
  	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  	Long64_t numberOfEntries = treeReader->GetEntries();
  	//numberOfEntries=1000;
  	int pos=7;
  	int contEntriesPostFilter=0;
  	cout<<"numberOfEntries: "<<numberOfEntries<<endl;
  	
  	string histJetEtaText, histJet1EtaText, histJet2EtaText, histJet3EtaText, histJet4EtaText, histJetCosThetaText, histJet1CosThetaText, histJet2CosThetaText, histJet3CosThetaText, histJet4CosThetaText, histSumJetPtText, histJetPtText, histJet1PtText, histJet2PtText, histJet3PtText, histJet4PtText, histJetB1MText, histJetB2MText, histMinJetMText, histJet1MText, histJet2MText, histJet3MText, histJet4MText, histSpherText, histAplanText, histNParticlesText, histTotalConstSizeText, histConstSizeB1Text, histConstSizeB2Text, histConstSizeB3Text, histConstSizeB4Text, histMinConstSizeText, histNEFlowTracksText, histNEFlowPhotonsText, histNEFlowNeutralHadronsText, histNEFlowObjectsText, histJetB1NChargedText, histJetB2NChargedText, histJetB3NChargedText, histJetB4NChargedText, histJetB1NNeutralsText, histJetB2NNeutralsText, histJetB3NNeutralsText, histJetB4NNeutralsText, histJetNObjectsText, histMinJetNObjectsText, histNJetsCompareAlgosText, histJetB1MAntiKt2JetsText, histJetB2MAntiKt2JetsText, histJetB1MAntiKt3JetsText, histJetB2MAntiKt3JetsText, histJetB1MAntiKt4JetsText, histJetB2MAntiKt4JetsText, histJetB1MAntiKt5JetsText, histJetB2MAntiKt5JetsText, histJetB1MAntiKt6JetsText, histJetB2MAntiKt6JetsText; 
	if(topology == 1) 
	{
  		histJetEtaText = jetAlgoText + "Jet Eta for events with HH to bb and bb (4b)";
  		histJet1EtaText = jetAlgoText + "JetB1 Eta for events with HH to bb and bb (4b)";
  		histJet2EtaText = jetAlgoText + "JetB2 Eta for events with HH to bb and bb (4b)";
  		histJet3EtaText = jetAlgoText + "JetB3 Eta for events with HH to bb and bb (4b)";
  		histJet4EtaText = jetAlgoText + "JetB4 Eta for events with HH to bb and bb (4b)";
  		
  		histJetCosThetaText = jetAlgoText + "Jet cosTheta for events with HH to bb and bb (4b)";
  		histJet1CosThetaText = jetAlgoText + "JetB1 cosTheta for events with HH to bb and bb (4b)";
  		histJet2CosThetaText = jetAlgoText + "JetB2 cosTheta for events with HH to bb and bb (4b)";
  		histJet3CosThetaText = jetAlgoText + "JetB3 cosTheta for events with HH to bb and bb (4b)";
  		histJet4CosThetaText = jetAlgoText + "JetB4 cosTheta for events with HH to bb and bb (4b)";
  		
  		histSumJetPtText = jetAlgoText + "Sum of jet Pt for events with HH to bb and bb (4b)";
  		histJetPtText = jetAlgoText + "Jet Pt for events with HH to bb and bb (4b)";
  		histJet4PtText = jetAlgoText + "Ind. jet Pt for events with HH to bb and bb (4b)";
  		
  		histJetB1MText = jetAlgoText + "b-tagged jet-pair 1 mass (HH)";
  		histJetB2MText = jetAlgoText + "b-tagged jet-pair 2 mass (HH)";
  		histMinJetMText = jetAlgoText + "Mass of the jet-pair comb. with the least inv. mass (HH)";
  		histJet1MText = jetAlgoText + "JetB1 mass (HH)";
  		histJet2MText = jetAlgoText + "JetB2 mass (HH)";
  		histJet3MText = jetAlgoText + "JetB3 mass (HH)";
  		histJet4MText = jetAlgoText + "JetB4 mass (HH)";
  		
  		histAplanText = jetAlgoText + "Aplanarity (HH)";
  		histSpherText = jetAlgoText + "Sphericity (HH)";
  		
  		histNParticlesText = jetAlgoText + "Number of particles (HH)";
  		histTotalConstSizeText = jetAlgoText + "Number of constituents in the event (HH)";
  		histConstSizeB1Text = jetAlgoText + "Number of constituents in jetB1 (HH)";
  		histConstSizeB2Text = jetAlgoText + "Number of constituents in jetB2 (HH)";
  		histConstSizeB3Text = jetAlgoText + "Number of constituents in jetB3 (HH)";
  		histConstSizeB4Text = jetAlgoText + "Number of constituents in jetB4 (HH)";
  		histMinConstSizeText = jetAlgoText + "Number of constituents in jet with min. number of consts. (HH)";
  		
  		histNEFlowTracksText = jetAlgoText + "Number of eflow tracks (HH)";
  		histNEFlowPhotonsText = jetAlgoText + "Number of eflow photons (HH)";
  		histNEFlowNeutralHadronsText = jetAlgoText + "Number of eflow neutral hadrons (HH)";
  		histNEFlowObjectsText = jetAlgoText + "Number of total eflow objects (i.e. tracks + photons + neutral hadrons) (HH)";
  		
  		histJetB1NChargedText = jetAlgoText + "Number of charged objects in jetB1 (HH)";
  		histJetB2NChargedText = jetAlgoText + "Number of charged objects in jetB2 (HH)";
  		histJetB3NChargedText = jetAlgoText + "Number of charged objects in jetB3 (HH)";
  		histJetB4NChargedText = jetAlgoText + "Number of charged objects in jetB4 (HH)";
  		histJetB1NNeutralsText = jetAlgoText + "Number of neutral objects in jetB1 (HH)";
  		histJetB2NNeutralsText = jetAlgoText + "Number of neutral objects in jetB2 (HH)";
  		histJetB3NNeutralsText = jetAlgoText + "Number of neutral objects in jetB3 (HH)";
  		histJetB4NNeutralsText = jetAlgoText + "Number of neutral objects in jetB4 (HH)";
  		histJetNObjectsText = jetAlgoText + "Total number of eflow objects (neutral and charged) in the event (HH)";
  		histMinJetNObjectsText = jetAlgoText + "Number of eflow objects in jet with min. number of eflow objects. (HH)";
  		
  		histNJetsCompareAlgosText = jetAlgoText + "nJets in durham vs. antiKt (HH)"; 
  		
  		histJetB1MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 2 antiKt (HH)";
  		histJetB2MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 2 antiKt (HH)";
  		histJetB1MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 3 antiKt (HH)";
  		histJetB2MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 3 antiKt (HH)";
  		histJetB1MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 4 antiKt (HH)";
  		histJetB2MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 4 antiKt (HH)";
  		histJetB1MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 5 antiKt (HH)";
  		histJetB2MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 5 antiKt (HH)";
  		histJetB1MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 6 antiKt (HH)";
  		histJetB2MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 6 antiKt (HH)";
  	}
  	if(topology == 2) 
	{
		histJetEtaText = jetAlgoText + "Jet Eta for events with qq (4b)";
  		histJet1EtaText = jetAlgoText + "JetB1 Eta for events with qq (4b)";
  		histJet2EtaText = jetAlgoText + "JetB2 Eta for events with qq (4b)";
  		histJet3EtaText = jetAlgoText + "JetB3 Eta for events with qq (4b)";
  		histJet4EtaText = jetAlgoText + "JetB4 Eta for events with qq (4b)";
  		
  		histJetCosThetaText = jetAlgoText + "Jet cosTheta for events with qq (4b)";
  		histJet1CosThetaText = jetAlgoText + "JetB1 cosTheta for events with qq (4b)";
  		histJet2CosThetaText = jetAlgoText + "JetB2 cosTheta for events with qq (4b)";
  		histJet3CosThetaText = jetAlgoText + "JetB3 cosTheta for events with qq (4b)";
  		histJet4CosThetaText = jetAlgoText + "JetB4 cosTheta for events with qq (4b)";
  		
  		histSumJetPtText = jetAlgoText + "Sum of jet Pt for events with qq (4b)";
		histJetPtText = jetAlgoText + "Jet Pt for events with qq (4b)";
		histJet4PtText = jetAlgoText + "Ind. jet Pt for events with qq (4b)";
		
		histJetB1MText = jetAlgoText + "b-tagged jet-pair 1 mass (qq)";
  		histJetB2MText = jetAlgoText + "b-tagged jet-pair 2 mass (qq)";
  		histMinJetMText = jetAlgoText + "Mass of the jet-pair comb. with the least inv. mass (qq)";
  		histJet1MText = jetAlgoText + "JetB1 mass (qq)";
  		histJet2MText = jetAlgoText + "JetB2 mass (qq)";
  		histJet3MText = jetAlgoText + "JetB3 mass (qq)";
  		histJet4MText = jetAlgoText + "JetB4 mass (qq)";
  		
  		histAplanText = jetAlgoText + "Aplanarity (qq)";
  		histSpherText = jetAlgoText + "Sphericity (qq)";
  		
  		histNParticlesText = jetAlgoText + "Number of particles (qq)";
  		histTotalConstSizeText = jetAlgoText + "Number of constituents in the event (qq)";
  		histConstSizeB1Text = jetAlgoText + "Number of constituents in jetB1 (qq)";
  		histConstSizeB2Text = jetAlgoText + "Number of constituents in jetB2 (qq)";
  		histConstSizeB3Text = jetAlgoText + "Number of constituents in jetB3 (qq)";
  		histConstSizeB4Text = jetAlgoText + "Number of constituents in jetB4 (qq)";
  		histMinConstSizeText = jetAlgoText + "Number of constituents in jet with min. number of consts. (qq)";
  		
  		histNEFlowTracksText = jetAlgoText + "Number of eflow tracks (qq)";
  		histNEFlowPhotonsText = jetAlgoText + "Number of eflow photons (qq)";
  		histNEFlowNeutralHadronsText = jetAlgoText + "Number of eflow neutral hadrons (qq)";
  		histNEFlowObjectsText = jetAlgoText + "Number of total eflow objects (i.e. tracks + photons + neutral hadrons) (qq)";
  		
  		histJetB1NChargedText = jetAlgoText + "Number of charged objects in jetB1 (qq)";
  		histJetB2NChargedText = jetAlgoText + "Number of charged objects in jetB2 (qq)";
  		histJetB3NChargedText = jetAlgoText + "Number of charged objects in jetB3 (qq)";
  		histJetB4NChargedText = jetAlgoText + "Number of charged objects in jetB4 (qq)";
  		histJetB1NNeutralsText = jetAlgoText + "Number of neutral objects in jetB1 (qq)";
  		histJetB2NNeutralsText = jetAlgoText + "Number of neutral objects in jetB2 (qq)";
  		histJetB3NNeutralsText = jetAlgoText + "Number of neutral objects in jetB3 (qq)";
  		histJetB4NNeutralsText = jetAlgoText + "Number of neutral objects in jetB4 (qq)";
  		histJetNObjectsText = jetAlgoText + "Total number of eflow objects in the event (qq)";
  		histMinJetNObjectsText = jetAlgoText + "Number of eflow objects in jet with min. number of eflow objects. (qq)";
  		
  		histNJetsCompareAlgosText = jetAlgoText + "nJets in durham vs. antiKt (qq)";
  		
  		histJetB1MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 2 antiKt (qq)";
  		histJetB2MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 2 antiKt (qq)";
  		histJetB1MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 3 antiKt (qq)";
  		histJetB2MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 3 antiKt (qq)";
  		histJetB1MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 4 antiKt (qq)";
  		histJetB2MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 4 antiKt (qq)";
  		histJetB1MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 5 antiKt (qq)";
  		histJetB2MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 5 antiKt (qq)";
  		histJetB1MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 6 antiKt (qq)";
  		histJetB2MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 6 antiKt (qq)";
  	}
  	if(topology == 3) 
	{
		histJetEtaText = jetAlgoText + "Jet Eta for events with ttbar (4b)";
  		histJet1EtaText = jetAlgoText + "JetB1 Eta for events with ttbar (4b)";
  		histJet2EtaText = jetAlgoText + "JetB2 Eta for events with ttbar (4b)";
  		histJet3EtaText = jetAlgoText + "JetB3 Eta for events with ttbar (4b)";
  		histJet4EtaText = jetAlgoText + "JetB4 Eta for events with ttbar (4b)";
  		
  		histJetCosThetaText = jetAlgoText + "Jet cosTheta for events with ttbar (4b)";
  		histJet1CosThetaText = jetAlgoText + "JetB1 cosTheta for events with ttbar (4b)";
  		histJet2CosThetaText = jetAlgoText + "JetB2 cosTheta for events with ttbar (4b)";
  		histJet3CosThetaText = jetAlgoText + "JetB3 cosTheta for events with ttbar (4b)";
  		histJet4CosThetaText = jetAlgoText + "JetB4 cosTheta for events with ttbar (4b)";
  		
  		histSumJetPtText = jetAlgoText + "Sum of jet Pt for events with ttbar (4b)";
		histJetPtText = jetAlgoText + "Jet Pt for events with ttbar (4b)";
		histJet4PtText = jetAlgoText + "Ind. jet Pt for events with ttbar (4b)";
		
		histJetB1MText = jetAlgoText + "b-tagged jet-pair 1 mass (ttbar)";
  		histJetB2MText = jetAlgoText + "b-tagged jet-pair 2 mass (ttbar)";
  		histMinJetMText = jetAlgoText + "Mass of the jet-pair comb. with the least inv. mass (ttbar)";
  		histJet1MText = jetAlgoText + "JetB1 mass (ttbar)";
  		histJet2MText = jetAlgoText + "JetB2 mass (ttbar)";
  		histJet3MText = jetAlgoText + "JetB3 mass (ttbar)";
  		histJet4MText = jetAlgoText + "JetB4 mass (ttbar)";
  		
  		histAplanText = jetAlgoText + "Aplanarity (ttbar)";
  		histSpherText = jetAlgoText + "Sphericity (ttbar)";
  		
  		histNParticlesText = jetAlgoText + "Number of particles (ttbar)";
  		histTotalConstSizeText = jetAlgoText + "Number of constituents in the event (ttbar)";
  		histConstSizeB1Text = jetAlgoText + "Number of constituents in jetB1 (ttbar)";
  		histConstSizeB2Text = jetAlgoText + "Number of constituents in jetB2 (ttbar)";
  		histConstSizeB3Text = jetAlgoText + "Number of constituents in jetB3 (ttbar)";
  		histConstSizeB4Text = jetAlgoText + "Number of constituents in jetB4 (ttbar)";
  		histMinConstSizeText = jetAlgoText + "Number of constituents in jet with min. number of consts. (ttbar)";
  		
		histNEFlowTracksText = jetAlgoText + "Number of eflow tracks (ttbar)";
  		histNEFlowPhotonsText = jetAlgoText + "Number of eflow photons (ttbar)";
  		histNEFlowNeutralHadronsText = jetAlgoText + "Number of eflow neutral hadrons (ttbar)";  
  		histNEFlowObjectsText = jetAlgoText + "Number of total eflow objects (i.e. tracks + photons + neutral hadrons) (ttbar)";	
  		
  		histJetB1NChargedText = jetAlgoText + "Number of charged objects in jetB1 (ttbar)";
  		histJetB2NChargedText = jetAlgoText + "Number of charged objects in jetB2 (ttbar)";
  		histJetB3NChargedText = jetAlgoText + "Number of charged objects in jetB3 (ttbar)";
  		histJetB4NChargedText = jetAlgoText + "Number of charged objects in jetB4 (ttbar)";
  		histJetB1NNeutralsText = jetAlgoText + "Number of neutral objects in jetB1 (ttbar)";
  		histJetB2NNeutralsText = jetAlgoText + "Number of neutral objects in jetB2 (ttbar)";
  		histJetB3NNeutralsText = jetAlgoText + "Number of neutral objects in jetB3 (ttbar)";
  		histJetB4NNeutralsText = jetAlgoText + "Number of neutral objects in jetB4 (ttbar)";
  		histJetNObjectsText = jetAlgoText + "Total number of eflow objects in the event (ttbar)";
  		histMinJetNObjectsText = jetAlgoText + "Number of eflow objects in jet with min. number of eflow objects. (ttbar)";
  		
  		histNJetsCompareAlgosText = jetAlgoText + "nJets in durham vs. antiKt (ttbar)";
  		
  		histJetB1MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 2 antiKt (ttbar)";
  		histJetB2MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 2 antiKt (ttbar)";
  		histJetB1MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 3 antiKt (ttbar)";
  		histJetB2MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 3 antiKt (ttbar)";
  		histJetB1MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 4 antiKt (ttbar)";
  		histJetB2MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 4 antiKt (ttbar)";
  		histJetB1MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 5 antiKt (ttbar)";
  		histJetB2MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 5 antiKt (ttbar)";
  		histJetB1MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 6 antiKt (ttbar)";
  		histJetB2MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 6 antiKt (ttbar)";
  	}
  	if(topology == 4) 
	{
		histJetEtaText = jetAlgoText + "Jet Eta for events with ZZ (4b)";
  		histJet1EtaText = jetAlgoText + "JetB1 Eta for events with ZZ (4b)";
  		histJet2EtaText = jetAlgoText + "JetB2 Eta for events with ZZ (4b)";
  		histJet3EtaText = jetAlgoText + "JetB3 Eta for events with ZZ (4b)";
  		histJet4EtaText = jetAlgoText + "JetB4 Eta for events with ZZ (4b)";
  		
  		histJetCosThetaText = jetAlgoText + "Jet cosTheta for events with ZZ (4b)";
  		histJet1CosThetaText = jetAlgoText + "JetB1 cosTheta for events with ZZ (4b)";
  		histJet2CosThetaText = jetAlgoText + "JetB2 cosTheta for events with ZZ (4b)";
  		histJet3CosThetaText = jetAlgoText + "JetB3 cosTheta for events with ZZ (4b)";
  		histJet4CosThetaText = jetAlgoText + "JetB4 cosTheta for events with ZZ (4b)";
  		
  		histSumJetPtText = jetAlgoText + "Sum of jet Pt for events with ZZ (4b)";
		histJetPtText = jetAlgoText + "Jet Pt for events with ZZ (4b)";
		histJet4PtText = jetAlgoText + "Ind. jet Pt for events with ZZ (4b)";
		
		histJetB1MText = jetAlgoText + "b-tagged jet-pair 1 mass (ZZ)";
  		histJetB2MText = jetAlgoText + "b-tagged jet-pair 2 mass (ZZ)";
  		histMinJetMText = jetAlgoText + "Mass of the jet-pair comb. with the least inv. mass (ZZ)";
  		histJet1MText = jetAlgoText + "JetB1 mass (ZZ)";
  		histJet2MText = jetAlgoText + "JetB2 mass (ZZ)";
  		histJet3MText = jetAlgoText + "JetB3 mass (ZZ)";
  		histJet4MText = jetAlgoText + "JetB4 mass (ZZ)";
  		
  		histAplanText = jetAlgoText + "Aplanarity (ZZ)";
  		histSpherText = jetAlgoText + "Sphericity (ZZ)";
  		
  		histNParticlesText = jetAlgoText + "Number of particles (ZZ)";
  		histTotalConstSizeText = jetAlgoText + "Number of constituents in the event (ZZ)";
  		histConstSizeB1Text = jetAlgoText + "Number of constituents in jetB1 (ZZ)";
  		histConstSizeB2Text = jetAlgoText + "Number of constituents in jetB2 (ZZ)";
  		histConstSizeB3Text = jetAlgoText + "Number of constituents in jetB3 (ZZ)";
  		histConstSizeB4Text = jetAlgoText + "Number of constituents in jetB4 (ZZ)";
  		histMinConstSizeText = jetAlgoText + "Number of constituents in jet with min. number of consts. (ZZ)";
  		
  		histNEFlowTracksText = jetAlgoText + "Number of eflow tracks (ZZ)";
  		histNEFlowPhotonsText = jetAlgoText + "Number of eflow photons (ZZ)";
  		histNEFlowNeutralHadronsText = jetAlgoText + "Number of eflow neutral hadrons (ZZ)";
  		histNEFlowObjectsText = jetAlgoText + "Number of total eflow objects (i.e. tracks + photons + neutral hadrons) (ZZ)";
  		
  		histJetB1NChargedText = jetAlgoText + "Number of charged objects in jetB1 (ZZ)";
  		histJetB2NChargedText = jetAlgoText + "Number of charged objects in jetB2 (ZZ)";
  		histJetB3NChargedText = jetAlgoText + "Number of charged objects in jetB3 (ZZ)";
  		histJetB4NChargedText = jetAlgoText + "Number of charged objects in jetB4 (ZZ)";
  		histJetB1NNeutralsText = jetAlgoText + "Number of neutral objects in jetB1 (ZZ)";
  		histJetB2NNeutralsText = jetAlgoText + "Number of neutral objects in jetB2 (ZZ)";
  		histJetB3NNeutralsText = jetAlgoText + "Number of neutral objects in jetB3 (ZZ)";
  		histJetB4NNeutralsText = jetAlgoText + "Number of neutral objects in jetB4 (ZZ)";
  		histJetNObjectsText = jetAlgoText + "Total number of eflow objects in the event (ZZ)";
  		histMinJetNObjectsText = jetAlgoText + "Number of eflow objects in jet with min. number of eflow objects. (ZZ)";
  		
  		histNJetsCompareAlgosText = jetAlgoText + "nJets in durham vs. antiKt (ZZ)";
  		
  		histJetB1MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 2 antiKt (ZZ)";
  		histJetB2MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 2 antiKt (ZZ)";
  		histJetB1MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 3 antiKt (ZZ)";
  		histJetB2MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 3 antiKt (ZZ)";
  		histJetB1MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 4 antiKt (ZZ)";
  		histJetB2MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 4 antiKt (ZZ)";
  		histJetB1MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 5 antiKt (ZZ)";
  		histJetB2MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 5 antiKt (ZZ)";
  		histJetB1MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 6 antiKt (ZZ)";
  		histJetB2MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 6 antiKt (ZZ)";
  	}
	if(topology == 5) 
	{
		histJetEtaText = jetAlgoText + "Jet Eta for events with WW (4b)";
  		histJet1EtaText = jetAlgoText + "JetB1 Eta for events with WW (4b)";
  		histJet2EtaText = jetAlgoText + "JetB2 Eta for events with WW (4b)";
  		histJet3EtaText = jetAlgoText + "JetB3 Eta for events with WW (4b)";
  		histJet4EtaText = jetAlgoText + "JetB4 Eta for events with WW (4b)";
  		
  		histJetCosThetaText = jetAlgoText + "Jet cosTheta for events with WW (4b)";
  		histJet1CosThetaText = jetAlgoText + "JetB1 cosTheta for events with WW (4b)";
  		histJet2CosThetaText = jetAlgoText + "JetB2 cosTheta for events with WW (4b)";
  		histJet3CosThetaText = jetAlgoText + "JetB3 cosTheta for events with WW (4b)";
  		histJet4CosThetaText = jetAlgoText + "JetB4 cosTheta for events with WW (4b)";
  		
  		histSumJetPtText = jetAlgoText + "Sum of jet Pt for events with WW (4b)";
		histJetPtText = jetAlgoText + "Jet Pt for events with WW (4b)";
		histJet4PtText = jetAlgoText + "Ind. jet Pt for events with WW (4b)";
		
		histJetB1MText = jetAlgoText + "b-tagged jet-pair 1 mass (WW)";
  		histJetB2MText = jetAlgoText + "b-tagged jet-pair 2 mass (WW)";
  		histMinJetMText = jetAlgoText + "Mass of the jet-pair comb. with the least inv. mass (WW)";
  		histJet1MText = jetAlgoText + "JetB1 mass (WW)";
  		histJet2MText = jetAlgoText + "JetB2 mass (WW)";
  		histJet3MText = jetAlgoText + "JetB3 mass (WW)";
  		histJet4MText = jetAlgoText + "JetB4 mass (WW)";
  		
  		histAplanText = jetAlgoText + "Aplanarity (WW)";
  		histSpherText = jetAlgoText + "Sphericity (WW)";
  		
  		histNParticlesText = jetAlgoText + "Number of particles (WW)";
  		histTotalConstSizeText = jetAlgoText + "Number of constituents in the event (WW)";
  		histConstSizeB1Text = jetAlgoText + "Number of constituents in jetB1 (WW)";
  		histConstSizeB2Text = jetAlgoText + "Number of constituents in jetB2 (WW)";
  		histConstSizeB3Text = jetAlgoText + "Number of constituents in jetB3 (WW)";
  		histConstSizeB4Text = jetAlgoText + "Number of constituents in jetB4 (WW)";
  		histMinConstSizeText = jetAlgoText + "Number of constituents in jet with min. number of consts. (WW)";
  		
  		histNEFlowTracksText = jetAlgoText + "Number of eflow tracks (WW)";
  		histNEFlowPhotonsText = jetAlgoText + "Number of eflow photons (WW)";
  		histNEFlowNeutralHadronsText = jetAlgoText + "Number of eflow neutral hadrons (WW)";
  		histNEFlowObjectsText = jetAlgoText + "Number of total eflow objects (i.e. tracks + photons + neutral hadrons) (WW)";
  		
  		histJetB1NChargedText = jetAlgoText + "Number of charged objects in jetB1 (WW)";
  		histJetB2NChargedText = jetAlgoText + "Number of charged objects in jetB2 (WW)";
  		histJetB3NChargedText = jetAlgoText + "Number of charged objects in jetB3 (WW)";
  		histJetB4NChargedText = jetAlgoText + "Number of charged objects in jetB4 (WW)";
  		histJetB1NNeutralsText = jetAlgoText + "Number of neutral objects in jetB1 (WW)";
  		histJetB2NNeutralsText = jetAlgoText + "Number of neutral objects in jetB2 (WW)";
  		histJetB3NNeutralsText = jetAlgoText + "Number of neutral objects in jetB3 (WW)";
  		histJetB4NNeutralsText = jetAlgoText + "Number of neutral objects in jetB4 (WW)";
  		histJetNObjectsText = jetAlgoText + "Total number of eflow objects in the event (WW)";
  		histMinJetNObjectsText = jetAlgoText + "Number of eflow objects in jet with min. number of eflow objects. (WW)";
  		
  		histNJetsCompareAlgosText = jetAlgoText + "nJets in durham vs. antiKt (WW)";
  		
  		histJetB1MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 2 antiKt (WW)";
  		histJetB2MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 2 antiKt (WW)";
  		histJetB1MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 3 antiKt (WW)";
  		histJetB2MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 3 antiKt (WW)";
  		histJetB1MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 4 antiKt (WW)";
  		histJetB2MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 4 antiKt (WW)";
  		histJetB1MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 5 antiKt (WW)";
  		histJetB2MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 5 antiKt (WW)";
  		histJetB1MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 6 antiKt (WW)";
  		histJetB2MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 6 antiKt (WW)";
  	}  	
  	
  	TH1 *histJetEta = new TH1F("jeteta", histJetEtaText.c_str(), 142.0, -4.0, 4.0);
  	TH1 *histJet1Eta = new TH1F("jeteta1", histJet1EtaText.c_str(), 142.0, -4.0, 4.0);
  	TH1 *histJet2Eta = new TH1F("jeteta2", histJet2EtaText.c_str(), 142.0, -4.0, 4.0);
  	TH1 *histJet3Eta = new TH1F("jeteta3", histJet3EtaText.c_str(), 142.0, -4.0, 4.0);
  	TH1 *histJet4Eta = new TH1F("jeteta4", histJet4EtaText.c_str(), 142.0, -4.0, 4.0);
  	
  	TH1 *histJetCosTheta = new TH1F("jetcostheta", histJetCosThetaText.c_str(), 142.0, -1.1, 1.1);
  	TH1 *histJet1CosTheta = new TH1F("jetcostheta1", histJet1CosThetaText.c_str(), 142.0, -1.1, 1.1);
  	TH1 *histJet2CosTheta = new TH1F("jetcostheta2", histJet2CosThetaText.c_str(), 142.0, -1.1, 1.1);
  	TH1 *histJet3CosTheta = new TH1F("jetcostheta3", histJet3CosThetaText.c_str(), 142.0, -1.1, 1.1);
  	TH1 *histJet4CosTheta = new TH1F("jetcostheta4", histJet4CosThetaText.c_str(), 142.0, -1.1, 1.1);
  	
  	TH1 *histSumJetPt = new TH1F("jetSumPt", histSumJetPtText.c_str(), 142.0, 0.0, 450.0);
  	TH1 *histJetPt = new TH1F("jetPt", histJetPtText.c_str(), 142.0, 0.0, 200.0);
  	TH1 *histJet1Pt = new TH1F("jetPt1", histJet1PtText.c_str(), 142.0, 0.0, 200.0);
  	TH1 *histJet2Pt = new TH1F("jetPt2", histJet2PtText.c_str(), 142.0, 0.0, 200.0);
  	TH1 *histJet3Pt = new TH1F("jetPt3", histJet3PtText.c_str(), 142.0, 0.0, 200.0);
  	TH1 *histJet4Pt = new TH1F("jetPt4", histJet4PtText.c_str(), 142.0, 0.0, 200.0);
  	histJet1Pt->SetLineColor(kBlue);
  	histJet2Pt->SetLineColor(kCyan);
  	histJet3Pt->SetLineColor(kGreen+2);
  	histJet4Pt->SetLineColor(kGreen);
  	
  	TH1 *histJetB1M = new TH1F("jetb1,", histJetB1MText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB2M = new TH1F("jetb2,", histJetB2MText.c_str(), 142.0, -1.0, 350);
  	TH1 *histMinJetM = new TH1F("minjetm", histMinJetMText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJet1M = new TH1F("jet1m,", histJet1MText.c_str(), 142.0, -1.0, 115);
  	TH1 *histJet2M = new TH1F("jet2m,", histJet2MText.c_str(), 142.0, -1.0, 115);
  	TH1 *histJet3M = new TH1F("jet3m,", histJet3MText.c_str(), 142.0, -1.0, 115);
  	TH1 *histJet4M = new TH1F("jet4m,", histJet4MText.c_str(), 142.0, -1.0, 115);
  	
  	TH1 *histSpher = new TH1F("sphericity", histSpherText.c_str(), 142.0, 0.0, 1.0);
	TH1 *histAplan = new TH1F("aplanarity", histAplanText.c_str(), 142.0, 0.0, 0.5);
	
	TH1 *histNParticles = new TH1F("nparticles", histNParticlesText.c_str(), 142.0, 0.0, 600);
	TH1 *histTotalConstSize = new TH1F("totalconstsize", histTotalConstSizeText.c_str(), 142.0, 0.0, 600);
	TH1 *histConstSizeB1 = new TH1F("constsizeb1", histConstSizeB1Text.c_str(), 142.0, 0.0, 200);
	TH1 *histConstSizeB2 = new TH1F("constsizeb2", histConstSizeB2Text.c_str(), 142.0, 0.0, 200);
	TH1 *histConstSizeB3 = new TH1F("constsizeb3", histConstSizeB3Text.c_str(), 142.0, 0.0, 200);
	TH1 *histConstSizeB4 = new TH1F("constsizeb4", histConstSizeB4Text.c_str(), 142.0, 0.0, 200);
	TH1 *histMinConstSize = new TH1F("minconstsize", histMinConstSizeText.c_str(), 142.0, 0.0, 200);
	
	TH1 *histNEFlowTracks = new TH1F("NEFlowTracks", histNEFlowTracksText.c_str(), 142.0, 0.0, 200);
	TH1 *histNEFlowPhotons = new TH1F("NEFlowPhotons", histNEFlowPhotonsText.c_str(), 142.0, 0.0, 200);
	TH1 *histNEFlowNeutralHadrons = new TH1F("NEFlowNeutralHadrons", histNEFlowNeutralHadronsText.c_str(), 142.0, 0.0, 200);
	TH1 *histNEFlowObjects = new TH1F("NEFlowObjects", histNEFlowObjectsText.c_str(), 142.0, 0.0, 200);
	
	TH1 *histJetB1NCharged = new TH1F("jetB1NCharged", histJetB1NChargedText.c_str(), 142.0, 0.0, 200);
	TH1 *histJetB2NCharged = new TH1F("jetB2NCharged", histJetB2NChargedText.c_str(), 142.0, 0.0, 200);
	TH1 *histJetB3NCharged = new TH1F("jetB3NCharged", histJetB3NChargedText.c_str(), 142.0, 0.0, 200);
	TH1 *histJetB4NCharged = new TH1F("jetB4NCharged", histJetB4NChargedText.c_str(), 142.0, 0.0, 200);
	TH1 *histJetB1NNeutrals = new TH1F("jetB1NNeutrals", histJetB1NNeutralsText.c_str(), 142.0, 0.0, 200);
	TH1 *histJetB2NNeutrals = new TH1F("jetB2NNeutrals", histJetB2NNeutralsText.c_str(), 142.0, 0.0, 200);
	TH1 *histJetB3NNeutrals = new TH1F("jetB3NNeutrals", histJetB3NNeutralsText.c_str(), 142.0, 0.0, 200);
	TH1 *histJetB4NNeutrals = new TH1F("jetB4NNeutrals", histJetB4NNeutralsText.c_str(), 142.0, 0.0, 200);
	TH1 *histJetNObjects = new TH1F("jetNObjects", histJetNObjectsText.c_str(), 142.0, 0.0, 200);
	TH1 *histMinJetNObjects = new TH1F("minJetNObjects", histMinJetNObjectsText.c_str(), 142.0, 0.0, 200);
	
	TH2F *histNJetsCompareAlgos = new TH2F("histnjjetscomparealgos", histNJetsCompareAlgosText.c_str(), 100, -1.0, 7.0, 100, -1.0, 7.0);
  	histNJetsCompareAlgos->GetXaxis()->SetTitle("Number of reco. jets with durham");
  	histNJetsCompareAlgos->GetYaxis()->SetTitle("Number of reco. jets with antiKt");
  	
  	TH1 *histJetB1MAntiKt2Jets = new TH1F("jetb12,", histJetB1MAntiKt2JetsText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB2MAntiKt2Jets = new TH1F("jetb22,", histJetB2MAntiKt2JetsText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB1MAntiKt3Jets = new TH1F("jetb13,", histJetB1MAntiKt3JetsText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB2MAntiKt3Jets = new TH1F("jetb23,", histJetB2MAntiKt3JetsText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB1MAntiKt4Jets = new TH1F("jetb14,", histJetB1MAntiKt4JetsText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB2MAntiKt4Jets = new TH1F("jetb24,", histJetB2MAntiKt4JetsText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB1MAntiKt5Jets = new TH1F("jetb15,", histJetB1MAntiKt5JetsText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB2MAntiKt5Jets = new TH1F("jetb25,", histJetB2MAntiKt5JetsText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB1MAntiKt6Jets = new TH1F("jetb16,", histJetB1MAntiKt6JetsText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB2MAntiKt6Jets = new TH1F("jetb26,", histJetB2MAntiKt6JetsText.c_str(), 142.0, -1.0, 350);
	
  	
  	TClonesArray *branchParticle;
  	TClonesArray *branchEvent;
  	TClonesArray *branchJetAntiKt;
  	TClonesArray *branchGenJetAntiKt;
  	TClonesArray *branchElectronAntiKt;
  	TClonesArray *branchMuonAntiKt;
  	TClonesArray *branchJetDurham;
  	TClonesArray *branchGenJetDurham;
  	TClonesArray *branchElectronDurham;
  	TClonesArray *branchMuonDurham;
  	TClonesArray *branchEFlowTrack;
  	TClonesArray *branchEFlowPhoton;
  	TClonesArray *branchEFlowNeutralHadron;
  	
  	
	branchParticle = treeReader->UseBranch("Particle");
	branchEvent = treeReader->UseBranch("Event");
	branchJetAntiKt = treeReader->UseBranch("JetAntiKt");
	//branchGenJetAntiKt = treeReader->UseBranch("GenJet");
	branchElectronAntiKt = treeReader->UseBranch("Electron");
	branchMuonAntiKt = treeReader->UseBranch("Muon");
	branchJetDurham = treeReader->UseBranch(jetAlgo.c_str());
	branchGenJetDurham = treeReader->UseBranch(genJetAlgo.c_str());
	//branchElectronDurham = treeReader->UseBranch("ElectronDurham");
	//branchMuonDurham = treeReader->UseBranch("MuonDurham");
	branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
	branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
	branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  	
	
	int trueFullTopEvents=0, trueSemiTopEvents=0, trueFullZEvents=0, trueSemiZEvents=0, trueFullWEvents=0, trueSemiWWEvents=0, trueFullWWEvents=0, trueSemiWEvents=0, trueFullLeptonEvents=0, trueSemiLeptonEvents=0, trueFullHadronEvents=0, trueSemiHadronEvents=0, trueFullNeutrinoEvents=0, trueSemiNeutrinoEvents=0, weirdDecays=0, contEtaCut=0, contCosThetaCut=0, trueLeptonEvents=0, trueHadronEvents=0, trueNeutrinoEvents=0;
	int trueLeptonNeutrinoEvents=0, trueLeptonHadronEvents=0, trueHadronNeutrinoEvents=0; 
	float pt1=0, pt2=0, pt3=0, pt4=0, ptSum=0, contMassB1Window=0, contMassB2Window=0, windowBLeft=0, windowBRight=0, windowB2Left=0, windowB2Right=0, contEventsWindow=0, contEvents=0, contEventsPostFilter=0;
	
	float aplanarity, invMassB1, invMassB2, minJetM, sphericity, cosThetaB1, cosThetaB2, cosThetaB3, cosThetaB4, sumPt, jetB1Pt, jetB2Pt, jetB3Pt, jetB4Pt, jetB1M, jetB2M, jetB3M, jetB4M, etaB1, etaB2, etaB3, etaB4, nParticles, totalConstSize, constSizeB1, constSizeB2, constSizeB3, constSizeB4, minConstSize, jetB1NCharged, jetB2NCharged, jetB3NCharged, jetB4NCharged, jetB1NNeutrals, jetB2NNeutrals, jetB3NNeutrals, jetB4NNeutrals, jetNObjects, minJetNObjects, invMassB1AntiKt, invMassB2AntiKt, invMassB1AntiKt2Jets, invMassB2AntiKt2Jets, invMassB1AntiKt3Jets, invMassB2AntiKt3Jets, invMassB1AntiKt4Jets, invMassB2AntiKt4Jets, invMassB1AntiKt5Jets, invMassB2AntiKt5Jets, invMassB1AntiKt6Jets, invMassB2AntiKt6Jets, nJetsAntiKt, invMassB11Best, invMassB21Best, invMassB12Best, invMassB22Best, invMassB13Best, invMassB23Best, invMassB14Best, invMassB24Best, invMassB15Best, invMassB25Best, invMassB16Best, invMassB26Best, invMassB17Best, invMassB27Best, invMassB18Best, invMassB28Best;
	
	int contLargeBMass=0, contSmallBMass=0;
	
	/////Filling out tree for BDT
	TreeTrain.Branch("aplanarity",&aplanarity,"aplanarity/F");
	TreeTrain.Branch("invMassB1",&invMassB1,"invMassB1/F");
	TreeTrain.Branch("invMassB2",&invMassB2,"invMassB2/F");
	TreeTrain.Branch("minJetM",&minJetM,"minJetM/F");  
	TreeTrain.Branch("sphericity",&sphericity,"sphericity/F");
	TreeTrain.Branch("cosThetaB1",&cosThetaB1,"cosThetaB1/F");
	TreeTrain.Branch("cosThetaB2",&cosThetaB2,"cosThetaB2/F");
	TreeTrain.Branch("cosThetaB3",&cosThetaB3,"cosThetaB3/F");
	TreeTrain.Branch("cosThetaB4",&cosThetaB4,"cosThetaB4/F");
	TreeTrain.Branch("sumPt",&sumPt,"sumPt/F");
	TreeTrain.Branch("jetB1Pt",&jetB1Pt,"jetB1Pt/F");
	TreeTrain.Branch("jetB2Pt",&jetB2Pt,"jetB2Pt/F");
	TreeTrain.Branch("jetB3Pt",&jetB3Pt,"jetB3Pt/F");
	TreeTrain.Branch("jetB4Pt",&jetB4Pt,"jetB4Pt/F");
	TreeTrain.Branch("jetB1M",&jetB1M,"jetB1M/F");
	TreeTrain.Branch("jetB2M",&jetB2M,"jetB2M/F");
	TreeTrain.Branch("jetB3M",&jetB3M,"jetB3M/F");
	TreeTrain.Branch("jetB4M",&jetB4M,"jetB4M/F");
	TreeTrain.Branch("nParticles",&nParticles,"nParticles/F");
	TreeTrain.Branch("constSizeB1",&constSizeB1,"constSizeB1/F");
	TreeTrain.Branch("constSizeB2",&constSizeB2,"constSizeB2/F");
	TreeTrain.Branch("constSizeB3",&constSizeB3,"constSizeB3/F");
	TreeTrain.Branch("constSizeB4",&constSizeB4,"constSizeB4/F");
	TreeTrain.Branch("minConstSize",&minConstSize,"minConstSize/F");
	TreeTrain.Branch("jetNObjects",&jetNObjects,"jetNObjects/F");
	TreeTrain.Branch("minJetNObjects",&minJetNObjects,"minJetNObjects/F");
	TreeTrain.Branch("invMassB1AntiKt",&invMassB1AntiKt,"invMassB1AntiKt/F");
	TreeTrain.Branch("invMassB2AntiKt",&invMassB2AntiKt,"invMassB2AntiKt/F");
	TreeTrain.Branch("nJetsAntiKt",&nJetsAntiKt,"nJetsAntiKt/F");
	TreeTrain.Branch("invMassB11Best",&invMassB11Best,"invMassB11Best/F");
	TreeTrain.Branch("invMassB21Best",&invMassB21Best,"invMassB21Best/F");
	TreeTrain.Branch("invMassB12Best",&invMassB12Best,"invMassB12Best/F");
	TreeTrain.Branch("invMassB22Best",&invMassB22Best,"invMassB22Best/F");
	TreeTrain.Branch("invMassB13Best",&invMassB13Best,"invMassB13Best/F");
	TreeTrain.Branch("invMassB23Best",&invMassB23Best,"invMassB23Best/F");
	TreeTrain.Branch("invMassB14Best",&invMassB14Best,"invMassB14Best/F");
	TreeTrain.Branch("invMassB24Best",&invMassB24Best,"invMassB24Best/F");
	TreeTrain.Branch("invMassB15Best",&invMassB15Best,"invMassB15Best/F");
	TreeTrain.Branch("invMassB25Best",&invMassB25Best,"invMassB25Best/F");
	TreeTrain.Branch("invMassB16Best",&invMassB16Best,"invMassB16Best/F");
	TreeTrain.Branch("invMassB26Best",&invMassB26Best,"invMassB26Best/F");
	TreeTrain.Branch("invMassB17Best",&invMassB17Best,"invMassB17Best/F");
	TreeTrain.Branch("invMassB27Best",&invMassB27Best,"invMassB27Best/F");
	TreeTrain.Branch("invMassB18Best",&invMassB18Best,"invMassB18Best/F");
	TreeTrain.Branch("invMassB28Best",&invMassB28Best,"invMassB28Best/F");

	
	TreeTest.Branch("aplanarity",&aplanarity,"aplanarity/F");
	TreeTest.Branch("invMassB1",&invMassB1,"invMassB1/F");
	TreeTest.Branch("invMassB2",&invMassB2,"invMassB2/F");
	TreeTest.Branch("minJetM",&minJetM,"minJetM/F");  
	TreeTest.Branch("sphericity",&sphericity,"sphericity/F");
	TreeTest.Branch("cosThetaB1",&cosThetaB1,"cosThetaB1/F");
	TreeTest.Branch("cosThetaB2",&cosThetaB2,"cosThetaB2/F");
	TreeTest.Branch("cosThetaB3",&cosThetaB3,"cosThetaB3/F");
	TreeTest.Branch("cosThetaB4",&cosThetaB4,"cosThetaB4/F");
	TreeTest.Branch("sumPt",&sumPt,"sumPt/F");
	TreeTest.Branch("jetB1Pt",&jetB1Pt,"jetB1Pt/F");
	TreeTest.Branch("jetB2Pt",&jetB2Pt,"jetB2Pt/F");
	TreeTest.Branch("jetB3Pt",&jetB3Pt,"jetB3Pt/F");
	TreeTest.Branch("jetB4Pt",&jetB4Pt,"jetB4Pt/F");
	TreeTest.Branch("jetB1M",&jetB1M,"jetB1M/F");
	TreeTest.Branch("jetB2M",&jetB2M,"jetB2M/F");
	TreeTest.Branch("jetB3M",&jetB3M,"jetB3M/F");
	TreeTest.Branch("jetB4M",&jetB4M,"jetB4M/F");
	TreeTest.Branch("nParticles",&nParticles,"nParticles/F");
	TreeTest.Branch("constSizeB1",&constSizeB1,"constSizeB1/F");
	TreeTest.Branch("constSizeB2",&constSizeB2,"constSizeB2/F");
	TreeTest.Branch("constSizeB3",&constSizeB3,"constSizeB3/F");
	TreeTest.Branch("constSizeB4",&constSizeB4,"constSizeB4/F");
	TreeTest.Branch("minConstSize",&minConstSize,"minConstSize/F");
	TreeTest.Branch("jetNObjects",&jetNObjects,"jetNObjects/F");
	TreeTest.Branch("minJetNObjects",&minJetNObjects,"minJetNObjects/F");
	TreeTest.Branch("invMassB1AntiKt",&invMassB1AntiKt,"invMassB1AntiKt/F");
	TreeTest.Branch("invMassB2AntiKt",&invMassB2AntiKt,"invMassB2AntiKt/F");
	TreeTest.Branch("nJetsAntiKt",&nJetsAntiKt,"nJetsAntiKt/F");
	TreeTest.Branch("invMassB11Best",&invMassB11Best,"invMassB11Best/F");
	TreeTest.Branch("invMassB21Best",&invMassB21Best,"invMassB21Best/F");
	TreeTest.Branch("invMassB12Best",&invMassB12Best,"invMassB12Best/F");
	TreeTest.Branch("invMassB22Best",&invMassB22Best,"invMassB22Best/F");
	TreeTest.Branch("invMassB13Best",&invMassB13Best,"invMassB13Best/F");
	TreeTest.Branch("invMassB23Best",&invMassB23Best,"invMassB23Best/F");
	TreeTest.Branch("invMassB14Best",&invMassB14Best,"invMassB14Best/F");
	TreeTest.Branch("invMassB24Best",&invMassB24Best,"invMassB24Best/F");
	TreeTest.Branch("invMassB15Best",&invMassB15Best,"invMassB15Best/F");
	TreeTest.Branch("invMassB25Best",&invMassB25Best,"invMassB25Best/F");
	TreeTest.Branch("invMassB16Best",&invMassB16Best,"invMassB16Best/F");
	TreeTest.Branch("invMassB26Best",&invMassB26Best,"invMassB26Best/F");
	TreeTest.Branch("invMassB17Best",&invMassB17Best,"invMassB17Best/F");
	TreeTest.Branch("invMassB27Best",&invMassB27Best,"invMassB27Best/F");
	TreeTest.Branch("invMassB18Best",&invMassB18Best,"invMassB18Best/F");
	TreeTest.Branch("invMassB28Best",&invMassB28Best,"invMassB28Best/F");

	
	for(int entry=0; entry<numberOfEntries; entry++) 
	{
		treeReader->ReadEntry(entry);
		
		nJetsAntiKt=0;
		if(branchJetAntiKt->GetEntries() > 0)  nJetsAntiKt =  branchJetAntiKt->GetEntries();
		int nJetsDurham=0;
		if(branchJetDurham->GetEntries() > 0)  nJetsDurham =  branchJetDurham->GetEntries();
		
		//if(nJetsAntiKt == 4) 
		if(nJetsDurham == 4)
		{
			
			int contBJetsAntiKt=0, contBJetsDurham=0;
			int contB2JetsAntiKt=0, contNBJetsDurham=0;
			TLorentzVector jetB1AntiKt, jetB2AntiKt, jetB3AntiKt, jetB4AntiKt;
			TLorentzVector jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham;
			int indexB1AntiKt=0, indexB2AntiKt=0, indexB3AntiKt=0, indexB4AntiKt=0;
			int indexB1Durham=0, indexB2Durham=0, indexB3Durham=0, indexB4Durham=0;
			totalConstSize=0;
			constSizeB1=0;
			constSizeB2=0; 
			constSizeB3=0; 
			constSizeB4=0;
			
			//jetBTagging(branchJetAntiKt, nJetsAntiKt, contBJetsAntiKt, contB2JetsAntiKt, jetB1AntiKt, jetB2AntiKt, jetB3AntiKt, jetB4AntiKt, indexB1AntiKt, indexB2AntiKt, indexB3AntiKt, indexB4AntiKt, contEtaCut, contCosThetaCut, pt1, pt2, pt3, pt4, ptSum, jetB1M, jetB2M, jetB3M, jetB4M, constSizeB1, constSizeB2, constSizeB3, constSizeB4, jetB1NCharged, jetB2NCharged, jetB3NCharged, jetB4NCharged, jetB1NNeutrals, jetB2NNeutrals, jetB3NNeutrals, jetB4NNeutrals);
			jetBTagging(branchJetDurham, nJetsDurham, contBJetsDurham, contNBJetsDurham, jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham, indexB1Durham, indexB2Durham, indexB3Durham, indexB4Durham, contEtaCut, contCosThetaCut, pt1, pt2, pt3, pt4, ptSum, jetB1M, jetB2M, jetB3M, jetB4M, constSizeB1, constSizeB2, constSizeB3, constSizeB4, jetB1NCharged, jetB2NCharged, jetB3NCharged, jetB4NCharged, jetB1NNeutrals, jetB2NNeutrals, jetB3NNeutrals, jetB4NNeutrals);
			
	      		//if(contBJetsAntiKt == 4)
	      		if(contBJetsDurham == 4 || (contBJetsDurham == 3 && contNBJetsDurham == 1))
	      		//if(contBJetsDurham == 4)
	      		{ 
		      		
		      		nParticles=0;
		      		contEntriesPostFilter++;
		      		//if(contEntriesPostFilter>1) return;
		      		
				if(branchParticle->GetEntries() > 0)  nParticles =  branchParticle->GetEntries();
				
				int contTopEvent=0, contZEvent=0, contWEvent=0, contLeptonEvent=0, contHadronEvent=0, contNeutrinoEvent=0, contG=0, contB=0, contWWEvent=0;
				
				//////Looking only for HH->bbbb for signal
				if(topology == 1)
				{
					for(int i=0;i<nParticles;i++)
		    			{
			  			GenParticle *particle = (GenParticle*) branchParticle->At(i); 
			  			int pid = abs(particle->PID);
			  			if(pid == 25 || pid == 36)
			  			{
					  		if(findDaughter(branchParticle, particle, 1) == 5 && findDaughter(branchParticle, particle, 2) == 5) contB++;
					  	}
		     			}
	     			}
	     			///////Looking only for HH->bbbb for signal
				
				///////Lepton veto for ttbar
				if(topology == 3)
				{
					for(int i=0;i<nParticles;i++)
			    		{
			  			GenParticle *particle = (GenParticle*) branchParticle->At(i); 
			  			int pid = abs(particle->PID);
				  		vector<int> finalDecays;
				  		if(pid == 6 && findMother(branchParticle, particle, 1) == 22 && findMother(branchParticle, particle, 2) == 22)
				  		{
					  		findFinalDecaysttbar(branchParticle, particle, finalDecays, contWEvent, contLeptonEvent, contHadronEvent, i);
					  		contTopEvent++;
					  	}
			     		}
			     		if(contTopEvent == 1) trueSemiTopEvents++;
			     		if(contTopEvent >= 2) trueFullTopEvents++;
			     		if(contWEvent == 1) trueSemiWEvents++;
			     		if(contWEvent >= 2) trueFullWEvents++;
			     		if(contLeptonEvent == 2) trueSemiLeptonEvents++;
			     		if(contLeptonEvent == 4) trueFullLeptonEvents++;
			     		if(contHadronEvent == 2) trueSemiHadronEvents++;
			     		if(contHadronEvent == 4) trueFullHadronEvents++;
	     			}
	     			///////Lepton veto for ttbar
	     			
	     			///////Lepton veto for ZZ
				if(topology == 4)
				{
					for(int i=0;i<nParticles;i++)
			    		{
			  			GenParticle *particle = (GenParticle*) branchParticle->At(i); 
			  			int pid = abs(particle->PID);
				  		vector<int> finalDecays;
				  		if((pid == 23 || pid == 32 || pid == 33) && findMother(branchParticle, particle, 1) == 11 && findMother(branchParticle, particle, 2) == 11)
				  		{
					  		findFinalDecaysZZ(branchParticle, particle, finalDecays, contNeutrinoEvent, contLeptonEvent, contHadronEvent, i);
					  		contZEvent++;
					  	}
			     		}
			     		if(contLeptonEvent>0) trueLeptonEvents++;
			     		if(contHadronEvent>0) trueHadronEvents++;
			     		if(contNeutrinoEvent>0) trueNeutrinoEvents++;
			     		
			     		if(contZEvent == 1) trueSemiZEvents++;
			     		if(contZEvent >= 2) trueFullZEvents++;
			     		if(contLeptonEvent == 2 && contHadronEvent == 2) trueLeptonHadronEvents++;
			     		if(contLeptonEvent == 2 && contNeutrinoEvent == 2) trueLeptonNeutrinoEvents++;
			     		if(contNeutrinoEvent == 2 && contHadronEvent == 2) trueHadronNeutrinoEvents++;
			     		if(contLeptonEvent == 4) trueFullLeptonEvents++;
			     		if(contHadronEvent == 4) trueFullHadronEvents++;
			     		if(contNeutrinoEvent == 4) trueFullNeutrinoEvents++;
	     			}
	     			///////Lepton veto for ZZ
	     			
	     			///////Lepton veto for WW
				if(topology == 5)
				{
					for(int i=0;i<nParticles;i++)
			    		{
			  			GenParticle *particle = (GenParticle*) branchParticle->At(i); 
			  			int pid = abs(particle->PID);
				  		vector<int> finalDecays;
				  		if((pid == 24) && findMother(branchParticle, particle, 1) == 11)
				  		{
					  		findFinalDecaysWW(branchParticle, particle, finalDecays, contNeutrinoEvent, contLeptonEvent, contHadronEvent, i);
					  		contWWEvent++;
					  	}
			     		}
			     		if(contLeptonEvent>0) trueLeptonEvents++;
			     		if(contHadronEvent>0) trueHadronEvents++;
			     		if(contNeutrinoEvent>0) trueNeutrinoEvents++;
			     		
			     		
			     		if(contWWEvent == 1) trueSemiWWEvents++;
			     		if(contWWEvent >= 2) trueFullWWEvents++;
			     		if(contLeptonEvent == 1 && contNeutrinoEvent == 1 && contHadronEvent == 2) trueLeptonHadronEvents++;
			     		if(contLeptonEvent == 2 && contNeutrinoEvent == 2) trueLeptonNeutrinoEvents++;
			     		if(contNeutrinoEvent == 2 && contHadronEvent == 2) trueHadronNeutrinoEvents++;
			     		if(contLeptonEvent == 4) trueFullLeptonEvents++;
			     		if(contHadronEvent == 4) trueFullHadronEvents++;
			     		if(contNeutrinoEvent == 4) trueFullNeutrinoEvents++;
			     		
			     		
	     			}
	     			///////Lepton veto for WW
	     			
	     			if(topology != 1)
	     			{
	     				contB=2;
	     			}
	     			
	     			if(contB>=2)
		     		{
		     			contEvents++;
		     			
		     			if(contLeptonEvent == 0) //lepton veto
		     			{
			     			contEventsPostFilter++;
			     			histNJetsCompareAlgos->Fill(nJetsDurham, nJetsAntiKt);
			     			
			     			//////inv. mass
			     			double jetPairB1Index1, jetPairB1Index2, jetPairB2Index1, jetPairB2Index2;
			     			TLorentzVector jetPairB1, jetPairB2, jetPairB1Sen, jetPairB2Sen, jetPairB1AntiKt, jetPairB2AntiKt;
			     			TLorentzVector jetPairB1Check, jetPairB2Check;
			     			findJetPairs(jetPairB1, jetPairB2, jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham, jetPairB1Index1, jetPairB1Index2, jetPairB2Index1, jetPairB2Index2);
			     			invMassB1 = jetPairB1.M();
				     		invMassB2 = jetPairB2.M();
				     		
				     		
			     			vector<int> bestGroup1Durham, bestGroup2Durham, bestGroup1AntiKt, bestGroup2AntiKt;
			     			vector<double> chiSquares, chiSquaresAntiKt;
			     			vector<vector<int>> group1s, group2s, group1sAntiKt, group2sAntiKt;
    						findBestCombination(branchJetDurham, nJetsDurham, bestGroup1Durham, bestGroup2Durham, chiSquares, group1s, group2s);
    						findBestCombination(branchJetAntiKt, nJetsAntiKt, bestGroup1AntiKt, bestGroup2AntiKt, chiSquaresAntiKt, group1sAntiKt, group2sAntiKt);
    						sortVectorTriple(chiSquares, group1s, group2s);
    						double senChiSquared=0;
    						queue<int> combinationsIndices;
    						int chiSquaresIndex=0;
    						for(const auto& chiSquare : chiSquares) ///To find first n unique combinations that minimize chiSquared 
    						{ 
    							if(chiSquare != senChiSquared) 
    							{
    								combinationsIndices.push(chiSquaresIndex);
    								senChiSquared=chiSquare;
    							}
    							chiSquaresIndex++;
    						}
    						
    						int combinationsIndicesSize = combinationsIndices.size();
    						for(int i=0; i<combinationsIndicesSize; i++)
    						{
						    	for (const auto& index : group1s[combinationsIndices.front()])
						    	{ 
						    		Jet *jet = (Jet*) branchJetDurham->At(index);
						    		jetPairB1Sen += jet->P4();
						    	}
						    	for (const auto& index : group2s[combinationsIndices.front()])
						    	{ 
						    		Jet *jet = (Jet*) branchJetDurham->At(index);
						    		jetPairB2Sen += jet->P4();
						    	}
						    	if(i==0)
						    	{
							    	invMassB11Best = jetPairB1Sen.M();
				     				invMassB21Best = jetPairB2Sen.M();
			     				}
			     				if(i==1)
						    	{
							    	invMassB12Best = jetPairB1Sen.M();
				     				invMassB22Best = jetPairB2Sen.M();
			     				}
			     				if(i==2)
						    	{
							    	invMassB13Best = jetPairB1Sen.M();
				     				invMassB23Best = jetPairB2Sen.M();
			     				}
			     				if(i==3)
						    	{
							    	invMassB14Best = jetPairB1Sen.M();
				     				invMassB24Best = jetPairB2Sen.M();
			     				}
			     				if(i==4)
						    	{
							    	invMassB15Best = jetPairB1Sen.M();
				     				invMassB25Best = jetPairB2Sen.M();
			     				}
			     				if(i==5)
						    	{
							    	invMassB16Best = jetPairB1Sen.M();
				     				invMassB26Best = jetPairB2Sen.M();
			     				}
			     				if(i==6)
						    	{
							    	invMassB17Best = jetPairB1Sen.M();
				     				invMassB27Best = jetPairB2Sen.M();
			     				}
			     				if(i==7)
						    	{
							    	invMassB18Best = jetPairB1Sen.M();
				     				invMassB28Best = jetPairB2Sen.M();
			     				}
						    	combinationsIndices.pop();
					    	}
					    	
					    	for (const auto& index : bestGroup1AntiKt)
					    	{ 
					    		Jet *jet = (Jet*) branchJetAntiKt->At(index);
					    		jetPairB1AntiKt += jet->P4();
					    	}
					    	for (const auto& index : bestGroup2AntiKt)
					    	{ 
					    		Jet *jet = (Jet*) branchJetAntiKt->At(index);
					    		jetPairB2AntiKt += jet->P4();
					    	}
					    	/*if(contEventsPostFilter < 2)
    						{ 
    							cout<<endl<<"queue elements: ";
    							while(!combinationsIndices.empty())
    							{ 
    								cout<<combinationsIndices.front()<<", ";
    								combinationsIndices.pop();
    							}
    							cout<<endl<<endl<<"group1s.size(): "<<group1s.size()<<endl;
    							for(int i=0; i<group1s.size(); i++)
    							{
    								cout<<endl<<endl<<"For "<<i<<" combination: "<<endl<<"group1: ";
    								for(const auto& index : group1s[i]) cout<<index<<", ";
    								cout<<endl<<"group2: ";	
    								for(const auto& index : group2s[i]) cout<<index<<", ";
    								cout<<endl<<"chiSquared: "<<chiSquares[i];
    							}
    						}*/	
					    	/*if(bestGroup1.size() != bestGroup2.size()) 
					    	{
					    		cout<<endl<<"3-jet-pair vs 1: "<<entry<<endl;
					    		cout<<"bestGroup1.size: "<<bestGroup1.size()<<"    : ";
					    		for (const auto& index : bestGroup1) cout<<index;
					    		cout<<endl<<"bestGroup2.size: "<<bestGroup2.size()<<"    : ";
					    		for (const auto& index : bestGroup2) cout<<index;
					    		cout<<endl;
					    		cout<<"invMassB1: "<<jetPairB1.M()<<endl;
					    		cout<<"invMassB2: "<<jetPairB2.M()<<endl;
					    		cout<<"invMassB1Check: "<<jetPairB1Check.M()<<"    indices: "<<jetPairB1Index1<<", "<<jetPairB1Index2<<endl;
					    		cout<<"invMassB2Check: "<<jetPairB2Check.M()<<"    indices: "<<jetPairB2Index1<<", "<<jetPairB2Index2<<endl;
					    	}*/
			     			invMassB1AntiKt = jetPairB1AntiKt.M();
			     			invMassB2AntiKt = jetPairB2AntiKt.M();
			     			histJetB1M->Fill(invMassB1, weight);
						histJetB2M->Fill(invMassB2, weight);
						if(nJetsAntiKt == 2)
						{
							histJetB1MAntiKt2Jets->Fill(invMassB1AntiKt, weight);
							histJetB2MAntiKt2Jets->Fill(invMassB2AntiKt, weight);
							invMassB1AntiKt2Jets = invMassB1AntiKt;
							invMassB2AntiKt2Jets = invMassB2AntiKt;
						}
						if(nJetsAntiKt == 3)
						{
							histJetB1MAntiKt3Jets->Fill(invMassB1AntiKt, weight);
							histJetB2MAntiKt3Jets->Fill(invMassB2AntiKt, weight);
							invMassB1AntiKt3Jets = invMassB1AntiKt;
							invMassB2AntiKt3Jets = invMassB2AntiKt;
						}
						if(nJetsAntiKt == 4)
						{
							histJetB1MAntiKt4Jets->Fill(invMassB1AntiKt, weight);
							histJetB2MAntiKt4Jets->Fill(invMassB2AntiKt, weight);
							invMassB1AntiKt4Jets = invMassB1AntiKt;
							invMassB2AntiKt4Jets = invMassB2AntiKt;
						}
						if(nJetsAntiKt == 5)
						{
							histJetB1MAntiKt5Jets->Fill(invMassB1AntiKt, weight);
							histJetB2MAntiKt5Jets->Fill(invMassB2AntiKt, weight);
							invMassB1AntiKt5Jets = invMassB1AntiKt;
							invMassB2AntiKt5Jets = invMassB2AntiKt;
						}
						if(nJetsAntiKt == 6)
						{
							histJetB1MAntiKt6Jets->Fill(invMassB1AntiKt, weight);
							histJetB2MAntiKt6Jets->Fill(invMassB2AntiKt, weight);
							invMassB1AntiKt6Jets = invMassB1AntiKt;
							invMassB2AntiKt6Jets = invMassB2AntiKt;
						}
						jetB1M = jetB1Durham.M();
						jetB2M = jetB2Durham.M();
						jetB3M = jetB3Durham.M();
						jetB4M = jetB4Durham.M();
						histJet1M->Fill(jetB1M, weight);
						histJet2M->Fill(jetB2M, weight);
						histJet3M->Fill(jetB3M, weight);
						histJet4M->Fill(jetB4M, weight);
						minJetM = findMinJetM(jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham);
						histMinJetM->Fill(minJetM, weight);
						
						/*if(invMassB > 124 && invMassB < 126) contLargeBMass++;
						if(invMassB > 124 && invMassB < 126 && contLargeBMass < 4) cout<<"Average b mass event: "<<entry<<endl;
						if(invMassB < 30) contSmallBMass++;
						if(invMassB < 30 && contSmallBMass < 4) cout<<"Small b mass event: "<<entry<<endl;*/
			     			//////inv. mass
			     			
			     			//////angles
			     			etaB1 = jetB1Durham.Eta();
			     			etaB2 = jetB2Durham.Eta();
			     			etaB3 = jetB3Durham.Eta();
			     			etaB4 = jetB4Durham.Eta();
			     			
						histJetEta->Fill(etaB1, weight);
			     			histJetEta->Fill(etaB2, weight);
			     			histJetEta->Fill(etaB3, weight);
			     			histJetEta->Fill(etaB4, weight);
			     			
			     			histJet1Eta->Fill(etaB1, weight);
			     			histJet2Eta->Fill(etaB2, weight);
			     			histJet3Eta->Fill(etaB3, weight);
			     			histJet4Eta->Fill(etaB4, weight);
			     			
			     			cosThetaB1 = findCosTheta(etaB1);
			     			cosThetaB2 = findCosTheta(etaB2);
			     			cosThetaB3 = findCosTheta(etaB3);
			     			cosThetaB4 = findCosTheta(etaB4);
			     			
			     			histJetCosTheta->Fill(cosThetaB1, weight);
			     			histJetCosTheta->Fill(cosThetaB2, weight);
			     			histJetCosTheta->Fill(cosThetaB3, weight);
			     			histJetCosTheta->Fill(cosThetaB4, weight);
			     			
			     			histJet1CosTheta->Fill(cosThetaB1, weight);
			     			histJet2CosTheta->Fill(cosThetaB2, weight);
			     			histJet3CosTheta->Fill(cosThetaB3, weight);
			     			histJet4CosTheta->Fill(cosThetaB4, weight);
			     			//////angles
						
						//////Pt
						sumPt = jetB1Durham.Pt() + jetB2Durham.Pt() + jetB3Durham.Pt() + jetB4Durham.Pt();
						histSumJetPt->Fill(sumPt, weight);
						
			     			jetB1Pt = jetB1Durham.Pt();
			     			jetB2Pt = jetB2Durham.Pt();
			     			jetB3Pt = jetB4Durham.Pt();
			     			jetB4Pt = jetB4Durham.Pt();
			     			histJetPt->Fill(jetB1Pt, weight);
			     			histJetPt->Fill(jetB2Pt, weight);
			     			histJetPt->Fill(jetB3Pt, weight);
			     			histJetPt->Fill(jetB4Pt, weight);
			     			
			     			histJet1Pt->Fill(jetB1Pt, weight);
			     			histJet2Pt->Fill(jetB2Pt, weight);
			     			histJet3Pt->Fill(jetB3Pt, weight);
			     			histJet4Pt->Fill(jetB4Pt, weight);
						//////Pt
						
						//////Event shape
						aplanarity = findEventShape(jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham, "aplanarity");
						sphericity = findEventShape(jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham, "sphericity");
						histAplan->Fill(aplanarity, weight);
						histSpher->Fill(sphericity, weight);
						//////Event shape
						
						//////Jet constituents
						totalConstSize = constSizeB1 + constSizeB2 + constSizeB3 + constSizeB4;
						histNParticles->Fill(nParticles, weight);
						histTotalConstSize->Fill(totalConstSize, weight);
						histConstSizeB1->Fill(constSizeB1, weight);
						histConstSizeB2->Fill(constSizeB2, weight);
						histConstSizeB3->Fill(constSizeB3, weight);
						histConstSizeB4->Fill(constSizeB4, weight);
						minConstSize = TMath::Min(constSizeB1, constSizeB2);
						minConstSize = TMath::Min(minConstSize, constSizeB3);
						minConstSize = TMath::Min(minConstSize, constSizeB4);
						histMinConstSize->Fill(minConstSize, weight);
						//////Jet constituents
						
						/////PFOs
						int nEFlowTracks=0, nEFlowPhotons=0, nEFlowNeutralHadrons=0, nEFlowObjects=0;
						nEFlowTracks =  branchEFlowTrack->GetEntries();
						nEFlowPhotons =  branchEFlowPhoton->GetEntries();
						nEFlowNeutralHadrons =  branchEFlowNeutralHadron->GetEntries();
						nEFlowObjects = nEFlowTracks + nEFlowPhotons + nEFlowNeutralHadrons;
						histNEFlowTracks->Fill(nEFlowTracks, weight);
						histNEFlowPhotons->Fill(nEFlowPhotons, weight); 
						histNEFlowNeutralHadrons->Fill(nEFlowNeutralHadrons, weight);
						histNEFlowObjects->Fill(nEFlowObjects, weight);
						
						float jetB1NObjects = jetB1NCharged+jetB1NNeutrals;
						float jetB2NObjects = jetB2NCharged+jetB2NNeutrals;
						float jetB3NObjects = jetB3NCharged+jetB3NNeutrals;
						float jetB4NObjects = jetB4NCharged+jetB4NNeutrals;
						jetNObjects = jetB1NObjects + jetB2NObjects + jetB3NObjects +jetB4NObjects;
						
						minJetNObjects = TMath::Min(jetB1NObjects, jetB2NObjects);
						minJetNObjects = TMath::Min(minJetNObjects, jetB3NObjects);
						minJetNObjects = TMath::Min(minJetNObjects, jetB4NObjects);
						
						histJetB1NCharged->Fill(jetB1NCharged, weight);
						histJetB2NCharged->Fill(jetB2NCharged, weight);
						histJetB3NCharged->Fill(jetB3NCharged, weight);
						histJetB4NCharged->Fill(jetB4NCharged, weight);
						
						histJetB1NNeutrals->Fill(jetB1NNeutrals, weight);
						histJetB2NNeutrals->Fill(jetB2NNeutrals, weight);
						histJetB3NNeutrals->Fill(jetB3NNeutrals, weight);
						histJetB4NNeutrals->Fill(jetB4NNeutrals, weight);
						
						histJetNObjects->Fill(jetNObjects, weight);
						histMinJetNObjects->Fill(minJetNObjects, weight);
						/////PFOs
						
						/////Filling out tree for BDT
						random_device rd;  // Will be used to obtain a seed for the random number engine
						mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
						uniform_int_distribution<> distrib(1, 4); // Uniform distribution between 1 and 4
						if(distrib(gen) != 1)
						{ 
							TreeTrain.Fill();
							//if(contEventsPostFilter<30) cout<<"contEventsPostFilter: "<<contEventsPostFilter<<"    NO"<<endl;
						}
						else
						{ 
							TreeTest.Fill();
							//if(contEventsPostFilter<30) cout<<"contEventsPostFilter: "<<contEventsPostFilter<<"    EN CUATRO"<<endl;
						}
						/////Filling out tree for BDT
						
						
						
							     			
		     			}
		     		}
	     			
	     			
	     		}
     		}
	}
	
	//////printing histograms
	if(topology == 1)
	{ 
		/*//////eta
		TCanvas *c11 = new TCanvas();
		TCanvas *c12 = new TCanvas();
		TCanvas *c13 = new TCanvas();
		TCanvas *c14 = new TCanvas();
		TCanvas *c15 = new TCanvas();			
		c11->cd();
    		histJetEta->Draw("HIST");
    		c12->cd();
    		histJet1Eta->Draw("HIST");				
    		c13->cd();
    		histJet2Eta->Draw("HIST");				
    		c14->cd();
    		histJet3Eta->Draw("HIST");				
    		c15->cd();
    		histJet4Eta->Draw("HIST");
    		//////eta*/
    		
    		/*//////cosTheta
		TCanvas *c16 = new TCanvas();
		TCanvas *c17 = new TCanvas();
		TCanvas *c18 = new TCanvas();
		TCanvas *c19 = new TCanvas();
		TCanvas *c110 = new TCanvas();			
		c16->cd();
    		histJetCosTheta->Draw("HIST");
    		c17->cd();
    		histJet1CosTheta->Draw("HIST");				
    		c18->cd();
    		histJet2CosTheta->Draw("HIST");				
    		c19->cd();
    		histJet3CosTheta->Draw("HIST");				
    		c110->cd();
    		histJet4CosTheta->Draw("HIST");
    		//////cosTheta*/
    		
    		/*/////Pt
    		c1001->cd();
    		histJet4Pt->Draw("HIST");
    		histJet1Pt->Draw("HIST same");
    		histJet2Pt->Draw("HIST same");
    		histJet3Pt->Draw("HIST same");
    		TCanvas *c111 = new TCanvas();
		TCanvas *c112 = new TCanvas();	
		c111->cd();
    		histSumJetPt->Draw("HIST");
    		c112->cd();
    		histJetPt->Draw("HIST");	
    		/////Pt*/
    		
    		///////inv. masses
    		/*TCanvas *c113 = new TCanvas();
		TCanvas *c114 = new TCanvas();
		TCanvas *c115 = new TCanvas();
		TCanvas *c116 = new TCanvas();
		TCanvas *c117 = new TCanvas();
		TCanvas *c118 = new TCanvas();
		TCanvas *c119 = new TCanvas();
		c113->cd();
		histJetB1M->Draw("HIST");
		c114->cd();
		histJetB2M->Draw("HIST");
		c115->cd();
		histMinJetM->Draw("HIST");
		c116->cd();
		histJet1M->Draw("HIST");
		c117->cd();
		histJet2M->Draw("HIST");
		c118->cd();
		histJet3M->Draw("HIST");
		c119->cd();
		histJet4M->Draw("HIST");
    		//////inv. masses*/
    		
    		/*/////Event shape
    		TCanvas *c120 = new TCanvas();
		TCanvas *c121 = new TCanvas();
		c120->cd();
		histAplan->Draw("HIST");
		c121->cd();
		histSpher->Draw("HIST");
    		/////Event shape*/
    		
    		/*/////Jet constituents
    		TCanvas *c122 = new TCanvas();
    		TCanvas *c123 = new TCanvas();
    		TCanvas *c124 = new TCanvas();
    		TCanvas *c125 = new TCanvas();
    		TCanvas *c126 = new TCanvas();
    		TCanvas *c127 = new TCanvas();
    		TCanvas *c128 = new TCanvas();
    		c122->cd();
    		histNParticles->Draw("HIST");
    		c123->cd();
    		histTotalConstSize->Draw("HIST");
    		c124->cd();
    		histConstSizeB1->Draw("HIST");
    		c125->cd();
    		histConstSizeB2->Draw("HIST");
    		c126->cd();
    		histConstSizeB3->Draw("HIST");
    		c127->cd();
    		histConstSizeB4->Draw("HIST");
    		c128->cd();
    		histMinConstSize->Draw("HIST");
    		/////*/
    		
    		/*/////PFOs
    		TCanvas *c129 = new TCanvas();
    		c129->cd();
    		histNEFlowTracks->Draw("HIST");
    		TCanvas *c130 = new TCanvas();
    		c130->cd();
    		histNEFlowPhotons->Draw("HIST");
    		TCanvas *c131 = new TCanvas();
    		c131->cd();
    		histNEFlowNeutralHadrons->Draw("HIST");
    		TCanvas *c132 = new TCanvas();
    		c132->cd();
    		histNEFlowObjects->Draw("HIST");
    		TCanvas *c133 = new TCanvas();
    		c133->cd();
    		histJetB1NCharged->Draw("HIST");
    		TCanvas *c134 = new TCanvas();
    		c134->cd();
    		histJetB2NCharged->Draw("HIST");
    		TCanvas *c135 = new TCanvas();
    		c135->cd();
    		histJetB3NCharged->Draw("HIST");
    		TCanvas *c136 = new TCanvas();
    		c136->cd();
    		histJetB4NCharged->Draw("HIST");
    		TCanvas *c137 = new TCanvas();
    		c137->cd();
    		histJetB1NNeutrals->Draw("HIST");
    		TCanvas *c138 = new TCanvas();
    		c138->cd();
    		histJetB2NNeutrals->Draw("HIST");
    		TCanvas *c139 = new TCanvas();
    		c139->cd();
    		histJetB3NNeutrals->Draw("HIST");
    		TCanvas *c140 = new TCanvas();
    		c140->cd();
    		histJetB4NNeutrals->Draw("HIST");
    		TCanvas *c141 = new TCanvas();
    		c141->cd();
    		histJetNObjects->Draw("HIST");
    		TCanvas *c142 = new TCanvas();
    		c142->cd();
    		histMinJetNObjects->Draw("HIST");
    		/////PFOs*/
    		
    		/*///////nJetsCompareAlgos
    		TCanvas *c143 = new TCanvas();
    		c143->cd();
    		histNJetsCompareAlgos->Draw("TEXT");
    		///////nJetsCompareAlgos*/
    		
    		/*///////invMass for antiKt pairs when 4 Durham jets
    		TCanvas *c144 = new TCanvas();
    		c144->cd();
    		histJetB1MAntiKt2Jets->Draw("HIST");
    		TCanvas *c145 = new TCanvas();
    		c145->cd();
    		histJetB2MAntiKt2Jets->Draw("HIST");
    		TCanvas *c146 = new TCanvas();
    		c146->cd();
    		histJetB1MAntiKt3Jets->Draw("HIST");
    		TCanvas *c147 = new TCanvas();
    		c147->cd();
    		histJetB2MAntiKt3Jets->Draw("HIST");
    		TCanvas *c148 = new TCanvas();
    		c148->cd();
    		histJetB1MAntiKt4Jets->Draw("HIST");
    		TCanvas *c149 = new TCanvas();
    		c149->cd();
    		histJetB2MAntiKt4Jets->Draw("HIST");
    		TCanvas *c150 = new TCanvas();
    		c150->cd();
    		histJetB1MAntiKt5Jets->Draw("HIST");
    		TCanvas *c151 = new TCanvas();
    		c151->cd();
    		histJetB2MAntiKt5Jets->Draw("HIST");
    		TCanvas *c152 = new TCanvas();
    		c152->cd();
    		histJetB1MAntiKt6Jets->Draw("HIST");
    		TCanvas *c153 = new TCanvas();
    		c153->cd();
    		histJetB2MAntiKt6Jets->Draw("HIST");
    		///////invMass for antiKt pairs when 4 Durham jets*/
    		
    		
	} 
	else if(topology == 2)
	{ 
		/*//////eta
		TCanvas *c21 = new TCanvas();
		TCanvas *c22 = new TCanvas();
		TCanvas *c23 = new TCanvas();
		TCanvas *c24 = new TCanvas();
		TCanvas *c25 = new TCanvas();			
		c21->cd();
    		histJetEta->Draw("HIST");
    		c22->cd();
    		histJet1Eta->Draw("HIST");				
    		c23->cd();
    		histJet2Eta->Draw("HIST");				
    		c24->cd();
    		histJet3Eta->Draw("HIST");				
    		c25->cd();
    		histJet4Eta->Draw("HIST");
    		//////eta*/
    		
    		/*//////cosTheta
		TCanvas *c26 = new TCanvas();
		TCanvas *c27 = new TCanvas();
		TCanvas *c28 = new TCanvas();
		TCanvas *c29 = new TCanvas();
		TCanvas *c210 = new TCanvas();			
		c26->cd();
    		histJetCosTheta->Draw("HIST");
    		c27->cd();
    		histJet1CosTheta->Draw("HIST");				
    		c28->cd();
    		histJet2CosTheta->Draw("HIST");				
    		c29->cd();
    		histJet3CosTheta->Draw("HIST");				
    		c210->cd();
    		histJet4CosTheta->Draw("HIST");
    		//////cosTheta*/
    		
    		/*/////Pt
    		c1001->cd();
    		histJet4Pt->Draw("HIST");
    		histJet1Pt->Draw("HIST same");
    		histJet2Pt->Draw("HIST same");
    		histJet3Pt->Draw("HIST same");
    		TCanvas *c211 = new TCanvas();
		TCanvas *c212 = new TCanvas();	
		c211->cd();
    		histSumJetPt->Draw("HIST");
    		c212->cd();
    		histJetPt->Draw("HIST");	
    		/////Pt
    		
    		///////inv. masses
    		TCanvas *c213 = new TCanvas();
		TCanvas *c214 = new TCanvas();
		TCanvas *c215 = new TCanvas();
		TCanvas *c216 = new TCanvas();
		TCanvas *c217 = new TCanvas();
		TCanvas *c218 = new TCanvas();
		TCanvas *c219 = new TCanvas();
		c213->cd();
		histJetB1M->Draw("HIST");
		c214->cd();
		histJetB2M->Draw("HIST");
		c215->cd();
		histMinJetM->Draw("HIST");
		c216->cd();
		histJet1M->Draw("HIST");
		c217->cd();
		histJet2M->Draw("HIST");
		c218->cd();
		histJet3M->Draw("HIST");
		c219->cd();
		histJet4M->Draw("HIST");
    		//////inv. masses
    		
    		/////Event shape
    		TCanvas *c220 = new TCanvas();
		TCanvas *c221 = new TCanvas();
		c220->cd();
		histAplan->Draw("HIST");
		c221->cd();
		histSpher->Draw("HIST");
    		/////Event shape*/
    		
    		/*/////Jet constituents
    		TCanvas *c222 = new TCanvas();
    		TCanvas *c223 = new TCanvas();
    		TCanvas *c224 = new TCanvas();
    		TCanvas *c225 = new TCanvas();
    		TCanvas *c226 = new TCanvas();
    		TCanvas *c227 = new TCanvas();
    		TCanvas *c228 = new TCanvas();
    		c222->cd();
    		histNParticles->Draw("HIST");
    		c223->cd();
    		histTotalConstSize->Draw("HIST");
    		c224->cd();
    		histConstSizeB1->Draw("HIST");
    		c225->cd();
    		histConstSizeB2->Draw("HIST");
    		c226->cd();
    		histConstSizeB3->Draw("HIST");
    		c227->cd();
    		histConstSizeB4->Draw("HIST");
    		c228->cd();
    		histMinConstSize->Draw("HIST");
    		/////*/
    		
    		/*/////PFOs
    		TCanvas *c229 = new TCanvas();
    		c229->cd();
    		histNEFlowTracks->Draw("HIST");
    		TCanvas *c230 = new TCanvas();
    		c230->cd();
    		histNEFlowPhotons->Draw("HIST");
    		TCanvas *c231 = new TCanvas();
    		c231->cd();
    		histNEFlowNeutralHadrons->Draw("HIST");
    		TCanvas *c232 = new TCanvas();
    		c232->cd();
    		histNEFlowObjects->Draw("HIST");
    		TCanvas *c233 = new TCanvas();
    		c233->cd();
    		histJetB1NCharged->Draw("HIST");
    		TCanvas *c234 = new TCanvas();
    		c234->cd();
    		histJetB2NCharged->Draw("HIST");
    		TCanvas *c235 = new TCanvas();
    		c235->cd();
    		histJetB3NCharged->Draw("HIST");
    		TCanvas *c236 = new TCanvas();
    		c236->cd();
    		histJetB4NCharged->Draw("HIST");
    		TCanvas *c237 = new TCanvas();
    		c237->cd();
    		histJetB1NNeutrals->Draw("HIST");
    		TCanvas *c238 = new TCanvas();
    		c238->cd();
    		histJetB2NNeutrals->Draw("HIST");
    		TCanvas *c239 = new TCanvas();
    		c239->cd();
    		histJetB3NNeutrals->Draw("HIST");
    		TCanvas *c240 = new TCanvas();
    		c240->cd();
    		histJetB4NNeutrals->Draw("HIST");
    		TCanvas *c241 = new TCanvas();
    		c241->cd();
    		histJetNObjects->Draw("HIST");
    		TCanvas *c242 = new TCanvas();
    		c242->cd();
    		histMinJetNObjects->Draw("HIST");
    		/////PFOs*/
    		
    		/*///////nJetsCompareAlgos
    		TCanvas *c243 = new TCanvas();
    		c243->cd();
    		histNJetsCompareAlgos->Draw("TEXT");
    		///////nJetsCompareAlgos*/
    		
    		/*///////invMass for antiKt pairs when 4 Durham jets
    		TCanvas *c244 = new TCanvas();
    		c244->cd();
    		histJetB1MAntiKt2Jets->Draw("HIST");
    		TCanvas *c245 = new TCanvas();
    		c245->cd();
    		histJetB2MAntiKt2Jets->Draw("HIST");
    		TCanvas *c246 = new TCanvas();
    		c246->cd();
    		histJetB1MAntiKt3Jets->Draw("HIST");
    		TCanvas *c247 = new TCanvas();
    		c247->cd();
    		histJetB2MAntiKt3Jets->Draw("HIST");
    		TCanvas *c248 = new TCanvas();
    		c248->cd();
    		histJetB1MAntiKt4Jets->Draw("HIST");
    		TCanvas *c249 = new TCanvas();
    		c249->cd();
    		histJetB2MAntiKt4Jets->Draw("HIST");
    		TCanvas *c250 = new TCanvas();
    		c250->cd();
    		histJetB1MAntiKt5Jets->Draw("HIST");
    		TCanvas *c251 = new TCanvas();
    		c251->cd();
    		histJetB2MAntiKt5Jets->Draw("HIST");
    		TCanvas *c252 = new TCanvas();
    		c252->cd();
    		histJetB1MAntiKt6Jets->Draw("HIST");
    		TCanvas *c253 = new TCanvas();
    		c253->cd();
    		histJetB2MAntiKt6Jets->Draw("HIST");
    		///////invMass for antiKt pairs when 4 Durham jets*/
	} 
	else if(topology == 3)
	{ 
		/*//////eta
		TCanvas *c31 = new TCanvas();
		TCanvas *c32 = new TCanvas();
		TCanvas *c33 = new TCanvas();
		TCanvas *c34 = new TCanvas();
		TCanvas *c35 = new TCanvas();			
		c31->cd();
    		histJetEta->Draw("HIST");
    		c32->cd();
    		histJet1Eta->Draw("HIST");				
    		c33->cd();
    		histJet2Eta->Draw("HIST");				
    		c34->cd();
    		histJet3Eta->Draw("HIST");				
    		c35->cd();
    		histJet4Eta->Draw("HIST");
    		//////eta*/
    		
    		/*//////cosTheta
		TCanvas *c36 = new TCanvas();
		TCanvas *c37 = new TCanvas();
		TCanvas *c38 = new TCanvas();
		TCanvas *c39 = new TCanvas();
		TCanvas *c310 = new TCanvas();			
		c36->cd();
    		histJetCosTheta->Draw("HIST");
    		c37->cd();
    		histJet1CosTheta->Draw("HIST");				
    		c38->cd();
    		histJet2CosTheta->Draw("HIST");				
    		c39->cd();
    		histJet3CosTheta->Draw("HIST");				
    		c310->cd();
    		histJet4CosTheta->Draw("HIST");
    		//////cosTheta*/
    		
    		/*/////Pt
    		c1001->cd();
    		histJet4Pt->Draw("HIST");
    		histJet1Pt->Draw("HIST same");
    		histJet2Pt->Draw("HIST same");
    		histJet3Pt->Draw("HIST same");
    		TCanvas *c311 = new TCanvas();
		TCanvas *c312 = new TCanvas();	
		c311->cd();
    		histSumJetPt->Draw("HIST");
    		c312->cd();
    		histJetPt->Draw("HIST");	
    		/////Pt
    		
    		///////inv. masses
    		TCanvas *c313 = new TCanvas();
		TCanvas *c314 = new TCanvas();
		TCanvas *c315 = new TCanvas();
		TCanvas *c316 = new TCanvas();
		TCanvas *c317 = new TCanvas();
		TCanvas *c318 = new TCanvas();
		TCanvas *c319 = new TCanvas();
		c313->cd();
		histJetB1M->Draw("HIST");
		c314->cd();
		histJetB2M->Draw("HIST");
		c315->cd();
		histMinJetM->Draw("HIST");
		c316->cd();
		histJet1M->Draw("HIST");
		c317->cd();
		histJet2M->Draw("HIST");
		c318->cd();
		histJet3M->Draw("HIST");
		c319->cd();
		histJet4M->Draw("HIST");
    		//////inv. masses
    		
    		/////Event shape
    		TCanvas *c320 = new TCanvas();
		TCanvas *c321 = new TCanvas();
		c320->cd();
		histAplan->Draw("HIST");
		c321->cd();
		histSpher->Draw("HIST");
    		/////Event shape*/
    		
    		/*/////Jet constituents
    		TCanvas *c322 = new TCanvas();
    		TCanvas *c323 = new TCanvas();
    		TCanvas *c324 = new TCanvas();
    		TCanvas *c325 = new TCanvas();
    		TCanvas *c326 = new TCanvas();
    		TCanvas *c327 = new TCanvas();
    		TCanvas *c328 = new TCanvas();
    		c322->cd();
    		histNParticles->Draw("HIST");
    		c323->cd();
    		histTotalConstSize->Draw("HIST");
    		c324->cd();
    		histConstSizeB1->Draw("HIST");
    		c325->cd();
    		histConstSizeB2->Draw("HIST");
    		c326->cd();
    		histConstSizeB3->Draw("HIST");
    		c327->cd();
    		histConstSizeB4->Draw("HIST");
    		c328->cd();
    		histMinConstSize->Draw("HIST");
    		/////*/
    		
    		/*/////PFOs
    		TCanvas *c329 = new TCanvas();
    		c329->cd();
    		histNEFlowTracks->Draw("HIST");
    		TCanvas *c330 = new TCanvas();
    		c330->cd();
    		histNEFlowPhotons->Draw("HIST");
    		TCanvas *c331 = new TCanvas();
    		c331->cd();
    		histNEFlowNeutralHadrons->Draw("HIST");
    		TCanvas *c332 = new TCanvas();
    		c332->cd();
    		histNEFlowObjects->Draw("HIST");
    		TCanvas *c333 = new TCanvas();
    		c333->cd();
    		histJetB1NCharged->Draw("HIST");
    		TCanvas *c334 = new TCanvas();
    		c334->cd();
    		histJetB2NCharged->Draw("HIST");
    		TCanvas *c335 = new TCanvas();
    		c335->cd();
    		histJetB3NCharged->Draw("HIST");
    		TCanvas *c336 = new TCanvas();
    		c336->cd();
    		histJetB4NCharged->Draw("HIST");
    		TCanvas *c337 = new TCanvas();
    		c337->cd();
    		histJetB1NNeutrals->Draw("HIST");
    		TCanvas *c338 = new TCanvas();
    		c338->cd();
    		histJetB2NNeutrals->Draw("HIST");
    		TCanvas *c339 = new TCanvas();
    		c339->cd();
    		histJetB3NNeutrals->Draw("HIST");
    		TCanvas *c340 = new TCanvas();
    		c340->cd();
    		histJetB4NNeutrals->Draw("HIST");
    		TCanvas *c341 = new TCanvas();
    		c341->cd();
    		histJetNObjects->Draw("HIST");
    		TCanvas *c342 = new TCanvas();
    		c342->cd();
    		histMinJetNObjects->Draw("HIST");
    		/////PFOs*/
    		
    		/*///////nJetsCompareAlgos
    		TCanvas *c343 = new TCanvas();
    		c343->cd();
    		histNJetsCompareAlgos->Draw("TEXT");
    		///////nJetsCompareAlgos*/
    		
    		/*///////invMass for antiKt pairs when 4 Durham jets
    		TCanvas *c344 = new TCanvas();
    		c344->cd();
    		histJetB1MAntiKt2Jets->Draw("HIST");
    		TCanvas *c345 = new TCanvas();
    		c345->cd();
    		histJetB2MAntiKt2Jets->Draw("HIST");
    		TCanvas *c346 = new TCanvas();
    		c346->cd();
    		histJetB1MAntiKt3Jets->Draw("HIST");
    		TCanvas *c347 = new TCanvas();
    		c347->cd();
    		histJetB2MAntiKt3Jets->Draw("HIST");
    		TCanvas *c348 = new TCanvas();
    		c348->cd();
    		histJetB1MAntiKt4Jets->Draw("HIST");
    		TCanvas *c349 = new TCanvas();
    		c349->cd();
    		histJetB2MAntiKt4Jets->Draw("HIST");
    		TCanvas *c350 = new TCanvas();
    		c350->cd();
    		histJetB1MAntiKt5Jets->Draw("HIST");
    		TCanvas *c351 = new TCanvas();
    		c351->cd();
    		histJetB2MAntiKt5Jets->Draw("HIST");
    		TCanvas *c352 = new TCanvas();
    		c352->cd();
    		histJetB1MAntiKt6Jets->Draw("HIST");
    		TCanvas *c353 = new TCanvas();
    		c353->cd();
    		histJetB2MAntiKt6Jets->Draw("HIST");
    		///////invMass for antiKt pairs when 4 Durham jets*/
	}
	else if(topology == 4)
	{ 
		/*//////eta
		TCanvas *c41 = new TCanvas();
		TCanvas *c42 = new TCanvas();
		TCanvas *c43 = new TCanvas();
		TCanvas *c44 = new TCanvas();
		TCanvas *c45 = new TCanvas();			
		c41->cd();
    		histJetEta->Draw("HIST");
    		c42->cd();
    		histJet1Eta->Draw("HIST");				
    		c43->cd();
    		histJet2Eta->Draw("HIST");				
    		c44->cd();
    		histJet3Eta->Draw("HIST");				
    		c45->cd();
    		histJet4Eta->Draw("HIST");
    		//////eta*/
    		
    		/*//////cosTheta
		TCanvas *c46 = new TCanvas();
		TCanvas *c47 = new TCanvas();
		TCanvas *c48 = new TCanvas();
		TCanvas *c49 = new TCanvas();
		TCanvas *c410 = new TCanvas();			
		c46->cd();
    		histJetCosTheta->Draw("HIST");
    		c47->cd();
    		histJet1CosTheta->Draw("HIST");				
    		c48->cd();
    		histJet2CosTheta->Draw("HIST");				
    		c49->cd();
    		histJet3CosTheta->Draw("HIST");				
    		c410->cd();
    		histJet4CosTheta->Draw("HIST");
    		//////cosTheta*/
    		
    		/*/////Pt
    		c1001->cd();
    		histJet4Pt->Draw("HIST");
    		histJet1Pt->Draw("HIST same");
    		histJet2Pt->Draw("HIST same");
    		histJet3Pt->Draw("HIST same");
    		TCanvas *c411 = new TCanvas();
		TCanvas *c412 = new TCanvas();	
		c411->cd();
    		histSumJetPt->Draw("HIST");
    		c412->cd();
    		histJetPt->Draw("HIST");	
    		/////Pt
    		
    		///////inv. masses
    		TCanvas *c413 = new TCanvas();
		TCanvas *c414 = new TCanvas();
		TCanvas *c415 = new TCanvas();
		TCanvas *c416 = new TCanvas();
		TCanvas *c417 = new TCanvas();
		TCanvas *c418 = new TCanvas();
		TCanvas *c419 = new TCanvas();
		c413->cd();
		histJetB1M->Draw("HIST");
		c414->cd();
		histJetB2M->Draw("HIST");
		c415->cd();
		histMinJetM->Draw("HIST");
		c416->cd();
		histJet1M->Draw("HIST");
		c417->cd();
		histJet2M->Draw("HIST");
		c418->cd();
		histJet3M->Draw("HIST");
		c419->cd();
		histJet4M->Draw("HIST");
    		//////inv. masses
    		
    		/////Event shape
    		TCanvas *c420 = new TCanvas();
		TCanvas *c421 = new TCanvas();
		c420->cd();
		histAplan->Draw("HIST");
		c421->cd();
		histSpher->Draw("HIST");
    		/////Event shape*/
    		
    		/*/////Jet constituents
    		TCanvas *c422 = new TCanvas();
    		TCanvas *c423 = new TCanvas();
    		TCanvas *c424 = new TCanvas();
    		TCanvas *c425 = new TCanvas();
    		TCanvas *c426 = new TCanvas();
    		TCanvas *c427 = new TCanvas();
    		TCanvas *c428 = new TCanvas();
    		c422->cd();
    		histNParticles->Draw("HIST");
    		c423->cd();
    		histTotalConstSize->Draw("HIST");
    		c424->cd();
    		histConstSizeB1->Draw("HIST");
    		c425->cd();
    		histConstSizeB2->Draw("HIST");
    		c426->cd();
    		histConstSizeB3->Draw("HIST");
    		c427->cd();
    		histConstSizeB4->Draw("HIST");
    		c428->cd();
    		histMinConstSize->Draw("HIST");
    		/////*/
    		
    		/*/////PFOs
    		TCanvas *c429 = new TCanvas();
    		c429->cd();
    		histNEFlowTracks->Draw("HIST");
    		TCanvas *c430 = new TCanvas();
    		c430->cd();
    		histNEFlowPhotons->Draw("HIST");
    		TCanvas *c431 = new TCanvas();
    		c431->cd();
    		histNEFlowNeutralHadrons->Draw("HIST");
    		TCanvas *c432 = new TCanvas();
    		c432->cd();
    		histNEFlowObjects->Draw("HIST");
    		TCanvas *c433 = new TCanvas();
    		c433->cd();
    		histJetB1NCharged->Draw("HIST");
    		TCanvas *c434 = new TCanvas();
    		c434->cd();
    		histJetB2NCharged->Draw("HIST");
    		TCanvas *c435 = new TCanvas();
    		c435->cd();
    		histJetB3NCharged->Draw("HIST");
    		TCanvas *c436 = new TCanvas();
    		c436->cd();
    		histJetB4NCharged->Draw("HIST");
    		TCanvas *c437 = new TCanvas();
    		c437->cd();
    		histJetB1NNeutrals->Draw("HIST");
    		TCanvas *c438 = new TCanvas();
    		c438->cd();
    		histJetB2NNeutrals->Draw("HIST");
    		TCanvas *c439 = new TCanvas();
    		c439->cd();
    		histJetB3NNeutrals->Draw("HIST");
    		TCanvas *c440 = new TCanvas();
    		c440->cd();
    		histJetB4NNeutrals->Draw("HIST");
    		TCanvas *c441 = new TCanvas();
    		c441->cd();
    		histJetNObjects->Draw("HIST");
    		TCanvas *c442 = new TCanvas();
    		c442->cd();
    		histMinJetNObjects->Draw("HIST");
    		/////PFOs*/
    		
    		/*///////nJetsCompareAlgos
    		TCanvas *c443 = new TCanvas();
    		c443->cd();
    		histNJetsCompareAlgos->Draw("TEXT");
    		///////nJetsCompareAlgos*/
    		
    		/*///////invMass for antiKt pairs when 4 Durham jets
    		TCanvas *c444 = new TCanvas();
    		c444->cd();
    		histJetB1MAntiKt2Jets->Draw("HIST");
    		TCanvas *c445 = new TCanvas();
    		c445->cd();
    		histJetB2MAntiKt2Jets->Draw("HIST");
    		TCanvas *c446 = new TCanvas();
    		c446->cd();
    		histJetB1MAntiKt3Jets->Draw("HIST");
    		TCanvas *c447 = new TCanvas();
    		c447->cd();
    		histJetB2MAntiKt3Jets->Draw("HIST");
    		TCanvas *c448 = new TCanvas();
    		c448->cd();
    		histJetB1MAntiKt4Jets->Draw("HIST");
    		TCanvas *c449 = new TCanvas();
    		c449->cd();
    		histJetB2MAntiKt4Jets->Draw("HIST");
    		TCanvas *c450 = new TCanvas();
    		c450->cd();
    		histJetB1MAntiKt5Jets->Draw("HIST");
    		TCanvas *c451 = new TCanvas();
    		c451->cd();
    		histJetB2MAntiKt5Jets->Draw("HIST");
    		TCanvas *c452 = new TCanvas();
    		c452->cd();
    		histJetB1MAntiKt6Jets->Draw("HIST");
    		TCanvas *c453 = new TCanvas();
    		c453->cd();
    		histJetB2MAntiKt6Jets->Draw("HIST");
    		///////invMass for antiKt pairs when 4 Durham jets*/
	}
	else if(topology == 5)
	{ 
		/*//////eta
		TCanvas *c51 = new TCanvas();
		TCanvas *c52 = new TCanvas();
		TCanvas *c53 = new TCanvas();
		TCanvas *c54 = new TCanvas();
		TCanvas *c55 = new TCanvas();			
		c51->cd();
    		histJetEta->Draw("HIST");
    		c52->cd();
    		histJet1Eta->Draw("HIST");				
    		c53->cd();
    		histJet2Eta->Draw("HIST");				
    		c54->cd();
    		histJet3Eta->Draw("HIST");				
    		c55->cd();
    		histJet4Eta->Draw("HIST");
    		//////eta*/
    		
    		/*//////cosTheta
		TCanvas *c56 = new TCanvas();
		TCanvas *c57 = new TCanvas();
		TCanvas *c58 = new TCanvas();
		TCanvas *c59 = new TCanvas();
		TCanvas *c510 = new TCanvas();			
		c56->cd();
    		histJetCosTheta->Draw("HIST");
    		c57->cd();
    		histJet1CosTheta->Draw("HIST");				
    		c58->cd();
    		histJet2CosTheta->Draw("HIST");				
    		c59->cd();
    		histJet3CosTheta->Draw("HIST");				
    		c510->cd();
    		histJet4CosTheta->Draw("HIST");
    		//////cosTheta*/
    		
    		/*/////Pt
    		c1001->cd();
    		histJet4Pt->Draw("HIST");
    		histJet1Pt->Draw("HIST same");
    		histJet2Pt->Draw("HIST same");
    		histJet3Pt->Draw("HIST same");
    		TCanvas *c511 = new TCanvas();
		TCanvas *c512 = new TCanvas();	
		c511->cd();
    		histSumJetPt->Draw("HIST");
    		c512->cd();
    		histJetPt->Draw("HIST");	
    		/////Pt
    		
    		///////inv. masses
    		TCanvas *c513 = new TCanvas();
		TCanvas *c514 = new TCanvas();
		TCanvas *c515 = new TCanvas();
		TCanvas *c516 = new TCanvas();
		TCanvas *c517 = new TCanvas();
		TCanvas *c518 = new TCanvas();
		TCanvas *c519 = new TCanvas();
		c513->cd();
		histJetB1M->Draw("HIST");
		c514->cd();
		histJetB2M->Draw("HIST");
		c515->cd();
		histMinJetM->Draw("HIST");
		c516->cd();
		histJet1M->Draw("HIST");
		c517->cd();
		histJet2M->Draw("HIST");
		c518->cd();
		histJet3M->Draw("HIST");
		c519->cd();
		histJet4M->Draw("HIST");
    		//////inv. masses
    		
    		/////Event shape
    		TCanvas *c520 = new TCanvas();
		TCanvas *c521 = new TCanvas();
		c520->cd();
		histAplan->Draw("HIST");
		c521->cd();
		histSpher->Draw("HIST");
    		/////Event shape*/
    		
    		/*/////Jet constituents
    		TCanvas *c522 = new TCanvas();
    		TCanvas *c523 = new TCanvas();
    		TCanvas *c524 = new TCanvas();
    		TCanvas *c525 = new TCanvas();
    		TCanvas *c526 = new TCanvas();
    		TCanvas *c527 = new TCanvas();
    		TCanvas *c528 = new TCanvas();
    		c522->cd();
    		histNParticles->Draw("HIST");
    		c523->cd();
    		histTotalConstSize->Draw("HIST");
    		c524->cd();
    		histConstSizeB1->Draw("HIST");
    		c525->cd();
    		histConstSizeB2->Draw("HIST");
    		c526->cd();
    		histConstSizeB3->Draw("HIST");
    		c527->cd();
    		histConstSizeB4->Draw("HIST");
    		c528->cd();
    		histMinConstSize->Draw("HIST");
    		/////*/
    		
    		/*/////PFOs
    		TCanvas *c529 = new TCanvas();
    		c529->cd();
    		histNEFlowTracks->Draw("HIST");
    		TCanvas *c530 = new TCanvas();
    		c530->cd();
    		histNEFlowPhotons->Draw("HIST");
    		TCanvas *c531 = new TCanvas();
    		c531->cd();
    		histNEFlowNeutralHadrons->Draw("HIST");
    		TCanvas *c532 = new TCanvas();
    		c532->cd();
    		histNEFlowObjects->Draw("HIST");
    		TCanvas *c533 = new TCanvas();
    		c533->cd();
    		histJetB1NCharged->Draw("HIST");
    		TCanvas *c534 = new TCanvas();
    		c534->cd();
    		histJetB2NCharged->Draw("HIST");
    		TCanvas *c535 = new TCanvas();
    		c535->cd();
    		histJetB3NCharged->Draw("HIST");
    		TCanvas *c536 = new TCanvas();
    		c536->cd();
    		histJetB4NCharged->Draw("HIST");
    		TCanvas *c537 = new TCanvas();
    		c537->cd();
    		histJetB1NNeutrals->Draw("HIST");
    		TCanvas *c538 = new TCanvas();
    		c538->cd();
    		histJetB2NNeutrals->Draw("HIST");
    		TCanvas *c539 = new TCanvas();
    		c539->cd();
    		histJetB3NNeutrals->Draw("HIST");
    		TCanvas *c540 = new TCanvas();
    		c540->cd();
    		histJetB4NNeutrals->Draw("HIST");
    		TCanvas *c541 = new TCanvas();
    		c541->cd();
    		histJetNObjects->Draw("HIST");
    		TCanvas *c542 = new TCanvas();
    		c542->cd();
    		histMinJetNObjects->Draw("HIST");
    		/////PFOs*/
    		
    		/*///////nJetsCompareAlgos
    		TCanvas *c543 = new TCanvas();
    		c543->cd();
    		histNJetsCompareAlgos->Draw("TEXT");
    		///////nJetsCompareAlgos*/
    		
    		/*///////invMass for antiKt pairs when 4 Durham jets
    		TCanvas *c544 = new TCanvas();
    		c544->cd();
    		histJetB1MAntiKt2Jets->Draw("HIST");
    		TCanvas *c545 = new TCanvas();
    		c545->cd();
    		histJetB2MAntiKt2Jets->Draw("HIST");
    		TCanvas *c546 = new TCanvas();
    		c546->cd();
    		histJetB1MAntiKt3Jets->Draw("HIST");
    		TCanvas *c547 = new TCanvas();
    		c547->cd();
    		histJetB2MAntiKt3Jets->Draw("HIST");
    		TCanvas *c548 = new TCanvas();
    		c548->cd();
    		histJetB1MAntiKt4Jets->Draw("HIST");
    		TCanvas *c549 = new TCanvas();
    		c549->cd();
    		histJetB2MAntiKt4Jets->Draw("HIST");
    		TCanvas *c550 = new TCanvas();
    		c550->cd();
    		histJetB1MAntiKt5Jets->Draw("HIST");
    		TCanvas *c551 = new TCanvas();
    		c551->cd();
    		histJetB2MAntiKt5Jets->Draw("HIST");
    		TCanvas *c552 = new TCanvas();
    		c552->cd();
    		histJetB1MAntiKt6Jets->Draw("HIST");
    		TCanvas *c553 = new TCanvas();
    		c553->cd();
    		histJetB2MAntiKt6Jets->Draw("HIST");
    		///////invMass for antiKt pairs when 4 Durham jets*/
	}
	if(topology == -999) 
	{
		TCanvas *c9991 = new TCanvas();
    		c9991->cd();
    		histNEFlowTracks->Draw("HIST");
	}
	/////////printing histograms
	
	if(topology == 1) cout<<"For HH: "<<endl;
	if(topology == 2) cout<<"For qq: "<<endl;
	if(topology == 3) cout<<"For ttbar: "<<endl;
	if(topology == 4) cout<<"For ZZ: "<<endl;
	if(topology == 5) cout<<"For WW: "<<endl;
	cout<<"Events that have 4b: "<<contEvents<<"    Weighted: "<<contEvents*weight<<endl;
	cout<<"Events that pass: "<<contEventsPostFilter<<"    Weighted: "<<contEventsPostFilter*weight<<endl<<endl;
	
	/*cout<<"trueSemiWWEvents: "<<trueSemiWWEvents<<endl;
	cout<<"trueFullWWEvents: "<<trueFullWWEvents<<endl;
	cout<<"trueLeptonHadronEvents: "<<trueLeptonHadronEvents<<endl;
	cout<<"trueLeptonNeutrinoEvents: "<<trueLeptonNeutrinoEvents<<endl;
	cout<<"trueHadronNeutrinoEvents: "<<trueHadronNeutrinoEvents<<endl;
	cout<<"trueFullLeptonEvents: "<<trueFullLeptonEvents<<endl;
	cout<<"trueFullHadronEvents: "<<trueFullHadronEvents<<endl;
	cout<<"trueFullNeutrinoEvents: "<<trueFullNeutrinoEvents<<endl;
	cout<<"trueHadronEvents: "<<trueHadronEvents<<endl;
	cout<<"trueLeptonEvents: "<<trueLeptonEvents<<endl;
	cout<<"trueNeutrinoEvents: "<<trueNeutrinoEvents<<endl;	
			     		
	cout<<endl<<endl;*/
			     		
}

void trueBPairMass(const char *inputFile, int topology, string size, float weight)
{
	
	if(topology == 1) cout<<"HHtrueBPairMass working!"<<endl;
	if(topology == 2) cout<<"qqtrueBPairMass working!"<<endl;
	if(topology == 3) cout<<"tttrueBPairMass working!"<<endl;
	if(topology == 4) cout<<"ZZtrueBPairMass working!"<<endl;
	if(topology == 5) cout<<"WWtrueBPairMass working!"<<endl;
	
	
	
	functionLinking();
	
	gSystem->Load("libDelphes");

  	TChain chain("Delphes");
  	chain.Add(inputFile);
  	
  	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  	Long64_t numberOfEntries = treeReader->GetEntries();
  	//numberOfEntries=20;
  	
  	string histBPairMText;
  	if(topology == 1) histBPairMText = "True b pair inv. mass (HH)"; 
  	if(topology == 2) histBPairMText = "True b pair inv. mass (qq)";
  	if(topology == 3) histBPairMText = "True b pair inv. mass (ttbar)";
  	if(topology == 4) histBPairMText = "True b pair inv. mass (ZZ)";
  	if(topology == 5) histBPairMText = "True b pair inv. mass (WW)";
  	TH1 *histBPairM = new TH1F("bpairm", histBPairMText.c_str(), 142.0, -1.0, 400);
  	
  	TClonesArray *branchParticle;
  	TClonesArray *branchEvent;
  	
  	if(size == "large")
  	{
	  	branchParticle = treeReader->UseBranch("Particle");
		branchEvent = treeReader->UseBranch("Event");
  	}
  	else if(size == "small")
  	{
  		branchParticle = treeReader->UseBranch("Particle");
		branchEvent = treeReader->UseBranch("Event");
	}
	
	GenParticle *b1;
	GenParticle *b2;
	
	for(int entry=0; entry<numberOfEntries; entry++) 
	{
		treeReader->ReadEntry(entry);
		//cout<<endl<<endl<<endl<<"Event "<<entry<<endl<<endl;
		int nParticles=0, contB=0;
		if(branchParticle->GetEntries() > 0)  nParticles =  branchParticle->GetEntries();
		for(int i=0; i<nParticles; i++)
		{
			GenParticle *particle = (GenParticle*) branchParticle->At(i); 
			int pid = abs(particle->PID);
			int m1 = findMother(branchParticle, particle, 1);
			int m2 = findMother(branchParticle, particle, 2);
			int d1 = findDaughter(branchParticle, particle, 1);
			int d2 = findDaughter(branchParticle, particle, 2);
			
			if(pid == 25 || pid == 36)
			{
				//cout<<"Index: "<<i<<"    pid: "<<pid<<endl<<"d1: "<<d1<<"    d2: "<<d2<<endl;
			}
			if(pid == 5 && m1 != 5)
			{
				//cout<<"Index: "<<i<<"    pt: "<<particle->PT<<endl;
				contB++;
				if(contB == 1) b1 = particle;
				else if(contB == 2)
				{ 
					b2 = particle;
					//break;
				}
			}
		}
		if(contB == 2)
		{
			TLorentzVector b1P4, b2P4, bPair;
			b1P4 = b1->P4();
			b2P4 = b2->P4();
			bPair = b1P4 + b2P4;
			float bPairM = bPair.M();
			histBPairM->Fill(bPairM, weight);
			//cout<<"b1 mass: "<<b1P4.M()<<"    b2 mass: "<<b2P4.M()<<endl<<"b pair mass: "<<bPairM<<endl;
		}
		//else cout<<"No entra"<<endl;
	}
	
	if(topology == 1)
	{ 
		TCanvas *c1 = new TCanvas();
		histBPairM->Draw("HIST");
	}
	if(topology == 2)
	{ 
		TCanvas *c2 = new TCanvas();
		histBPairM->Draw("HIST");
	}
	if(topology == 3)
	{ 
		TCanvas *c3 = new TCanvas();
		histBPairM->Draw("HIST");
	}
	if(topology == 4)
	{ 
		TCanvas *c4 = new TCanvas();
		histBPairM->Draw("HIST");
	}
	if(topology == 5)
	{ 
		TCanvas *c5 = new TCanvas();
		histBPairM->Draw("HIST");
	}
}

void FSRGammaGammaHHbbbbAnalysis(const char *inputFileHH, const char *inputFileqq, const char *inputFilett, const char *inputFileZZ, const char *inputFileWW, string size)
{
	cout<<"New macro working"<<endl;
	
	double weightHH=0.001225, weightqq=0.0349, weightttbar=0.503, weightZZ=0.8167, weightWW=0.5149;
	
	string jetAlgoText = "(durham rtd_cut=10) ";
	string jetAlgo = "Jet10";
  	string genJetAlgo = "GenJet10";
  	string rtdCut = "10";
  	string jetAlgoOutputTreeSTrain = "outputTreeSHHbbbbESpreadDurham"+rtdCut+"Train.root";
  	string jetAlgoOutputTreeSTest = "outputTreeSHHbbbbESpreadDurham"+rtdCut+"Test.root";
  	string jetAlgoOutputTreeBqqTrain = "outputTreeBqqHHbbbbESpreadDurham"+rtdCut+"Train.root";
  	string jetAlgoOutputTreeBqqTest = "outputTreeBqqHHbbbbESpreadDurham"+rtdCut+"Test.root";
  	string jetAlgoOutputTreeBttTrain = "outputTreeBttHHbbbbESpreadDurham"+rtdCut+"Train.root";
  	string jetAlgoOutputTreeBttTest = "outputTreeBttHHbbbbESpreadDurham"+rtdCut+"Test.root";
  	string jetAlgoOutputTreeBZZTrain = "outputTreeBZZHHbbbbESpreadDurham"+rtdCut+"Train.root";
  	string jetAlgoOutputTreeBZZTest = "outputTreeBZZHHbbbbESpreadDurham"+rtdCut+"Test.root";
  	string jetAlgoOutputTreeBWWTrain = "outputTreeBWWHHbbbbESpreadDurham"+rtdCut+"Train.root";
  	string jetAlgoOutputTreeBWWTest = "outputTreeBWWHHbbbbESpreadDurham"+rtdCut+"Test.root";
  	
  	/*string jetAlgoText = "(antiKt R=0.5) ";
  	string jetAlgo = "JetAntiKt";
  	string genJetAlgo = "GenJetAntiKt";
  	string jetAlgoOutputTreeSTrain = "outputTreeSHHbbbbESpreadAntiKtTrain.root";
  	string jetAlgoOutputTreeSTest = "outputTreeSHHbbbbESpreadAntiKtTest.root";
  	string jetAlgoOutputTreeBqqTrain = "outputTreeBqqHHbbbbESpreadAntiKtTrain.root";
  	string jetAlgoOutputTreeBqqTest = "outputTreeBqqHHbbbbESpreadAntiKtTest.root";
  	string jetAlgoOutputTreeBttTrain = "outputTreeBttHHbbbbESpreadAntiKtTrain.root";
  	string jetAlgoOutputTreeBttTest = "outputTreeBttHHbbbbESpreadAntiKtTest.root";
  	string jetAlgoOutputTreeBZZTrain = "outputTreeBZZHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBZZTest = "outputTreeBZZHHbbbbESpreadAntiKtTest.root";
  	string jetAlgoOutputTreeBWWTrain = "outputTreeBWWHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBWWTest = "outputTreeBWWHHbbbbESpreadAntiKtTest.root";*/
  	
  	cout<<"jetAlgo: "<<jetAlgoText<<endl;
  	
  	////////creation of File for TMVA for signal HH
  	TFile *outputTreeSTrain = new TFile(jetAlgoOutputTreeSTrain.c_str(), "recreate");
  	TTree TreeSTrain("TreeSTrain","a simple Tree with simple variables (Train)");
  	TFile *outputTreeSTest = new TFile(jetAlgoOutputTreeSTest.c_str(), "recreate");
  	TTree TreeSTest("TreeSTest","a simple Tree with simple variables (Test)");
  	
  	analysis(inputFileHH, 1, size, weightHH, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest);
  	
  	outputTreeSTrain->cd();
   	TreeSTrain.Write();
   	
   	outputTreeSTest->cd();
   	TreeSTest.Write();
   	
   	outputTreeSTrain->Close();
	outputTreeSTest->Close();
  	////////creation of File for TMVA for signal HH
  	
	/*////////creation of File for TMVA for back qq  	
  	TFile *outputTreeBqqTrain = new TFile(jetAlgoOutputTreeBqqTrain.c_str(), "recreate");
  	TTree TreeBqqTrain("TreeBqqTrain","a bqqimple Tree with bqqimple variables (Train)");
  	TFile *outputTreeBqqTest = new TFile(jetAlgoOutputTreeBqqTest.c_str(), "recreate");
  	TTree TreeBqqTest("TreeBqqTest","a bqqimple Tree with bqqimple variables (Test)");
  	
  	analysis(inputFileqq, 2, size, weightqq, jetAlgoText, jetAlgo, genJetAlgo, TreeBqqTrain, TreeBqqTest);
  	
  	outputTreeBqqTrain->cd();
   	TreeBqqTrain.Write();
   	
   	outputTreeBqqTest->cd();
   	TreeBqqTest.Write();
  	
  	outputTreeBqqTrain->Close();
	outputTreeBqqTest->Close();
  	////////creation of File for TMVA for back qq*/ 	
  	
  	/*////////creation of File for TMVA for back ttbar  
	TFile *outputTreeBttTrain = new TFile(jetAlgoOutputTreeBttTrain.c_str(), "recreate");
  	TTree TreeBttTrain("TreeBttTrain","a bttimple Tree with bttimple variables (Train)");
  	TFile *outputTreeBttTest = new TFile(jetAlgoOutputTreeBttTest.c_str(), "recreate");
  	TTree TreeBttTest("TreeBttTest","a bttimple Tree with bttimple variables (Test)");
  	
  	analysis(inputFilett, 3, size, weightttbar, jetAlgoText, jetAlgo, genJetAlgo, TreeBttTrain, TreeBttTest);
  	
  	outputTreeBttTrain->cd();
   	TreeBttTrain.Write();
   	
   	outputTreeBttTest->cd();
   	TreeBttTest.Write();
	
	outputTreeBttTrain->Close();
	outputTreeBttTest->Close();
  	////////creation of File for TMVA for back ttbar*/
  	
  	////////creation of File for TMVA for back ZZ  	
  	TFile *outputTreeBZZTrain = new TFile(jetAlgoOutputTreeBZZTrain.c_str(), "recreate");
  	TTree TreeBZZTrain("TreeBZZTrain","a bZZimple Tree with bZZimple variables (Train)");
  	TFile *outputTreeBZZTest = new TFile(jetAlgoOutputTreeBZZTest.c_str(), "recreate");
  	TTree TreeBZZTest("TreeBZZTest","a bZZimple Tree with bZZimple variables (Test)");
  	
  	analysis(inputFileZZ, 4, size, weightZZ, jetAlgoText, jetAlgo, genJetAlgo, TreeBZZTrain, TreeBZZTest);
  	
  	outputTreeBZZTrain->cd();
   	TreeBZZTrain.Write();
   	
   	outputTreeBZZTest->cd();
   	TreeBZZTest.Write();
  	
  	outputTreeBZZTrain->Close();
	outputTreeBZZTest->Close();
  	////////creation of File for TMVA for back ZZ
  	
  	/*///////creation of File for TMVA for back WW  	
  	TFile *outputTreeBWWTrain = new TFile(jetAlgoOutputTreeBWWTrain.c_str(), "recreate");
  	TTree TreeBWWTrain("TreeBWWTrain","a bWWimple Tree with bWWimple variables (Train)");
  	TFile *outputTreeBWWTest = new TFile(jetAlgoOutputTreeBWWTest.c_str(), "recreate");
  	TTree TreeBWWTest("TreeBWWTest","a bWWimple Tree with bWWimple variables (Test)");
  	
  	analysis(inputFileWW, 5, size, weightWW, jetAlgoText, jetAlgo, genJetAlgo, TreeBWWTrain, TreeBWWTest);
  	
  	outputTreeBWWTrain->cd();
   	TreeBWWTrain.Write();
   	
   	outputTreeBWWTest->cd();
   	TreeBWWTest.Write();
  	
  	outputTreeBWWTrain->Close();
	outputTreeBWWTest->Close();
  	////////creation of File for TMVA for back WW*/

  	
	/*//////Pt
	TCanvas *c1001 = new TCanvas();
	TCanvas *c1002 = new TCanvas();
	TCanvas *c1003 = new TCanvas();
	TCanvas *c1004 = new TCanvas();
	//////Pt*/
	
	
	/////Just to plot; no file creation
	/*TTree TreeSTrain("TreeSTrain","a simple Tree with simple variables (Train)");
  	TTree TreeSTest("TreeSTest","a simple Tree with simple variables (Test)");
	analysis(inputFileHH, 1, size, weightHH, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest);*/
	/*TTree TreeBqqTrain("TreeBqqTrain","a bqqimple Tree with simple variables (Train)");
  	TTree TreeBqqTest("TreeBqqTest","a bqqimple Tree with simple variables (Test)");
	analysis(inputFileqq, 2, size, weightqq, jetAlgoText, jetAlgo, genJetAlgo, TreeBqqTrain, TreeBqqTest);
	TTree TreeBttTrain("TreeBttTrain","a bttimple Tree with bttimple variables (Train)");
	TTree TreeBttTest("TreeBttTest","a bttimple Tree with bttimple variables (Test)");
	analysis(inputFilett, 3, size, weightttbar, jetAlgoText, jetAlgo, genJetAlgo, TreeBttTrain, TreeBttTest);
	TTree TreeBZZTrain("TreeBZZTrain","a bZZimple Tree with bZZimple variables (Train)");
	TTree TreeBZZTest("TreeBZZTest","a bZZimple Tree with bZZimple variables (Test)");
	analysis(inputFileZZ, 4, size, weightZZ, jetAlgoText, jetAlgo, genJetAlgo, TreeBZZTrain, TreeBZZTest);
	TTree TreeBWWTrain("TreeBWWTrain","a bWWimple Tree with bWWimple variables (Train)");
	TTree TreeBWWTest("TreeBWWTest","a bWWimple Tree with bWWimple variables (Test)");
	analysis(inputFileWW, 4, size, weightWW, jetAlgoText, jetAlgo, genJetAlgo, TreeBWWTrain, TreeBWWTest);
	TTree TreeDummy("dummy","dummy");
	analysis("hh", -999, size, 0.001225, jetAlgoText, jetAlgo, genJetAlgo, TreeDummy, TreeDummy); /////Just to give time for the hists. to load; delete if not plotting.*/
	/////Just to plot; no file creation
	
	/*trueBPairMass(inputFileHH, 1, size, 0.001225);
	trueBPairMass(inputFileqq, 2, size, 0.3062);
	trueBPairMass(inputFilett, 3, size, 0.503);*/

}














