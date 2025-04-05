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
#include <stdexcept>
#include <sys/stat.h>
#include <sstream>

///kinematic fit
#include "epConstrainHH.h"
#include "LorentzVectorWithErrors.h"
#include "pxyConstrainHH.h"
#include "eqmConstrainHH.h"
///////kinematic fit

#endif

using namespace std;

void functionLinking()
{	
	cout<<5+9-7<<endl;
	cout<<"functionLinking working"<<endl;
}

/////kinematic fit
Double_t funPoly4(Double_t* x, Double_t* par) {

  Double_t xx=x[0];
  Double_t yy=x[1];
  
  Double_t result = par[0]*exp(par[1]+par[2]*xx+par[3]*yy+par[4]*pow(xx,2)+par[5]*pow(yy,2)+par[6]*xx*yy+par[7]*pow(xx,3)+par[8]*pow(yy,3)+par[9]*pow(xx,2)*yy+par[10]*xx*pow(yy,2)+par[11]*pow(xx,4)+par[12]*pow(yy,4)+par[13]*pow(xx,3)*yy+par[14]*xx*pow(yy,3)+par[15]*pow(xx*yy,2));
  cout << " funPoly4 " << " x= " << xx << " " << yy 
       << " result= " << result << endl;
  return result;
}

Double_t funPoly2(Double_t* x, Double_t* par) {
  
  Double_t result = par[0]+par[1]*x[0]+par[2]*x[1]+par[3]*pow(x[0],2)+par[4]*pow(x[1],2)+par[5]*x[0]*x[1];
  cout << " funPoly2 " << " x= " << x[0] << " " << x[1] 
       << " par= "  << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5]
       << " result= " << result << endl;
  return result;
}
///////kinematic fit

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

//Function that recursively goes over the daughters of a particle until it finds a quark, a lepton or a b. It also finds if the top decays into a W. Used only for ZH.
void findFinalDecaysZH(TClonesArray *branchParticle, GenParticle *particle, vector<int>& finalDecays, int& contNeutrinoEvent, int& contLeptonEvent, int& contHadronEvent, int index) 
{	
  	int pid = abs(particle->PID);
  	finalDecays.push_back(pid);
  	int m1 = abs(findMother(branchParticle, particle, 1));
  	//cout<<endl<<endl;
  	//cout<<"Index: "<<index<<endl<<"Pid func: "<<pid<<endl<<"m1 func: "<<m1<<endl;
  	
  	
  	if((pid == 1 || pid == 2 || pid == 3 || pid == 4 || pid == 5) && (m1 == 23 || m1 == 25 || m1 == 36)) //lookinf for u, d, s, c quarks from Z or H
  	{
  		contHadronEvent++;
  		//cout<<"Entra hadron"<<endl;
  		return;
  	}
  	if((pid == 11 || pid == 13 || pid == 15 || pid == 17) && (m1 == 23 || m1 == 25 || m1 == 36))  //looking for leptons from Z or H
  	{
  		contLeptonEvent++;
  		//cout<<"Entra lepton"<<endl;
  		return;
  	}
  	if((pid == 12 || pid == 14 || pid == 16 || pid == 18) && (m1 == 23 || m1 == 25 || m1 == 36))  //looking for neutrino from Z or H
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
  		findFinalDecaysZH(branchParticle, d1, finalDecays, contNeutrinoEvent, contLeptonEvent, contHadronEvent, d1Index);
  	}
  	int d2Index = findDaughterIndex(branchParticle, particle, 2);
  	if(d2Index != -999 && d1Index != d2Index)
  	{
  		GenParticle *d2 = (GenParticle*) branchParticle->At(d2Index);
  		//cout<<"Entra d2"<<endl;
  		findFinalDecaysZH(branchParticle, d2, finalDecays, contNeutrinoEvent, contLeptonEvent, contHadronEvent, d2Index);
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

/////function that receives eta and returns theta
double findTheta(double eta)
{
	double theta = pow(TMath::E(), -eta);
	theta = 2*TMath::ATan(theta);
	
	return theta;
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

//////Function that finds the jet pair combinations that minimizes chiSquared for H mass. 
void findJetPairs(TLorentzVector& jetPairB1, TLorentzVector& jetPairB2, TLorentzVector jetB1, TLorentzVector jetB2, TLorentzVector jetB3, TLorentzVector jetB4, double& jetPairB1Index1, double& jetPairB1Index2, double& jetPairB2Index1, double& jetPairB2Index2, double& minJetChiS)
{
	double errorJetPairMass=1; ///////Should eventually measure this (? -- CHECK)
	jetPairB1=jetB1+jetB2;
	jetPairB2=jetB3+jetB4;
	double jetChiS12=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	double jetChiS12ZZ=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	jetPairB1=jetB1+jetB3;
	jetPairB2=jetB2+jetB4;
	double jetChiS13=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	double jetChiS13ZZ=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	jetPairB1=jetB1+jetB4;
	jetPairB2=jetB2+jetB3;
	double jetChiS14=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	minJetChiS=TMath::Min(TMath::Min(jetChiS12, jetChiS13), jetChiS14);
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

//////Function that finds the jet pair combinations that minimizes chiSquared for Z mass. 
void findJetPairsZZ(TLorentzVector& jetPairB1, TLorentzVector& jetPairB2, TLorentzVector jetB1, TLorentzVector jetB2, TLorentzVector jetB3, TLorentzVector jetB4, double& jetPairB1Index1, double& jetPairB1Index2, double& jetPairB2Index1, double& jetPairB2Index2, bool& flagZZMass, float& minChiSquaredZZMass, double distanceZMass, int& contEventsZWindow10, int& contEventsZWindow01, int& contEventsZWindow05, int& contEventsZWindow1, int& contEventsZWindow15, int& contEventsZWindow2, int& contEventsZWindow5, float& distanceZ1MinChiSquaredZZMass, float& distanceZ2MinChiSquaredZZMass)
{
	double errorJetPairMass=1; ///////Should eventually measure this (? -- CHECK)
	jetPairB1=jetB1+jetB2;
	jetPairB2=jetB3+jetB4;
	//if(abs(jetPairB1.M()-90)<distanceZMass && abs(jetPairB2.M()-90)<distanceZMass) flagZZMass=true; 
	double jetChiS12=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	jetPairB1=jetB1+jetB3;
	jetPairB2=jetB2+jetB4;
	//if(abs(jetPairB1.M()-90)<distanceZMass && abs(jetPairB2.M()-90)<distanceZMass) flagZZMass=true;
	double jetChiS13=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	jetPairB1=jetB1+jetB4;
	jetPairB2=jetB2+jetB3;
	//if(abs(jetPairB1.M()-90)<distanceZMass && abs(jetPairB2.M()-90)<distanceZMass) flagZZMass=true;
	double jetChiS14=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	minChiSquaredZZMass = TMath::Min(TMath::Min(jetChiS12, jetChiS13), jetChiS14);
	if(minChiSquaredZZMass==jetChiS12)
	{
		jetPairB1=jetB1+jetB2;
		jetPairB2=jetB3+jetB4;
		jetPairB1Index1=0;
		jetPairB1Index2=1;
		jetPairB2Index1=2;
		jetPairB2Index2=3;
	}
	else if(minChiSquaredZZMass==jetChiS13)
	{
		jetPairB1=jetB1+jetB3;
		jetPairB2=jetB2+jetB4;
		jetPairB1Index1=0;
		jetPairB1Index2=2;
		jetPairB2Index1=1;
		jetPairB2Index2=3;
	}
	else if(minChiSquaredZZMass==jetChiS14)
	{
		jetPairB1=jetB1+jetB4;
		jetPairB2=jetB2+jetB3;
		jetPairB1Index1=0;
		jetPairB1Index2=3;
		jetPairB2Index1=1;
		jetPairB2Index2=2;
	}
	
	if(abs(jetPairB1.M()-90)<distanceZMass && abs(jetPairB2.M()-90)<distanceZMass) flagZZMass=true; 
	if(abs(jetPairB1.M()-90)<0.1 && abs(jetPairB2.M()-90)<0.1) contEventsZWindow01++;
	if(abs(jetPairB1.M()-90)<0.5 && abs(jetPairB2.M()-90)<0.5) contEventsZWindow05++;
	if(abs(jetPairB1.M()-90)<1 && abs(jetPairB2.M()-90)<1) contEventsZWindow1++;
	if(abs(jetPairB1.M()-90)<1.5 && abs(jetPairB2.M()-90)<1.5) contEventsZWindow15++;
	if(abs(jetPairB1.M()-90)<2 && abs(jetPairB2.M()-90)<2) contEventsZWindow2++;
	if(abs(jetPairB1.M()-90)<5 && abs(jetPairB2.M()-90)<5) contEventsZWindow5++;
	if(abs(jetPairB1.M()-90)<10 && abs(jetPairB2.M()-90)<10) contEventsZWindow10++;
	distanceZ1MinChiSquaredZZMass=abs(jetPairB1.M()-90);
	distanceZ2MinChiSquaredZZMass=abs(jetPairB2.M()-90);
	
	
	return;
}
///////

//////Function that finds the jet pair combinations that minimizes chiSquared for Z mass. 
void findJetPairsZH(TLorentzVector& jetPairB1, TLorentzVector& jetPairB2, TLorentzVector jetB1, TLorentzVector jetB2, TLorentzVector jetB3, TLorentzVector jetB4, double& jetPairB1Index1, double& jetPairB1Index2, double& jetPairB2Index1, double& jetPairB2Index2, bool& flagZHMass, float& minChiSquaredZHMass, double distanceZMass, int& contEventsZWindow10, int& contEventsZWindow01, int& contEventsZWindow05, int& contEventsZWindow1, int& contEventsZWindow15, int& contEventsZWindow2, int& contEventsZWindow5, float& distanceZ1MinChiSquaredZHMass, float& distanceZ2MinChiSquaredZHMass, float& invMassZHZ, float& invMassZHH, int contEventsPostFilter)
{
	double errorJetPairMass=1; ///////Should eventually measure this (? -- CHECK)
	jetPairB1=jetB1+jetB2;
	jetPairB2=jetB3+jetB4;
	//if(abs(jetPairB1.M()-90)<distanceZMass && abs(jetPairB2.M()-90)<distanceZMass) flagZHMass=true; 
	float jetChiS12=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	float jetChiS21=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	if(contEventsPostFilter<20) 
	{
		cout<<"jetPairB1 (1,2): Mass: "<<jetPairB1.M()<<endl;
		cout<<"jetPairB2 (3,4): Mass: "<<jetPairB2.M()<<endl;
		cout<<"jetChiS12 (B1 125, B2 90): "<<jetChiS12<<endl;
		cout<<"jetChiS21 (B1 90, B2 125): "<<jetChiS21<<endl<<endl;
	}
	jetPairB1=jetB1+jetB3;
	jetPairB2=jetB2+jetB4;
	//if(abs(jetPairB1.M()-90)<distanceZMass && abs(jetPairB2.M()-90)<distanceZMass) flagZHMass=true;
	float jetChiS13=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	float jetChiS31=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	if(contEventsPostFilter<20) 
	{
		cout<<"jetPairB1 (1,3): Mass: "<<jetPairB1.M()<<endl;
		cout<<"jetPairB2 (2,4): Mass: "<<jetPairB2.M()<<endl;
		cout<<"jetChiS13 (B1 125, B2 90): "<<jetChiS13<<endl;
		cout<<"jetChiS31 (B1 90, B2 125): "<<jetChiS31<<endl<<endl;
	}
	jetPairB1=jetB1+jetB4;
	jetPairB2=jetB2+jetB3;
	//if(abs(jetPairB1.M()-90)<distanceZMass && abs(jetPairB2.M()-90)<distanceZMass) flagZHMass=true;
	float jetChiS14=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	float jetChiS41=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	if(contEventsPostFilter<20) 
	{
		cout<<"jetPairB1 (1,4): Mass: "<<jetPairB1.M()<<endl;
		cout<<"jetPairB2 (2,3): Mass: "<<jetPairB2.M()<<endl;
		cout<<"jetChiS14 (B1 125, B2 90): "<<jetChiS14<<endl;
		cout<<"jetChiS41 (B1 90, B2 125): "<<jetChiS41<<endl<<endl;
	}
	//minChiSquaredZHMass = TMath::Min(TMath::Min(TMath::Min(TMath::Min(TMath::Min(jetChiS21, jetChiS31), jetChiS41), jetChiS12), jetChiS13), jetChiS14);
	minChiSquaredZHMass = std::min({jetChiS21, jetChiS31, jetChiS41, jetChiS12, jetChiS13, jetChiS14});
	if(contEventsPostFilter<20) cout<<"minChiSquaredZHMass: "<<minChiSquaredZHMass<<"    ";
	if(minChiSquaredZHMass==jetChiS12)
	{
		jetPairB1=jetB1+jetB2;
		jetPairB2=jetB3+jetB4;
		jetPairB1Index1=0;
		jetPairB1Index2=1;
		jetPairB2Index1=2;
		jetPairB2Index2=3;
		invMassZHZ = jetPairB2.M();
		invMassZHH = jetPairB1.M();
		if(contEventsPostFilter<20) cout<<"entra 12"<<endl<<endl<<endl;
	}
	if(minChiSquaredZHMass==jetChiS21)
	{
		jetPairB1=jetB1+jetB2;
		jetPairB2=jetB3+jetB4;
		jetPairB1Index1=0;
		jetPairB1Index2=1;
		jetPairB2Index1=2;
		jetPairB2Index2=3;
		invMassZHZ = jetPairB1.M();
		invMassZHH = jetPairB2.M();
		if(contEventsPostFilter<20) cout<<"entra 21"<<endl<<endl<<endl;
	}
	else if(minChiSquaredZHMass==jetChiS13)
	{
		jetPairB1=jetB1+jetB3;
		jetPairB2=jetB2+jetB4;
		jetPairB1Index1=0;
		jetPairB1Index2=2;
		jetPairB2Index1=1;
		jetPairB2Index2=3;
		invMassZHZ = jetPairB2.M();
		invMassZHH = jetPairB1.M();
		if(contEventsPostFilter<20) cout<<"entra 13"<<endl<<endl<<endl;
	}
	else if(minChiSquaredZHMass==jetChiS31)
	{
		jetPairB1=jetB1+jetB3;
		jetPairB2=jetB2+jetB4;
		jetPairB1Index1=0;
		jetPairB1Index2=2;
		jetPairB2Index1=1;
		jetPairB2Index2=3;
		invMassZHZ = jetPairB1.M();
		invMassZHH = jetPairB2.M();
		if(contEventsPostFilter<20) cout<<"entra 31"<<endl<<endl<<endl;
	}
	else if(minChiSquaredZHMass==jetChiS14)
	{
		jetPairB1=jetB1+jetB4;
		jetPairB2=jetB2+jetB3;
		jetPairB1Index1=0;
		jetPairB1Index2=3;
		jetPairB2Index1=1;
		jetPairB2Index2=2;
		invMassZHZ = jetPairB2.M();
		invMassZHH = jetPairB1.M();
		if(contEventsPostFilter<20) cout<<"entra 14"<<endl<<endl<<endl;
	}
	else if(minChiSquaredZHMass==jetChiS41)
	{
		jetPairB1=jetB1+jetB4;
		jetPairB2=jetB2+jetB3;
		jetPairB1Index1=0;
		jetPairB1Index2=3;
		jetPairB2Index1=1;
		jetPairB2Index2=2;
		invMassZHZ = jetPairB1.M();
		invMassZHH = jetPairB2.M();
		if(contEventsPostFilter<20) cout<<"entra 41"<<endl<<endl<<endl;
	}
	else 
	{
		if(contEventsPostFilter<20) cout<<"no entra lpm"<<endl<<endl<<endl;
	}
	if(abs(jetPairB1.M()-90)<distanceZMass && abs(jetPairB2.M()-90)<distanceZMass) flagZHMass=true; 
	if(abs(jetPairB1.M()-90)<0.1 && abs(jetPairB2.M()-90)<0.1) contEventsZWindow01++;
	if(abs(jetPairB1.M()-90)<0.5 && abs(jetPairB2.M()-90)<0.5) contEventsZWindow05++;
	if(abs(jetPairB1.M()-90)<1 && abs(jetPairB2.M()-90)<1) contEventsZWindow1++;
	if(abs(jetPairB1.M()-90)<1.5 && abs(jetPairB2.M()-90)<1.5) contEventsZWindow15++;
	if(abs(jetPairB1.M()-90)<2 && abs(jetPairB2.M()-90)<2) contEventsZWindow2++;
	if(abs(jetPairB1.M()-90)<5 && abs(jetPairB2.M()-90)<5) contEventsZWindow5++;
	if(abs(jetPairB1.M()-90)<10 && abs(jetPairB2.M()-90)<10) contEventsZWindow10++;
	distanceZ1MinChiSquaredZHMass=abs(invMassZHZ-90);
	distanceZ2MinChiSquaredZHMass=abs(invMassZHH-125);
	
	
	return;
}
///////

////////Function for ZH that given two jet pairs, uses chiSquare min to say which pair recos. the Z and which the H
bool tagZHJetPairMasses(TLorentzVector& jetPairB1, TLorentzVector& jetPairB2)
{
	double errorJetPairMass=1;
	double chiSquared1=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	double chiSquared2=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
	if(chiSquared2 < chiSquared1) return true;
	else return false;
}

//////// Function to calculate chi-squared for Z (90 GeV) and H (125 GeV)
double calculateChiSquaredZH(TClonesArray *branchJet, const vector<int>& group1, const vector<int>& group2) {
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
    double chiSquared1=pow((jetPairB1.M()-125), 2)/errorJetPairMass + pow((jetPairB2.M()-90), 2)/errorJetPairMass;
	double chiSquared2=pow((jetPairB1.M()-90), 2)/errorJetPairMass + pow((jetPairB2.M()-125), 2)/errorJetPairMass;
    return TMath::Min(chiSquared1, chiSquared2);
}
///////

///////Function to try all jet pair combinations for ZH mass reco
void findBestCombinationZH(TClonesArray *branchJet, int n, vector<int>& bestGroup1, vector<int>& bestGroup2, vector<double>& chiSquares, vector<vector<int>>& group1s, vector<vector<int>>& group2s) 
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

			double chiSquared = calculateChiSquaredZH(branchJet, group1, group2);
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

// Function to calculate the magnitude of a vector
float Magnitude(TLorentzVector vec) 
{
    return vec.P();
}
//////

//////Function that returns the maxThrust of the events
float findThrust(vector<TLorentzVector>& momenta, TLorentzVector& thrustAxis)
{
	float maxThrust = 0.0;
	
	for(float theta = 0.0; theta < TMath::Pi(); theta += 0.1) 
        {
        	for(float phi = 0.0; phi < 2.0 * TMath::Pi(); phi += 0.1) 
       		{
                	TLorentzVector axis(TMath::Sin(theta) * TMath::Cos(phi), 
                                    TMath::Sin(theta) * TMath::Sin(phi), 
                                    TMath::Cos(theta), 0.0);
                
                	float sumP = 0.0;
                	float sumP_parallel = 0.0;
		        for (auto& p : momenta) 
		        {
		            sumP += p.P();
		            sumP_parallel += TMath::Abs(p.Dot(axis) / Magnitude(axis));
		        }

		        float thrust = sumP_parallel / sumP;
		        if (thrust > maxThrust) 
		        {
		            maxThrust = thrust;
		            thrustAxis = axis;
		        }
            	}
	}
	return maxThrust;
}
//////////

/*//////kinematic fit
////////function that uses eqmConstrainHH to do the jet pairing given four jets
void eqmConstrainHHPairing(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector jet3, TLorentzVector jet4, TLorentzVector& jetPair1, TLorentzVector& jetPair2, double& chi2ndf, vector<double> BfractionalEnergyErrors, vector<double> BbetaXErrors, vector<double> BbetaYErrors, vector<double> BbetaZErrors, double Ecm, bool enableExtraTries, double nSigVar, TLorentzVector& jetB1Fiteqm, TLorentzVector& jetB2Fiteqm, TLorentzVector& jetB3Fiteqm, TLorentzVector& jetB4Fiteqm)
{
	TLorentzVector jet1Copy12, jet2Copy12, jet3Copy12, jet4Copy12;
	TLorentzVector jet1Copy13, jet2Copy13, jet3Copy13, jet4Copy13;
	TLorentzVector jet1Copy14, jet2Copy14, jet3Copy14, jet4Copy14;

	jet1Copy12 = jet1;
	jet2Copy12 = jet2;
	jet3Copy12 = jet3;
	jet4Copy12 = jet4;
	jet1Copy13 = jet1;
	jet2Copy13 = jet2;
	jet3Copy13 = jet3;
	jet4Copy13 = jet4;
	jet1Copy14 = jet1;
	jet2Copy14 = jet2;
	jet3Copy14 = jet3;
	jet4Copy14 = jet4;

	vector<LorentzVectorWithErrors>* inputLVWE12 = new vector<LorentzVectorWithErrors>(4);
	vector<LorentzVectorWithErrors>* inputLVWE13 = new vector<LorentzVectorWithErrors>(4);
	vector<LorentzVectorWithErrors>* inputLVWE14 = new vector<LorentzVectorWithErrors>(4);
	inputLVWE12->at(0)=LorentzVectorWithErrors(jet1Copy12,BfractionalEnergyErrors[0]*jet1Copy12.Energy(),BbetaXErrors[0],BbetaYErrors[0],BbetaZErrors[0]);
	inputLVWE12->at(1)=LorentzVectorWithErrors(jet2Copy12,BfractionalEnergyErrors[1]*jet2Copy12.Energy(),BbetaXErrors[1],BbetaYErrors[1],BbetaZErrors[1]);
	inputLVWE12->at(2)=LorentzVectorWithErrors(jet3Copy12,BfractionalEnergyErrors[2]*jet3Copy12.Energy(),BbetaXErrors[2],BbetaYErrors[2],BbetaZErrors[2]);
	inputLVWE12->at(3)=LorentzVectorWithErrors(jet4Copy12,BfractionalEnergyErrors[3]*jet4Copy12.Energy(),BbetaXErrors[3],BbetaYErrors[3],BbetaZErrors[3]);
	inputLVWE13->at(0)=LorentzVectorWithErrors(jet1Copy13,BfractionalEnergyErrors[0]*jet1Copy13.Energy(),BbetaXErrors[0],BbetaYErrors[0],BbetaZErrors[0]);
	inputLVWE13->at(1)=LorentzVectorWithErrors(jet3Copy13,BfractionalEnergyErrors[2]*jet3Copy13.Energy(),BbetaXErrors[2],BbetaYErrors[2],BbetaZErrors[2]);
	inputLVWE13->at(2)=LorentzVectorWithErrors(jet2Copy13,BfractionalEnergyErrors[1]*jet2Copy13.Energy(),BbetaXErrors[1],BbetaYErrors[1],BbetaZErrors[1]);
	inputLVWE13->at(3)=LorentzVectorWithErrors(jet4Copy13,BfractionalEnergyErrors[3]*jet4Copy13.Energy(),BbetaXErrors[3],BbetaYErrors[3],BbetaZErrors[3]);
	inputLVWE14->at(0)=LorentzVectorWithErrors(jet1Copy14,BfractionalEnergyErrors[0]*jet1Copy14.Energy(),BbetaXErrors[0],BbetaYErrors[0],BbetaZErrors[0]);
	inputLVWE14->at(1)=LorentzVectorWithErrors(jet4Copy14,BfractionalEnergyErrors[3]*jet4Copy14.Energy(),BbetaXErrors[3],BbetaYErrors[3],BbetaZErrors[3]);
	inputLVWE14->at(2)=LorentzVectorWithErrors(jet2Copy14,BfractionalEnergyErrors[1]*jet2Copy14.Energy(),BbetaXErrors[1],BbetaYErrors[1],BbetaZErrors[1]);
	inputLVWE14->at(3)=LorentzVectorWithErrors(jet3Copy14,BfractionalEnergyErrors[2]*jet3Copy14.Energy(),BbetaXErrors[2],BbetaYErrors[2],BbetaZErrors[2]);

	cout << " signal eqmConstrainHH12 " << endl;
	eqmConstrainHH eqmCHH12(inputLVWE12,Ecm,enableExtraTries,nSigVar);
	vector<LorentzVectorWithErrors>* outputLVeqm12=eqmCHH12.NumericalMinimization(); 
	cout << " signal eqmConstrainHH13 " << endl;
	eqmConstrainHH eqmCHH13(inputLVWE13,Ecm,enableExtraTries,nSigVar);
	vector<LorentzVectorWithErrors>* outputLVeqm13=eqmCHH13.NumericalMinimization(); 
	cout << " signal eqmConstrainHH14 " << endl;
	eqmConstrainHH eqmCHH14(inputLVWE14,Ecm,enableExtraTries,nSigVar);
	vector<LorentzVectorWithErrors>* outputLVeqm14=eqmCHH14.NumericalMinimization(); 


	double chi2ndfeqm12, chi2ndfeqm13, chi2ndfeqm14;
	TLorentzVector jetB1Fiteqm12, jetB2Fiteqm12, jetPairB1Fiteqm12, jetB3Fiteqm12, jetB4Fiteqm12, jetPairB2Fiteqm12;
	TLorentzVector jetB1Fiteqm13, jetB2Fiteqm13, jetPairB1Fiteqm13, jetB3Fiteqm13, jetB4Fiteqm13, jetPairB2Fiteqm13;
	TLorentzVector jetB1Fiteqm14, jetB2Fiteqm14, jetPairB1Fiteqm14, jetB3Fiteqm14, jetB4Fiteqm14, jetPairB2Fiteqm14;

	jetB1Fiteqm12=outputLVeqm12->at(0).getLV();
	jetB2Fiteqm12=outputLVeqm12->at(1).getLV();
	jetB3Fiteqm12=outputLVeqm12->at(2).getLV();
	jetB4Fiteqm12=outputLVeqm12->at(3).getLV();
	chi2ndfeqm12=outputLVeqm12->at(0).getChi2perNDF();
	jetB1Fiteqm13=outputLVeqm13->at(0).getLV();
	jetB2Fiteqm13=outputLVeqm13->at(1).getLV();
	jetB3Fiteqm13=outputLVeqm13->at(2).getLV();
	jetB4Fiteqm13=outputLVeqm13->at(3).getLV();
	chi2ndfeqm13=outputLVeqm13->at(0).getChi2perNDF();
	jetB1Fiteqm14=outputLVeqm14->at(0).getLV();
	jetB2Fiteqm14=outputLVeqm14->at(1).getLV();
	jetB3Fiteqm14=outputLVeqm14->at(2).getLV();
	jetB4Fiteqm14=outputLVeqm14->at(3).getLV();
	chi2ndfeqm14=outputLVeqm14->at(0).getChi2perNDF();

	if(chi2ndfeqm12 < chi2ndfeqm13 && chi2ndfeqm12 < chi2ndfeqm14)
	{
		jetPair1 = jetB1Fiteqm12 + jetB2Fiteqm12;
		jetPair2 = jetB3Fiteqm12 + jetB4Fiteqm12;
		jetB1Fiteqm = jetB1Fiteqm12;
		jetB2Fiteqm = jetB2Fiteqm12;
		jetB3Fiteqm = jetB3Fiteqm12;
		jetB4Fiteqm = jetB4Fiteqm12;
		chi2ndf = chi2ndfeqm12;
	}
	if(chi2ndfeqm13 < chi2ndfeqm12 && chi2ndfeqm13 < chi2ndfeqm14)
	{
		jetPair1 = jetB1Fiteqm13 + jetB2Fiteqm13;
		jetPair2 = jetB3Fiteqm13 + jetB4Fiteqm13;
		jetB1Fiteqm = jetB1Fiteqm13;
		jetB2Fiteqm = jetB2Fiteqm13;
		jetB3Fiteqm = jetB3Fiteqm13;
		jetB4Fiteqm = jetB4Fiteqm13;
		chi2ndf = chi2ndfeqm13;
	}
	if(chi2ndfeqm14 < chi2ndfeqm12 && chi2ndfeqm14 < chi2ndfeqm13)
	{
		jetPair1 = jetB1Fiteqm14 + jetB2Fiteqm14;
		jetPair2 = jetB3Fiteqm14 + jetB4Fiteqm14;
		jetB1Fiteqm = jetB1Fiteqm14;
		jetB2Fiteqm = jetB2Fiteqm14;
		jetB3Fiteqm = jetB3Fiteqm14;
		jetB4Fiteqm = jetB4Fiteqm14;
		chi2ndf = chi2ndfeqm14;
	}

	delete inputLVWE12;
	delete inputLVWE13;
	delete inputLVWE14;
	delete outputLVeqm12;
	delete outputLVeqm13;
	delete outputLVeqm14;
}
///////kinematic fit*/

/*//////kinematic fit pythia
////////function that uses eqmConstrainHH to do the jet pairing given four jets
void eqmConstrainHHPairingPythia(vector<TLorentzVector> jetBoffset, vector<TLorentzVector> jetsOrderedPT, TLorentzVector& jetPair1, TLorentzVector& jetPair2, double& chi2ndf, vector<double> BfractionalEnergyErrors, vector<double> BbetaXErrors, vector<double> BbetaYErrors, vector<double> BbetaZErrors, double Ecm, bool enableExtraTries, double nSigVar, TLorentzVector& jetB1Fiteqm, TLorentzVector& jetB2Fiteqm, TLorentzVector& jetB3Fiteqm, TLorentzVector& jetB4Fiteqm)
{
	TLorentzVector jet1Copy12, jet2Copy12, jet3Copy12, jet4Copy12;
	TLorentzVector jet1Copy13, jet2Copy13, jet3Copy13, jet4Copy13;
	TLorentzVector jet1Copy14, jet2Copy14, jet3Copy14, jet4Copy14;

	jet1Copy12 = jetsOrderedPT[0];
	jet2Copy12 = jetsOrderedPT[1];
	jet3Copy12 = jetsOrderedPT[2];
	jet4Copy12 = jetsOrderedPT[3];
	jet1Copy13 = jetsOrderedPT[0];
	jet2Copy13 = jetsOrderedPT[1];
	jet3Copy13 = jetsOrderedPT[2];
	jet4Copy13 = jetsOrderedPT[3];
	jet1Copy14 = jetsOrderedPT[0];
	jet2Copy14 = jetsOrderedPT[1];
	jet3Copy14 = jetsOrderedPT[2];
	jet4Copy14 = jetsOrderedPT[3];

	jet1Copy12Offset = jetBoffset[0];
	jet2Copy12Offset = jetBoffset[1];
	jet3Copy12Offset = jetBoffset[2];
	jet4Copy12Offset = jetBoffset[3];
	jet1Copy13Offset = jetBoffset[0];
	jet2Copy13Offset = jetBoffset[1];
	jet3Copy13Offset = jetBoffset[2];
	jet4Copy13Offset = jetBoffset[3];
	jet1Copy14Offset = jetBoffset[0];
	jet2Copy14Offset = jetBoffset[1];
	jet3Copy14Offset = jetBoffset[2];
	jet4Copy14Offset = jetBoffset[3];

	vector<LorentzVectorWithErrors>* inputLVWE12 = new vector<LorentzVectorWithErrors>(4);
	vector<LorentzVectorWithErrors>* inputLVWE13 = new vector<LorentzVectorWithErrors>(4);
	vector<LorentzVectorWithErrors>* inputLVWE14 = new vector<LorentzVectorWithErrors>(4);
	inputLVWE12->at(0)=LorentzVectorWithErrors(jet1Copy12Offset,BfractionalEnergyErrors[0]*jet1Copy12.Energy(),BbetaXErrors[0],BbetaYErrors[0],BbetaZErrors[0]);
	inputLVWE12->at(1)=LorentzVectorWithErrors(jet2Copy12Offset,BfractionalEnergyErrors[1]*jet2Copy12.Energy(),BbetaXErrors[1],BbetaYErrors[1],BbetaZErrors[1]);
	inputLVWE12->at(2)=LorentzVectorWithErrors(jet3Copy12Offset,BfractionalEnergyErrors[2]*jet3Copy12.Energy(),BbetaXErrors[2],BbetaYErrors[2],BbetaZErrors[2]);
	inputLVWE12->at(3)=LorentzVectorWithErrors(jet4Copy12Offset,BfractionalEnergyErrors[3]*jet4Copy12.Energy(),BbetaXErrors[3],BbetaYErrors[3],BbetaZErrors[3]);
	inputLVWE13->at(0)=LorentzVectorWithErrors(jet1Copy13Offset,BfractionalEnergyErrors[0]*jet1Copy13.Energy(),BbetaXErrors[0],BbetaYErrors[0],BbetaZErrors[0]);
	inputLVWE13->at(1)=LorentzVectorWithErrors(jet3Copy13Offset,BfractionalEnergyErrors[2]*jet3Copy13.Energy(),BbetaXErrors[2],BbetaYErrors[2],BbetaZErrors[2]);
	inputLVWE13->at(2)=LorentzVectorWithErrors(jet2Copy13Offset,BfractionalEnergyErrors[1]*jet2Copy13.Energy(),BbetaXErrors[1],BbetaYErrors[1],BbetaZErrors[1]);
	inputLVWE13->at(3)=LorentzVectorWithErrors(jet4Copy13Offset,BfractionalEnergyErrors[3]*jet4Copy13.Energy(),BbetaXErrors[3],BbetaYErrors[3],BbetaZErrors[3]);
	inputLVWE14->at(0)=LorentzVectorWithErrors(jet1Copy14Offset,BfractionalEnergyErrors[0]*jet1Copy14.Energy(),BbetaXErrors[0],BbetaYErrors[0],BbetaZErrors[0]);
	inputLVWE14->at(1)=LorentzVectorWithErrors(jet4Copy14Offset,BfractionalEnergyErrors[3]*jet4Copy14.Energy(),BbetaXErrors[3],BbetaYErrors[3],BbetaZErrors[3]);
	inputLVWE14->at(2)=LorentzVectorWithErrors(jet2Copy14Offset,BfractionalEnergyErrors[1]*jet2Copy14.Energy(),BbetaXErrors[1],BbetaYErrors[1],BbetaZErrors[1]);
	inputLVWE14->at(3)=LorentzVectorWithErrors(jet3Copy14Offset,BfractionalEnergyErrors[2]*jet3Copy14.Energy(),BbetaXErrors[2],BbetaYErrors[2],BbetaZErrors[2]);

	//cout << " signal eqmConstrainHH12 " << endl;
	eqmConstrainHH eqmCHH12(inputLVWE12,Ecm,enableExtraTries,nSigVar);
	vector<LorentzVectorWithErrors>* outputLVeqm12=eqmCHH12.NumericalMinimization(); 
	//cout << " signal eqmConstrainHH13 " << endl;
	eqmConstrainHH eqmCHH13(inputLVWE13,Ecm,enableExtraTries,nSigVar);
	vector<LorentzVectorWithErrors>* outputLVeqm13=eqmCHH13.NumericalMinimization(); 
	//cout << " signal eqmConstrainHH14 " << endl;
	eqmConstrainHH eqmCHH14(inputLVWE14,Ecm,enableExtraTries,nSigVar);
	vector<LorentzVectorWithErrors>* outputLVeqm14=eqmCHH14.NumericalMinimization(); 


	double chi2ndfeqm12, chi2ndfeqm13, chi2ndfeqm14;
	TLorentzVector jetB1Fiteqm12, jetB2Fiteqm12, jetPairB1Fiteqm12, jetB3Fiteqm12, jetB4Fiteqm12, jetPairB2Fiteqm12;
	TLorentzVector jetB1Fiteqm13, jetB2Fiteqm13, jetPairB1Fiteqm13, jetB3Fiteqm13, jetB4Fiteqm13, jetPairB2Fiteqm13;
	TLorentzVector jetB1Fiteqm14, jetB2Fiteqm14, jetPairB1Fiteqm14, jetB3Fiteqm14, jetB4Fiteqm14, jetPairB2Fiteqm14;

	jetB1Fiteqm12=outputLVeqm12->at(0).getLV();
	jetB2Fiteqm12=outputLVeqm12->at(1).getLV();
	jetB3Fiteqm12=outputLVeqm12->at(2).getLV();
	jetB4Fiteqm12=outputLVeqm12->at(3).getLV();
	chi2ndfeqm12=outputLVeqm12->at(0).getChi2perNDF();
	jetB1Fiteqm13=outputLVeqm13->at(0).getLV();
	jetB2Fiteqm13=outputLVeqm13->at(1).getLV();
	jetB3Fiteqm13=outputLVeqm13->at(2).getLV();
	jetB4Fiteqm13=outputLVeqm13->at(3).getLV();
	chi2ndfeqm13=outputLVeqm13->at(0).getChi2perNDF();
	jetB1Fiteqm14=outputLVeqm14->at(0).getLV();
	jetB2Fiteqm14=outputLVeqm14->at(1).getLV();
	jetB3Fiteqm14=outputLVeqm14->at(2).getLV();
	jetB4Fiteqm14=outputLVeqm14->at(3).getLV();
	chi2ndfeqm14=outputLVeqm14->at(0).getChi2perNDF();

	if(chi2ndfeqm12 < chi2ndfeqm13 && chi2ndfeqm12 < chi2ndfeqm14)
	{
		jetPair1 = jetB1Fiteqm12 + jetB2Fiteqm12;
		jetPair2 = jetB3Fiteqm12 + jetB4Fiteqm12;
		jetB1Fiteqm = jetB1Fiteqm12;
		jetB2Fiteqm = jetB2Fiteqm12;
		jetB3Fiteqm = jetB3Fiteqm12;
		jetB4Fiteqm = jetB4Fiteqm12;
		chi2ndf = chi2ndfeqm12;
	}
	if(chi2ndfeqm13 < chi2ndfeqm12 && chi2ndfeqm13 < chi2ndfeqm14)
	{
		jetPair1 = jetB1Fiteqm13 + jetB2Fiteqm13;
		jetPair2 = jetB3Fiteqm13 + jetB4Fiteqm13;
		jetB1Fiteqm = jetB1Fiteqm13;
		jetB2Fiteqm = jetB2Fiteqm13;
		jetB3Fiteqm = jetB3Fiteqm13;
		jetB4Fiteqm = jetB4Fiteqm13;
		chi2ndf = chi2ndfeqm13;
	}
	if(chi2ndfeqm14 < chi2ndfeqm12 && chi2ndfeqm14 < chi2ndfeqm13)
	{
		jetPair1 = jetB1Fiteqm14 + jetB2Fiteqm14;
		jetPair2 = jetB3Fiteqm14 + jetB4Fiteqm14;
		jetB1Fiteqm = jetB1Fiteqm14;
		jetB2Fiteqm = jetB2Fiteqm14;
		jetB3Fiteqm = jetB3Fiteqm14;
		jetB4Fiteqm = jetB4Fiteqm14;
		chi2ndf = chi2ndfeqm14;
	}

	delete inputLVWE12;
	delete inputLVWE13;
	delete inputLVWE14;
	delete outputLVeqm12;
	delete outputLVeqm13;
	delete outputLVeqm14;
}
///////kinematic fit pythia*/

/*gamma gamma -> xx analysis. For topology: 
1) HH
2) qqbar
3) ttbar
4) ZZ
5) WW
6) eGamma->qqX
7) eGamma->qqqqX
8) eGamma->qqHX
9) ZH
10) pebb
11) pebbqq
12) peqqH
13) pett*/ 
void analysis(vector<const char *> inputFiles, int topology, float weight, string jetAlgoText, string jetAlgo, string genJetAlgo, TTree& TreeTrain, TTree& TreeTest, TTree& TreeMerge, TTree& TreeFull, string fileFunction, string preselection)
{
	if(topology == 1) cout<<"HHAnalysis working!"<<endl;
	else if(topology == 2) cout<<"qqAnalysis working!"<<endl;
	else if(topology == 3) cout<<"ttAnalysis working!"<<endl;
	else if(topology == 4) cout<<"ZZAnalysis working!"<<endl;
	else if(topology == 5) cout<<"WWAnalysis working!"<<endl;
	else if(topology == 6) cout<<"qqXAnalysis working!"<<endl;
	else if(topology == 7) cout<<"qqqqXAnalysis working!"<<endl;
	else if(topology == 8) cout<<"qqHXAnalysis working!"<<endl;
	else if(topology == 9) cout<<"ZHAnalysis working!"<<endl;
	else if(topology == 10) cout<<"pebbAnalysis working!"<<endl;
	else if(topology == 11) cout<<"pebbqqAnalysis working!"<<endl;
	else if(topology == 12) cout<<"peqqHAnalysis working!"<<endl;
	else if(topology == 13) cout<<"pettAnalysis working!"<<endl;
	//else if(topology == -999) return;
	
	
	functionLinking();
	
	gSystem->Load("libDelphes");

  	TChain chain("Delphes");
	for(int file=0;file<inputFiles.size();file++)
	{
		const char *inputFileSen = inputFiles[file];
		chain.Add(inputFileSen);
	}
  	
  	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  	Long64_t numberOfEntries = treeReader->GetEntries();
	//numberOfEntries=10000;
  	//int pos=7;
  	int contEntriesPostFilter=0;
  	cout<<"numberOfEntries: "<<numberOfEntries<<endl;

	///////kinematic fit
	/* int nBins2dMass=380;
	double min2dMass=0.;
	double max2dMass=380.; */
	
	int nBins2dMass=180;
	double min2dMass=-5.0;
	double max2dMass=250.;
	
	/*int nBins2dMass=36;
	double min2dMass=114.;
	double max2dMass=150.; */
	
	/*  int nBins2dMass=80;
	double min2dMass=110.;
	double max2dMass=150.; */
	///////kinematic fit
  	
  	string histJetEtaText, histJet1EtaText, histJet2EtaText, histJet3EtaText, histJet4EtaText, histJetCosThetaText, histJet1CosThetaText, histJet2CosThetaText, histJet3CosThetaText, histJet4CosThetaText, histSumJetPtText, histJetPtText, histJet1PtText, histJet2PtText, histJet3PtText, histJet4PtText, histJetB1MText, histJetB2MText, histMinJetMText, histJet1MText, histJet2MText, histJet3MText, histJet4MText, histSpherText, histAplanText, histNParticlesText, histTotalConstSizeText, histConstSizeB1Text, histConstSizeB2Text, histConstSizeB3Text, histConstSizeB4Text, histMinConstSizeText, histNEFlowTracksText, histNEFlowPhotonsText, histNEFlowNeutralHadronsText, histNEFlowObjectsText, histJetB1NChargedText, histJetB2NChargedText, histJetB3NChargedText, histJetB4NChargedText, histJetB1NNeutralsText, histJetB2NNeutralsText, histJetB3NNeutralsText, histJetB4NNeutralsText, histJetNObjectsText, histMinJetNObjectsText, histNJetsCompareAlgosText, histJetB1MAntiKt2JetsText, histJetB2MAntiKt2JetsText, histJetB1MAntiKt3JetsText, histJetB2MAntiKt3JetsText, histJetB1MAntiKt4JetsText, histJetB2MAntiKt4JetsText, histJetB1MAntiKt5JetsText, histJetB2MAntiKt5JetsText, histJetB1MAntiKt6JetsText, histJetB2MAntiKt6JetsText, histNJetsDurham0Text, histNJetsDurham5Text, histNJetsDurham10Text, histNJetsDurham15Text, histNJetsDurham20Text, histNJetsDurham25Text, histNJetsDurham30Text, histMinChiSquaredZZMassText, histInvMassZZ1Text, histInvMassZZ2Text, histDistanceZ1MinChiSquaredZZMassText, histDistanceZ2MinChiSquaredZZMassText, histMinChiSquaredZHMassText, histInvMassZH1Text, histInvMassZH2Text, histDistanceZ1MinChiSquaredZHMassText, histDistanceZ2MinChiSquaredZHMassText, histExclYmerge12Text, histExclYmerge23Text, histExclYmerge34Text, histExclYmerge45Text, histExclYmerge56Text, histJetB1M1BestText, histJetB2M1BestText, histDiHMText, histInvMassB1FitepText, histInvMassB2FitepText, histChi2ndfepText, histInvMassB1FitpxyText, histInvMassB2FitpxyText, histChi2ndfpxyText, histInvMassB1FiteqmText, histInvMassB2FiteqmText, histChi2ndfeqmText, histInvMassB1FitBestText, histInvMassB2FitBestText, histChi2ndfBestText, histMinJetChiSText, histInvMass4JetsText, histMinDeltaRJetsText, histDeltaRJetPairsText, histEAssymetryText, histPTAssymetryText, histNSubJetsTrimmedB1Text, histNSubJetsTrimmedB2Text, histNSubJetsTrimmedB3Text, histNSubJetsTrimmedB4Text, histNSubJetsPrunnedB1Text, histNSubJetsPrunnedB2Text, histNSubJetsPrunnedB3Text, histNSubJetsPrunnedB4Text, histNSubJetsSoftDroppedB1Text, histNSubJetsSoftDroppedB2Text, histNSubJetsSoftDroppedB3Text, histNSubJetsSoftDroppedB4Text, histJetPairMassDeltaText; 
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
  		histJetB2MText = jetAlgoText + "b-tagged jet-pair 2 mass (HH";
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
  		
  		histNJetsDurham0Text = jetAlgoText + "nJets in durham with rtdCut = 0 (HH)"; 
  		histNJetsDurham5Text = jetAlgoText + "nJets in durham with rtdCut = 5 (HH)"; 
  		histNJetsDurham10Text = jetAlgoText + "nJets in durham with rtdCut = 10 (HH)"; 
  		histNJetsDurham15Text = jetAlgoText + "nJets in durham with rtdCut = 15 (HH)"; 
  		histNJetsDurham20Text = jetAlgoText + "nJets in durham with rtdCut = 20 (HH)"; 
  		histNJetsDurham25Text = jetAlgoText + "nJets in durham with rtdCut = 25 (HH)"; 
  		histNJetsDurham30Text = jetAlgoText + "nJets in durham with rtdCut = 30 (HH)"; 
  		
  		histMinChiSquaredZZMassText = jetAlgoText + "minChiSquared for the three possible jet pair comb. for ZZ mass (HH)"; 
  		histInvMassZZ1Text =  jetAlgoText + "Mass of jetPair1 of comb. with min. chiSquared for ZZ (HH)";
  		histInvMassZZ2Text =  jetAlgoText + "Mass of jetPair2 of comb. with min. chiSquared for ZZ (HH)";
  		histDistanceZ1MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair1 of comb. with min. chiSquared for ZZ (HH)";
  		histDistanceZ2MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair2 of comb. with min. chiSquared for ZZ (HH)";

		histMinChiSquaredZHMassText = jetAlgoText + "minChiSquared for the three possible jet pair comb. for ZH mass (HH)"; 
  		histInvMassZH1Text =  jetAlgoText + "Mass of jetPair1 of comb. with min. chiSquared for ZH (HH)";
  		histInvMassZH2Text =  jetAlgoText + "Mass of jetPair2 of comb. with min. chiSquared for ZH (HH)";
  		histDistanceZ1MinChiSquaredZHMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair1 of comb. with min. chiSquared for ZH (HH)";
  		histDistanceZ2MinChiSquaredZHMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair2 of comb. with min. chiSquared for ZH (HH)";
  		
  		histExclYmerge12Text = jetAlgoText + "y value when merging from 1 to 2 jets (HH)"; 
  		histExclYmerge23Text = jetAlgoText + "y value when merging from 2 to 3 jets (HH)"; 
  		histExclYmerge34Text = jetAlgoText + "y value when merging from 3 to 4 jets (HH)"; 
  		histExclYmerge45Text = jetAlgoText + "y value when merging from 4 to 5 jets (HH)"; 
  		histExclYmerge56Text = jetAlgoText + "y value when merging from 5 to 6 jets (HH)"; 

		histJetB1M1BestText = jetAlgoText + "b-tagged jet-pair 1 mass for best comb. for HH mass (HH";
  		histJetB2M1BestText = jetAlgoText + "b-tagged jet-pair 2 mass for best comb. for HH mass (HH)";

		histDiHMText = jetAlgoText + "Di-Higgs Inv. Mass";

		histInvMassB1FitepText = jetAlgoText + "fitted (ep) b-tagged jet-pair 1 mass (HH)";
  		histInvMassB2FitepText = jetAlgoText + "fitted (ep) b-tagged jet-pair 2 mass (HH)";
		histChi2ndfepText = jetAlgoText + "chi2ndfep (ep) (HH)";
		histInvMassB1FitpxyText = jetAlgoText + "fitted (pxy) b-tagged jet-pair 1 mass (HH)";
  		histInvMassB2FitpxyText = jetAlgoText + "fitted (pxy) b-tagged jet-pair 2 mass (HH)";
		histChi2ndfpxyText = jetAlgoText + "chi2ndfpxy (pxy) (HH)";
		histInvMassB1FiteqmText = jetAlgoText + "fitted (eqm, paired with eqm) b-tagged jet-pair 1 mass (HH)";
  		histInvMassB2FiteqmText = jetAlgoText + "fitted (eqm, paired with eqm) b-tagged jet-pair 2 mass (HH)";
		histChi2ndfeqmText = jetAlgoText + "chi2ndfeqm (eqm, paired with eqm) (HH)";
		histInvMassB1FitBestText = jetAlgoText + "fitted (least chi2nf) b-tagged jet-pair 1 mass (HH)";
  		histInvMassB2FitBestText = jetAlgoText + "fitted (least chi2nf) b-tagged jet-pair 2 mass (HH)";
		histChi2ndfBestText = jetAlgoText + "chi2ndfeqm (least chi2nf) (HH)";

		histMinJetChiSText = jetAlgoText + "min. jet chi2nd (HH)";

		histInvMass4JetsText = jetAlgoText + "4-jet Inv. Mass (HH)";
		histMinDeltaRJetsText = jetAlgoText + "min. deltaR between jets (HH)";
		histDeltaRJetPairsText = jetAlgoText + "deltaR between jet-pairs (HH)";
		histEAssymetryText = jetAlgoText + "E assymetry (HH)";
		histPTAssymetryText = jetAlgoText + "PT assymetry (HH)";

		histNSubJetsTrimmedB1Text = jetAlgoText + "nSubJetsTrimmed for jet B1 (HH)";
		histNSubJetsTrimmedB2Text = jetAlgoText + "nSubJetsTrimmed for jet B2 (HH)";
		histNSubJetsTrimmedB3Text = jetAlgoText + "nSubJetsTrimmed for jet B3 (HH)";
		histNSubJetsTrimmedB4Text = jetAlgoText + "nSubJetsTrimmed for jet B4 (HH)";
		histNSubJetsPrunnedB1Text = jetAlgoText + "nSubJetsPrunned for jet B1 (HH)";
		histNSubJetsPrunnedB2Text = jetAlgoText + "nSubJetsPrunned for jet B2 (HH)";
		histNSubJetsPrunnedB3Text = jetAlgoText + "nSubJetsPrunned for jet B3 (HH)";
		histNSubJetsPrunnedB4Text = jetAlgoText + "nSubJetsPrunned for jet B4 (HH)";
		histNSubJetsSoftDroppedB1Text = jetAlgoText + "nSubJetsSoftDropped for jet B1 (HH)";
		histNSubJetsSoftDroppedB2Text = jetAlgoText + "nSubJetsSoftDropped for jet B2 (HH)";
		histNSubJetsSoftDroppedB3Text = jetAlgoText + "nSubJetsSoftDropped for jet B3 (HH)";
		histNSubJetsSoftDroppedB4Text = jetAlgoText + "nSubJetsSoftDropped for jet B4 (HH)";

		histJetPairMassDeltaText = jetAlgoText + "JetPairMassDelta (HH)";
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
  		
  		histNJetsDurham0Text = jetAlgoText + "nJets in durham with rtdCut = 0 (qq)"; 
  		histNJetsDurham5Text = jetAlgoText + "nJets in durham with rtdCut = 5 (qq)"; 
  		histNJetsDurham10Text = jetAlgoText + "nJets in durham with rtdCut = 10 (qq)"; 
  		histNJetsDurham15Text = jetAlgoText + "nJets in durham with rtdCut = 15 (qq)"; 
  		histNJetsDurham20Text = jetAlgoText + "nJets in durham with rtdCut = 20 (qq)"; 
  		histNJetsDurham25Text = jetAlgoText + "nJets in durham with rtdCut = 25 (qq)"; 
  		histNJetsDurham30Text = jetAlgoText + "nJets in durham with rtdCut = 30 (qq)";
  		
  		histMinChiSquaredZZMassText = jetAlgoText + "minChiSquared for the three possible jet pair comb. for ZZ mass (qq)";
  		histInvMassZZ1Text =  jetAlgoText + "Mass of jetPair1 of comb. with min. chiSquared for ZZ (qq)";
  		histInvMassZZ2Text =  jetAlgoText + "Mass of jetPair2 of comb. with min. chiSquared for ZZ (qq)"; 
  		histDistanceZ1MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair1 of comb. with min. chiSquared for ZZ (qq)";
  		histDistanceZ2MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair2 of comb. with min. chiSquared for ZZ (qq)";
  		
  		histExclYmerge12Text = jetAlgoText + "y value when merging from 1 to 2 jets (qq)"; 
  		histExclYmerge23Text = jetAlgoText + "y value when merging from 2 to 3 jets (qq)"; 
  		histExclYmerge34Text = jetAlgoText + "y value when merging from 3 to 4 jets (qq)"; 
  		histExclYmerge45Text = jetAlgoText + "y value when merging from 4 to 5 jets (qq)"; 
  		histExclYmerge56Text = jetAlgoText + "y value when merging from 5 to 6 jets (qq)"; 
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
  		
  		histNJetsDurham0Text = jetAlgoText + "nJets in durham with rtdCut = 0 (ttbar)"; 
  		histNJetsDurham5Text = jetAlgoText + "nJets in durham with rtdCut = 5 (ttbar)"; 
  		histNJetsDurham10Text = jetAlgoText + "nJets in durham with rtdCut = 10 (ttbar)"; 
  		histNJetsDurham15Text = jetAlgoText + "nJets in durham with rtdCut = 15 (ttbar)"; 
  		histNJetsDurham20Text = jetAlgoText + "nJets in durham with rtdCut = 20 (ttbar)"; 
  		histNJetsDurham25Text = jetAlgoText + "nJets in durham with rtdCut = 25 (ttbar)"; 
  		histNJetsDurham30Text = jetAlgoText + "nJets in durham with rtdCut = 30 (ttbar)";
  		
  		histMinChiSquaredZZMassText = jetAlgoText + "minChiSquared for the three possible jet pair comb. for ZZ mass (ttbar)"; 
  		histInvMassZZ1Text =  jetAlgoText + "Mass of jetPair1 of comb. with min. chiSquared for ZZ (ttbar)";
  		histInvMassZZ2Text =  jetAlgoText + "Mass of jetPair2 of comb. with min. chiSquared for ZZ (ttbar)";
  		histDistanceZ1MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair1 of comb. with min. chiSquared for ZZ (ttbar)";
  		histDistanceZ2MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair2 of comb. with min. chiSquared for ZZ (ttbar)";
  		
  		histExclYmerge12Text = jetAlgoText + "y value when merging from 1 to 2 jets (ttbar)"; 
  		histExclYmerge23Text = jetAlgoText + "y value when merging from 2 to 3 jets (ttbar)"; 
  		histExclYmerge34Text = jetAlgoText + "y value when merging from 3 to 4 jets (ttbar)"; 
  		histExclYmerge45Text = jetAlgoText + "y value when merging from 4 to 5 jets (ttbar)"; 
  		histExclYmerge56Text = jetAlgoText + "y value when merging from 5 to 6 jets (ttbar)"; 
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
  		
  		histNJetsDurham0Text = jetAlgoText + "nJets in durham with rtdCut = 0 (ZZ)"; 
  		histNJetsDurham5Text = jetAlgoText + "nJets in durham with rtdCut = 5 (ZZ)"; 
  		histNJetsDurham10Text = jetAlgoText + "nJets in durham with rtdCut = 10 (ZZ)"; 
  		histNJetsDurham15Text = jetAlgoText + "nJets in durham with rtdCut = 15 (ZZ)"; 
  		histNJetsDurham20Text = jetAlgoText + "nJets in durham with rtdCut = 20 (ZZ)"; 
  		histNJetsDurham25Text = jetAlgoText + "nJets in durham with rtdCut = 25 (ZZ)"; 
  		histNJetsDurham30Text = jetAlgoText + "nJets in durham with rtdCut = 30 (ZZ)";
  		
  		histMinChiSquaredZZMassText = jetAlgoText + "minChiSquared for the three possible jet pair comb. for ZZ mass (ZZ)"; 
  		histInvMassZZ1Text =  jetAlgoText + "Mass of jetPair1 of comb. with min. chiSquared for ZZ (ZZ)";
  		histInvMassZZ2Text =  jetAlgoText + "Mass of jetPair2 of comb. with min. chiSquared for ZZ (ZZ)";
  		histDistanceZ1MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair1 of comb. with min. chiSquared for ZZ (ZZ)";
  		histDistanceZ2MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair2 of comb. with min. chiSquared for ZZ (ZZ)";
  		
  		histExclYmerge12Text = jetAlgoText + "y value when merging from 1 to 2 jets (ZZ)"; 
  		histExclYmerge23Text = jetAlgoText + "y value when merging from 2 to 3 jets (ZZ)"; 
  		histExclYmerge34Text = jetAlgoText + "y value when merging from 3 to 4 jets (ZZ)"; 
  		histExclYmerge45Text = jetAlgoText + "y value when merging from 4 to 5 jets (ZZ)"; 
  		histExclYmerge56Text = jetAlgoText + "y value when merging from 5 to 6 jets (ZZ)"; 
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
  		
  		histNJetsDurham0Text = jetAlgoText + "nJets in durham with rtdCut = 0 (WW)"; 
  		histNJetsDurham5Text = jetAlgoText + "nJets in durham with rtdCut = 5 (WW)"; 
  		histNJetsDurham10Text = jetAlgoText + "nJets in durham with rtdCut = 10 (WW)"; 
  		histNJetsDurham15Text = jetAlgoText + "nJets in durham with rtdCut = 15 (WW)"; 
  		histNJetsDurham20Text = jetAlgoText + "nJets in durham with rtdCut = 20 (WW)"; 
  		histNJetsDurham25Text = jetAlgoText + "nJets in durham with rtdCut = 25 (WW)"; 
  		histNJetsDurham30Text = jetAlgoText + "nJets in durham with rtdCut = 30 (WW)";
  		
  		histMinChiSquaredZZMassText = jetAlgoText + "minChiSquared for the three possible jet pair comb. for ZZ mass (WW)";
  		histInvMassZZ1Text =  jetAlgoText + "Mass of jetPair1 of comb. with min. chiSquared for ZZ (WW)";
  		histInvMassZZ2Text =  jetAlgoText + "Mass of jetPair2 of comb. with min. chiSquared for ZZ (WW)";
  		histDistanceZ1MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair1 of comb. with min. chiSquared for ZZ (WW)";
  		histDistanceZ2MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair2 of comb. with min. chiSquared for ZZ (WW)";
  		
  		histExclYmerge12Text = jetAlgoText + "y value when merging from 1 to 2 jets (WW)"; 
  		histExclYmerge23Text = jetAlgoText + "y value when merging from 2 to 3 jets (WW)"; 
  		histExclYmerge34Text = jetAlgoText + "y value when merging from 3 to 4 jets (WW)"; 
  		histExclYmerge45Text = jetAlgoText + "y value when merging from 4 to 5 jets (WW)"; 
  		histExclYmerge56Text = jetAlgoText + "y value when merging from 5 to 6 jets (WW)"; 
  	} 
	if(topology == 9) 
	{
		histJetEtaText = jetAlgoText + "Jet Eta for events with ZH (4b)";
  		histJet1EtaText = jetAlgoText + "JetB1 Eta for events with ZH (4b)";
  		histJet2EtaText = jetAlgoText + "JetB2 Eta for events with ZH (4b)";
  		histJet3EtaText = jetAlgoText + "JetB3 Eta for events with ZH (4b)";
  		histJet4EtaText = jetAlgoText + "JetB4 Eta for events with ZH (4b)";
  		
  		histJetCosThetaText = jetAlgoText + "Jet cosTheta for events with ZH (4b)";
  		histJet1CosThetaText = jetAlgoText + "JetB1 cosTheta for events with ZH (4b)";
  		histJet2CosThetaText = jetAlgoText + "JetB2 cosTheta for events with ZH (4b)";
  		histJet3CosThetaText = jetAlgoText + "JetB3 cosTheta for events with ZH (4b)";
  		histJet4CosThetaText = jetAlgoText + "JetB4 cosTheta for events with ZH (4b)";
  		
  		histSumJetPtText = jetAlgoText + "Sum of jet Pt for events with ZH (4b)";
		histJetPtText = jetAlgoText + "Jet Pt for events with ZH (4b)";
		histJet4PtText = jetAlgoText + "Ind. jet Pt for events with ZH (4b)";
		
		histJetB1MText = jetAlgoText + "b-tagged jet-pair 1 mass (ZH)";
  		histJetB2MText = jetAlgoText + "b-tagged jet-pair 2 mass (ZH)";
  		histMinJetMText = jetAlgoText + "Mass of the jet-pair comb. with the least inv. mass (ZH)";
  		histJet1MText = jetAlgoText + "JetB1 mass (ZH)";
  		histJet2MText = jetAlgoText + "JetB2 mass (ZH)";
  		histJet3MText = jetAlgoText + "JetB3 mass (ZH)";
  		histJet4MText = jetAlgoText + "JetB4 mass (ZH)";
  		
  		histAplanText = jetAlgoText + "Aplanarity (ZH)";
  		histSpherText = jetAlgoText + "Sphericity (ZH)";
  		
  		histNParticlesText = jetAlgoText + "Number of particles (ZH)";
  		histTotalConstSizeText = jetAlgoText + "Number of constituents in the event (ZH)";
  		histConstSizeB1Text = jetAlgoText + "Number of constituents in jetB1 (ZH)";
  		histConstSizeB2Text = jetAlgoText + "Number of constituents in jetB2 (ZH)";
  		histConstSizeB3Text = jetAlgoText + "Number of constituents in jetB3 (ZH)";
  		histConstSizeB4Text = jetAlgoText + "Number of constituents in jetB4 (ZH)";
  		histMinConstSizeText = jetAlgoText + "Number of constituents in jet with min. number of consts. (ZH)";
  		
  		histNEFlowTracksText = jetAlgoText + "Number of eflow tracks (ZH)";
  		histNEFlowPhotonsText = jetAlgoText + "Number of eflow photons (ZH)";
  		histNEFlowNeutralHadronsText = jetAlgoText + "Number of eflow neutral hadrons (ZH)";
  		histNEFlowObjectsText = jetAlgoText + "Number of total eflow objects (i.e. tracks + photons + neutral hadrons) (ZH)";
  		
  		histJetB1NChargedText = jetAlgoText + "Number of charged objects in jetB1 (ZH)";
  		histJetB2NChargedText = jetAlgoText + "Number of charged objects in jetB2 (ZH)";
  		histJetB3NChargedText = jetAlgoText + "Number of charged objects in jetB3 (ZH)";
  		histJetB4NChargedText = jetAlgoText + "Number of charged objects in jetB4 (ZH)";
  		histJetB1NNeutralsText = jetAlgoText + "Number of neutral objects in jetB1 (ZH)";
  		histJetB2NNeutralsText = jetAlgoText + "Number of neutral objects in jetB2 (ZH)";
  		histJetB3NNeutralsText = jetAlgoText + "Number of neutral objects in jetB3 (ZH)";
  		histJetB4NNeutralsText = jetAlgoText + "Number of neutral objects in jetB4 (ZH)";
  		histJetNObjectsText = jetAlgoText + "Total number of eflow objects in the event (ZH)";
  		histMinJetNObjectsText = jetAlgoText + "Number of eflow objects in jet with min. number of eflow objects. (ZH)";
  		
  		histNJetsCompareAlgosText = jetAlgoText + "nJets in durham vs. antiKt (ZH)";
  		
  		histJetB1MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 2 antiKt (ZH)";
  		histJetB2MAntiKt2JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 2 antiKt (ZH)";
  		histJetB1MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 3 antiKt (ZH)";
  		histJetB2MAntiKt3JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 3 antiKt (ZH)";
  		histJetB1MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 4 antiKt (ZH)";
  		histJetB2MAntiKt4JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 4 antiKt (ZH)";
  		histJetB1MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 5 antiKt (ZH)";
  		histJetB2MAntiKt5JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 5 antiKt (ZH)";
  		histJetB1MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 1 mass -- 4 durham, 6 antiKt (ZH)";
  		histJetB2MAntiKt6JetsText = jetAlgoText + "b-tagged jet-pair 2 mass -- 4 durham, 6 antiKt (ZH)";
  		
  		histNJetsDurham0Text = jetAlgoText + "nJets in durham with rtdCut = 0 (ZH)"; 
  		histNJetsDurham5Text = jetAlgoText + "nJets in durham with rtdCut = 5 (ZH)"; 
  		histNJetsDurham10Text = jetAlgoText + "nJets in durham with rtdCut = 10 (ZH)"; 
  		histNJetsDurham15Text = jetAlgoText + "nJets in durham with rtdCut = 15 (ZH)"; 
  		histNJetsDurham20Text = jetAlgoText + "nJets in durham with rtdCut = 20 (ZH)"; 
  		histNJetsDurham25Text = jetAlgoText + "nJets in durham with rtdCut = 25 (ZH)"; 
  		histNJetsDurham30Text = jetAlgoText + "nJets in durham with rtdCut = 30 (ZH)";
  		
  		histMinChiSquaredZZMassText = jetAlgoText + "minChiSquared for the three possible jet pair comb. for ZZ mass (ZH)"; 
  		histInvMassZZ1Text =  jetAlgoText + "Mass of jetPair1 of comb. with min. chiSquared for ZZ (ZH)";
  		histInvMassZZ2Text =  jetAlgoText + "Mass of jetPair2 of comb. with min. chiSquared for ZZ (ZH)";
  		histDistanceZ1MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair1 of comb. with min. chiSquared for ZZ (ZH)";
  		histDistanceZ2MinChiSquaredZZMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair2 of comb. with min. chiSquared for ZZ (ZH)";

		histMinChiSquaredZHMassText = jetAlgoText + "minChiSquared for the three possible jet pair comb. for ZH mass (ZH)"; 
  		histInvMassZH1Text =  jetAlgoText + "Mass of jetPair1 of comb. with min. chiSquared for ZH (ZH)";
  		histInvMassZH2Text =  jetAlgoText + "Mass of jetPair2 of comb. with min. chiSquared for ZH (ZH)";
  		histDistanceZ1MinChiSquaredZHMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair1 of comb. with min. chiSquared for ZH (ZH)";
  		histDistanceZ2MinChiSquaredZHMassText = jetAlgoText + "Abs. delta between Z mass and mass of jetPair2 of comb. with min. chiSquared for ZH (ZH)";
  		
  		histExclYmerge12Text = jetAlgoText + "y value when merging from 1 to 2 jets (ZH)"; 
  		histExclYmerge23Text = jetAlgoText + "y value when merging from 2 to 3 jets (ZH)"; 
  		histExclYmerge34Text = jetAlgoText + "y value when merging from 3 to 4 jets (ZH)"; 
  		histExclYmerge45Text = jetAlgoText + "y value when merging from 4 to 5 jets (ZH)"; 
  		histExclYmerge56Text = jetAlgoText + "y value when merging from 5 to 6 jets (ZH)";

		histJetB1M1BestText = jetAlgoText + "b-tagged jet-pair 1 mass for best comb. for HH mass (ZH)";
  		histJetB2M1BestText = jetAlgoText + "b-tagged jet-pair 2 mass for best comb. for HH mass (ZH)"; 
  	} 

	if(topology == 12) 
	{
		histAplanText = jetAlgoText + "Aplanarity (qqH)";

		histJetB1MText = jetAlgoText + "b-tagged jet-pair 1 mass (qqH)";
  		histJetB2MText = jetAlgoText + "b-tagged jet-pair 2 mass (qqH)";

		histJetB1M1BestText = jetAlgoText + "b-tagged jet-pair 1 mass for best comb. for HH mass (qqH)";
  		histJetB2M1BestText = jetAlgoText + "b-tagged jet-pair 2 mass for best comb. for HH mass (qqH)"; 
		
		histInvMass4JetsText = jetAlgoText + "4-jet Inv. Mass (qqH)";
		histMinDeltaRJetsText = jetAlgoText + "min. deltaR between jets (qqH)";
		histDeltaRJetPairsText = jetAlgoText + "deltaR between jet-pairs (qqH)";
		histEAssymetryText = jetAlgoText + "E assymetry (qqH)";
		histPTAssymetryText = jetAlgoText + "PT assymetry (qqH)";

		histNSubJetsTrimmedB1Text = jetAlgoText + "nSubJetsTrimmed for jet B1 (qqH)";
		histNSubJetsTrimmedB2Text = jetAlgoText + "nSubJetsTrimmed for jet B2 (qqH)";
		histNSubJetsTrimmedB3Text = jetAlgoText + "nSubJetsTrimmed for jet B3 (qqH)";
		histNSubJetsTrimmedB4Text = jetAlgoText + "nSubJetsTrimmed for jet B4 (qqH)";
		histNSubJetsPrunnedB1Text = jetAlgoText + "nSubJetsPrunned for jet B1 (qqH)";
		histNSubJetsPrunnedB2Text = jetAlgoText + "nSubJetsPrunned for jet B2 (qqH)";
		histNSubJetsPrunnedB3Text = jetAlgoText + "nSubJetsPrunned for jet B3 (qqH)";
		histNSubJetsPrunnedB4Text = jetAlgoText + "nSubJetsPrunned for jet B4 (qqH)";
		histNSubJetsSoftDroppedB1Text = jetAlgoText + "nSubJetsSoftDropped for jet B1 (qqH)";
		histNSubJetsSoftDroppedB2Text = jetAlgoText + "nSubJetsSoftDropped for jet B2 (qqH)";
		histNSubJetsSoftDroppedB3Text = jetAlgoText + "nSubJetsSoftDropped for jet B3 (qqH)";
		histNSubJetsSoftDroppedB4Text = jetAlgoText + "nSubJetsSoftDropped for jet B4 (qqH)";

		histJetPairMassDeltaText = jetAlgoText + "JetPairMassDelta (qqH)";
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
  	
  	TH1 *histNJetsDurham0 = new TH1F("nJetsDurham0,", histNJetsDurham0Text.c_str(), 100.0, -1.0, 6);
  	TH1 *histNJetsDurham5 = new TH1F("nJetsDurham5,", histNJetsDurham5Text.c_str(), 100.0, -1.0, 6);
  	TH1 *histNJetsDurham10 = new TH1F("nJetsDurham10,", histNJetsDurham10Text.c_str(), 100.0, -1.0, 6);
  	TH1 *histNJetsDurham15 = new TH1F("nJetsDurham15,", histNJetsDurham15Text.c_str(), 100.0, -1.0, 6);
  	TH1 *histNJetsDurham20 = new TH1F("nJetsDurham20,", histNJetsDurham20Text.c_str(), 100.0, -1.0, 6);
  	TH1 *histNJetsDurham25 = new TH1F("nJetsDurham25,", histNJetsDurham25Text.c_str(), 100.0, -1.0, 6);
  	TH1 *histNJetsDurham30 = new TH1F("nJetsDurham30,", histNJetsDurham30Text.c_str(), 100.0, -1.0, 6);
  	
  	TH1 *histMinChiSquaredZZMass = new TH1F("MinChiSquaredZZMass,", histMinChiSquaredZZMassText.c_str(), 142.0, -100.0, 6000);
  	TH1 *histInvMassZZ1 = new TH1F("InvMassZZ1,", histInvMassZZ1Text.c_str(), 142.0, -1.0, 350);
  	TH1 *histInvMassZZ2 = new TH1F("InvMassZZ2,", histInvMassZZ2Text.c_str(), 142.0, -1.0, 350);
  	TH1 *histDistanceZ1MinChiSquaredZZMass = new TH1F("distanceZ1MinChiSquaredZZMass,", histDistanceZ1MinChiSquaredZZMassText.c_str(), 142.0, -1.0, 100);
  	TH1 *histDistanceZ2MinChiSquaredZZMass = new TH1F("distanceZ2MinChiSquaredZZMass,", histDistanceZ2MinChiSquaredZZMassText.c_str(), 142.0, -1.0, 100);

	TH1 *histMinChiSquaredZHMass = new TH1F("MinChiSquaredZHMass,", histMinChiSquaredZHMassText.c_str(), 142.0, -100.0, 6000);
  	TH1 *histInvMassZH1 = new TH1F("InvMassZH1,", histInvMassZH1Text.c_str(), 142.0, -1.0, 350);
  	TH1 *histInvMassZH2 = new TH1F("InvMassZH2,", histInvMassZH2Text.c_str(), 142.0, -1.0, 350);
  	TH1 *histDistanceZ1MinChiSquaredZHMass = new TH1F("distanceZ1MinChiSquaredZHMass,", histDistanceZ1MinChiSquaredZHMassText.c_str(), 142.0, -1.0, 100);
  	TH1 *histDistanceZ2MinChiSquaredZHMass = new TH1F("distanceZ2MinChiSquaredZHMass,", histDistanceZ2MinChiSquaredZHMassText.c_str(), 142.0, -1.0, 100);
  	
  	TH1 *histExclYmerge12 = new TH1F("ExclYmerge12,", histExclYmerge12Text.c_str(), 142.0, -0.1, 1);
  	TH1 *histExclYmerge23 = new TH1F("ExclYmerge23,", histExclYmerge23Text.c_str(), 142.0, -0.05, 0.35);
  	TH1 *histExclYmerge34 = new TH1F("ExclYmerge34,", histExclYmerge34Text.c_str(), 142.0, -0.05, 0.20);
  	TH1 *histExclYmerge45 = new TH1F("ExclYmerge45,", histExclYmerge45Text.c_str(), 142.0, -0.005, 0.04);
  	TH1 *histExclYmerge56 = new TH1F("ExclYmerge56,", histExclYmerge56Text.c_str(), 142.0, -0.001, 0.015);

	TH1 *histJetB1M1Best = new TH1F("jetb11Best,", histJetB1M1BestText.c_str(), 142.0, -1.0, 350);
  	TH1 *histJetB2M1Best = new TH1F("jetb21Best,", histJetB2M1BestText.c_str(), 142.0, -1.0, 350);

	TH2F *histNJets2D = new TH2F("histnjetsd2", "(rtd_min=10) b-tagged vs. non-b-taged jets (bbbb)", 100, -1.0, 7.0, 100, -1.0, 7.0);
	histNJets2D->GetXaxis()->SetTitle("b-tagged jets");
  	histNJets2D->GetYaxis()->SetTitle("non-b-tagged jets");

	TH1 *histDiHM = new TH1F("histDiHM,", histDiHMText.c_str(), 142.0, -1.0, 500);

	TH1 *histInvMassB1Fitep = new TH1F("jetb1pairfitep,", histInvMassB1FitepText.c_str(), 142.0, -10.0, 350);
  	TH1 *histInvMassB2Fitep = new TH1F("jetb2pairfitep,", histInvMassB2FitepText.c_str(), 142.0, -10.0, 350);
	TH1 *histChi2ndfep = new TH1F("histChi2ndfep,", histChi2ndfepText.c_str(), 142.0, -5.0, 50);
	TH1 *histInvMassB1Fitpxy = new TH1F("jetb1pairfitpxy,", histInvMassB1FitpxyText.c_str(), 142.0, -10.0, 350);
  	TH1 *histInvMassB2Fitpxy = new TH1F("jetb2pairfitpxy,", histInvMassB2FitpxyText.c_str(), 142.0, -10.0, 350);
	TH1 *histChi2ndfpxy = new TH1F("histChi2ndfpxy,", histChi2ndfpxyText.c_str(), 142.0, -5.0, 50);
	TH1 *histInvMassB1Fiteqm = new TH1F("jetb1pairfiteqm,", histInvMassB1FiteqmText.c_str(), 142.0, -10.0, 350);
  	TH1 *histInvMassB2Fiteqm = new TH1F("jetb2pairfiteqm,", histInvMassB2FiteqmText.c_str(), 142.0, -10.0, 350);
  	TH1 *histChi2ndfeqm = new TH1F("histChi2ndfeqm,", histChi2ndfeqmText.c_str(), 142.0, -5.0, 50);
  	TH1 *histInvMassB1FitBest = new TH1F("InvMassB1FitBest,", histInvMassB1FitBestText.c_str(), 142.0, -10.0, 350);
	TH1 *histInvMassB2FitBest = new TH1F("InvMassB2FitBest,", histInvMassB2FitBestText.c_str(), 142.0, -10.0, 350);
	TH1 *histChi2ndfBest = new TH1F("Chi2ndfBest,", histChi2ndfBestText.c_str(), 142.0, -5.0, 50);

	TH1 *histMinJetChiS = new TH1F("histMinJetChiS,", histMinJetChiSText.c_str(), 142.0, -5.0, 50);

	TH1 *histInvMass4Jets = new TH1F("histInvMass4Jets,", histInvMass4JetsText.c_str(), 142.0, -10.0, 420);
  	TH1 *histMinDeltaRJets = new TH1F("histMinDeltaRJets,", histMinDeltaRJetsText.c_str(), 142.0, -0.1, 3.14);
  	TH1 *histDeltaRJetPairs = new TH1F("histDeltaRJetPairs,", histDeltaRJetPairsText.c_str(), 142.0, -0.1, 10);
  	TH1 *histEAssymetry = new TH1F("histEAssymetry,", histEAssymetryText.c_str(), 142.0, -0.1, 1.1);
  	TH1 *histPTAssymetry = new TH1F("histPTAssymetry,", histPTAssymetryText.c_str(), 142.0, -0.1, 1.1);

	TH1 *histNSubJetsTrimmedB1 = new TH1F("histNSubJetsTrimmedB1,", histNSubJetsTrimmedB1Text.c_str(), 142.0, -1, 10);
	TH1 *histNSubJetsTrimmedB2 = new TH1F("histNSubJetsTrimmedB2,", histNSubJetsTrimmedB2Text.c_str(), 142.0, -1, 10);
	TH1 *histNSubJetsTrimmedB3 = new TH1F("histNSubJetsTrimmedB3,", histNSubJetsTrimmedB3Text.c_str(), 142.0, -1, 10);
	TH1 *histNSubJetsTrimmedB4 = new TH1F("histNSubJetsTrimmedB4,", histNSubJetsTrimmedB4Text.c_str(), 142.0, -1, 10);
  	TH1 *histNSubJetsPrunnedB1 = new TH1F("histNSubJetsPrunnedB1,", histNSubJetsPrunnedB1Text.c_str(), 142.0, -1, 10);
  	TH1 *histNSubJetsPrunnedB2 = new TH1F("histNSubJetsPrunnedB2,", histNSubJetsPrunnedB2Text.c_str(), 142.0, -1, 10);
  	TH1 *histNSubJetsPrunnedB3 = new TH1F("histNSubJetsPrunnedB3,", histNSubJetsPrunnedB3Text.c_str(), 142.0, -1, 10);
  	TH1 *histNSubJetsPrunnedB4 = new TH1F("histNSubJetsPrunnedB4,", histNSubJetsPrunnedB4Text.c_str(), 142.0, -1, 10);
  	TH1 *histNSubJetsSoftDroppedB1 = new TH1F("histNSubJetsSoftDroppedB1,", histNSubJetsSoftDroppedB1Text.c_str(), 142.0, -1, 10);
  	TH1 *histNSubJetsSoftDroppedB2 = new TH1F("histNSubJetsSoftDroppedB2,", histNSubJetsSoftDroppedB2Text.c_str(), 142.0, -1, 10);
  	TH1 *histNSubJetsSoftDroppedB3 = new TH1F("histNSubJetsSoftDroppedB3,", histNSubJetsSoftDroppedB3Text.c_str(), 142.0, -1, 10);
  	TH1 *histNSubJetsSoftDroppedB4 = new TH1F("histNSubJetsSoftDroppedB4,", histNSubJetsSoftDroppedB4Text.c_str(), 142.0, -1, 10);

  	TH1 *histJetPairMassDelta = new TH1F("histJetPairMassDelta,", histJetPairMassDeltaText.c_str(), 142.0, -100.0, 4000);
  	
  	TClonesArray *branchParticle;
  	TClonesArray *branchEvent;
  	TClonesArray *branchJetAntiKt;
  	TClonesArray *branchGenJetAntiKt;
  	TClonesArray *branchElectron;
  	TClonesArray *branchMuon;
  	TClonesArray *branchJetDurham;
  	TClonesArray *branchGenJetDurham;
  	TClonesArray *branchElectronDurham;
  	TClonesArray *branchMuonDurham;
  	TClonesArray *branchEFlowTrack;
  	TClonesArray *branchEFlowPhoton;
  	TClonesArray *branchEFlowNeutralHadron;
  	TClonesArray *branchJetDurham0;
  	TClonesArray *branchJetDurham5;
  	TClonesArray *branchJetDurham10;
  	TClonesArray *branchJetDurham15;
  	TClonesArray *branchJetDurham20;
  	TClonesArray *branchJetDurham25;
  	TClonesArray *branchJetDurham30;
	TClonesArray *branchMissingET;
  	
  	
	branchParticle = treeReader->UseBranch("Particle");
	branchEvent = treeReader->UseBranch("Event");
	branchJetAntiKt = treeReader->UseBranch("JetAntiKt");
	//branchGenJetAntiKt = treeReader->UseBranch("GenJet");
	branchElectron = treeReader->UseBranch("Electron");
	branchMuon = treeReader->UseBranch("Muon");
	branchJetDurham = treeReader->UseBranch(jetAlgo.c_str());
	branchGenJetDurham = treeReader->UseBranch(genJetAlgo.c_str());
	//branchElectronDurham = treeReader->UseBranch("ElectronDurham");
	//branchMuonDurham = treeReader->UseBranch("MuonDurham");
	branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
	branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
	branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
	branchJetDurham0 = treeReader->UseBranch("Jet0");
	branchJetDurham5 = treeReader->UseBranch("Jet5");
	branchJetDurham10 = treeReader->UseBranch("Jet10");
	branchJetDurham15 = treeReader->UseBranch("Jet15");
	branchJetDurham20 = treeReader->UseBranch("Jet20");
	branchJetDurham25 = treeReader->UseBranch("Jet25");
	branchJetDurham30 = treeReader->UseBranch("Jet30");
	branchMissingET = treeReader->UseBranch("MissingET");
  	
	
	int trueFullTopEvents=0, trueSemiTopEvents=0, trueFullZEvents=0, trueSemiZEvents=0, trueFullZHEvents=0, trueSemiZHEvents=0, trueFullWEvents=0, trueSemiWWEvents=0, trueFullWWEvents=0, trueSemiWEvents=0, trueFullLeptonEvents=0, trueSemiLeptonEvents=0, trueFullHadronEvents=0, trueSemiHadronEvents=0, trueFullNeutrinoEvents=0, trueSemiNeutrinoEvents=0, weirdDecays=0, contEtaCut=0, contCosThetaCut=0, trueLeptonEvents=0, trueHadronEvents=0, trueNeutrinoEvents=0;
	int trueLeptonNeutrinoEvents=0, trueLeptonHadronEvents=0, trueHadronNeutrinoEvents=0; 
	float pt1=0, pt2=0, pt3=0, pt4=0, ptSum=0, contMassB1Window=0, contMassB2Window=0, windowBLeft=0, windowBRight=0, windowB2Left=0, windowB2Right=0, contEventsWindow=0, contEvents=0, contEventsPostFilter=0, contEvents4Jets=0;
	double distanceZMass;
	double contLeptonEventsZH=0;
	int contEventsZWindow10=0, contEventsZWindow01=0, contEventsZWindow05=0, contEventsZWindow1=0, contEventsZWindow15=0, contEventsZWindow2=0, contEventsZWindow5=0;
	
	float aplanarity, invMassB1, invMassB2, minJetM, sphericity, cosThetaB1, cosThetaB2, cosThetaB3, cosThetaB4, sumPt, jetB1Pt, jetB2Pt, jetB3Pt, jetB4Pt, jetB1M, jetB2M, jetB3M, jetB4M, etaB1, etaB2, etaB3, etaB4, nParticles, totalConstSize, constSizeB1, constSizeB2, constSizeB3, constSizeB4, minConstSize, jetB1NCharged, jetB2NCharged, jetB3NCharged, jetB4NCharged, jetB1NNeutrals, jetB2NNeutrals, jetB3NNeutrals, jetB4NNeutrals, jetNObjects, minJetNObjects, invMassB1AntiKt, invMassB2AntiKt, invMassB1AntiKt2Jets, invMassB2AntiKt2Jets, invMassB1AntiKt3Jets, invMassB2AntiKt3Jets, invMassB1AntiKt4Jets, invMassB2AntiKt4Jets, invMassB1AntiKt5Jets, invMassB2AntiKt5Jets, invMassB1AntiKt6Jets, invMassB2AntiKt6Jets, nJetsAntiKt, invMassB11Best, invMassB21Best, invMassB12Best, invMassB22Best, invMassB13Best, invMassB23Best, invMassB14Best, invMassB24Best, invMassB15Best, invMassB25Best, invMassB16Best, invMassB26Best, invMassB17Best, invMassB27Best, invMassB18Best, invMassB28Best, entryIndex, nJetsDurham0, nJetsDurham5, nJetsDurham10, nJetsDurham15, nJetsDurham20, nJetsDurham25, nJetsDurham30, minChiSquaredZZMass, distanceZ1MinChiSquaredZZMass, distanceZ2MinChiSquaredZZMass, exclYmerge12, exclYmerge23, exclYmerge34, exclYmerge45, exclYmerge56, invMassZZ1, invMassZZ2, thrust, boostB1, boostB2, boostB3, boostB4, boostSystem, missingET;
	float minChiSquaredZHMass, distanceZ1MinChiSquaredZHMass, distanceZ2MinChiSquaredZHMass, invMassZHZ, invMassZHH;
	int contLargeBMass=0, contSmallBMass=0;
	double minJetChiS=999, minJetChiSep=999, minJetChiSpxy=999, minJetChiSeqm=999, minJetChiSBest=999;
	float invMass4Jets=0, minDeltaRJets=0, deltaRJetPairs=0, eAssymetry=0, pTAssymetry=0, jetPairMassDelta=0, nSubJetsTrimmedB1=0, nSubJetsTrimmedB2=0, nSubJetsTrimmedB3=0, nSubJetsTrimmedB4=0, nSubJetsPrunnedB1=0,  nSubJetsPrunnedB2=0,  nSubJetsPrunnedB3=0,  nSubJetsPrunnedB4=0, nSubJetsSoftDroppedB1=0, nSubJetsSoftDroppedB2=0, nSubJetsSoftDroppedB3=0, nSubJetsSoftDroppedB4=0;
	
	////kinematic fit
	float chi2ndfep, chi2ndfpxy, chi2ndfeqm, chi2ndfeqmback, chi2ndfBest;
	TLorentzVector jetB1Fitep, jetB2Fitep, jetPairB1Fitep, jetB3Fitep, jetB4Fitep, jetPairB2Fitep;
	TLorentzVector jetB1Fitpxy, jetB2Fitpxy, jetPairB1Fitpxy, jetB3Fitpxy, jetB4Fitpxy, jetPairB2Fitpxy;
	TLorentzVector jetB1Fiteqm, jetB2Fiteqm, jetPairB1Fiteqm, jetB3Fiteqm, jetB4Fiteqm, jetPairB2Fiteqm;
	TLorentzVector jetB1Fitini, jetB2Fitini, jetPairB1Fitini, jetB3Fitini, jetB4Fitini, jetPairB2Fitini;
	float invMassB1Fitep, invMassB2Fitep, invMassB1Fitpxy, invMassB2Fitpxy, invMassB1Fiteqm, invMassB2Fiteqm, invMassB1FitBest, invMassB2FitBest;
	float leadingJetPairChiSquaredReco, leadingJetPairChiSquaredTrue, nonLeadingJetPairChiSquaredReco, nonLeadingJetPairChiSquaredTrue, jetPairsChiSquaredReco, jetPairsChiSquaredTrue;
	float jetPairMassDeltaFit;
	//////kinematic fit
	
	if(fileFunction == "generate")
	{
		/////Filling out tree for BDT
		TreeTrain.Branch("entryIndex",&entryIndex,"entryIndex/F");
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
		//TreeTrain.Branch("nJetsDurham0",&nJetsDurham0,"nJetsDurham0/F");
		//TreeTrain.Branch("nJetsDurham5",&nJetsDurham5,"nJetsDurham5/F");
		//TreeTrain.Branch("nJetsDurham10",&nJetsDurham10,"nJetsDurham10/F");
		TreeTrain.Branch("nJetsDurham15",&nJetsDurham15,"nJetsDurham15/F");
		TreeTrain.Branch("nJetsDurham20",&nJetsDurham20,"nJetsDurham20/F");
		TreeTrain.Branch("nJetsDurham25",&nJetsDurham25,"nJetsDurham25/F");
		TreeTrain.Branch("nJetsDurham30",&nJetsDurham30,"nJetsDurham30/F");
		TreeTrain.Branch("distanceZ1MinChiSquaredZZMass",&distanceZ1MinChiSquaredZZMass,"distanceZ1MinChiSquaredZZMass/F");
		TreeTrain.Branch("distanceZ2MinChiSquaredZZMass",&distanceZ2MinChiSquaredZZMass,"distanceZ2MinChiSquaredZZMass/F");
		TreeTrain.Branch("exclYmerge12",&exclYmerge12,"exclYmerge12/F");
		TreeTrain.Branch("exclYmerge23",&exclYmerge23,"exclYmerge23/F");
		TreeTrain.Branch("exclYmerge34",&exclYmerge34,"exclYmerge34/F");
		TreeTrain.Branch("exclYmerge45",&exclYmerge45,"exclYmerge45/F");
		TreeTrain.Branch("exclYmerge56",&exclYmerge56,"exclYmerge56/F");
		TreeTrain.Branch("invMassZZ1",&invMassZZ1,"invMassZZ1/F");
		TreeTrain.Branch("invMassZZ2",&invMassZZ2,"invMassZZ2/F");
		TreeTrain.Branch("thrust",&thrust,"thrust/F");
		TreeTrain.Branch("boostB1",&boostB1,"boostB1/F");
		TreeTrain.Branch("boostB2",&boostB2,"boostB2/F");
		TreeTrain.Branch("boostB3",&boostB3,"boostB3/F");
		TreeTrain.Branch("boostB4",&boostB4,"boostB4/F");
		TreeTrain.Branch("boostSystem",&boostSystem,"boostSystem/F");
		TreeTrain.Branch("missingET",&missingET,"missingET/F");
		TreeTrain.Branch("invMass4Jets",&invMass4Jets,"invMass4Jets/F");
		TreeTrain.Branch("deltaRJetPairs",&deltaRJetPairs,"deltaRJetPairs/F");
		TreeTrain.Branch("eAssymetry",&eAssymetry,"eAssymetry/F");
		TreeTrain.Branch("pTAssymetry",&pTAssymetry,"pTAssymetry/F");
		TreeTrain.Branch("jetPairMassDelta",&jetPairMassDelta,"jetPairMassDelta/F");
		/////kinematic fit
		TreeTrain.Branch("invMassB1FitBest",&invMassB1FitBest,"invMassB1FitBest/F");
		TreeTrain.Branch("invMassB2FitBest",&invMassB2FitBest,"invMassB2FitBest/F");
		TreeTrain.Branch("chi2ndfBest",&chi2ndfBest,"chi2ndfBest/F");
		TreeTrain.Branch("jetPairMassDeltaFit",&jetPairMassDeltaFit,"jetPairMassDeltaFit/F");
		///////kinematic fit

		TreeTest.Branch("entryIndex",&entryIndex,"entryIndex/F");
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
		//TreeTest.Branch("nJetsDurham0",&nJetsDurham0,"nJetsDurham0/F");
		//TreeTest.Branch("nJetsDurham5",&nJetsDurham5,"nJetsDurham5/F");
		//TreeTest.Branch("nJetsDurham10",&nJetsDurham10,"nJetsDurham10/F");
		TreeTest.Branch("nJetsDurham15",&nJetsDurham15,"nJetsDurham15/F");
		TreeTest.Branch("nJetsDurham20",&nJetsDurham20,"nJetsDurham20/F");
		TreeTest.Branch("nJetsDurham25",&nJetsDurham25,"nJetsDurham25/F");
		TreeTest.Branch("nJetsDurham30",&nJetsDurham30,"nJetsDurham30/F");
		TreeTest.Branch("distanceZ1MinChiSquaredZZMass",&distanceZ1MinChiSquaredZZMass,"distanceZ1MinChiSquaredZZMass/F");
		TreeTest.Branch("distanceZ2MinChiSquaredZZMass",&distanceZ2MinChiSquaredZZMass,"distanceZ2MinChiSquaredZZMass/F");
		TreeTest.Branch("exclYmerge12",&exclYmerge12,"exclYmerge12/F");
		TreeTest.Branch("exclYmerge23",&exclYmerge23,"exclYmerge23/F");
		TreeTest.Branch("exclYmerge34",&exclYmerge34,"exclYmerge34/F");
		TreeTest.Branch("exclYmerge45",&exclYmerge45,"exclYmerge45/F");
		TreeTest.Branch("exclYmerge56",&exclYmerge56,"exclYmerge56/F");
		TreeTest.Branch("invMassZZ1",&invMassZZ1,"invMassZZ1/F");
		TreeTest.Branch("invMassZZ2",&invMassZZ2,"invMassZZ2/F");
		TreeTest.Branch("thrust",&thrust,"thrust/F");
		TreeTest.Branch("boostB1",&boostB1,"boostB1/F");
		TreeTest.Branch("boostB2",&boostB2,"boostB2/F");
		TreeTest.Branch("boostB3",&boostB3,"boostB3/F");
		TreeTest.Branch("boostB4",&boostB4,"boostB4/F");
		TreeTest.Branch("boostSystem",&boostSystem,"boostSystem/F");
		TreeTest.Branch("missingET",&missingET,"missingET/F");
		TreeTest.Branch("invMass4Jets",&invMass4Jets,"invMass4Jets/F");
		TreeTest.Branch("deltaRJetPairs",&deltaRJetPairs,"deltaRJetPairs/F");
		TreeTest.Branch("eAssymetry",&eAssymetry,"eAssymetry/F");
		TreeTest.Branch("pTAssymetry",&pTAssymetry,"pTAssymetry/F");
		TreeTest.Branch("jetPairMassDelta",&jetPairMassDelta,"jetPairMassDelta/F");
		/////kinematic fit
		TreeTest.Branch("invMassB1FitBest",&invMassB1FitBest,"invMassB1FitBest/F");
		TreeTest.Branch("invMassB2FitBest",&invMassB2FitBest,"invMassB2FitBest/F");
		TreeTest.Branch("chi2ndfBest",&chi2ndfBest,"chi2ndfBest/F");
		TreeTest.Branch("jetPairMassDeltaFit",&jetPairMassDeltaFit,"jetPairMassDeltaFit/F");
		///////kinematic fit

		TreeFull.Branch("entryIndex",&entryIndex,"entryIndex/F");
		TreeFull.Branch("aplanarity",&aplanarity,"aplanarity/F");
		TreeFull.Branch("invMassB1",&invMassB1,"invMassB1/F");
		TreeFull.Branch("invMassB2",&invMassB2,"invMassB2/F");
		TreeFull.Branch("minJetM",&minJetM,"minJetM/F");  
		TreeFull.Branch("sphericity",&sphericity,"sphericity/F");
		TreeFull.Branch("cosThetaB1",&cosThetaB1,"cosThetaB1/F");
		TreeFull.Branch("cosThetaB2",&cosThetaB2,"cosThetaB2/F");
		TreeFull.Branch("cosThetaB3",&cosThetaB3,"cosThetaB3/F");
		TreeFull.Branch("cosThetaB4",&cosThetaB4,"cosThetaB4/F");
		TreeFull.Branch("sumPt",&sumPt,"sumPt/F");
		TreeFull.Branch("jetB1Pt",&jetB1Pt,"jetB1Pt/F");
		TreeFull.Branch("jetB2Pt",&jetB2Pt,"jetB2Pt/F");
		TreeFull.Branch("jetB3Pt",&jetB3Pt,"jetB3Pt/F");
		TreeFull.Branch("jetB4Pt",&jetB4Pt,"jetB4Pt/F");
		TreeFull.Branch("jetB1M",&jetB1M,"jetB1M/F");
		TreeFull.Branch("jetB2M",&jetB2M,"jetB2M/F");
		TreeFull.Branch("jetB3M",&jetB3M,"jetB3M/F");
		TreeFull.Branch("jetB4M",&jetB4M,"jetB4M/F");
		TreeFull.Branch("nParticles",&nParticles,"nParticles/F");
		TreeFull.Branch("constSizeB1",&constSizeB1,"constSizeB1/F");
		TreeFull.Branch("constSizeB2",&constSizeB2,"constSizeB2/F");
		TreeFull.Branch("constSizeB3",&constSizeB3,"constSizeB3/F");
		TreeFull.Branch("constSizeB4",&constSizeB4,"constSizeB4/F");
		TreeFull.Branch("minConstSize",&minConstSize,"minConstSize/F");
		TreeFull.Branch("jetNObjects",&jetNObjects,"jetNObjects/F");
		TreeFull.Branch("minJetNObjects",&minJetNObjects,"minJetNObjects/F");
		TreeFull.Branch("invMassB1AntiKt",&invMassB1AntiKt,"invMassB1AntiKt/F");
		TreeFull.Branch("invMassB2AntiKt",&invMassB2AntiKt,"invMassB2AntiKt/F");
		TreeFull.Branch("nJetsAntiKt",&nJetsAntiKt,"nJetsAntiKt/F");
		TreeFull.Branch("invMassB11Best",&invMassB11Best,"invMassB11Best/F");
		TreeFull.Branch("invMassB21Best",&invMassB21Best,"invMassB21Best/F");
		TreeFull.Branch("invMassB12Best",&invMassB12Best,"invMassB12Best/F");
		TreeFull.Branch("invMassB22Best",&invMassB22Best,"invMassB22Best/F");
		TreeFull.Branch("invMassB13Best",&invMassB13Best,"invMassB13Best/F");
		TreeFull.Branch("invMassB23Best",&invMassB23Best,"invMassB23Best/F");
		TreeFull.Branch("invMassB14Best",&invMassB14Best,"invMassB14Best/F");
		TreeFull.Branch("invMassB24Best",&invMassB24Best,"invMassB24Best/F");
		TreeFull.Branch("invMassB15Best",&invMassB15Best,"invMassB15Best/F");
		TreeFull.Branch("invMassB25Best",&invMassB25Best,"invMassB25Best/F");
		TreeFull.Branch("invMassB16Best",&invMassB16Best,"invMassB16Best/F");
		TreeFull.Branch("invMassB26Best",&invMassB26Best,"invMassB26Best/F");
		TreeFull.Branch("invMassB17Best",&invMassB17Best,"invMassB17Best/F");
		TreeFull.Branch("invMassB27Best",&invMassB27Best,"invMassB27Best/F");
		TreeFull.Branch("invMassB18Best",&invMassB18Best,"invMassB18Best/F");
		TreeFull.Branch("invMassB28Best",&invMassB28Best,"invMassB28Best/F");
		//TreeFull.Branch("nJetsDurham0",&nJetsDurham0,"nJetsDurham0/F");
		//TreeFull.Branch("nJetsDurham5",&nJetsDurham5,"nJetsDurham5/F");
		//TreeFull.Branch("nJetsDurham10",&nJetsDurham10,"nJetsDurham10/F");
		TreeFull.Branch("nJetsDurham15",&nJetsDurham15,"nJetsDurham15/F");
		TreeFull.Branch("nJetsDurham20",&nJetsDurham20,"nJetsDurham20/F");
		TreeFull.Branch("nJetsDurham25",&nJetsDurham25,"nJetsDurham25/F");
		TreeFull.Branch("nJetsDurham30",&nJetsDurham30,"nJetsDurham30/F");
		TreeFull.Branch("distanceZ1MinChiSquaredZZMass",&distanceZ1MinChiSquaredZZMass,"distanceZ1MinChiSquaredZZMass/F");
		TreeFull.Branch("distanceZ2MinChiSquaredZZMass",&distanceZ2MinChiSquaredZZMass,"distanceZ2MinChiSquaredZZMass/F");
		TreeFull.Branch("exclYmerge12",&exclYmerge12,"exclYmerge12/F");
		TreeFull.Branch("exclYmerge23",&exclYmerge23,"exclYmerge23/F");
		TreeFull.Branch("exclYmerge34",&exclYmerge34,"exclYmerge34/F");
		TreeFull.Branch("exclYmerge45",&exclYmerge45,"exclYmerge45/F");
		TreeFull.Branch("exclYmerge56",&exclYmerge56,"exclYmerge56/F");
		TreeFull.Branch("invMassZZ1",&invMassZZ1,"invMassZZ1/F");
		TreeFull.Branch("invMassZZ2",&invMassZZ2,"invMassZZ2/F");
		TreeFull.Branch("thrust",&thrust,"thrust/F");
		TreeFull.Branch("boostB1",&boostB1,"boostB1/F");
		TreeFull.Branch("boostB2",&boostB2,"boostB2/F");
		TreeFull.Branch("boostB3",&boostB3,"boostB3/F");
		TreeFull.Branch("boostB4",&boostB4,"boostB4/F");
		TreeFull.Branch("boostSystem",&boostSystem,"boostSystem/F");
		TreeFull.Branch("missingET",&missingET,"missingET/F");
		TreeFull.Branch("invMass4Jets",&invMass4Jets,"invMass4Jets/F");
		TreeFull.Branch("deltaRJetPairs",&deltaRJetPairs,"deltaRJetPairs/F");
		TreeFull.Branch("eAssymetry",&eAssymetry,"eAssymetry/F");
		TreeFull.Branch("pTAssymetry",&pTAssymetry,"pTAssymetry/F");
		TreeFull.Branch("jetPairMassDelta",&jetPairMassDelta,"jetPairMassDelta/F");
		/////kinematic fit
		TreeFull.Branch("invMassB1FitBest",&invMassB1FitBest,"invMassB1FitBest/F");
		TreeFull.Branch("invMassB2FitBest",&invMassB2FitBest,"invMassB2FitBest/F");
		TreeFull.Branch("chi2ndfBest",&chi2ndfBest,"chi2ndfBest/F");
		TreeFull.Branch("jetPairMassDeltaFit",&jetPairMassDeltaFit,"jetPairMassDeltaFit/F");
		///////kinematic fit
	}
	
	TreeMerge.Branch("entryIndex",&entryIndex,"entryIndex/F");	
	//TreeMerge.Branch("nJetsDurham0",&nJetsDurham0,"nJetsDurham0/F");
	//TreeMerge.Branch("nJetsDurham5",&nJetsDurham5,"nJetsDurham5/F");
	//TreeMerge.Branch("nJetsDurham10",&nJetsDurham10,"nJetsDurham10/F");
	TreeMerge.Branch("nJetsDurham15",&nJetsDurham15,"nJetsDurham15/F");
	TreeMerge.Branch("nJetsDurham20",&nJetsDurham20,"nJetsDurham20/F");
	TreeMerge.Branch("nJetsDurham25",&nJetsDurham25,"nJetsDurham25/F");
	TreeMerge.Branch("nJetsDurham30",&nJetsDurham30,"nJetsDurham30/F");
	TreeMerge.Branch("distanceZ1MinChiSquaredZZMass",&distanceZ1MinChiSquaredZZMass,"distanceZ1MinChiSquaredZZMass/F");
	TreeMerge.Branch("distanceZ2MinChiSquaredZZMass",&distanceZ2MinChiSquaredZZMass,"distanceZ2MinChiSquaredZZMass/F");
	TreeMerge.Branch("exclYmerge12",&exclYmerge12,"exclYmerge12/F");
	TreeMerge.Branch("exclYmerge23",&exclYmerge23,"exclYmerge23/F");
	TreeMerge.Branch("exclYmerge34",&exclYmerge34,"exclYmerge34/F");
	TreeMerge.Branch("exclYmerge45",&exclYmerge45,"exclYmerge45/F");
	TreeMerge.Branch("exclYmerge56",&exclYmerge56,"exclYmerge56/F");

	//////kinematic fit
	double Ecm=380.;
	bool enableExtraTries=false;  //   if true then fit is performed several times with different intial betaX, betaY, betaZ valeus
	double nSigVar=3.;   // if enableExtraTries=true, controls spread in initial beta values w.r.t. values given by jetB1,jetB2,jetNB1,jetNB2

	//////pythia errors
	vector<double> BfractionalEnergyOffset;
	BfractionalEnergyOffset.push_back(-0.043);
	BfractionalEnergyOffset.push_back(-0.141);
	BfractionalEnergyOffset.push_back(-0.158);
	BfractionalEnergyOffset.push_back(-0.255);
	
	vector<double> BfractionalEnergyErrors;
	BfractionalEnergyErrors.push_back(0.117);
	BfractionalEnergyErrors.push_back(0.244);
	BfractionalEnergyErrors.push_back(0.280);
	BfractionalEnergyErrors.push_back(0.359);

	vector<double> BbetaXErrors;
	BbetaXErrors.push_back(0.00658);
	BbetaXErrors.push_back(0.00811);
	BbetaXErrors.push_back(0.01320);
	BbetaXErrors.push_back(0.02061);

	vector<double> BbetaYErrors;
	BbetaYErrors.push_back(0.00658);
	BbetaYErrors.push_back(0.00838);
	BbetaYErrors.push_back(0.01297);
	BbetaYErrors.push_back(0.02037);

	vector<double> BbetaZErrors;
	BbetaZErrors.push_back(0.00532);
	BbetaZErrors.push_back(0.00740);
	BbetaZErrors.push_back(0.01124);
	BbetaZErrors.push_back(0.01812);
	//////pythia errors

	/*double BfractionalEnergyError1=0.0990;
	double BfractionalEnergyError2=0.1329;
	double BfractionalEnergyError3=0.1617;
	double BfractionalEnergyError4=0.1662;*/
	/*double BfractionalEnergyError1=0.08244;
	double BfractionalEnergyError2=0.1015;
	double BfractionalEnergyError3=0.114;
	double BfractionalEnergyError4=0.1208;*/
	/*double BfractionalEnergyError2=0.1;
	double BfractionalEnergyError3=0.15;
	double BfractionalEnergyError4=0.20;
	vector<double> BfractionalEnergyErrors = {BfractionalEnergyError1, BfractionalEnergyError2, BfractionalEnergyError3, BfractionalEnergyError4};*/
	//double BbetaError=0.001;
	/*double BbetaXError1=0.1293;
	double BbetaXError2=0.1721;
	double BbetaXError3=0.1972;
	double BbetaXError4=0.2072;*/
	/*double BbetaXError1=0.04055;
	double BbetaXError2=0.04525;
	double BbetaXError3=0.0532;
	double BbetaXError4=0.06258;*/
	/*double BbetaXError1=0.015;
	double BbetaXError2=0.015;
	double BbetaXError3=0.015;
	double BbetaXError4=0.015;*/
	/*double BbetaXError1=0.001;
	double BbetaXError2=0.001;
	double BbetaXError3=0.001;
	double BbetaXError4=0.001;
	vector<double> BbetaXErrors = {BbetaXError1, BbetaXError2, BbetaXError3, BbetaXError4};*/
	/*double BbetaYError1=0.1286;
	double BbetaYError2=0.1725;
	double BbetaYError3=0.1974;
	double BbetaYError4=0.2063;*/
	/*double BbetaYError1=0.04048;
	double BbetaYError2=0.04518;
	double BbetaYError3=0.05313;
	double BbetaYError4=0.06216;*/
	/*double BbetaYError1=0.015;
	double BbetaYError2=0.015;
	double BbetaYError3=0.015;
	double BbetaYError4=0.015;*/
	/*double BbetaYError1=0.001;
	double BbetaYError2=0.001;
	double BbetaYError3=0.001;
	double BbetaYError4=0.001;
	vector<double> BbetaYErrors = {BbetaYError1, BbetaYError2, BbetaYError3, BbetaYError4};*/
	/*double BbetaZError1=0.1663;
	double BbetaZError2=0.1925;
	double BbetaZError3=0.1958;
	double BbetaZError4=0.1715;*/
	/*double BbetaZError1=0.04656;
	double BbetaZError2=0.04733;
	double BbetaZError3=0.05042;
	double BbetaZError4=0.04766;*/
	/*double BbetaZError1=0.015;
	double BbetaZError2=0.015;
	double BbetaZError3=0.015;
	double BbetaZError4=0.015;*/
	/*double BbetaZError1=0.001;
	double BbetaZError2=0.001;
	double BbetaZError3=0.001;
	double BbetaZError4=0.001;
	vector<double> BbetaZErrors = {BbetaZError1, BbetaZError2, BbetaZError3, BbetaZError4};*/

	//cout << " BfractionalEnergyError= " << BfractionalEnergyError << " BbetaError= " << BbetaError << endl;
	
	bool testFittedMass=false;
	double massMaxB=134e10;
	double massMinB=-120;
	
	cout << " testFittedMass= " <<  testFittedMass << " massMaxB= " << massMaxB << " massMinB= " << massMinB << endl;
				
	//  double scaleMaxEntries=0.1;
	//  cout << " scaleMaxEntries= " << scaleMaxEntries << endl;

	//  maxEntries *= scaleMaxEntries;
	/////kinematic fit
	
	for(int entry=0; entry<numberOfEntries; entry++) 
	{
		treeReader->ReadEntry(entry);
		entryIndex=entry;

		if(entry % 250000 == 0) cout<<"Entry: "<<entry<<endl;
		
		nJetsAntiKt=0;
		if(branchJetAntiKt->GetEntries() > 0)  nJetsAntiKt =  branchJetAntiKt->GetEntries();
		int nJetsDurham=0;
		if(branchJetDurham->GetEntries() > 0)  nJetsDurham =  branchJetDurham->GetEntries();
		nJetsDurham0=0;
		if(branchJetDurham0->GetEntries() > 0)  nJetsDurham0 =  branchJetDurham0->GetEntries();
		nJetsDurham5=0;
		if(branchJetDurham5->GetEntries() > 0)  nJetsDurham5 =  branchJetDurham5->GetEntries();
		nJetsDurham10=0;
		if(branchJetDurham10->GetEntries() > 0)  nJetsDurham10 =  branchJetDurham10->GetEntries();
		nJetsDurham15=0;
		if(branchJetDurham15->GetEntries() > 0)  nJetsDurham15 =  branchJetDurham15->GetEntries();
		nJetsDurham20=0;
		if(branchJetDurham20->GetEntries() > 0)  nJetsDurham20 =  branchJetDurham20->GetEntries();
		nJetsDurham25=0;
		if(branchJetDurham25->GetEntries() > 0)  nJetsDurham25 =  branchJetDurham25->GetEntries();
		nJetsDurham30=0;
		if(branchJetDurham30->GetEntries() > 0)  nJetsDurham30 =  branchJetDurham30->GetEntries();
		
		//if(nJetsAntiKt == 4) 
		if(nJetsDurham == 4)
		//if(1==1)
		{
			contEvents4Jets++;
			/*random_device rd;  // Will be used to obtain a seed for the random number engine
			mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
			uniform_int_distribution<> distrib(1, 10); // Uniform distribution between 1 and 10
			if(distrib(gen) != 1) continue;*/
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
			
				bool flagPreselection=false;
	      		if(preselection == "4BSplit")
	      		{
	      			if(contBJetsDurham == 4) flagPreselection=true;
	      		}
	      		else if(preselection == "34BSplit")
	      		{
	      			if(contBJetsDurham == 4 || (contBJetsDurham == 3 && contNBJetsDurham == 1)) flagPreselection=true;
	      		}
	      		else throw std::runtime_error("ERROR: unknown preselection!");
	      		
	      		//if(contBJetsAntiKt == 4)
	      		if(flagPreselection==true)
	      		//if(1==1)
	      		{ 
		      		
		      		nParticles=0;
		      		contEntriesPostFilter++;
		      		//if(contEntriesPostFilter>1) return;
		      		
				if(branchParticle->GetEntries() > 0)  nParticles =  branchParticle->GetEntries();
				
				int contTopEvent=0, contZEvent=0, contWEvent=0, contLeptonEvent=0, contHadronEvent=0, contNeutrinoEvent=0, contG=0, contB=0, contWWEvent=0, contZHEvent=0;
				
				//////Looking only for HH->bbbb for signal
				if(topology == 1)
				{
					int loopCont=nParticles;
    				if(nParticles>25) loopCont=25;
					for(int i=0;i<loopCont;i++)
		    		{
			  			GenParticle *particle = (GenParticle*) branchParticle->At(i); 
			  			int pid = abs(particle->PID);
			  			int m1 = findMother(branchParticle, particle, 1);
      					if(pid == 5 && (m1 == 25 || m1 == 36)) contB++;
		     		}
	     		}
	     		///////Looking only for HH->bbbb for signal
				
				/*///////Lepton veto for ttbar
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

					//contLeptonEvent=0;

					//for(int i=0; i<nParticles; i++)
					//{
					//	GenParticle *particle = (GenParticle*) branchParticle->At(i);
					//	int pid = abs(particle->PID);
					//	int m1Pid = findMother(branchParticle, particle, 1);
					//	//int m2Pid = findMother(branchParticle, particle, 2);
					//	if((pid == 11 || pid == 13 || pid == 15) && (m1Pid == 6))
					//	{
					//		double electronCosTheta = abs(findCosTheta(particle->Eta));
					//		if(electronCosTheta < 0.95) contLeptonEvent++;
					//	} 
					//}

					//int nElectrons = branchElectron->GetEntries();
					//int nMuons = branchMuon->GetEntries();
					//contLeptonEvent += nElectrons + nMuons;

					//contLeptonEvent=0;
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

				///////Lepton veto for qqX, qqqqX, or qqHX
				if(topology == 6 || topology == 7 || topology == 8) 
				{
					int contEEGamma=0;
					float electronCosTheta=0;
					for(int i=0; i<13; i++)
					{
						GenParticle *particle = (GenParticle*) branchParticle->At(i);
						int pid = abs(particle->PID);
						int m1Pid = findMother(branchParticle, particle, 1);
						int m2Pid = findMother(branchParticle, particle, 2);
						if(contEEGamma == 1 && pid == 11 && (m1Pid == 22 || m1Pid == 11) && (m2Pid == 22 || m2Pid == 11)) electronCosTheta = findCosTheta(particle->Eta);
						if(pid == 11  && (m1Pid == 22 || m1Pid == 11) && (m2Pid == 22 || m2Pid == 11)) contEEGamma++;
					}
					electronCosTheta = abs(electronCosTheta);
					if(contEEGamma >= 2 && electronCosTheta < 0.95) contLeptonEvent++;
				}
				///////Lepton veto for qqX, qqqqX, or qqHX

				///////Lepton veto for ZH
				if(topology == 9)
				{
					for(int i=0;i<nParticles;i++)
			    	{
			  			GenParticle *particle = (GenParticle*) branchParticle->At(i); 
			  			int pid = abs(particle->PID);
				  		vector<int> finalDecays;
				  		if(((pid == 23 || pid == 32 || pid == 33) || (pid == 25 || pid == 36)) && abs(findMother(branchParticle, particle, 1)) == 22 && abs(findMother(branchParticle, particle, 2)) == 22)
				  		{
					  		findFinalDecaysZH(branchParticle, particle, finalDecays, contNeutrinoEvent, contLeptonEvent, contHadronEvent, i);
					  		contZHEvent++;
					  	}
			    	}
					if(contLeptonEvent>0) trueLeptonEvents++;
					if(contHadronEvent>0) trueHadronEvents++;
					if(contNeutrinoEvent>0) trueNeutrinoEvents++;
					
					if(contZHEvent == 1) trueSemiZHEvents++;
					if(contZHEvent >= 2) trueFullZHEvents++;
					if(contLeptonEvent == 2 && contHadronEvent == 2) trueLeptonHadronEvents++;
					if(contLeptonEvent == 2 && contNeutrinoEvent == 2) trueLeptonNeutrinoEvents++;
					if(contNeutrinoEvent == 2 && contHadronEvent == 2) trueHadronNeutrinoEvents++;
					if(contLeptonEvent == 4) trueFullLeptonEvents++;
					if(contHadronEvent == 4) trueFullHadronEvents++;
					if(contNeutrinoEvent == 4) trueFullNeutrinoEvents++;
	     		}
	     		///////Lepton veto for ZH

				///////Lepton veto for pebb
				if(topology == 10)
				{
					int contpebb=0;
					float electronCosTheta=0;
					for(int i=6; i<25; i++)
					{
						GenParticle *particle = (GenParticle*) branchParticle->At(i);
						int pid = abs(particle->PID);
						int m1Pid = findMother(branchParticle, particle, 1);
						//int m2Pid = findMother(branchParticle, particle, 2);
						if((pid == 11 || pid == 13 || pid == 15) && (m1Pid == 5))
						{
							electronCosTheta = abs(findCosTheta(particle->Eta));
							if(electronCosTheta < 0.95) contLeptonEvent++;
						} 
					}
				}
				//////////Lepton veto for pebb

				///////Lepton veto for pebbqq
				if(topology == 11)
				{
					int contpebbqq=0;
					float electronCosTheta=0;
					for(int i=4; i<25; i++)
					{
						GenParticle *particle = (GenParticle*) branchParticle->At(i);
						int pid = abs(particle->PID);
						int m1Pid = findMother(branchParticle, particle, 1);
						//int m2Pid = findMother(branchParticle, particle, 2);
						if((pid == 11 || pid == 13 || pid == 15) && (m1Pid >= 1 && m1Pid <= 8))
						{
							electronCosTheta = abs(findCosTheta(particle->Eta));
							if(electronCosTheta < 0.95) contLeptonEvent++;
						} 
					}
				}
				//////////Lepton veto for pebbqq

				///////Lepton veto for peqqH
				if(topology == 12)
				{
					int contpeqqH=0;
					float electronCosTheta=0;
					for(int i=4; i<35; i++)
					{
						GenParticle *particle = (GenParticle*) branchParticle->At(i);
						int pid = abs(particle->PID);
						int m1Pid = findMother(branchParticle, particle, 1);
						//int m2Pid = findMother(branchParticle, particle, 2);
						if((pid == 11 || pid == 13 || pid == 15) && ((m1Pid >= 1 && m1Pid <= 8) || (m1Pid == 25)))
						{
							electronCosTheta = abs(findCosTheta(particle->Eta));
							if(electronCosTheta < 0.95) contLeptonEvent++;
						} 
					}
				}
				//////////Lepton veto for peqqH


				/////Lepton veto for pett
				////////CHECK TERMINAR
				/////Lepton veto for pett*/
				

				/////General lepton veto
				contLeptonEvent=0;
				int nElectrons = branchElectron->GetEntries();
				int nMuons = branchMuon->GetEntries();
				contLeptonEvent += nElectrons + nMuons;
				//////

			
					
	     			
	     			if(topology != 1)
	     			{
	     				contB=4;
	     			}
	     			
	     			if(contB>=4)
		     		//if(1==1)
		     		{
		     			contEvents++;
		     			histNJets2D->Fill(contBJetsDurham,contNBJetsDurham);
		     			if(contLeptonEvent == 0) //lepton veto
		     			//if(1==1)
		     			{
			     			contEventsPostFilter++;
			     			histNJetsCompareAlgos->Fill(nJetsDurham, nJetsAntiKt);
			     			histNJetsDurham0->Fill(nJetsDurham0, weight);
			     			histNJetsDurham5->Fill(nJetsDurham5, weight);
			     			histNJetsDurham10->Fill(nJetsDurham10, weight);
			     			histNJetsDurham15->Fill(nJetsDurham15, weight);
			     			histNJetsDurham20->Fill(nJetsDurham20, weight);
			     			histNJetsDurham25->Fill(nJetsDurham25, weight);
			     			histNJetsDurham30->Fill(nJetsDurham30, weight);

							///////kinematic fit pythia
							TLorentzVector jetB1DurhamCopy, jetB2DurhamCopy, jetB3DurhamCopy, jetB4DurhamCopy;
							jetB1DurhamCopy = jetB1Durham;
							jetB2DurhamCopy = jetB2Durham;
							jetB3DurhamCopy = jetB3Durham;
							jetB4DurhamCopy = jetB4Durham;
							vector<TLorentzVector> jetBoffset;
							vector<TLorentzVector> jetsOrderedPT = {jetB1DurhamCopy, jetB2DurhamCopy, jetB3DurhamCopy, jetB4DurhamCopy};
							int indexBOffset=0;
							for(auto tlvb : jetsOrderedPT) {
								jetBoffset.push_back(tlvb);
								double e0=tlvb.E();
								double ef=e0*(1-BfractionalEnergyOffset[indexBOffset]);
								double pxf=tlvb.Px()*ef/e0;
								double pyf=tlvb.Py()*ef/e0;
								double pzf=tlvb.Pz()*ef/e0;
								jetBoffset[indexBOffset].SetPxPyPzE(pxf,pyf,pzf,ef);
								indexBOffset++;
							}
							///////////kinematic fit pythia

							/*///////kinematic fit Eqm
							eqmConstrainHHPairing(jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham, jetPairB1Fiteqm, jetPairB2Fiteqm, chi2ndfeqm, BfractionalEnergyErrors, BbetaXErrors, BbetaYErrors, BbetaZErrors, Ecm, enableExtraTries, nSigVar, jetB1Fiteqm, jetB2Fiteqm, jetB3Fiteqm, jetB4Fiteqm);
							invMassB1Fiteqm = jetPairB1Fiteqm.M();
				     		invMassB2Fiteqm = jetPairB2Fiteqm.M();
							histInvMassB1Fiteqm->Fill(invMassB1Fiteqm, weight);
							histInvMassB2Fiteqm->Fill(invMassB2Fiteqm, weight);
							//histChi2ndfeqm->Fill(chi2ndfeqm, weight);
							//histInvMass2DFiteqm->Fill(invMassB1Fiteqm, invMassB2Fiteqm, weight);
							minJetChiSeqm = pow((jetPairB1Fiteqm.M()-125), 2)/1 + pow((jetPairB2Fiteqm.M()-125), 2)/1;
							histChi2ndfeqm->Fill(minJetChiSeqm, weight);
							///////////kinematic fit Eqm*/

							/*///////kinematic fit Eqm pythia
							eqmConstrainHHPairingPythia(jetBoffset, jetsOrderedPT, jetPairB1Fiteqm, jetPairB2Fiteqm, chi2ndfeqm, BfractionalEnergyErrors, BbetaXErrors, BbetaYErrors, BbetaZErrors, Ecm, enableExtraTries, nSigVar, jetB1Fiteqm, jetB2Fiteqm, jetB3Fiteqm, jetB4Fiteqm);
							invMassB1Fiteqm = jetPairB1Fiteqm.M();
				     		invMassB2Fiteqm = jetPairB2Fiteqm.M();
							histInvMassB1Fiteqm->Fill(invMassB1Fiteqm, weight);
							histInvMassB2Fiteqm->Fill(invMassB2Fiteqm, weight);
							//histChi2ndfeqm->Fill(chi2ndfeqm, weight);
							//histInvMass2DFiteqm->Fill(invMassB1Fiteqm, invMassB2Fiteqm, weight);
							minJetChiSeqm = pow((jetPairB1Fiteqm.M()-125), 2)/1 + pow((jetPairB2Fiteqm.M()-125), 2)/1;
							histChi2ndfeqm->Fill(minJetChiSeqm, weight);
							///////////kinematic fit Eqm pythia*/
			     			
							/*////////kinematic fit
							TLorentzVector jetB1DurhamCopy, jetB2DurhamCopy, jetB3DurhamCopy, jetB4DurhamCopy;
							jetB1DurhamCopy = jetB1Durham;
							jetB2DurhamCopy = jetB2Durham;
							jetB3DurhamCopy = jetB3Durham;
							jetB4DurhamCopy = jetB4Durham;
							vector<LorentzVectorWithErrors>* inputLVWE=new vector<LorentzVectorWithErrors>(4);
							inputLVWE->at(0)=LorentzVectorWithErrors(jetB1DurhamCopy,BfractionalEnergyErrors[0]*jetB1DurhamCopy.Energy(),BbetaXErrors[0],BbetaYErrors[0],BbetaZErrors[0]);
		      				inputLVWE->at(1)=LorentzVectorWithErrors(jetB2DurhamCopy,BfractionalEnergyErrors[1]*jetB2DurhamCopy.Energy(),BbetaXErrors[1],BbetaYErrors[1],BbetaZErrors[1]);
							inputLVWE->at(2)=LorentzVectorWithErrors(jetB3DurhamCopy,BfractionalEnergyErrors[2]*jetB3DurhamCopy.Energy(),BbetaXErrors[2],BbetaYErrors[2],BbetaZErrors[2]);
							inputLVWE->at(3)=LorentzVectorWithErrors(jetB4DurhamCopy,BfractionalEnergyErrors[3]*jetB4DurhamCopy.Energy(),BbetaXErrors[3],BbetaYErrors[3],BbetaZErrors[3]);

							cout << " signal epConstrainHH " << endl;
							epConstrainHH epCHH(inputLVWE,Ecm,enableExtraTries,nSigVar, false);
							vector<LorentzVectorWithErrors>* outputLVep=epCHH.NumericalMinimization();

							cout << " signal pxyConstrainHH " << endl;
							pxyConstrainHH pxyCHH(inputLVWE,Ecm,enableExtraTries,nSigVar);
							vector<LorentzVectorWithErrors>* outputLVpxy=pxyCHH.NumericalMinimization();

							jetB1Fitep=outputLVep->at(0).getLV();
							jetB2Fitep=outputLVep->at(1).getLV();
							jetB3Fitep=outputLVep->at(2).getLV();
							jetB4Fitep=outputLVep->at(3).getLV();
							chi2ndfep=outputLVep->at(0).getChi2perNDF();

							jetB1Fitpxy=outputLVpxy->at(0).getLV();
							jetB2Fitpxy=outputLVpxy->at(1).getLV();
							jetB3Fitpxy=outputLVpxy->at(2).getLV();
							jetB4Fitpxy=outputLVpxy->at(3).getLV();
							chi2ndfpxy=outputLVpxy->at(0).getChi2perNDF();

							delete inputLVWE;
							delete outputLVep;
							delete outputLVpxy;
							//////////kinematic fit*/

							///////////kinematic fit pythia
							vector<LorentzVectorWithErrors>* inputLVWE=new vector<LorentzVectorWithErrors>(4);
							inputLVWE->at(0)=LorentzVectorWithErrors(jetBoffset[0],BfractionalEnergyErrors[0]*jetsOrderedPT[0].Energy(),BbetaXErrors[0],BbetaYErrors[0],BbetaZErrors[0]);
		      				inputLVWE->at(1)=LorentzVectorWithErrors(jetBoffset[1],BfractionalEnergyErrors[1]*jetsOrderedPT[1].Energy(),BbetaXErrors[1],BbetaYErrors[1],BbetaZErrors[1]);
							inputLVWE->at(2)=LorentzVectorWithErrors(jetBoffset[2],BfractionalEnergyErrors[2]*jetsOrderedPT[2].Energy(),BbetaXErrors[2],BbetaYErrors[2],BbetaZErrors[2]);
							inputLVWE->at(3)=LorentzVectorWithErrors(jetBoffset[3],BfractionalEnergyErrors[3]*jetsOrderedPT[3].Energy(),BbetaXErrors[3],BbetaYErrors[3],BbetaZErrors[3]);

							cout << " signal epConstrainHH " << endl;
							epConstrainHH epCHH(inputLVWE,Ecm,enableExtraTries,nSigVar, false);
							vector<LorentzVectorWithErrors>* outputLVep=epCHH.NumericalMinimization();

							cout << " signal pxyConstrainHH " << endl;
							pxyConstrainHH pxyCHH(inputLVWE,Ecm,enableExtraTries,nSigVar);
							vector<LorentzVectorWithErrors>* outputLVpxy=pxyCHH.NumericalMinimization();

							jetB1Fitep=outputLVep->at(0).getLV();
							jetB2Fitep=outputLVep->at(1).getLV();
							jetB3Fitep=outputLVep->at(2).getLV();
							jetB4Fitep=outputLVep->at(3).getLV();
							chi2ndfep=outputLVep->at(0).getChi2perNDF();

							jetB1Fitpxy=outputLVpxy->at(0).getLV();
							jetB2Fitpxy=outputLVpxy->at(1).getLV();
							jetB3Fitpxy=outputLVpxy->at(2).getLV();
							jetB4Fitpxy=outputLVpxy->at(3).getLV();
							chi2ndfpxy=outputLVpxy->at(0).getChi2perNDF();

							delete inputLVWE;
							delete outputLVep;
							delete outputLVpxy;
							///////////kinematic fit pythia
							
							/////kinematic fitted inv. mass
							double jetPairB1Index1Fitep, jetPairB1Index2Fitep, jetPairB2Index1Fitep, jetPairB2Index2Fitep;
			     			TLorentzVector jetPairB1Fitep, jetPairB2Fitep, jetPairB1SenFitep, jetPairB2SenFitep, jetPairB1AntiKtFitep, jetPairB2AntiKtFitep;
			     			TLorentzVector jetPairB1CheckFitep, jetPairB2CheckFitep;
			     			findJetPairs(jetPairB1Fitep, jetPairB2Fitep, jetB1Fitep, jetB2Fitep, jetB3Fitep, jetB4Fitep, jetPairB1Index1Fitep, jetPairB1Index2Fitep, jetPairB2Index1Fitep, jetPairB2Index2Fitep, minJetChiSep);
			     			invMassB1Fitep = jetPairB1Fitep.M();
				     		invMassB2Fitep = jetPairB2Fitep.M();
							histInvMassB1Fitep->Fill(invMassB1Fitep, weight);
							histInvMassB2Fitep->Fill(invMassB2Fitep, weight);
							histChi2ndfep->Fill(chi2ndfep, weight);
							//histChi2ndfep->Fill(minJetChiSep, weight);
							//histInvMass2DFitep->Fill(invMassB1Fitep, invMassB2Fitep, weight);

							double jetPairB1Index1Fitpxy, jetPairB1Index2Fitpxy, jetPairB2Index1Fitpxy, jetPairB2Index2Fitpxy;
			     			TLorentzVector jetPairB1Fitpxy, jetPairB2Fitpxy, jetPairB1SenFitpxy, jetPairB2SenFitpxy, jetPairB1AntiKtFitpxy, jetPairB2AntiKtFitpxy;
			     			TLorentzVector jetPairB1CheckFitpxy, jetPairB2CheckFitpxy;
			     			findJetPairs(jetPairB1Fitpxy, jetPairB2Fitpxy, jetB1Fitpxy, jetB2Fitpxy, jetB3Fitpxy, jetB4Fitpxy, jetPairB1Index1Fitpxy, jetPairB1Index2Fitpxy, jetPairB2Index1Fitpxy, jetPairB2Index2Fitpxy, minJetChiSpxy);
			     			invMassB1Fitpxy = jetPairB1Fitpxy.M();
				     		invMassB2Fitpxy = jetPairB2Fitpxy.M();
							histInvMassB1Fitpxy->Fill(invMassB1Fitpxy, weight);
							histInvMassB2Fitpxy->Fill(invMassB2Fitpxy, weight);
							histChi2ndfpxy->Fill(chi2ndfpxy, weight);
							//histChi2ndfpxy->Fill(minJetChiSpxy, weight);
							//histInvMass2DFitpxy->Fill(invMassB1Fitpxy, invMassB2Fitpxy, weight);

							// if(chi2ndfep < chi2ndfpxy && chi2ndfep < chi2ndfeqm)
							// {
							//     invMassB1FitBest = invMassB1Fitep;
							//     invMassB2FitBest = invMassB2Fitep;
							//     chi2ndfBest = chi2ndfep;
							// }
							// else if(chi2ndfpxy < chi2ndfep && chi2ndfpxy < chi2ndfeqm)
							// {
							//     invMassB1FitBest = invMassB1Fitpxy;
							//     invMassB2FitBest = invMassB2Fitpxy;
							//     chi2ndfBest = chi2ndfpxy;
							// }
							// else if(chi2ndfeqm < chi2ndfep && chi2ndfeqm < chi2ndfpxy)
							// {
							//     invMassB1FitBest = invMassB1Fiteqm;
							//     invMassB2FitBest = invMassB2Fiteqm;
							//     chi2ndfBest = chi2ndfeqm;
							// }
							// histInvMassB1FitBest->Fill(invMassB1FitBest, weight);
							// histInvMassB2FitBest->Fill(invMassB2FitBest, weight);
							// histChi2ndfBest->Fill(chi2ndfBest, weight);
							// histInvMass2DFitBest->Fill(invMassB1FitBest, invMassB2FitBest, weight);

								//////kinematic fit pythia
							if(chi2ndfep < 0.4)
							{
								invMassB1FitBest = invMassB1Fitep;
								invMassB2FitBest = invMassB2Fitep;
								chi2ndfBest = chi2ndfep;
							}
							else
							{
								invMassB1FitBest = invMassB1Fitpxy;
								invMassB2FitBest = invMassB2Fitpxy;
								chi2ndfBest = chi2ndfpxy;
							}
							histInvMassB1FitBest->Fill(invMassB1FitBest, weight);
							histInvMassB2FitBest->Fill(invMassB2FitBest, weight);
							histChi2ndfBest->Fill(chi2ndfBest, weight);
							jetPairMassDeltaFit = pow(invMassB1FitBest - 125.2, 2) + pow(invMassB2FitBest - 125.2, 2);
								///////kinematic fit pythia
							
							//////kinematic fitted inv. mass

			     			//////inv. mass
			     			double jetPairB1Index1, jetPairB1Index2, jetPairB2Index1, jetPairB2Index2;
			     			TLorentzVector jetPairB1, jetPairB2, jetPairB1Sen, jetPairB2Sen, jetPairB1AntiKt, jetPairB2AntiKt;
			     			TLorentzVector jetPairB1Check, jetPairB2Check;
			     			findJetPairs(jetPairB1, jetPairB2, jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham, jetPairB1Index1, jetPairB1Index2, jetPairB2Index1, jetPairB2Index2, minJetChiS);
			     			invMassB1 = jetPairB1.M();
				     		invMassB2 = jetPairB2.M();
							invMass4Jets = (jetB1Durham + jetB2Durham + jetB3Durham + jetB4Durham).M();
							histInvMass4Jets->Fill(invMass4Jets, weight);
							jetPairMassDelta = pow(invMassB1 - 125.2, 2) + pow(invMassB2 - 125.2, 2);
							histJetPairMassDelta->Fill(jetPairMassDelta, weight);
				     		
			     			bool flagZZMass=false;
			     			TLorentzVector jetPairB1ZZ, jetPairB2ZZ;
			     			double jetPairB1ZZIndex1, jetPairB1ZZIndex2, jetPairB2ZZIndex1, jetPairB2ZZIndex2;
			     			distanceZMass=0.05;
			     			findJetPairsZZ(jetPairB1ZZ, jetPairB2ZZ, jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham, jetPairB1ZZIndex1, jetPairB1ZZIndex2, jetPairB2ZZIndex1, jetPairB2ZZIndex2, flagZZMass, minChiSquaredZZMass, distanceZMass, contEventsZWindow10, contEventsZWindow01, contEventsZWindow05, contEventsZWindow1, contEventsZWindow15, contEventsZWindow2, contEventsZWindow5, distanceZ1MinChiSquaredZZMass, distanceZ2MinChiSquaredZZMass);
			     			histMinChiSquaredZZMass->Fill(minChiSquaredZZMass, weight);
			     			invMassZZ1 = jetPairB1ZZ.M();
			     			invMassZZ2 = jetPairB2ZZ.M();
			     			histInvMassZZ1->Fill(invMassZZ1, weight);
			     			histInvMassZZ2->Fill(invMassZZ2, weight);
			     			histDistanceZ1MinChiSquaredZZMass->Fill(distanceZ1MinChiSquaredZZMass, weight);
			     			histDistanceZ2MinChiSquaredZZMass->Fill(distanceZ2MinChiSquaredZZMass, weight);

							/*//////kinematic fit
							if(minJetChiSep < minJetChiSpxy && minJetChiSep < minJetChiSeqm)
							{
								invMassB1FitBest = invMassB1Fitep;
								invMassB2FitBest = invMassB2Fitep;
								minJetChiSBest = minJetChiSep;
							}
							else if(minJetChiSpxy < minJetChiSep && minJetChiSpxy < minJetChiSeqm)
							{
								invMassB1FitBest = invMassB1Fitpxy;
								invMassB2FitBest = invMassB2Fitpxy;
								minJetChiSBest = minJetChiSpxy;
							}
							else if(minJetChiSeqm < minJetChiSep && minJetChiSeqm < minJetChiSpxy)
							{
								invMassB1FitBest = invMassB1Fiteqm;
								invMassB2FitBest = invMassB2Fiteqm;
								minJetChiSBest = minJetChiSeqm;
							}
							histInvMassB1FitBest->Fill(invMassB1FitBest, weight);
							histInvMassB2FitBest->Fill(invMassB2FitBest, weight);
							histChi2ndfBest->Fill(minJetChiSBest, weight);
							//invMassB1 = invMassB1FitBest;
							//invMassB2 = invMassB2FitBest;
							/////CHECK CHANGE FOR FITTED ANALYSIS KINEMATIC FIT
							//////kinematic fit*/

							/*bool flagZHMass=false;
			     			TLorentzVector jetPairB1ZH, jetPairB2ZH;
			     			double jetPairB1ZHIndex1, jetPairB1ZHIndex2, jetPairB2ZHIndex1, jetPairB2ZHIndex2;
			     			distanceZMass=0.05;
			     			findJetPairsZH(jetPairB1ZH, jetPairB2ZH, jetB1Durham, jetB2Durham, jetB3Durham, jetB4Durham, jetPairB1ZHIndex1, jetPairB1ZHIndex2, jetPairB2ZHIndex1, jetPairB2ZHIndex2, flagZHMass, minChiSquaredZHMass, distanceZMass, contEventsZWindow10, contEventsZWindow01, contEventsZWindow05, contEventsZWindow1, contEventsZWindow15, contEventsZWindow2, contEventsZWindow5, distanceZ1MinChiSquaredZHMass, distanceZ2MinChiSquaredZHMass, invMassZHZ, invMassZHH, contEventsPostFilter);
			     			histMinChiSquaredZHMass->Fill(minChiSquaredZHMass, weight);
			     			//histInvMassZH1->Fill(invMassZHZ, weight);
			     			//histInvMassZH2->Fill(invMassZHH, weight);
			     			histDistanceZ1MinChiSquaredZHMass->Fill(distanceZ1MinChiSquaredZHMass, weight);
			     			histDistanceZ2MinChiSquaredZHMass->Fill(distanceZ2MinChiSquaredZHMass, weight);
			     			//if(flagZHMass==true) break;*/
				     		
			     			vector<int> bestGroup1Durham, bestGroup2Durham, bestGroup1AntiKt, bestGroup2AntiKt;
			     			vector<double> chiSquares, chiSquaresAntiKt;
			     			vector<vector<int>> group1s, group2s, group1sAntiKt, group2sAntiKt;
    						findBestCombination(branchJetDurham, nJetsDurham, bestGroup1Durham, bestGroup2Durham, chiSquares, group1s, group2s);
    						findBestCombination(branchJetAntiKt, nJetsAntiKt, bestGroup1AntiKt, bestGroup2AntiKt, chiSquaresAntiKt, group1sAntiKt, group2sAntiKt);

							//findBestCombinationZH(branchJetDurham, nJetsDurham, bestGroup1Durham, bestGroup2Durham, chiSquares, group1s, group2s);
							sortVectorTriple(chiSquares, group1s, group2s);
    						double senChiSquared=0;
    						queue<int> combinationsIndices;
    						int chiSquaresIndex=0;
							double recoMassZ=0, recoMassH=0;
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
							    	
									/*bool flagZ1H2 = tagZHJetPairMasses(jetPairB1Sen, jetPairB2Sen);
									if(flagZ1H2)
									{
										recoMassZ = jetPairB1Sen.M();
										recoMassH = jetPairB2Sen.M();
									}
									else
									{
										recoMassH = jetPairB1Sen.M();
										recoMassZ = jetPairB2Sen.M();
									}*/
									invMassB11Best = jetPairB1Sen.M();
				     				invMassB21Best = jetPairB2Sen.M();
									histJetB1M1Best->Fill(invMassB11Best, weight);
									histJetB2M1Best->Fill(invMassB21Best, weight);
									/*histInvMassZH1->Fill(recoMassZ, weight);
			     					histInvMassZH2->Fill(recoMassH, weight);*/
									TLorentzVector diJets = jetPairB1Sen + jetPairB2Sen;
									double invMassDiH = diJets.M();
									histDiHM->Fill(invMassDiH, weight);
									
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
			     			histMinJetChiS->Fill(minJetChiS, weight);
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

						deltaRJetPairs = jetPairB1.DeltaR(jetPairB2);
						histDeltaRJetPairs->Fill(deltaRJetPairs, weight);
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
						
						///////////////////y_nm for jets
						if(nJetsDurham>0)
						{
							for(int i=0;i<1;i++)
							{
								Jet *jet = (Jet*) branchJetDurham->At(i);
								exclYmerge12 = jet->ExclYmerge12;
								exclYmerge23 = jet->ExclYmerge23;
								exclYmerge34 = jet->ExclYmerge34;
								exclYmerge45 = jet->ExclYmerge45;
								exclYmerge56 = jet->ExclYmerge56;
								histExclYmerge12->Fill(exclYmerge12, weight);
								histExclYmerge23->Fill(exclYmerge23, weight);
								histExclYmerge34->Fill(exclYmerge34, weight);
								histExclYmerge45->Fill(exclYmerge45, weight);
								histExclYmerge56->Fill(exclYmerge56, weight);
	    					}
    					}
    						
						/*std::vector<fastjet::PseudoJet> eflowObjects;
						for (int j=0; j<branchEFlowTrack->GetEntries(); j++) 
						{
							TLorentzVector *eflowTrack = (TLorentzVector*)branchEFlowTrack->At(j);
							eflowObjects.emplace_back(eflowTrack->Px(), eflowTrack->Py(), eflowTrack->Pz(), eflowTrack->E());
							
						}
						for (int j=0; j<branchEFlowPhoton->GetEntries(); j++) 
						{
							TLorentzVector *eflowPhoton = (TLorentzVector*)branchEFlowPhoton->At(j);
						    	eflowObjects.emplace_back(eflowPhoton->Px(), eflowPhoton->Py(), eflowPhoton->Pz(), eflowPhoton->E());
						}
						for (int j=0; j<branchEFlowNeutralHadron->GetEntries(); j++) 
						{
							TLorentzVector *eflowNeutralHadron = (TLorentzVector*)branchEFlowNeutralHadron->At(j);
						    	eflowObjects.emplace_back(eflowNeutralHadron->Px(), eflowNeutralHadron->Py(), eflowNeutralHadron->Pz(), eflowNeutralHadron->E());
						}
						fastjet::JetDefinition jet_def(fastjet::ee_kt_algorithm);
						fastjet::ClusterSequence cs(eflowObjects, jet_def);
						std::vector<fastjet::PseudoJet> exclusive_jets = fastjet::sorted_by_pt(cs.exclusive_jets_up_to(4));
						//Apply minimum PT cut
						std::vector<fastjet::PseudoJet> final_jets;
						for (const auto& jet : exclusive_jets) {
						    if (jet.pt() > 10.0) {
							final_jets.push_back(jet);
						    }
						}
						//Apply minimum PT cut
						double y_34 = cs.exclusive_ymerge(3); // 3-jets to 4-jets
        					double y_45 = cs.exclusive_ymerge(4); // 4-jets to 5-jets
        					if(contEventsPostFilter<10) cout<<"Entry " << entry << ": y_34 = " << y_34 << ", y_45 = " << y_45 <<endl;*/
						//////////////////////y_nm for jets
						
						/////Thrust
						vector<TLorentzVector> momenta;
						for (int i = 0; i < branchEFlowTrack->GetEntries(); ++i) 
						{
						    	Track *track = (Track*) branchEFlowTrack->At(i);
						    	TLorentzVector track4V = track->P4();
						    	momenta.push_back(track4V);
						}
						for (int i = 0; i < branchEFlowPhoton->GetEntries(); ++i) 
						{
						    	Tower *photon = (Tower*) branchEFlowPhoton->At(i);
						    	TLorentzVector photon4V = photon->P4();
						    	momenta.push_back(photon4V);
						}
						for (int i = 0; i < branchEFlowNeutralHadron->GetEntries(); ++i) 
						{
						    	Tower *neutralHadron = (Tower*) branchEFlowNeutralHadron->At(i);
						    	TLorentzVector neutralHadron4V = neutralHadron->P4();
						    	momenta.push_back(neutralHadron4V);
						}
						TLorentzVector thrustAxis;
						thrust = findThrust(momenta, thrustAxis);
						////////Thrust

						//////// Boost
						TLorentzVector totalSystem = jetB1Durham + jetB2Durham + jetB3Durham + jetB4Durham;
						boostB1 = jetB1Durham.Pz()/jetB1Durham.E();
						boostB2 = jetB2Durham.Pz()/jetB2Durham.E();
						boostB3 = jetB3Durham.Pz()/jetB3Durham.E();
						boostB4 = jetB4Durham.Pz()/jetB4Durham.E();
						boostSystem = totalSystem.Pz()/totalSystem.E();
						/////////// Boost

						///////Missing ET
						MissingET *missinget = (MissingET*) branchMissingET->At(0); 
						missingET = missinget->MET;
						/////Missing ET

						///////Energy and PT Assymetries
						/*eAssymetry = (abs(jetPairB1.E() - jetPairB2.E()))/(jetPairB1.E() + jetPairB2.E());
						pTAssymetry = (abs(jetPairB1.Pt() - jetPairB2.Pt()))/(jetPairB1.Pt() + jetPairB2.Pt());
						histEAssymetry->Fill(eAssymetry, weight);
						histPTAssymetry->Fill(pTAssymetry, weight);*/

						float E1 = jetPairB1.E();
						float E2 = jetPairB2.E();
						if((E1 + E2) > 0) {
							eAssymetry = fabs(E1 - E2) / (E1 + E2);
							histEAssymetry->Fill(eAssymetry);
						}
						float PT1 = jetPairB1.Pt();
						float PT2 = jetPairB2.Pt();
						if((PT1 + PT2) > 0) {
							pTAssymetry = fabs(PT1 - PT2) / (PT1 + PT2);
							histPTAssymetry->Fill(pTAssymetry);
						}
						///////Energy and PT Assymetries

						/////////N-subjetiness
						for(int i = 0; i < branchJetDurham->GetEntries(); ++i)
						{
							Jet *jet = (Jet*) branchJetDurham->At(i);
							if(i==0)
							{
								nSubJetsTrimmedB1 = jet->NSubJetsTrimmed;
								nSubJetsPrunnedB1 = jet->NSubJetsPruned;
								nSubJetsSoftDroppedB1 = jet->NSubJetsSoftDropped;
							} else if(i==1)
							{
								nSubJetsTrimmedB2 = jet->NSubJetsTrimmed;
								nSubJetsPrunnedB2 = jet->NSubJetsPruned;
								nSubJetsSoftDroppedB2 = jet->NSubJetsSoftDropped;
							}
							else if(i==2)
							{
								nSubJetsTrimmedB3 = jet->NSubJetsTrimmed;
								nSubJetsPrunnedB3 = jet->NSubJetsPruned;
								nSubJetsSoftDroppedB3 = jet->NSubJetsSoftDropped;
							}
							else if(i==3)
							{
								nSubJetsTrimmedB4 = jet->NSubJetsTrimmed;
								nSubJetsPrunnedB4 = jet->NSubJetsPruned;
								nSubJetsSoftDroppedB4 = jet->NSubJetsSoftDropped;
							}
						}
						histNSubJetsTrimmedB1->Fill(nSubJetsTrimmedB1);
						histNSubJetsTrimmedB2->Fill(nSubJetsTrimmedB2);
						histNSubJetsTrimmedB3->Fill(nSubJetsTrimmedB3);
						histNSubJetsTrimmedB4->Fill(nSubJetsTrimmedB4);
						histNSubJetsPrunnedB1->Fill(nSubJetsPrunnedB1);
						histNSubJetsPrunnedB2->Fill(nSubJetsPrunnedB2);
						histNSubJetsPrunnedB3->Fill(nSubJetsPrunnedB3);
						histNSubJetsPrunnedB4->Fill(nSubJetsPrunnedB4);
						histNSubJetsSoftDroppedB1->Fill(nSubJetsSoftDroppedB1);
						histNSubJetsSoftDroppedB2->Fill(nSubJetsSoftDroppedB2);
						histNSubJetsSoftDroppedB3->Fill(nSubJetsSoftDroppedB3);
						histNSubJetsSoftDroppedB4->Fill(nSubJetsSoftDroppedB4);
						//////N-subjetiness
						
						
						/////Filling out tree for BDT
						if(fileFunction == "generate")
						{
							TreeFull.Fill();
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
						}
						else if(fileFunction == "merge")
						{
							if(contEventsPostFilter<2) cout<<"Merging files!"<<endl;
							TreeMerge.Fill();
						}
						else if(fileFunction == "plot")
						{ 
							if(contEventsPostFilter<2) cout<<"Only plotting -- no file generation!"<<endl;
						}
						else throw std::runtime_error("ERROR: unknown fileFunction!");
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
		TCanvas *c114 = new TCanvas();*/
		/*TCanvas *c115 = new TCanvas();
		TCanvas *c116 = new TCanvas();
		TCanvas *c117 = new TCanvas();
		TCanvas *c118 = new TCanvas();
		TCanvas *c119 = new TCanvas();*/
		/*c113->cd();
		histJetB1M->Draw("HIST");
		c114->cd();
		histJetB2M->Draw("HIST");*/
		/*c115->cd();
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
		histSpher->Draw("HIST")
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
    		
    		/*////nJetsDurham
    		TCanvas *c154 = new TCanvas();
    		c154->cd();
    		histNJetsDurham0->Draw("HIST");
    		TCanvas *c155 = new TCanvas();
    		c155->cd();
    		histNJetsDurham5->Draw("HIST");
    		TCanvas *c156 = new TCanvas();
    		c156->cd();
    		histNJetsDurham10->Draw("HIST");
    		TCanvas *c157 = new TCanvas();
    		c157->cd();
    		histNJetsDurham15->Draw("HIST");
    		TCanvas *c158 = new TCanvas();
    		c158->cd();
    		histNJetsDurham20->Draw("HIST");
    		TCanvas *c159 = new TCanvas();
    		c159->cd();
    		histNJetsDurham25->Draw("HIST");
    		TCanvas *c160 = new TCanvas();
    		c160->cd();
    		histNJetsDurham30->Draw("HIST");
    		////nJetsDurham*/
    		
    		/*///////minChiSquaredZZMass
    		TCanvas *c161 = new TCanvas();
    		c161->cd();
    		histMinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c162 = new TCanvas();
    		c162->cd();
    		histInvMassZZ1->Draw("HIST");
    		TLine* distanceZMassLeft = new TLine(90-distanceZMass, histInvMassZZ1->GetMinimum(), 90-distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassLeft->SetLineColor(kRed);
	        distanceZMassLeft->SetLineWidth(1);
	        TLine* distanceZMassRight = new TLine(90+distanceZMass, histInvMassZZ1->GetMinimum(), 90+distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassRight->SetLineColor(kRed);
	        distanceZMassRight->SetLineWidth(1);
	        distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
    		TCanvas *c163 = new TCanvas();
    		c163->cd();
    		histInvMassZZ2->Draw("HIST");
    		distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
	        TCanvas *c164 = new TCanvas();
    		c164->cd();
    		histDistanceZ1MinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c165 = new TCanvas();
    		c165->cd();
    		histDistanceZ2MinChiSquaredZZMass->Draw("HIST");
    		/////minChiSquaredZZMass*/
    		
    		/*////////y_nm
    		TCanvas *c166 = new TCanvas();
    		c166->cd();
    		histExclYmerge12->Draw("HIST");
    		TCanvas *c167 = new TCanvas();
    		c167->cd();
    		histExclYmerge23->Draw("HIST");
    		TCanvas *c168 = new TCanvas();
    		c168->cd();
    		histExclYmerge34->Draw("HIST");
    		TCanvas *c169 = new TCanvas();
    		c169->cd();
    		histExclYmerge45->Draw("HIST");
    		TCanvas *c170 = new TCanvas();
    		c170->cd();
    		histExclYmerge56->Draw("HIST");
    		///////y_nm*/

		/*/////////inv. mass for best comb
		TCanvas *c171 = new TCanvas();
		TCanvas *c172 = new TCanvas();
		c171->cd();
		histJetB1M1Best->Draw("HIST");
		c172->cd();
		histJetB2M1Best->Draw("HIST");
		/////////inv. mass for best comb*/

		/*////jet dist.
		TCanvas *c173 = new TCanvas();
		c173->cd();
    	histNJets2D->Draw("TEXT");
		/////jet dist.*/

		///////minChiSquaredZHMass
    	/*TCanvas *c174 = new TCanvas();
    	c174->cd();
		histMinChiSquaredZHMass->Draw("HIST");
		TCanvas *c175 = new TCanvas();
		c175->cd();
		histInvMassZH1->Draw("HIST");
		TLine* distanceZMassLeft = new TLine(90-distanceZMass, histInvMassZZ1->GetMinimum(), 90-distanceZMass, histInvMassZZ1->GetMaximum());
		distanceZMassLeft->SetLineColor(kRed);
		distanceZMassLeft->SetLineWidth(1);
		TLine* distanceZMassRight = new TLine(90+distanceZMass, histInvMassZZ1->GetMinimum(), 90+distanceZMass, histInvMassZZ1->GetMaximum());
		distanceZMassRight->SetLineColor(kRed);
		distanceZMassRight->SetLineWidth(1);
		distanceZMassLeft->Draw("same");
		distanceZMassRight->Draw("same");
		TCanvas *c176 = new TCanvas();
		c176->cd();
		histInvMassZH2->Draw("HIST");
		distanceZMassLeft->Draw("same");
		distanceZMassRight->Draw("same");
		TCanvas *c177 = new TCanvas();
		c177->cd();
		histDistanceZ1MinChiSquaredZZMass->Draw("HIST");
		TCanvas *c178 = new TCanvas();
		c178->cd();
		histDistanceZ2MinChiSquaredZZMass->Draw("HIST");*/
		/////minChiSquaredZHMass

		/*/////////di-Higgs mass
		TCanvas *c179 = new TCanvas();
		c179->cd();
		histDiHM->Draw("HIST");
		//////////di-Higgs Mass*/

		/*/////kinematic fitted inv. mass
		TCanvas *canvasB1Masses = new TCanvas("canvasB1Masses", "B1 Masses", 800, 800);
		canvasB1Masses->Divide(3, 2); // Dividing into 3 columns and 2 rows for 5 histograms
		canvasB1Masses->cd(1);
		histJetB1M->Draw("HIST");
		canvasB1Masses->cd(2);
		histInvMassB1Fitep->Draw("HIST");
		canvasB1Masses->cd(3);
		histInvMassB1Fitpxy->Draw("HIST");
		//canvasB1Masses->cd(4);
		//histInvMassB1Fiteqm->Draw("HIST");
		canvasB1Masses->cd(5);
		histInvMassB1FitBest->Draw("HIST");*/

		/*TCanvas *canvasB2Masses = new TCanvas("canvasB2Masses", "B2 Masses", 800, 800);
		canvasB2Masses->Divide(3, 2); // Dividing into 3 columns and 2 rows for 5 histograms
		canvasB2Masses->cd(1);
		histJetB2M->Draw("HIST");
		canvasB2Masses->cd(2);
		histInvMassB2Fitep->Draw("HIST");
		canvasB2Masses->cd(3);
		histInvMassB2Fitpxy->Draw("HIST");
		canvasB2Masses->cd(4);
		histInvMassB2Fiteqm->Draw("HIST");
		canvasB2Masses->cd(5);
		histInvMassB2FitBest->Draw("HIST");*/

		/*TCanvas *canvasChi2ndf = new TCanvas("canvasChi2ndf", "Chi2ndf", 800, 800);
		canvasChi2ndf->Divide(3, 2);
		//canvasChi2ndf->cd(1);
		//histMinJetChiS->Draw("HIST");
		canvasChi2ndf->cd(2);
		histChi2ndfep->Draw("HIST");
		canvasChi2ndf->cd(3);
		histChi2ndfpxy->Draw("HIST");
		//canvasChi2ndf->cd(4);
		//histChi2ndfeqm->Draw("HIST");
		canvasChi2ndf->cd(5);
		histChi2ndfBest->Draw("HIST");
		//////////kinematic fitted inv. mass*/

		/*//////qqH study
		TCanvas *c80canvasCombined = new TCanvas("c80canvasCombined", "Combined Plots", 1200, 800);
		c80canvasCombined->Divide(3, 2); // Dividing into 3 columns and 2 rows for 6 histograms
		c80canvasCombined->cd(1);
		histInvMass4Jets->Draw("HIST");
		c80canvasCombined->cd(2);
		histDeltaRJetPairs->Draw("HIST");
		c80canvasCombined->cd(3);
		histEAssymetry->Draw("HIST");
		c80canvasCombined->cd(4);
		histPTAssymetry->Draw("HIST");
		c80canvasCombined->cd(5);
		histJetPairMassDelta->Draw("HIST");
		/////*/

		/*//////N-subjetiness
		TCanvas *c81canvasNSubjetiness = new TCanvas("c81canvasNSubjetiness", "N-Subjetiness Plots", 1200, 900);
		c81canvasNSubjetiness->Divide(3, 4); // Dividing into 3 columns and 4 rows for 12 histograms
		
		c81canvasNSubjetiness->cd(1);
		histNSubJetsPrunnedB1->Draw("HIST");
		c81canvasNSubjetiness->cd(2);
		histNSubJetsPrunnedB2->Draw("HIST");
		c81canvasNSubjetiness->cd(3);
		histNSubJetsPrunnedB3->Draw("HIST");
		c81canvasNSubjetiness->cd(4);
		histNSubJetsPrunnedB4->Draw("HIST");
		c81canvasNSubjetiness->cd(5);
		histNSubJetsSoftDroppedB1->Draw("HIST");
		c81canvasNSubjetiness->cd(6);
		histNSubJetsSoftDroppedB2->Draw("HIST");
		c81canvasNSubjetiness->cd(7);
		histNSubJetsSoftDroppedB3->Draw("HIST");
		c81canvasNSubjetiness->cd(8);
		histNSubJetsSoftDroppedB4->Draw("HIST");
		c81canvasNSubjetiness->cd(9);
		histNSubJetsSoftDroppedB1->Draw("HIST");
		c81canvasNSubjetiness->cd(10);
		histNSubJetsSoftDroppedB2->Draw("HIST");
		c81canvasNSubjetiness->cd(11);
		histNSubJetsSoftDroppedB3->Draw("HIST");
		c81canvasNSubjetiness->cd(12);
		histNSubJetsSoftDroppedB4->Draw("HIST");
		///////N-subjetiness*/


    		
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
    		
    		/*////nJetsDurham
    		TCanvas *c254 = new TCanvas();
    		c254->cd();
    		histNJetsDurham0->Draw("HIST");
    		TCanvas *c255 = new TCanvas();
    		c255->cd();
    		histNJetsDurham5->Draw("HIST");
    		TCanvas *c256 = new TCanvas();
    		c256->cd();
    		histNJetsDurham10->Draw("HIST");
    		TCanvas *c257 = new TCanvas();
    		c257->cd();
    		histNJetsDurham15->Draw("HIST");
    		TCanvas *c258 = new TCanvas();
    		c258->cd();
    		histNJetsDurham20->Draw("HIST");
    		TCanvas *c259 = new TCanvas();
    		c259->cd();
    		histNJetsDurham25->Draw("HIST");
    		TCanvas *c260 = new TCanvas();
    		c260->cd();
    		histNJetsDurham30->Draw("HIST");
    		////nJetsDurham*/
    		
    		/*///////minChiSquaredZZMass
    		TCanvas *c261 = new TCanvas();
    		c261->cd();
    		histMinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c262 = new TCanvas();
    		c262->cd();
    		histInvMassZZ1->Draw("HIST");
    		TLine* distanceZMassLeft = new TLine(90-distanceZMass, histInvMassZZ1->GetMinimum(), 90-distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassLeft->SetLineColor(kRed);
	        distanceZMassLeft->SetLineWidth(1);
	        TLine* distanceZMassRight = new TLine(90+distanceZMass, histInvMassZZ1->GetMinimum(), 90+distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassRight->SetLineColor(kRed);
	        distanceZMassRight->SetLineWidth(1);
	        distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
    		TCanvas *c263 = new TCanvas();
    		c263->cd();
    		histInvMassZZ2->Draw("HIST");
    		distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
	        TCanvas *c264 = new TCanvas();
    		c264->cd();
    		histDistanceZ1MinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c265 = new TCanvas();
    		c265->cd();
    		histDistanceZ2MinChiSquaredZZMass->Draw("HIST");
    		/////minChiSquaredZZMass*/
    		
    		/*////////y_nm
    		TCanvas *c266 = new TCanvas();
    		c266->cd();
    		histExclYmerge12->Draw("HIST");
    		TCanvas *c267 = new TCanvas();
    		c267->cd();
    		histExclYmerge23->Draw("HIST");
    		TCanvas *c268 = new TCanvas();
    		c268->cd();
    		histExclYmerge34->Draw("HIST");
    		TCanvas *c269 = new TCanvas();
    		c269->cd();
    		histExclYmerge45->Draw("HIST");
    		TCanvas *c270 = new TCanvas();
    		c270->cd();
    		histExclYmerge56->Draw("HIST");
    		///////y_nm*/
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
    		
    		/*////nJetsDurham
    		TCanvas *c354 = new TCanvas();
    		c354->cd();
    		histNJetsDurham0->Draw("HIST");
    		TCanvas *c355 = new TCanvas();
    		c355->cd();
    		histNJetsDurham5->Draw("HIST");
    		TCanvas *c356 = new TCanvas();
    		c356->cd();
    		histNJetsDurham10->Draw("HIST");
    		TCanvas *c357 = new TCanvas();
    		c357->cd();
    		histNJetsDurham15->Draw("HIST");
    		TCanvas *c358 = new TCanvas();
    		c358->cd();
    		histNJetsDurham20->Draw("HIST");
    		TCanvas *c359 = new TCanvas();
    		c359->cd();
    		histNJetsDurham25->Draw("HIST");
    		TCanvas *c360 = new TCanvas();
    		c360->cd();
    		histNJetsDurham30->Draw("HIST");
    		////nJetsDurham*/
    		
    		/*///////minChiSquaredZZMass
    		TCanvas *c361 = new TCanvas();
    		c361->cd();
    		histMinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c362 = new TCanvas();
    		c362->cd();
    		histInvMassZZ1->Draw("HIST");
    		TLine* distanceZMassLeft = new TLine(90-distanceZMass, histInvMassZZ1->GetMinimum(), 90-distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassLeft->SetLineColor(kRed);
	        distanceZMassLeft->SetLineWidth(1);
	        TLine* distanceZMassRight = new TLine(90+distanceZMass, histInvMassZZ1->GetMinimum(), 90+distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassRight->SetLineColor(kRed);
	        distanceZMassRight->SetLineWidth(1);
	        distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
    		TCanvas *c363 = new TCanvas();
    		c363->cd();
    		histInvMassZZ2->Draw("HIST");
    		distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
	        TCanvas *c364 = new TCanvas();
    		c364->cd();
    		histDistanceZ1MinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c365 = new TCanvas();
    		c365->cd();
    		histDistanceZ2MinChiSquaredZZMass->Draw("HIST");
    		/////minChiSquaredZZMass*/
    		
    		/*////////y_nm
    		TCanvas *c266 = new TCanvas();
    		c266->cd();
    		histExclYmerge12->Draw("HIST");
    		TCanvas *c267 = new TCanvas();
    		c267->cd();
    		histExclYmerge23->Draw("HIST");
    		TCanvas *c268 = new TCanvas();
    		c268->cd();
    		histExclYmerge34->Draw("HIST");
    		TCanvas *c269 = new TCanvas();
    		c269->cd();
    		histExclYmerge45->Draw("HIST");
    		TCanvas *c270 = new TCanvas();
    		c270->cd();
    		histExclYmerge56->Draw("HIST");
    		///////y_nm*/
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
    		
    		/*////nJetsDurham
    		TCanvas *c454 = new TCanvas();
    		c454->cd();
    		histNJetsDurham0->Draw("HIST");
    		TCanvas *c455 = new TCanvas();
    		c455->cd();
    		histNJetsDurham5->Draw("HIST");
    		TCanvas *c456 = new TCanvas();
    		c456->cd();
    		histNJetsDurham10->Draw("HIST");
    		TCanvas *c457 = new TCanvas();
    		c457->cd();
    		histNJetsDurham15->Draw("HIST");
    		TCanvas *c458 = new TCanvas();
    		c458->cd();
    		histNJetsDurham20->Draw("HIST");
    		TCanvas *c459 = new TCanvas();
    		c459->cd();
    		histNJetsDurham25->Draw("HIST");
    		TCanvas *c460 = new TCanvas();
    		c460->cd();
    		histNJetsDurham30->Draw("HIST");
    		////nJetsDurham*/
    		
    		/*///////minChiSquaredZZMass
    		TCanvas *c461 = new TCanvas();
    		c461->cd();
    		histMinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c462 = new TCanvas();
    		c462->cd();
    		histInvMassZZ1->Draw("HIST");
    		TLine* distanceZMassLeft = new TLine(90-distanceZMass, histInvMassZZ1->GetMinimum(), 90-distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassLeft->SetLineColor(kRed);
	        distanceZMassLeft->SetLineWidth(1);
	        TLine* distanceZMassRight = new TLine(90+distanceZMass, histInvMassZZ1->GetMinimum(), 90+distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassRight->SetLineColor(kRed);
	        distanceZMassRight->SetLineWidth(1);
	        distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
    		TCanvas *c463 = new TCanvas();
    		c463->cd();
    		histInvMassZZ2->Draw("HIST");
    		distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
	        TCanvas *c464 = new TCanvas();
    		c464->cd();
    		histDistanceZ1MinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c465 = new TCanvas();
    		c465->cd();
    		histDistanceZ2MinChiSquaredZZMass->Draw("HIST");
    		/////minChiSquaredZZMass*/
    		
    		/*/////////y_nm
    		TCanvas *c466 = new TCanvas();
    		c466->cd();
    		histExclYmerge12->Draw("HIST");
    		TCanvas *c467 = new TCanvas();
    		c467->cd();
    		histExclYmerge23->Draw("HIST");
    		TCanvas *c468 = new TCanvas();
    		c468->cd();
    		histExclYmerge34->Draw("HIST");
    		TCanvas *c469 = new TCanvas();
    		c469->cd();
    		histExclYmerge45->Draw("HIST");
    		TCanvas *c470 = new TCanvas();
    		c470->cd();
    		histExclYmerge56->Draw("HIST");
    		///////y_nm*/
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
    		
    		/*////nJetsDurham
    		TCanvas *c554 = new TCanvas();
    		c554->cd();
    		histNJetsDurham0->Draw("HIST");
    		TCanvas *c555 = new TCanvas();
    		c555->cd();
    		histNJetsDurham5->Draw("HIST");
    		TCanvas *c556 = new TCanvas();
    		c556->cd();
    		histNJetsDurham10->Draw("HIST");
    		TCanvas *c557 = new TCanvas();
    		c557->cd();
    		histNJetsDurham15->Draw("HIST");
    		TCanvas *c558 = new TCanvas();
    		c558->cd();
    		histNJetsDurham20->Draw("HIST");
    		TCanvas *c559 = new TCanvas();
    		c559->cd();
    		histNJetsDurham25->Draw("HIST");
    		TCanvas *c560 = new TCanvas();
    		c560->cd();
    		histNJetsDurham30->Draw("HIST");
    		////nJetsDurham*/
    		
    		/*///////minChiSquaredZZMass
    		TCanvas *c561 = new TCanvas();
    		c561->cd();
    		histMinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c562 = new TCanvas();
    		c562->cd();
    		histInvMassZZ1->Draw("HIST");
    		TLine* distanceZMassLeft = new TLine(90-distanceZMass, histInvMassZZ1->GetMinimum(), 90-distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassLeft->SetLineColor(kRed);
	        distanceZMassLeft->SetLineWidth(1);
	        TLine* distanceZMassRight = new TLine(90+distanceZMass, histInvMassZZ1->GetMinimum(), 90+distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassRight->SetLineColor(kRed);
	        distanceZMassRight->SetLineWidth(1);
	        distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
    		TCanvas *c563 = new TCanvas();
    		c563->cd();
    		histInvMassZZ2->Draw("HIST");
    		distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
	        TCanvas *c564 = new TCanvas();
    		c564->cd();
    		histDistanceZ1MinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c565 = new TCanvas();
    		c565->cd();
    		histDistanceZ2MinChiSquaredZZMass->Draw("HIST");
    		/////minChiSquaredZZMass*/
    		
    		/*////////y_nm
    		TCanvas *c566 = new TCanvas();
    		c566->cd();
    		histExclYmerge12->Draw("HIST");
    		TCanvas *c567 = new TCanvas();
    		c567->cd();
    		histExclYmerge23->Draw("HIST");
    		TCanvas *c568 = new TCanvas();
    		c568->cd();
    		histExclYmerge34->Draw("HIST");
    		TCanvas *c569 = new TCanvas();
    		c569->cd();
    		histExclYmerge45->Draw("HIST");
    		TCanvas *c570 = new TCanvas();
    		c570->cd();
    		histExclYmerge56->Draw("HIST");
    		///////y_nm*/
	}
	else if(topology == 9)
	{ 
		/*//////eta
		TCanvas *c91 = new TCanvas();
		TCanvas *c92 = new TCanvas();
		TCanvas *c93 = new TCanvas();
		TCanvas *c94 = new TCanvas();
		TCanvas *c95 = new TCanvas();			
		c91->cd();
    		histJetEta->Draw("HIST");
    		c92->cd();
    		histJet1Eta->Draw("HIST");				
    		c93->cd();
    		histJet2Eta->Draw("HIST");				
    		c94->cd();
    		histJet3Eta->Draw("HIST");				
    		c95->cd();
    		histJet4Eta->Draw("HIST");
    		//////eta*/
    		
    		/*//////cosTheta
		TCanvas *c96 = new TCanvas();
		TCanvas *c97 = new TCanvas();
		TCanvas *c98 = new TCanvas();
		TCanvas *c99 = new TCanvas();
		TCanvas *c910 = new TCanvas();			
		c96->cd();
    		histJetCosTheta->Draw("HIST");
    		c97->cd();
    		histJet1CosTheta->Draw("HIST");				
    		c98->cd();
    		histJet2CosTheta->Draw("HIST");				
    		c99->cd();
    		histJet3CosTheta->Draw("HIST");				
    		c910->cd();
    		histJet4CosTheta->Draw("HIST");
    		//////cosTheta*/
    		
    		/*/////Pt
    		c1001->cd();
    		histJet4Pt->Draw("HIST");
    		histJet1Pt->Draw("HIST same");
    		histJet2Pt->Draw("HIST same");
    		histJet3Pt->Draw("HIST same");
    		TCanvas *c911 = new TCanvas();
		TCanvas *c912 = new TCanvas();	
		c911->cd();
    		histSumJetPt->Draw("HIST");
    		c912->cd();
    		histJetPt->Draw("HIST");	
    		/////Pt
    		
    		///////inv. masses
    		TCanvas *c913 = new TCanvas();
		TCanvas *c914 = new TCanvas();
		TCanvas *c915 = new TCanvas();
		TCanvas *c916 = new TCanvas();
		TCanvas *c917 = new TCanvas();
		TCanvas *c918 = new TCanvas();
		TCanvas *c919 = new TCanvas();
		c913->cd();
		histJetB1M->Draw("HIST");
		c914->cd();
		histJetB2M->Draw("HIST");
		c915->cd();
		histMinJetM->Draw("HIST");
		c916->cd();
		histJet1M->Draw("HIST");
		c917->cd();
		histJet2M->Draw("HIST");
		c918->cd();
		histJet3M->Draw("HIST");
		c919->cd();
		histJet4M->Draw("HIST");
    		//////inv. masses
    		
    		/////Event shape
    		TCanvas *c920 = new TCanvas();
		TCanvas *c921 = new TCanvas();
		c920->cd();
		histAplan->Draw("HIST");
		c921->cd();
		histSpher->Draw("HIST");
    		/////Event shape*/
    		
    		/*/////Jet constituents
    		TCanvas *c922 = new TCanvas();
    		TCanvas *c923 = new TCanvas();
    		TCanvas *c924 = new TCanvas();
    		TCanvas *c925 = new TCanvas();
    		TCanvas *c926 = new TCanvas();
    		TCanvas *c927 = new TCanvas();
    		TCanvas *c928 = new TCanvas();
    		c922->cd();
    		histNParticles->Draw("HIST");
    		c923->cd();
    		histTotalConstSize->Draw("HIST");
    		c924->cd();
    		histConstSizeB1->Draw("HIST");
    		c925->cd();
    		histConstSizeB2->Draw("HIST");
    		c926->cd();
    		histConstSizeB3->Draw("HIST");
    		c927->cd();
    		histConstSizeB4->Draw("HIST");
    		c928->cd();
    		histMinConstSize->Draw("HIST");
    		/////*/
    		
    		/*/////PFOs
    		TCanvas *c929 = new TCanvas();
    		c929->cd();
    		histNEFlowTracks->Draw("HIST");
    		TCanvas *c930 = new TCanvas();
    		c930->cd();
    		histNEFlowPhotons->Draw("HIST");
    		TCanvas *c931 = new TCanvas();
    		c931->cd();
    		histNEFlowNeutralHadrons->Draw("HIST");
    		TCanvas *c932 = new TCanvas();
    		c932->cd();
    		histNEFlowObjects->Draw("HIST");
    		TCanvas *c933 = new TCanvas();
    		c933->cd();
    		histJetB1NCharged->Draw("HIST");
    		TCanvas *c934 = new TCanvas();
    		c934->cd();
    		histJetB2NCharged->Draw("HIST");
    		TCanvas *c935 = new TCanvas();
    		c935->cd();
    		histJetB3NCharged->Draw("HIST");
    		TCanvas *c936 = new TCanvas();
    		c936->cd();
    		histJetB4NCharged->Draw("HIST");
    		TCanvas *c937 = new TCanvas();
    		c937->cd();
    		histJetB1NNeutrals->Draw("HIST");
    		TCanvas *c938 = new TCanvas();
    		c938->cd();
    		histJetB2NNeutrals->Draw("HIST");
    		TCanvas *c939 = new TCanvas();
    		c939->cd();
    		histJetB3NNeutrals->Draw("HIST");
    		TCanvas *c940 = new TCanvas();
    		c940->cd();
    		histJetB4NNeutrals->Draw("HIST");
    		TCanvas *c941 = new TCanvas();
    		c941->cd();
    		histJetNObjects->Draw("HIST");
    		TCanvas *c942 = new TCanvas();
    		c942->cd();
    		histMinJetNObjects->Draw("HIST");
    		/////PFOs*/
    		
    		/*///////nJetsCompareAlgos
    		TCanvas *c943 = new TCanvas();
    		c943->cd();
    		histNJetsCompareAlgos->Draw("TEXT");
    		///////nJetsCompareAlgos*/
    		
    		/*///////invMass for antiKt pairs when 4 Durham jets
    		TCanvas *c944 = new TCanvas();
    		c944->cd();
    		histJetB1MAntiKt2Jets->Draw("HIST");
    		TCanvas *c945 = new TCanvas();
    		c945->cd();
    		histJetB2MAntiKt2Jets->Draw("HIST");
    		TCanvas *c946 = new TCanvas();
    		c946->cd();
    		histJetB1MAntiKt3Jets->Draw("HIST");
    		TCanvas *c947 = new TCanvas();
    		c947->cd();
    		histJetB2MAntiKt3Jets->Draw("HIST");
    		TCanvas *c948 = new TCanvas();
    		c948->cd();
    		histJetB1MAntiKt4Jets->Draw("HIST");
    		TCanvas *c949 = new TCanvas();
    		c949->cd();
    		histJetB2MAntiKt4Jets->Draw("HIST");
    		TCanvas *c950 = new TCanvas();
    		c950->cd();
    		histJetB1MAntiKt5Jets->Draw("HIST");
    		TCanvas *c951 = new TCanvas();
    		c951->cd();
    		histJetB2MAntiKt5Jets->Draw("HIST");
    		TCanvas *c952 = new TCanvas();
    		c952->cd();
    		histJetB1MAntiKt6Jets->Draw("HIST");
    		TCanvas *c953 = new TCanvas();
    		c953->cd();
    		histJetB2MAntiKt6Jets->Draw("HIST");
    		///////invMass for antiKt pairs when 4 Durham jets*/
    		
    		/*////nJetsDurham
    		TCanvas *c954 = new TCanvas();
    		c954->cd();
    		histNJetsDurham0->Draw("HIST");
    		TCanvas *c955 = new TCanvas();
    		c955->cd();
    		histNJetsDurham5->Draw("HIST");
    		TCanvas *c956 = new TCanvas();
    		c956->cd();
    		histNJetsDurham10->Draw("HIST");
    		TCanvas *c957 = new TCanvas();
    		c957->cd();
    		histNJetsDurham15->Draw("HIST");
    		TCanvas *c958 = new TCanvas();
    		c958->cd();
    		histNJetsDurham20->Draw("HIST");
    		TCanvas *c959 = new TCanvas();
    		c959->cd();
    		histNJetsDurham25->Draw("HIST");
    		TCanvas *c960 = new TCanvas();
    		c960->cd();
    		histNJetsDurham30->Draw("HIST");
    		////nJetsDurham*/
    		
    		/*///////minChiSquaredZZMass
    		TCanvas *c961 = new TCanvas();
    		c961->cd();
    		histMinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c962 = new TCanvas();
    		c962->cd();
    		histInvMassZZ1->Draw("HIST");
    		TLine* distanceZMassLeft = new TLine(90-distanceZMass, histInvMassZZ1->GetMinimum(), 90-distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassLeft->SetLineColor(kRed);
	        distanceZMassLeft->SetLineWidth(1);
	        TLine* distanceZMassRight = new TLine(90+distanceZMass, histInvMassZZ1->GetMinimum(), 90+distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassRight->SetLineColor(kRed);
	        distanceZMassRight->SetLineWidth(1);
	        distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
    		TCanvas *c963 = new TCanvas();
    		c963->cd();
    		histInvMassZZ2->Draw("HIST");
    		distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
	        TCanvas *c964 = new TCanvas();
    		c964->cd();
    		histDistanceZ1MinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c965 = new TCanvas();
    		c965->cd();
    		histDistanceZ2MinChiSquaredZZMass->Draw("HIST");
    		/////minChiSquaredZZMass*/

    		/*////////y_nm
    		TCanvas *c966 = new TCanvas();
    		c966->cd();
    		histExclYmerge12->Draw("HIST");
    		TCanvas *c967 = new TCanvas();
    		c967->cd();
    		histExclYmerge23->Draw("HIST");
    		TCanvas *c968 = new TCanvas();
    		c968->cd();
    		histExclYmerge34->Draw("HIST");
    		TCanvas *c969 = new TCanvas();
    		c969->cd();
    		histExclYmerge45->Draw("HIST");
    		TCanvas *c970 = new TCanvas();
    		c970->cd();
    		histExclYmerge56->Draw("HIST");
    		///////y_nm*/

			/*/////////inv. mass for best comb
			TCanvas *c971 = new TCanvas();
			TCanvas *c972 = new TCanvas();
			c971->cd();
			histJetB1M1Best->Draw("HIST");
			c972->cd();
			histJetB2M1Best->Draw("HIST");
			/////////inv. mass for best comb*/

			///////minChiSquaredZHMass
    		/*TCanvas *c973 = new TCanvas();
    		c973->cd();
    		histMinChiSquaredZHMass->Draw("HIST");
    		TCanvas *c974 = new TCanvas();
    		c974->cd();
    		histInvMassZH1->Draw("HIST");
    		TLine* distanceZMassLeft = new TLine(90-distanceZMass, histInvMassZZ1->GetMinimum(), 90-distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassLeft->SetLineColor(kRed);
	        distanceZMassLeft->SetLineWidth(1);
	        TLine* distanceZMassRight = new TLine(90+distanceZMass, histInvMassZZ1->GetMinimum(), 90+distanceZMass, histInvMassZZ1->GetMaximum());
	        distanceZMassRight->SetLineColor(kRed);
	        distanceZMassRight->SetLineWidth(1);
	        distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
    		TCanvas *c975 = new TCanvas();
    		c975->cd();
    		histInvMassZH2->Draw("HIST");
    		distanceZMassLeft->Draw("same");
	        distanceZMassRight->Draw("same");
	        TCanvas *c976 = new TCanvas();
    		c976->cd();
    		histDistanceZ1MinChiSquaredZZMass->Draw("HIST");
    		TCanvas *c977 = new TCanvas();
    		c977->cd();
    		histDistanceZ2MinChiSquaredZZMass->Draw("HIST");*/
    		/////minChiSquaredZHMass
	}
	if(topology == 12)
	{
		/*/////Event shape
		TCanvas *c0121 = new TCanvas();
		c0121->cd();
		histAplan->Draw("HIST");
		//////Event shape*/

		/*///////////inv. masses
		TCanvas *c0122 = new TCanvas();
		c0122->cd();
		histJetB1M->Draw("HIST");
		TCanvas *c0123 = new TCanvas();
		c0123->cd();
		histJetB2M->Draw("HIST");
		TCanvas *c0124 = new TCanvas();
		c0124->cd();
		histJetB1M1Best->Draw("HIST");
		TCanvas *c0125 = new TCanvas();
		c0125->cd();
		histJetB2M1Best->Draw("HIST");
		//////Inv. masses*/

		/*/////qqH study
		TCanvas *c0126 = new TCanvas("c0126", "qqH Study Plots", 1200, 800);
		c0126->Divide(3, 2);

		c0126->cd(1);
		histInvMass4Jets->Draw("HIST");

		c0126->cd(2);
		histDeltaRJetPairs->Draw("HIST");

		c0126->cd(3);
		histEAssymetry->Draw("HIST");

		c0126->cd(4);
		histPTAssymetry->Draw("HIST");

		c0126->cd(5);
		histJetPairMassDelta->Draw("HIST");
		/////qqH study*/
		
		/*/////N-subjetiness
		TCanvas *c0127canvasNSubjetiness = new TCanvas("c0127canvasNSubjetiness", "N-Subjetiness Plots", 1200, 900);
		c0127canvasNSubjetiness->Divide(3, 4); // Dividing into 3 columns and 4 rows for 12 histograms
		
		c0127canvasNSubjetiness->cd(1);
		histNSubJetsPrunnedB1->Draw("HIST");
		c0127canvasNSubjetiness->cd(2);
		histNSubJetsPrunnedB2->Draw("HIST");
		c0127canvasNSubjetiness->cd(3);
		histNSubJetsPrunnedB3->Draw("HIST");
		c0127canvasNSubjetiness->cd(4);
		histNSubJetsPrunnedB4->Draw("HIST");
		c0127canvasNSubjetiness->cd(5);
		histNSubJetsSoftDroppedB1->Draw("HIST");
		c0127canvasNSubjetiness->cd(6);
		histNSubJetsSoftDroppedB2->Draw("HIST");
		c0127canvasNSubjetiness->cd(7);
		histNSubJetsSoftDroppedB3->Draw("HIST");
		c0127canvasNSubjetiness->cd(8);
		histNSubJetsSoftDroppedB4->Draw("HIST");
		c0127canvasNSubjetiness->cd(9);
		histNSubJetsSoftDroppedB1->Draw("HIST");
		c0127canvasNSubjetiness->cd(10);
		histNSubJetsSoftDroppedB2->Draw("HIST");
		c0127canvasNSubjetiness->cd(11);
		histNSubJetsSoftDroppedB3->Draw("HIST");
		c0127canvasNSubjetiness->cd(12);
		histNSubJetsSoftDroppedB4->Draw("HIST");
		///////N-subjetiness*/
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
	if(topology == 6) cout<<"For qqX: "<<endl;
	if(topology == 7) cout<<"For qqqqX: "<<endl;
	if(topology == 8) cout<<"For qqHX: "<<endl;
	if(topology == 9) cout<<"For ZH: "<<endl;
	if(topology == 10) cout<<"For pebb: "<<endl;
	if(topology == 11) cout<<"For pebbqq: "<<endl;
	if(topology == 12) cout<<"For peqqH: "<<endl;
	if(topology == 13) cout<<"For pett: "<<endl;
	cout<<"Events that have 4b: "<<contEvents<<"    Weighted: "<<contEvents*weight<<endl;
	cout<<"Events that pass: "<<contEventsPostFilter<<"    Weighted: "<<contEventsPostFilter*weight<<endl<<endl;

	//cout<<"Events 4 jets: "<<contEvents4Jets<<endl<<endl<<endl;
	/*cout<<"Events in 0.1 Z mass window: "<<contEventsZWindow01<<endl;
	cout<<"Events in 0.5 Z mass window: "<<contEventsZWindow05<<endl;
	cout<<"Events in 1.5 Z mass window: "<<contEventsZWindow15<<endl;
	cout<<"Events in 2 Z mass window: "<<contEventsZWindow2<<endl;
	cout<<"Events in 5 Z mass window: "<<contEventsZWindow5<<endl;
	cout<<"Events in 10 Z mass window: "<<contEventsZWindow10<<endl<<endl<<endl;*/
	
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

/*void trueBPairMass(const char *inputFile, int topology, float weight)
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
  	branchParticle = treeReader->UseBranch("Particle");
	branchEvent = treeReader->UseBranch("Event");
  	
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
}*/

//function that fills out to sets, one for trainign and one for teting, with the entry index of the events in the trees for training and testing.
void generateSetsMerge(int topology, set<int>& setTrain, set<int>& setTest, string rtdCut)
{
	TFile *fileTrain, *fileTest;
	TTree *treeTrain, *treeTest;
	if(topology == 1)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeSHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeSHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeSTrain", treeTrain);
		fileTest->GetObject("TreeSTest", treeTest);
	}
	if(topology == 2)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBqqTrain", treeTrain);
		fileTest->GetObject("TreeBqqTest", treeTest);
	}
	if(topology == 3)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBttHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBttHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBttTrain", treeTrain);
		fileTest->GetObject("TreeBttTest", treeTest);
	}
	if(topology == 4)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBZZTrain", treeTrain);
		fileTest->GetObject("TreeBZZTest", treeTest);
	}
	if(topology == 5)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBWWTrain", treeTrain);
		fileTest->GetObject("TreeBWWTest", treeTest);
	}
	if(topology == 6)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBqqXTrain", treeTrain);
		fileTest->GetObject("TreeBqqXTest", treeTest);
	}
	if(topology == 7)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBqqqqXTrain", treeTrain);
		fileTest->GetObject("TreeBqqqqXTest", treeTest);
	}
	if(topology == 8)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBqqHXTrain", treeTrain);
		fileTest->GetObject("TreeBqqHXTest", treeTest);
	}
	if(topology == 9)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBZHHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBZHHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBZHTrain", treeTrain);
		fileTest->GetObject("TreeBZHTest", treeTest);
	}
	if(topology == 10)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBpebbHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBpebbHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBpebbTrain", treeTrain);
		fileTest->GetObject("TreeBpebbTest", treeTest);
	}
	if(topology == 11)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBpebbqqHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBpebbqqHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBpebbqqTrain", treeTrain);
		fileTest->GetObject("TreeBpebbqqTest", treeTest);
	}
	if(topology == 12)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBpeqqHHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBpeqqHHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBpeqqHTrain", treeTrain);
		fileTest->GetObject("TreeBpeqqHTest", treeTest);
	}
	if(topology == 13)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBpettHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBpettHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBpettTrain", treeTrain);
		fileTest->GetObject("TreeBpettTest", treeTest);
	}
	
	
	cout<<"train file: "<<fileTrain->GetName()<<endl;
	cout<<"test file: "<<fileTest->GetName()<<endl;
	Float_t entryIndexTrain, entryIndexTest;
	treeTrain->SetBranchAddress("entryIndex", &entryIndexTrain);
	treeTest->SetBranchAddress("entryIndex", &entryIndexTest);
	
	Long64_t nEntries = treeTrain->GetEntries();
	for (Long64_t i = 0; i < nEntries; ++i) {
	 	treeTrain->GetEntry(i);
	 	setTrain.insert(entryIndexTrain);
	 }
	nEntries = treeTest->GetEntries();
	for (Long64_t i = 0; i < nEntries; ++i) {
	 	treeTest->GetEntry(i);
	 	setTest.insert(entryIndexTest);
	 }
	 
	/*for(int i=0; i<100; i++)
	{
		if(setTrain.find(i) != setTrain.end()) cout<<"Entry "<<i<<" is in training set"<<endl;
		else if(setTest.find(i) != setTest.end()) cout<<"Entry "<<i<<" is in testing set"<<endl; 
		else cout<<"entry "<<i<<" does not pass the filter"<<endl;
	}*/
}

////function that receives 43trees, two original for training and 2 for testing and the new one with all the new vars., and merges the respective pairs. It returns the two original trees but now will all the info. from the 3.
void mergeTrees(int topology, set<int> setTrain, set<int> setTest, TTree& TreeTrainNew, TTree& TreeTestNew, TTree& TreeMerge, string rtdCut)
{
	TFile *fileTrain, *fileTest;
	TTree *TreeTrain, *TreeTest;
	if(topology == 1)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeSHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeSHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeSTrain", TreeTrain);
		fileTest->GetObject("TreeSTest", TreeTest);
	}
	if(topology == 2)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBqqTrain", TreeTrain);
		fileTest->GetObject("TreeBqqTest", TreeTest);
	}
	if(topology == 3)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBttHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBttHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBttTrain", TreeTrain);
		fileTest->GetObject("TreeBttTest", TreeTest);
	}
	if(topology == 4)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBZZTrain", TreeTrain);
		fileTest->GetObject("TreeBZZTest", TreeTest);
	}
	if(topology == 5)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBWWTrain", TreeTrain);
		fileTest->GetObject("TreeBWWTest", TreeTest);
	}
	if(topology == 6)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBqqXTrain", TreeTrain);
		fileTest->GetObject("TreeBqqXTest", TreeTest);
	}
	if(topology == 7)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBqqqqXTrain", TreeTrain);
		fileTest->GetObject("TreeBqqqqXTest", TreeTest);
	}
	if(topology == 8)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBqqHXTrain", TreeTrain);
		fileTest->GetObject("TreeBqqHXTest", TreeTest);
	}
	if(topology == 9)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBZHHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBZHHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBZHTrain", TreeTrain);
		fileTest->GetObject("TreeBZHTest", TreeTest);
	}

	if(topology == 10)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBpebbHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBpebbHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBpebbTrain", TreeTrain);
		fileTest->GetObject("TreeBpebbTest", TreeTest);
	}

	if(topology == 11)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBpebbqqHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBpebbqqHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBpebbqqTrain", TreeTrain);
		fileTest->GetObject("TreeBpebbqqTest", TreeTest);
	}

	if(topology == 12)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBpeqqHHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBpeqqHHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBpeqqHTrain", TreeTrain);
		fileTest->GetObject("TreeBpeqqHTest", TreeTest);
	}

	if(topology == 13)
	{
		string fileTrainText = "analysis/SampleOG/outputTreeBpettHHbbbbESpreadDurham"+rtdCut+"Train.root";
		string fileTestText = "analysis/SampleOG/outputTreeBpettHHbbbbESpreadDurham"+rtdCut+"Test.root";
		fileTrain = TFile::Open(fileTrainText.c_str());
		fileTest = TFile::Open(fileTestText.c_str());
		fileTrain->GetObject("TreeBpettTrain", TreeTrain);
		fileTest->GetObject("TreeBpettTest", TreeTest);
	}
	
	float aplanarity, invMassB1, invMassB2, minJetM, sphericity, cosThetaB1, cosThetaB2, cosThetaB3, cosThetaB4, sumPt, jetB1Pt, jetB2Pt, jetB3Pt, jetB4Pt, jetB1M, jetB2M, jetB3M, jetB4M, etaB1, etaB2, etaB3, etaB4, nParticles, totalConstSize, constSizeB1, constSizeB2, constSizeB3, constSizeB4, minConstSize, jetB1NCharged, jetB2NCharged, jetB3NCharged, jetB4NCharged, jetB1NNeutrals, jetB2NNeutrals, jetB3NNeutrals, jetB4NNeutrals, jetNObjects, minJetNObjects, invMassB1AntiKt, invMassB2AntiKt, invMassB1AntiKt2Jets, invMassB2AntiKt2Jets, invMassB1AntiKt3Jets, invMassB2AntiKt3Jets, invMassB1AntiKt4Jets, invMassB2AntiKt4Jets, invMassB1AntiKt5Jets, invMassB2AntiKt5Jets, invMassB1AntiKt6Jets, invMassB2AntiKt6Jets, nJetsAntiKt, invMassB11Best, invMassB21Best, invMassB12Best, invMassB22Best, invMassB13Best, invMassB23Best, invMassB14Best, invMassB24Best, invMassB15Best, invMassB25Best, invMassB16Best, invMassB26Best, invMassB17Best, invMassB27Best, invMassB18Best, invMassB28Best, entryIndex;
	float entryIndexMerge, nJetsDurham0, nJetsDurham5, nJetsDurham10, nJetsDurham15, nJetsDurham20, nJetsDurham25, nJetsDurham30, minChiSquaredZZMass, distanceZ1MinChiSquaredZZMass, distanceZ2MinChiSquaredZZMass, exclYmerge12, exclYmerge23, exclYmerge34, exclYmerge45, exclYmerge56;
	float entryIndexMergeHolder, nJetsDurham0Holder, nJetsDurham5Holder, nJetsDurham10Holder, nJetsDurham15Holder, nJetsDurham20Holder, nJetsDurham25Holder, nJetsDurham30Holder, minChiSquaredZZMassHolder, distanceZ1MinChiSquaredZZMassHolder, distanceZ2MinChiSquaredZZMassHolder, exclYmerge12Holder, exclYmerge23Holder, exclYmerge34Holder, exclYmerge45Holder, exclYmerge56Holder;
	
	
	TreeTrainNew.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeTrainNew.Branch("aplanarity",&aplanarity,"aplanarity/F");
	TreeTrainNew.Branch("invMassB1",&invMassB1,"invMassB1/F");
	TreeTrainNew.Branch("invMassB2",&invMassB2,"invMassB2/F");
	TreeTrainNew.Branch("minJetM",&minJetM,"minJetM/F");  
	TreeTrainNew.Branch("sphericity",&sphericity,"sphericity/F");
	TreeTrainNew.Branch("cosThetaB1",&cosThetaB1,"cosThetaB1/F");
	TreeTrainNew.Branch("cosThetaB2",&cosThetaB2,"cosThetaB2/F");
	TreeTrainNew.Branch("cosThetaB3",&cosThetaB3,"cosThetaB3/F");
	TreeTrainNew.Branch("cosThetaB4",&cosThetaB4,"cosThetaB4/F");
	TreeTrainNew.Branch("sumPt",&sumPt,"sumPt/F");
	TreeTrainNew.Branch("jetB1Pt",&jetB1Pt,"jetB1Pt/F");
	TreeTrainNew.Branch("jetB2Pt",&jetB2Pt,"jetB2Pt/F");
	TreeTrainNew.Branch("jetB3Pt",&jetB3Pt,"jetB3Pt/F");
	TreeTrainNew.Branch("jetB4Pt",&jetB4Pt,"jetB4Pt/F");
	TreeTrainNew.Branch("jetB1M",&jetB1M,"jetB1M/F");
	TreeTrainNew.Branch("jetB2M",&jetB2M,"jetB2M/F");
	TreeTrainNew.Branch("jetB3M",&jetB3M,"jetB3M/F");
	TreeTrainNew.Branch("jetB4M",&jetB4M,"jetB4M/F");
	TreeTrainNew.Branch("nParticles",&nParticles,"nParticles/F");
	TreeTrainNew.Branch("constSizeB1",&constSizeB1,"constSizeB1/F");
	TreeTrainNew.Branch("constSizeB2",&constSizeB2,"constSizeB2/F");
	TreeTrainNew.Branch("constSizeB3",&constSizeB3,"constSizeB3/F");
	TreeTrainNew.Branch("constSizeB4",&constSizeB4,"constSizeB4/F");
	TreeTrainNew.Branch("minConstSize",&minConstSize,"minConstSize/F");
	TreeTrainNew.Branch("jetNObjects",&jetNObjects,"jetNObjects/F");
	TreeTrainNew.Branch("minJetNObjects",&minJetNObjects,"minJetNObjects/F");
	TreeTrainNew.Branch("invMassB1AntiKt",&invMassB1AntiKt,"invMassB1AntiKt/F");
	TreeTrainNew.Branch("invMassB2AntiKt",&invMassB2AntiKt,"invMassB2AntiKt/F");
	TreeTrainNew.Branch("nJetsAntiKt",&nJetsAntiKt,"nJetsAntiKt/F");
	TreeTrainNew.Branch("invMassB11Best",&invMassB11Best,"invMassB11Best/F");
	TreeTrainNew.Branch("invMassB21Best",&invMassB21Best,"invMassB21Best/F");
	TreeTrainNew.Branch("invMassB12Best",&invMassB12Best,"invMassB12Best/F");
	TreeTrainNew.Branch("invMassB22Best",&invMassB22Best,"invMassB22Best/F");
	TreeTrainNew.Branch("invMassB13Best",&invMassB13Best,"invMassB13Best/F");
	TreeTrainNew.Branch("invMassB23Best",&invMassB23Best,"invMassB23Best/F");
	TreeTrainNew.Branch("invMassB14Best",&invMassB14Best,"invMassB14Best/F");
	TreeTrainNew.Branch("invMassB24Best",&invMassB24Best,"invMassB24Best/F");
	TreeTrainNew.Branch("invMassB15Best",&invMassB15Best,"invMassB15Best/F");
	TreeTrainNew.Branch("invMassB25Best",&invMassB25Best,"invMassB25Best/F");
	TreeTrainNew.Branch("invMassB16Best",&invMassB16Best,"invMassB16Best/F");
	TreeTrainNew.Branch("invMassB26Best",&invMassB26Best,"invMassB26Best/F");
	TreeTrainNew.Branch("invMassB17Best",&invMassB17Best,"invMassB17Best/F");
	TreeTrainNew.Branch("invMassB27Best",&invMassB27Best,"invMassB27Best/F");
	TreeTrainNew.Branch("invMassB18Best",&invMassB18Best,"invMassB18Best/F");
	TreeTrainNew.Branch("invMassB28Best",&invMassB28Best,"invMassB28Best/F");
	//TreeTrainNew.Branch("nJetsDurham0",&nJetsDurham0,"nJetsDurham0/F");
	//TreeTrainNew.Branch("nJetsDurham5",&nJetsDurham5,"nJetsDurham5/F");
	//TreeTrainNew.Branch("nJetsDurham10",&nJetsDurham10,"nJetsDurham10/F");
	TreeTrainNew.Branch("nJetsDurham15",&nJetsDurham15,"nJetsDurham15/F");
	TreeTrainNew.Branch("nJetsDurham20",&nJetsDurham20,"nJetsDurham20/F");
	TreeTrainNew.Branch("nJetsDurham25",&nJetsDurham25,"nJetsDurham25/F");
	TreeTrainNew.Branch("nJetsDurham30",&nJetsDurham30,"nJetsDurham30/F");
	TreeTrainNew.Branch("distanceZ1MinChiSquaredZZMass",&distanceZ1MinChiSquaredZZMass,"distanceZ1MinChiSquaredZZMass/F");
	TreeTrainNew.Branch("distanceZ2MinChiSquaredZZMass",&distanceZ2MinChiSquaredZZMass,"distanceZ2MinChiSquaredZZMass/F");
	TreeTrainNew.Branch("exclYmerge12",&exclYmerge12,"exclYmerge12/F");
	TreeTrainNew.Branch("exclYmerge23",&exclYmerge23,"exclYmerge23/F");
	TreeTrainNew.Branch("exclYmerge34",&exclYmerge34,"exclYmerge34/F");
	TreeTrainNew.Branch("exclYmerge45",&exclYmerge45,"exclYmerge45/F");
	TreeTrainNew.Branch("exclYmerge56",&exclYmerge56,"exclYmerge56/F");
	
	TreeTestNew.Branch("entryIndex",&entryIndex,"entryIndex/F");
	TreeTestNew.Branch("aplanarity",&aplanarity,"aplanarity/F");
	TreeTestNew.Branch("invMassB1",&invMassB1,"invMassB1/F");
	TreeTestNew.Branch("invMassB2",&invMassB2,"invMassB2/F");
	TreeTestNew.Branch("minJetM",&minJetM,"minJetM/F");  
	TreeTestNew.Branch("sphericity",&sphericity,"sphericity/F");
	TreeTestNew.Branch("cosThetaB1",&cosThetaB1,"cosThetaB1/F");
	TreeTestNew.Branch("cosThetaB2",&cosThetaB2,"cosThetaB2/F");
	TreeTestNew.Branch("cosThetaB3",&cosThetaB3,"cosThetaB3/F");
	TreeTestNew.Branch("cosThetaB4",&cosThetaB4,"cosThetaB4/F");
	TreeTestNew.Branch("sumPt",&sumPt,"sumPt/F");
	TreeTestNew.Branch("jetB1Pt",&jetB1Pt,"jetB1Pt/F");
	TreeTestNew.Branch("jetB2Pt",&jetB2Pt,"jetB2Pt/F");
	TreeTestNew.Branch("jetB3Pt",&jetB3Pt,"jetB3Pt/F");
	TreeTestNew.Branch("jetB4Pt",&jetB4Pt,"jetB4Pt/F");
	TreeTestNew.Branch("jetB1M",&jetB1M,"jetB1M/F");
	TreeTestNew.Branch("jetB2M",&jetB2M,"jetB2M/F");
	TreeTestNew.Branch("jetB3M",&jetB3M,"jetB3M/F");
	TreeTestNew.Branch("jetB4M",&jetB4M,"jetB4M/F");
	TreeTestNew.Branch("nParticles",&nParticles,"nParticles/F");
	TreeTestNew.Branch("constSizeB1",&constSizeB1,"constSizeB1/F");
	TreeTestNew.Branch("constSizeB2",&constSizeB2,"constSizeB2/F");
	TreeTestNew.Branch("constSizeB3",&constSizeB3,"constSizeB3/F");
	TreeTestNew.Branch("constSizeB4",&constSizeB4,"constSizeB4/F");
	TreeTestNew.Branch("minConstSize",&minConstSize,"minConstSize/F");
	TreeTestNew.Branch("jetNObjects",&jetNObjects,"jetNObjects/F");
	TreeTestNew.Branch("minJetNObjects",&minJetNObjects,"minJetNObjects/F");
	TreeTestNew.Branch("invMassB1AntiKt",&invMassB1AntiKt,"invMassB1AntiKt/F");
	TreeTestNew.Branch("invMassB2AntiKt",&invMassB2AntiKt,"invMassB2AntiKt/F");
	TreeTestNew.Branch("nJetsAntiKt",&nJetsAntiKt,"nJetsAntiKt/F");
	TreeTestNew.Branch("invMassB11Best",&invMassB11Best,"invMassB11Best/F");
	TreeTestNew.Branch("invMassB21Best",&invMassB21Best,"invMassB21Best/F");
	TreeTestNew.Branch("invMassB12Best",&invMassB12Best,"invMassB12Best/F");
	TreeTestNew.Branch("invMassB22Best",&invMassB22Best,"invMassB22Best/F");
	TreeTestNew.Branch("invMassB13Best",&invMassB13Best,"invMassB13Best/F");
	TreeTestNew.Branch("invMassB23Best",&invMassB23Best,"invMassB23Best/F");
	TreeTestNew.Branch("invMassB14Best",&invMassB14Best,"invMassB14Best/F");
	TreeTestNew.Branch("invMassB24Best",&invMassB24Best,"invMassB24Best/F");
	TreeTestNew.Branch("invMassB15Best",&invMassB15Best,"invMassB15Best/F");
	TreeTestNew.Branch("invMassB25Best",&invMassB25Best,"invMassB25Best/F");
	TreeTestNew.Branch("invMassB16Best",&invMassB16Best,"invMassB16Best/F");
	TreeTestNew.Branch("invMassB26Best",&invMassB26Best,"invMassB26Best/F");
	TreeTestNew.Branch("invMassB17Best",&invMassB17Best,"invMassB17Best/F");
	TreeTestNew.Branch("invMassB27Best",&invMassB27Best,"invMassB27Best/F");
	TreeTestNew.Branch("invMassB18Best",&invMassB18Best,"invMassB18Best/F");
	TreeTestNew.Branch("invMassB28Best",&invMassB28Best,"invMassB28Best/F");
	//TreeTestNew.Branch("nJetsDurham0",&nJetsDurham0,"nJetsDurham0/F");
	//TreeTestNew.Branch("nJetsDurham5",&nJetsDurham5,"nJetsDurham5/F");
	//TreeTestNew.Branch("nJetsDurham10",&nJetsDurham10,"nJetsDurham10/F");
	TreeTestNew.Branch("nJetsDurham15",&nJetsDurham15,"nJetsDurham15/F");
	TreeTestNew.Branch("nJetsDurham20",&nJetsDurham20,"nJetsDurham20/F");
	TreeTestNew.Branch("nJetsDurham25",&nJetsDurham25,"nJetsDurham25/F");
	TreeTestNew.Branch("nJetsDurham30",&nJetsDurham30,"nJetsDurham30/F");
	TreeTestNew.Branch("distanceZ1MinChiSquaredZZMass",&distanceZ1MinChiSquaredZZMass,"distanceZ1MinChiSquaredZZMass/F");
	TreeTestNew.Branch("distanceZ2MinChiSquaredZZMass",&distanceZ2MinChiSquaredZZMass,"distanceZ2MinChiSquaredZZMass/F");
	TreeTestNew.Branch("exclYmerge12",&exclYmerge12,"exclYmerge12/F");
	TreeTestNew.Branch("exclYmerge23",&exclYmerge23,"exclYmerge23/F");
	TreeTestNew.Branch("exclYmerge34",&exclYmerge34,"exclYmerge34/F");
	TreeTestNew.Branch("exclYmerge45",&exclYmerge45,"exclYmerge45/F");
	TreeTestNew.Branch("exclYmerge56",&exclYmerge56,"exclYmerge56/F");
	
	TreeTrain->SetBranchAddress("entryIndex", &entryIndex);
	TreeTrain->SetBranchAddress("aplanarity", &aplanarity);
	TreeTrain->SetBranchAddress("invMassB1", &invMassB1);
	TreeTrain->SetBranchAddress("invMassB2", &invMassB2);
	TreeTrain->SetBranchAddress("minJetM", &minJetM);  
	TreeTrain->SetBranchAddress("sphericity", &sphericity);
	TreeTrain->SetBranchAddress("cosThetaB1", &cosThetaB1);
	TreeTrain->SetBranchAddress("cosThetaB2", &cosThetaB2);
	TreeTrain->SetBranchAddress("cosThetaB3", &cosThetaB3);
	TreeTrain->SetBranchAddress("cosThetaB4", &cosThetaB4);
	TreeTrain->SetBranchAddress("sumPt", &sumPt);
	TreeTrain->SetBranchAddress("jetB1Pt", &jetB1Pt);
	TreeTrain->SetBranchAddress("jetB2Pt", &jetB2Pt);
	TreeTrain->SetBranchAddress("jetB3Pt", &jetB3Pt);
	TreeTrain->SetBranchAddress("jetB4Pt", &jetB4Pt);
	TreeTrain->SetBranchAddress("jetB1M", &jetB1M);
	TreeTrain->SetBranchAddress("jetB2M", &jetB2M);
	TreeTrain->SetBranchAddress("jetB3M", &jetB3M);
	TreeTrain->SetBranchAddress("jetB4M", &jetB4M);
	TreeTrain->SetBranchAddress("nParticles", &nParticles);
	TreeTrain->SetBranchAddress("constSizeB1", &constSizeB1);
	TreeTrain->SetBranchAddress("constSizeB2", &constSizeB2);
	TreeTrain->SetBranchAddress("constSizeB3", &constSizeB3);
	TreeTrain->SetBranchAddress("constSizeB4", &constSizeB4);
	TreeTrain->SetBranchAddress("minConstSize", &minConstSize);
	TreeTrain->SetBranchAddress("jetNObjects", &jetNObjects);
	TreeTrain->SetBranchAddress("minJetNObjects", &minJetNObjects);
	TreeTrain->SetBranchAddress("invMassB1AntiKt", &invMassB1AntiKt);
	TreeTrain->SetBranchAddress("invMassB2AntiKt", &invMassB2AntiKt);
	TreeTrain->SetBranchAddress("nJetsAntiKt", &nJetsAntiKt);
	TreeTrain->SetBranchAddress("invMassB11Best", &invMassB11Best);
	TreeTrain->SetBranchAddress("invMassB21Best", &invMassB21Best);
	TreeTrain->SetBranchAddress("invMassB12Best", &invMassB12Best);
	TreeTrain->SetBranchAddress("invMassB22Best", &invMassB22Best);
	TreeTrain->SetBranchAddress("invMassB13Best", &invMassB13Best);
	TreeTrain->SetBranchAddress("invMassB23Best", &invMassB23Best);
	TreeTrain->SetBranchAddress("invMassB14Best", &invMassB14Best);
	TreeTrain->SetBranchAddress("invMassB24Best", &invMassB24Best);
	TreeTrain->SetBranchAddress("invMassB15Best", &invMassB15Best);
	TreeTrain->SetBranchAddress("invMassB25Best", &invMassB25Best);
	TreeTrain->SetBranchAddress("invMassB16Best", &invMassB16Best);
	TreeTrain->SetBranchAddress("invMassB26Best", &invMassB26Best);
	TreeTrain->SetBranchAddress("invMassB17Best", &invMassB17Best);
	TreeTrain->SetBranchAddress("invMassB27Best", &invMassB27Best);
	TreeTrain->SetBranchAddress("invMassB18Best", &invMassB18Best);
	TreeTrain->SetBranchAddress("invMassB28Best", &invMassB28Best);
	
	TreeTest->SetBranchAddress("entryIndex", &entryIndex);
	TreeTest->SetBranchAddress("aplanarity", &aplanarity);
	TreeTest->SetBranchAddress("invMassB1", &invMassB1);
	TreeTest->SetBranchAddress("invMassB2", &invMassB2);
	TreeTest->SetBranchAddress("minJetM", &minJetM);  
	TreeTest->SetBranchAddress("sphericity", &sphericity);
	TreeTest->SetBranchAddress("cosThetaB1", &cosThetaB1);
	TreeTest->SetBranchAddress("cosThetaB2", &cosThetaB2);
	TreeTest->SetBranchAddress("cosThetaB3", &cosThetaB3);
	TreeTest->SetBranchAddress("cosThetaB4", &cosThetaB4);
	TreeTest->SetBranchAddress("sumPt", &sumPt);
	TreeTest->SetBranchAddress("jetB1Pt", &jetB1Pt);
	TreeTest->SetBranchAddress("jetB2Pt", &jetB2Pt);
	TreeTest->SetBranchAddress("jetB3Pt", &jetB3Pt);
	TreeTest->SetBranchAddress("jetB4Pt", &jetB4Pt);
	TreeTest->SetBranchAddress("jetB1M", &jetB1M);
	TreeTest->SetBranchAddress("jetB2M", &jetB2M);
	TreeTest->SetBranchAddress("jetB3M", &jetB3M);
	TreeTest->SetBranchAddress("jetB4M", &jetB4M);
	TreeTest->SetBranchAddress("nParticles", &nParticles);
	TreeTest->SetBranchAddress("constSizeB1", &constSizeB1);
	TreeTest->SetBranchAddress("constSizeB2", &constSizeB2);
	TreeTest->SetBranchAddress("constSizeB3", &constSizeB3);
	TreeTest->SetBranchAddress("constSizeB4", &constSizeB4);
	TreeTest->SetBranchAddress("minConstSize", &minConstSize);
	TreeTest->SetBranchAddress("jetNObjects", &jetNObjects);
	TreeTest->SetBranchAddress("minJetNObjects", &minJetNObjects);
	TreeTest->SetBranchAddress("invMassB1AntiKt", &invMassB1AntiKt);
	TreeTest->SetBranchAddress("invMassB2AntiKt", &invMassB2AntiKt);
	TreeTest->SetBranchAddress("nJetsAntiKt", &nJetsAntiKt);
	TreeTest->SetBranchAddress("invMassB11Best", &invMassB11Best);
	TreeTest->SetBranchAddress("invMassB21Best", &invMassB21Best);
	TreeTest->SetBranchAddress("invMassB12Best", &invMassB12Best);
	TreeTest->SetBranchAddress("invMassB22Best", &invMassB22Best);
	TreeTest->SetBranchAddress("invMassB13Best", &invMassB13Best);
	TreeTest->SetBranchAddress("invMassB23Best", &invMassB23Best);
	TreeTest->SetBranchAddress("invMassB14Best", &invMassB14Best);
	TreeTest->SetBranchAddress("invMassB24Best", &invMassB24Best);
	TreeTest->SetBranchAddress("invMassB15Best", &invMassB15Best);
	TreeTest->SetBranchAddress("invMassB25Best", &invMassB25Best);
	TreeTest->SetBranchAddress("invMassB16Best", &invMassB16Best);
	TreeTest->SetBranchAddress("invMassB26Best", &invMassB26Best);
	TreeTest->SetBranchAddress("invMassB17Best", &invMassB17Best);
	TreeTest->SetBranchAddress("invMassB27Best", &invMassB27Best);
	TreeTest->SetBranchAddress("invMassB18Best", &invMassB18Best);
	TreeTest->SetBranchAddress("invMassB28Best", &invMassB28Best);
	
	TreeMerge.SetBranchAddress("entryIndex", &entryIndexMergeHolder);	
	//TreeMerge.SetBranchAddress("nJetsDurham0", &nJetsDurham0Holder);
	//TreeMerge.SetBranchAddress("nJetsDurham5", &nJetsDurham5Holder);
	//TreeMerge.SetBranchAddress("nJetsDurham10", &nJetsDurham10Holder);
	TreeMerge.SetBranchAddress("nJetsDurham15", &nJetsDurham15Holder);
	TreeMerge.SetBranchAddress("nJetsDurham20", &nJetsDurham20Holder);
	TreeMerge.SetBranchAddress("nJetsDurham25", &nJetsDurham25Holder);
	TreeMerge.SetBranchAddress("nJetsDurham30", &nJetsDurham30Holder);
	TreeMerge.SetBranchAddress("distanceZ1MinChiSquaredZZMass", &distanceZ1MinChiSquaredZZMassHolder);
	TreeMerge.SetBranchAddress("distanceZ2MinChiSquaredZZMass", &distanceZ2MinChiSquaredZZMassHolder);
	TreeMerge.SetBranchAddress("exclYmerge12", &exclYmerge12Holder);
	TreeMerge.SetBranchAddress("exclYmerge23", &exclYmerge23Holder);
	TreeMerge.SetBranchAddress("exclYmerge34", &exclYmerge34Holder);
	TreeMerge.SetBranchAddress("exclYmerge45", &exclYmerge45Holder);
	TreeMerge.SetBranchAddress("exclYmerge56", &exclYmerge56Holder);
	
	
	int entries = TreeMerge.GetEntries();
	//entries=50;
	int contEntriesTrain=0, contEntriesTest=0;
	for(int entry=0; entry<entries; entry++)
	{
		TreeMerge.GetEntry(entry);
		//cout<<"entryIndexMerge: "<<entryIndexMerge<<endl;
		entryIndexMerge = entryIndexMergeHolder;
		//nJetsDurham0 = nJetsDurham0Holder;
		//nJetsDurham5 = nJetsDurham5Holder;
		//nJetsDurham10 = nJetsDurham10Holder;
		nJetsDurham15 = nJetsDurham15Holder;
		nJetsDurham20 = nJetsDurham20Holder;
		nJetsDurham25 = nJetsDurham25Holder;
		nJetsDurham30 = nJetsDurham30Holder;
		distanceZ1MinChiSquaredZZMass = distanceZ1MinChiSquaredZZMassHolder;
		distanceZ2MinChiSquaredZZMass = distanceZ2MinChiSquaredZZMassHolder;
		exclYmerge12 = exclYmerge12Holder;
		exclYmerge23 = exclYmerge23Holder;
		exclYmerge34 = exclYmerge34Holder;
		exclYmerge45 = exclYmerge45Holder;
		exclYmerge56 = exclYmerge56Holder;
		
		//if(entry<10) cout<<"entryIndexMerge: "<<entryIndexMerge<<endl<<"distanceZ1MinChiSquaredZZMass: "<<distanceZ1MinChiSquaredZZMass<<endl<<"exclYmerge34Holder: "<<exclYmerge34Holder<<endl;

		if(setTrain.find(entryIndexMerge) != setTrain.end())
		{
			TreeTrain->GetEntry(contEntriesTrain);
			TreeTrainNew.Fill();
			//if(entry<10) cout<<"Entra train"<<endl<<endl;
			contEntriesTrain++;
		}
		else if(setTest.find(entryIndexMerge) != setTest.end())
		{
			TreeTest->GetEntry(contEntriesTest);
			TreeTestNew.Fill();
			//if(entry<10) cout<<"Entra test"<<endl<<endl;
			
			contEntriesTest++;
		}
		//else if(entry<10) cout<<"No entra ninguno"<<endl<<endl;
	}
}

// Function to check if a file exists
bool fileExists(const string& filename) 
{
	struct stat buffer;
	return (stat(filename.c_str(), &buffer) == 0);
}

void FSRGammaGammaHHbbbbAnalysis()
{
	string fileFunction = "generate";
	string rtdCut = "10";
	string preselection = "34BSplit";
	int sampleIndex = 0;
	string sampleName = "Sample" + to_string(sampleIndex);
	sampleName = "SampleN";
	
	cout<<"New macro working"<<endl;
	
	//double weightHH=0.001225;
	//double weightHH=0.0015552;
	double weightHH=0.0015552*1.155;
	double weightqq=0.0349, weightttbar=0.503;
	//double weightZZ=0.8167; 
	double weightZZ=0.17088*1.155; 
	double weightWW=0.5149, weightqqX=0.04347826, weightqqqqX=0.04, weightqqHX=0.001;
	//double weightZH=0.00207445;
	double weightZH=0.00207445*1.155;
	double weightpebb=0.7536;
	double weightpebbqq=0.1522;
	double weightpeqqH=0.1237;
	double weightpett=0.0570;
	
	string jetAlgoText = "(durham rtd_cut="+rtdCut+") ";
	string jetAlgo = "Jet"+rtdCut;
  	string genJetAlgo = "GenJet"+rtdCut;
  	//string rtdCut = "10";
  	string jetAlgoOutputTreeSTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeSHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeSTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeSHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	string jetAlgoOutputTreeBqqTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBqqTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	string jetAlgoOutputTreeBttTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBttHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBttTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBttHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	string jetAlgoOutputTreeBZZTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBZZTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
  	string jetAlgoOutputTreeBWWTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBWWTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBqqXTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBqqXTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBqqqqXTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBqqqqXTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBqqHXTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBqqHXTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBZHTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBZHHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBZHTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBZHHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBpebbTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpebbHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
  	string jetAlgoOutputTreeBpebbTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpebbHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBpebbqqTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpebbqqHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	string jetAlgoOutputTreeBpebbqqTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpebbqqHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBpeqqHTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpeqqHHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	string jetAlgoOutputTreeBpeqqHTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpeqqHHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";
	string jetAlgoOutputTreeBpettTrain = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpettHHbbbbESpreadDurham"+rtdCut+preselection+"Train"+sampleName+".root";
	string jetAlgoOutputTreeBpettTest = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpettHHbbbbESpreadDurham"+rtdCut+preselection+"Test"+sampleName+".root";

	string jetAlgoOutputTreeS = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeSHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBqq = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBtt = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBttHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBZZ = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBZZHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBWW = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBWWHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBqqX = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqXHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBqqqqX = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqqqXHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBqqHX = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBqqHXHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBZH = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBZHHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBpebb = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpebbHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBpebbqq = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpebbqqHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBpeqqH = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpeqqHHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";
	string jetAlgoOutputTreeBpett = "analysis/kinematicFitAnalysisFiles/newVarspe/fittedPythia/BTag85SiD2024XCC/outputTreeBpettHHbbbbESpreadDurham"+rtdCut+preselection+sampleName+".root";

  	
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
  	string jetAlgoOutputTreeBWWTest = "outputTreeBWWHHbbbbESpreadAntiKtTest.root";
	string jetAlgoOutputTreeBqqXTrain = "outputTreeBqqXHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBqqXTest = "outputTreeBqqXHHbbbbESpreadAntiKtTest.root";
	string jetAlgoOutputTreeBqqqqXTrain = "outputTreeBqqqqXHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBqqqqXTest = "outputTreeBqqqqXHHbbbbESpreadAntiKtTest.root";
	string jetAlgoOutputTreeBqqHXTrain = "outputTreeBqqHXHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBqqHXTest = "outputTreeBqqHXHHbbbbESpreadAntiKtTest.root";
	string jetAlgoOutputTreeBZHTrain = "outputTreeBZHHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBZHTest = "outputTreeBZHHHbbbbESpreadAntiKtTest.root";
	string jetAlgoOutputTreeBpebbTrain = "outputTreeBpebbHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBpebbTest = "outputTreeBpebbHHbbbbESpreadAntiKtTest.root";
	string jetAlgoOutputTreeBpebbqqTrain = "outputTreeBpebbqqHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBpebbqqTest = "outputTreeBpebbqqHHbbbbESpreadAntiKtTest.root";
	string jetAlgoOutputTreeBpeqqHTrain = "outputTreeBpeqqHHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBpeqqHTest = "outputTreeBpeqqHHHbbbbESpreadAntiKtTest.root";
	string jetAlgoOutputTreeBpettTrain = "outputTreeBpettHHbbbbESpreadAntiKtTrain.root"
  	string jetAlgoOutputTreeBpettTest = "outputTreeBpettHHbbbbESpreadAntiKtTest.root";
	*/
  	
  	cout<<"jetAlgo: "<<jetAlgoText<<endl;
  	if(fileFunction == "generate" || fileFunction == "merge")
  	{
	  	/*////////creation of File for TMVA for signal HH
		vector<const char *> inputFilesHH;
	  	//const char *inputFileHH = "analysis/FilesPostDelphes/GammaGammaHHESpreadAll.root";
		//const char *inputFileHH = "analysis/FilesPostDelphes/GammaGammaHHESpreadAllILDDSiDi.root";
		const char *inputFileHH = "analysis/FilesPostDelphes/GammaGammaHHESpreadAllSiD2024XCC.root";
	  	//const char *inputFileHH = "analysis/FilesPostDelphes/GammaGammaHH380All.root";
	  	inputFilesHH.push_back(inputFileHH);
	  	// Check if file exists and increment sampleIndex if necessary
	  	sampleIndex=0;
	  	TFile *outputTreeSTrain = new TFile(jetAlgoOutputTreeSTrain.c_str(), "recreate");
	  	TTree TreeSTrain("TreeSTrain","a simple Tree with simple variables (Train)");
	  	TFile *outputTreeSTest = new TFile(jetAlgoOutputTreeSTest.c_str(), "recreate");
	  	TTree TreeSTest("TreeSTest","a simple Tree with simple variables (Test)");
	  	TTree TreeSMerge("TreeSMerge","a simple Tree with simple variables (merge)");
	  	set<int> setTrainHH, setTestHH;
		TFile *outputTreeS = new TFile(jetAlgoOutputTreeS.c_str(), "recreate");
	  	TTree TreeS("TreeS","a simple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(1, setTrainHH, setTestHH, rtdCut); ///arguments: topology, set for training events, set for testing events
		//cout<<"size before analysis: "<<TreeS.GetEntries();
	  	analysis(inputFilesHH, 1, weightHH, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(1, setTrainHH, setTestHH, TreeSTrain, TreeSTest, TreeSMerge, rtdCut);
	  	
	  	outputTreeSTrain->cd();
	   	TreeSTrain.Write();
	   	
	   	outputTreeSTest->cd();
	   	TreeSTest.Write();

		outputTreeS->cd();
	   	TreeS.Write();
	   	
	   	outputTreeSTrain->Close();
		outputTreeSTest->Close();
		outputTreeS->Close();
	  	////////creation of File for TMVA for signal HH*/
	  	
		/*////////creation of File for TMVA for back qq  
		vector<const char *> inputFilesqq;
		//const char *inputFileqq = "analysis/FilesPostDelphes/GammaGammabbbbqqESpreadAll.root";
		//const char *inputFileqq = "analysis/FilesPostDelphes/GammaGammabbbbqqESpreadAllILDDSiDi.root";
		const char *inputFileqq = "analysis/FilesPostDelphes/GammaGammabbbbqqESpreadAllSiD2024XCC.root";
		//const char *inputFileqq = "analysis/FilesPostDelphes/GammaGammabbbbqq380All.root";
		inputFilesqq.push_back(inputFileqq);
		// Check if file exists and increment sampleIndex if necessary
		sampleIndex=0;
	  	TFile *outputTreeBqqTrain = new TFile(jetAlgoOutputTreeBqqTrain.c_str(), "recreate");
	  	TTree TreeBqqTrain("TreeBqqTrain","a bqqimple Tree with bqqimple variables (Train)");
	  	TFile *outputTreeBqqTest = new TFile(jetAlgoOutputTreeBqqTest.c_str(), "recreate");
	  	TTree TreeBqqTest("TreeBqqTest","a bqqimple Tree with bqqimple variables (Test)");
	  	TTree TreeBqqMerge("TreeBqqMerge","a bqqimple Tree with bqqimple variables (merge)");
	  	set<int> setTrainqq, setTestqq;
		TFile *outputTreeBqq = new TFile(jetAlgoOutputTreeBqq.c_str(), "recreate");
	  	TTree TreeBqq("TreeBqq","a bqqimple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(2, setTrainqq, setTestqq, rtdCut); ///arguments: topology, set for training events, set for testing events
	  	analysis(inputFilesqq, 2, weightqq, jetAlgoText, jetAlgo, genJetAlgo, TreeBqqTrain, TreeBqqTest, TreeBqqMerge, TreeBqq, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(2, setTrainqq, setTestqq, TreeBqqTrain, TreeBqqTest, TreeBqqMerge, rtdCut);
	  	
	  	outputTreeBqqTrain->cd();
	   	TreeBqqTrain.Write();
	   	
	   	outputTreeBqqTest->cd();
	   	TreeBqqTest.Write();

		outputTreeBqq->cd();
	   	TreeBqq.Write();
	  	
	  	outputTreeBqqTrain->Close();
		outputTreeBqqTest->Close();
		outputTreeBqq->Close();
	  	////////creation of File for TMVA for back qq*/
	  	
	  	/*////////creation of File for TMVA for back ttbar 
		vector<const char *> inputFilestt; 
	  	//const char *inputFilett = "analysis/FilesPostDelphes/GammaGammattAll.root";
		//const char *inputFilett = "analysis/FilesPostDelphes/GammaGammattAllILDDSiDi.root";
		const char *inputFilett = "analysis/FilesPostDelphes/GammaGammattAllSiD2024XCC.root";
	  	//const char *inputFilett = "analysis/FilesPostDelphes/GammaGammattSmall.root";
		inputFilestt.push_back(inputFilett);
	  	sampleIndex=0;
		TFile *outputTreeBttTrain = new TFile(jetAlgoOutputTreeBttTrain.c_str(), "recreate");
	  	TTree TreeBttTrain("TreeBttTrain","a bttimple Tree with bttimple variables (Train)");
	  	TFile *outputTreeBttTest = new TFile(jetAlgoOutputTreeBttTest.c_str(), "recreate");
	  	TTree TreeBttTest("TreeBttTest","a bttimple Tree with bttimple variables (Test)");
	  	TTree TreeBttMerge("TreeBttMerge","a bttimple Tree with bttimple variables (merge)");
	  	set<int> setTrainttbar, setTestttbar;
		TFile *outputTreeBtt = new TFile(jetAlgoOutputTreeBtt.c_str(), "recreate");
	  	TTree TreeBtt("TreeBtt","a bttimple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(3, setTrainttbar, setTestttbar, rtdCut); ///arguments: topology, set for training events, set for testing events
	  	analysis(inputFilestt, 3, weightttbar, jetAlgoText, jetAlgo, genJetAlgo, TreeBttTrain, TreeBttTest, TreeBttMerge, TreeBtt, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(3, setTrainttbar, setTestttbar, TreeBttTrain, TreeBttTest, TreeBttMerge, rtdCut);
	  	
	  	outputTreeBttTrain->cd();
	   	TreeBttTrain.Write();
	   	
	   	outputTreeBttTest->cd();
	   	TreeBttTest.Write();

		outputTreeBtt->cd();
	   	TreeBtt.Write();
		
		outputTreeBttTrain->Close();
		outputTreeBttTest->Close();
		outputTreeBtt->Close();
	  	////////creation of File for TMVA for back ttbar*/
	  	
	  	/*////////creation of File for TMVA for back ZZ  	
		vector<const char *> inputFilesZZ; 
		//const char *inputFileZZ = "analysis/FilesPostDelphes/GammaGammaZZESpreadAll.root";
		//const char *inputFileZZ = "analysis/FilesPostDelphes/GammaGammaZZESpreadAllILDDSiDi.root";
		const char *inputFileZZ = "analysis/FilesPostDelphes/GammaGammaZZESpreadAllSiD2024XCC.root";
		//const char *inputFileZZ = "analysis/FilesPostDelphes/GammaGammaZZ380All.root";
		inputFilesZZ.push_back(inputFileZZ);
		sampleIndex=0;
	  	TFile *outputTreeBZZTrain = new TFile(jetAlgoOutputTreeBZZTrain.c_str(), "recreate");
	  	TTree TreeBZZTrain("TreeBZZTrain","a bZZimple Tree with bZZimple variables (Train)");
	  	TFile *outputTreeBZZTest = new TFile(jetAlgoOutputTreeBZZTest.c_str(), "recreate");
	  	TTree TreeBZZTest("TreeBZZTest","a bZZimple Tree with bZZimple variables (Test)");
	  	TTree TreeBZZMerge("TreeBZZMerge","a bZZimple Tree with bZZimple variables (merge)");
	  	set<int> setTrainZZ, setTestZZ;
		TFile *outputTreeBZZ = new TFile(jetAlgoOutputTreeBZZ.c_str(), "recreate");
	  	TTree TreeBZZ("TreeBZZ","a bZZimple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(4, setTrainZZ, setTestZZ, rtdCut); ///arguments: topology, set for training events, set for testing events
	  	analysis(inputFilesZZ, 4, weightZZ, jetAlgoText, jetAlgo, genJetAlgo, TreeBZZTrain, TreeBZZTest, TreeBZZMerge, TreeBZZ, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(4, setTrainZZ, setTestZZ, TreeBZZTrain, TreeBZZTest, TreeBZZMerge, rtdCut);
	  	
	  	outputTreeBZZTrain->cd();
	   	TreeBZZTrain.Write();
	   	
	   	outputTreeBZZTest->cd();
	   	TreeBZZTest.Write();

		outputTreeBZZ->cd();
	   	TreeBZZ.Write();
	  	
	  	outputTreeBZZTrain->Close();
		outputTreeBZZTest->Close();
		outputTreeBZZ->Close();
	  	////////creation of File for TMVA for back ZZ
	  	
	  	///////creation of File for TMVA for back WW  
		vector<const char *> inputFilesWW; 
	  	//const char *inputFileWW = "analysis/FilesPostDelphes/GammaGammaWWESpreadAll.root";
		//const char *inputFileWW = "analysis/FilesPostDelphes/GammaGammaWWESpreadAllILDDSiDi.root";
		const char *inputFileWW = "analysis/FilesPostDelphes/GammaGammaWWESpreadAllSiD2024XCC.root";
	  	//const char *inputFileWW = "analysis/FilesPostDelphes/GammaGammaWW380All.root";
		inputFilesWW.push_back(inputFileWW);
	  	sampleIndex=0;
	  	TFile *outputTreeBWWTrain = new TFile(jetAlgoOutputTreeBWWTrain.c_str(), "recreate");
	  	TTree TreeBWWTrain("TreeBWWTrain","a bWWimple Tree with bWWimple variables (Train)");
	  	TFile *outputTreeBWWTest = new TFile(jetAlgoOutputTreeBWWTest.c_str(), "recreate");
	  	TTree TreeBWWTest("TreeBWWTest","a bWWimple Tree with bWWimple variables (Test)");
	  	TTree TreeBWWMerge("TreeBWWMerge","a bWWimple Tree with bWWimple variables (merge)");
	  	set<int> setTrainWW, setTestWW;
		TFile *outputTreeBWW = new TFile(jetAlgoOutputTreeBWW.c_str(), "recreate");
	  	TTree TreeBWW("TreeBWW","a bWWimple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(5, setTrainWW, setTestWW, rtdCut); ///arguments: topology, set for training events, set for testing events
	  	analysis(inputFilesWW, 5, weightWW, jetAlgoText, jetAlgo, genJetAlgo, TreeBWWTrain, TreeBWWTest, TreeBWWMerge, TreeBWW, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(5, setTrainWW, setTestWW, TreeBWWTrain, TreeBWWTest, TreeBWWMerge, rtdCut);
	  	
	  	outputTreeBWWTrain->cd();
	   	TreeBWWTrain.Write();
	   	
	   	outputTreeBWWTest->cd();
	   	TreeBWWTest.Write();

		outputTreeBWW->cd();
	   	TreeBWW.Write();
	  	
	  	outputTreeBWWTrain->Close();
		outputTreeBWWTest->Close();
		outputTreeBWW->Close();
	  	////////creation of File for TMVA for back WW*/

		/*///////creation of File for TMVA for back qqX
		vector<const char *> inputFilesqqX; 
		//const char *inputFileqqX = "analysis/FilesPostDelphes/eGammaqqXAll.root";
		//const char *inputFileqqX = "analysis/FilesPostDelphes/eGammaqqXAllILDDSiDi.root";
		const char *inputFileqqX = "analysis/FilesPostDelphes/eGammaqqXAllSiD2024XCC.root";
	  	//const char *inputFileqqX = "analysis/FilesPostDelphes/GammaGammaZZ380All.root";
		inputFilesqqX.push_back(inputFileqqX);
	  	sampleIndex=0;
		TFile *outputTreeBqqXTrain = new TFile(jetAlgoOutputTreeBqqXTrain.c_str(), "recreate");
		TTree TreeBqqXTrain("TreeBqqXTrain","a bqqXimple Tree with bqqXimple variables (Train)");
		TFile *outputTreeBqqXTest = new TFile(jetAlgoOutputTreeBqqXTest.c_str(), "recreate");
		TTree TreeBqqXTest("TreeBqqXTest","a bqqXimple Tree with bqqXimple variables (Test)");
		TTree TreeBqqXMerge("TreeBqqXMerge","a bqqXimple Tree with bqqXimple variables (merge)");
		set<int> setTrainqqX, setTestqqX;
		TFile *outputTreeBqqX = new TFile(jetAlgoOutputTreeBqqX.c_str(), "recreate");
	  	TTree TreeBqqX("TreeBqqX","a bqqXimple Tree with simple variables (full)");
		
		if(fileFunction=="merge") generateSetsMerge(6, setTrainqqX, setTestqqX, rtdCut); ///arguments: topology, set for training events, set for testing events
		analysis(inputFilesqqX, 6, weightqqX, jetAlgoText, jetAlgo, genJetAlgo, TreeBqqXTrain, TreeBqqXTest, TreeBqqXMerge, TreeBqqX, fileFunction, preselection);
		if(fileFunction=="merge") mergeTrees(6, setTrainqqX, setTestqqX, TreeBqqXTrain, TreeBqqXTest, TreeBqqXMerge, rtdCut);
		
		outputTreeBqqXTrain->cd();
		TreeBqqXTrain.Write();
		
		outputTreeBqqXTest->cd();
		TreeBqqXTest.Write();

		outputTreeBqqX->cd();
	   	TreeBqqX.Write();
		
		outputTreeBqqXTrain->Close();
		outputTreeBqqXTest->Close();
		outputTreeBqqX->Close();
	  	////////creation of File for TMVA for back qqX*/

		///////creation of File for TMVA for back qqqqX
		vector<const char *> inputFilesqqqqX; 
		//const char *inputFileqqqqX = "analysis/FilesPostDelphes/eGammaqqqqXAll.root";
		//const char *inputFileqqqqX = "analysis/FilesPostDelphes/eGammaqqqqXAllILDDSiDi.root";
	  	const char *inputFileqqqqX = "analysis/FilesPostDelphes/eGammaqqqqXAllSiD2024XCC.root";
	  	//const char *inputFileqqqqX = "analysis/FilesPostDelphes/GammaGammaZZ380All.root";
		inputFilesqqqqX.push_back(inputFileqqqqX);
	  	sampleIndex=0;
		TFile *outputTreeBqqqqXTrain = new TFile(jetAlgoOutputTreeBqqqqXTrain.c_str(), "recreate");
		TTree TreeBqqqqXTrain("TreeBqqqqXTrain","a bqqqqXimple Tree with bqqqqXimple variables (Train)");
		TFile *outputTreeBqqqqXTest = new TFile(jetAlgoOutputTreeBqqqqXTest.c_str(), "recreate");
		TTree TreeBqqqqXTest("TreeBqqqqXTest","a bqqqqXimple Tree with bqqqqXimple variables (Test)");
		TTree TreeBqqqqXMerge("TreeBqqqqXMerge","a bqqqqXimple Tree with bqqqqXimple variables (merge)");
		set<int> setTrainqqqqX, setTestqqqqX;
		TFile *outputTreeBqqqqX = new TFile(jetAlgoOutputTreeBqqqqX.c_str(), "recreate");
	  	TTree TreeBqqqqX("TreeBqqqqX","a bqqqqXimple Tree with bqqqqXimple variables (full)");
		
		if(fileFunction=="merge") generateSetsMerge(7, setTrainqqqqX, setTestqqqqX, rtdCut); ///arguments: topology, set for training events, set for testing events
		analysis(inputFilesqqqqX, 7, weightqqqqX, jetAlgoText, jetAlgo, genJetAlgo, TreeBqqqqXTrain, TreeBqqqqXTest, TreeBqqqqXMerge, TreeBqqqqX, fileFunction, preselection);
		if(fileFunction=="merge") mergeTrees(7, setTrainqqqqX, setTestqqqqX, TreeBqqqqXTrain, TreeBqqqqXTest, TreeBqqqqXMerge, rtdCut);
		
		outputTreeBqqqqXTrain->cd();
		TreeBqqqqXTrain.Write();
		
		outputTreeBqqqqXTest->cd();
		TreeBqqqqXTest.Write();

		outputTreeBqqqqX->cd();
	   	TreeBqqqqX.Write();
		
		outputTreeBqqqqXTrain->Close();
		outputTreeBqqqqXTest->Close();
		outputTreeBqqqqX->Close();
		////////creation of File for TMVA for back qqqqX
		
		/*///////creation of File for TMVA for back qqHX
		vector<const char *> inputFilesqqHX; 
		//const char *inputFileqqHX = "analysis/FilesPostDelphes/eGammaqqHXAll.root";
		//const char *inputFileqqHX = "analysis/FilesPostDelphes/eGammaqqHXAllILDDSiDi.root";
		const char *inputFileqqHX1 = "analysis/FilesPostDelphes/eGammaqqHXAllSiD2024XCC1.root";
		const char *inputFileqqHX2 = "analysis/FilesPostDelphes/eGammaqqHXAllSiD2024XCC2.root";
	  	//const char *inputFileqqHX = "analysis/FilesPostDelphes/GammaGammaZZ380All.root";
		inputFilesqqHX.push_back(inputFileqqHX1);
		inputFilesqqHX.push_back(inputFileqqHX2);
	  	sampleIndex=0;
		TFile *outputTreeBqqHXTrain = new TFile(jetAlgoOutputTreeBqqHXTrain.c_str(), "recreate");
		TTree TreeBqqHXTrain("TreeBqqHXTrain","a bqqHXimple Tree with bqqHXimple variables (Train)");
		TFile *outputTreeBqqHXTest = new TFile(jetAlgoOutputTreeBqqHXTest.c_str(), "recreate");
		TTree TreeBqqHXTest("TreeBqqHXTest","a bqqHXimple Tree with bqqHXimple variables (Test)");
		TTree TreeBqqHXMerge("TreeBqqHXMerge","a bqqHXimple Tree with bqqHXimple variables (merge)");
		set<int> setTrainqqHX, setTestqqHX;
		TFile *outputTreeBqqHX = new TFile(jetAlgoOutputTreeBqqHX.c_str(), "recreate");
	  	TTree TreeBqqHX("TreeBqqHX","a bqqHXimple Tree with simple variables (full)");
		
		if(fileFunction=="merge") generateSetsMerge(8, setTrainqqHX, setTestqqHX, rtdCut); ///arguments: topology, set for training events, set for testing events
		analysis(inputFilesqqHX, 8, weightqqHX, jetAlgoText, jetAlgo, genJetAlgo, TreeBqqHXTrain, TreeBqqHXTest, TreeBqqHXMerge, TreeBqqHX, fileFunction, preselection);
		if(fileFunction=="merge") mergeTrees(8, setTrainqqHX, setTestqqHX, TreeBqqHXTrain, TreeBqqHXTest, TreeBqqHXMerge, rtdCut);
		
		outputTreeBqqHXTrain->cd();
		TreeBqqHXTrain.Write();
		
		outputTreeBqqHXTest->cd();
		TreeBqqHXTest.Write();

		outputTreeBqqHX->cd();
	   	TreeBqqHX.Write();
		
		outputTreeBqqHXTrain->Close();
		outputTreeBqqHXTest->Close();
		outputTreeBqqHX->Close();
	  	////////creation of File for TMVA for back qqHX*/

		/*////////creation of File for TMVA for back ZH 	
		vector<const char *> inputFilesZH; 
		//const char *inputFileZH = "analysis/FilesPostDelphes/GammaGammaZHESpreadAll.root";
		//const char *inputFileZH = "analysis/FilesPostDelphes/GammaGammaZHESpreadAllILDDSiDi.root";
		const char *inputFileZH = "analysis/FilesPostDelphes/GammaGammaZHESpreadAllSiD2024XCC.root";
		//const char *inputFileZH = "analysis/FilesPostDelphes/GammaGammaZH380AllILDDSiDi.root";
		inputFilesZH.push_back(inputFileZH);
		sampleIndex=0;
	  	TFile *outputTreeBZHTrain = new TFile(jetAlgoOutputTreeBZHTrain.c_str(), "recreate");
	  	TTree TreeBZHTrain("TreeBZHTrain","a bZHimple Tree with bZHimple variables (Train)");
	  	TFile *outputTreeBZHTest = new TFile(jetAlgoOutputTreeBZHTest.c_str(), "recreate");
	  	TTree TreeBZHTest("TreeBZHTest","a bZHimple Tree with bZHimple variables (Test)");
	  	TTree TreeBZHMerge("TreeBZHMerge","a bZHimple Tree with bZHimple variables (merge)");
	  	set<int> setTrainZH, setTestZH;
		TFile *outputTreeBZH = new TFile(jetAlgoOutputTreeBZH.c_str(), "recreate");
	  	TTree TreeBZH("TreeBZH","a bZHimple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(9, setTrainZH, setTestZH, rtdCut); ///arguments: topology, set for training events, set for testing events
	  	analysis(inputFilesZH, 9, weightZH, jetAlgoText, jetAlgo, genJetAlgo, TreeBZHTrain, TreeBZHTest, TreeBZHMerge, TreeBZH, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(9, setTrainZH, setTestZH, TreeBZHTrain, TreeBZHTest, TreeBZHMerge, rtdCut);
	  	
	  	outputTreeBZHTrain->cd();
	   	TreeBZHTrain.Write();
	   	
	   	outputTreeBZHTest->cd();
	   	TreeBZHTest.Write();

		outputTreeBZH->cd();
	   	TreeBZH.Write();
	  	
	  	outputTreeBZHTrain->Close();
		outputTreeBZHTest->Close();
		outputTreeBZH->Close();
	  	////////creation of File for TMVA for back ZH*/

		/*////////creation of File for TMVA for back pebb 
		vector<const char *> inputFilespebb; 
		//const char *inputFilepebb = "analysis/FilesPostDelphes/pebbESpreadAllILDDSiDi.root";
		const char *inputFilepebb = "analysis/FilesPostDelphes/pebbESpreadAllSiD2024XCC.root";
		//const char *inputFilepebb = "analysis/FilesPostDelphes/pebb380AllILDDSiDi.root";
		inputFilespebb.push_back(inputFilepebb);
		sampleIndex=0;
	  	TFile *outputTreeBpebbTrain = new TFile(jetAlgoOutputTreeBpebbTrain.c_str(), "recreate");
	  	TTree TreeBpebbTrain("TreeBpebbTrain","a bpebbimple Tree with bpebbimple variables (Train)");
	  	TFile *outputTreeBpebbTest = new TFile(jetAlgoOutputTreeBpebbTest.c_str(), "recreate");
	  	TTree TreeBpebbTest("TreeBpebbTest","a bpebbimple Tree with bpebbimple variables (Test)");
	  	TTree TreeBpebbMerge("TreeBpebbMerge","a bpebbimple Tree with bpebbimple variables (merge)");
	  	set<int> setTrainpebb, setTestpebb;
		TFile *outputTreeBpebb = new TFile(jetAlgoOutputTreeBpebb.c_str(), "recreate");
	  	TTree TreeBpebb("TreeBpebb","a bpebbimple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(10, setTrainpebb, setTestpebb, rtdCut); ///arguments: topology, set for training events, set for testing events
	  	analysis(inputFilespebb, 10, weightpebb, jetAlgoText, jetAlgo, genJetAlgo, TreeBpebbTrain, TreeBpebbTest, TreeBpebbMerge, TreeBpebb, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(10, setTrainpebb, setTestpebb, TreeBpebbTrain, TreeBpebbTest, TreeBpebbMerge, rtdCut);
	  	
	  	outputTreeBpebbTrain->cd();
	   	TreeBpebbTrain.Write();
	   	
	   	outputTreeBpebbTest->cd();
	   	TreeBpebbTest.Write();

		outputTreeBpebb->cd();
	   	TreeBpebb.Write();
	  	
	  	outputTreeBpebbTrain->Close();
		outputTreeBpebbTest->Close();
		outputTreeBpebb->Close();
	  	////////creation of File for TMVA for back pebb

		////////creation of File for TMVA for back pebbqq 	
		vector<const char *> inputFilespebbqq; 
		//const char *inputFilepebbqq = "analysis/FilesPostDelphes/pebbqqESpreadAll.root";
		//const char *inputFilepebbqq = "analysis/FilesPostDelphes/pebbqqESpreadAllILDDSiDi.root";
		const char *inputFilepebbqq = "analysis/FilesPostDelphes/pebbqqESpreadAllSiD2024XCC.root";
		//const char *inputFilepebbqq = "analysis/FilesPostDelphes/pebbqq380AllILDDSiDi.root";
		inputFilespebbqq.push_back(inputFilepebbqq);
		sampleIndex=0;
	  	TFile *outputTreeBpebbqqTrain = new TFile(jetAlgoOutputTreeBpebbqqTrain.c_str(), "recreate");
	  	TTree TreeBpebbqqTrain("TreeBpebbqqTrain","a bpebbqqimple Tree with bpebbqqimple variables (Train)");
	  	TFile *outputTreeBpebbqqTest = new TFile(jetAlgoOutputTreeBpebbqqTest.c_str(), "recreate");
	  	TTree TreeBpebbqqTest("TreeBpebbqqTest","a bpebbqqimple Tree with bpebbqqimple variables (Test)");
	  	TTree TreeBpebbqqMerge("TreeBpebbqqMerge","a bpebbqqimple Tree with bpebbqqimple variables (merge)");
	  	set<int> setTrainpebbqq, setTestpebbqq;
		TFile *outputTreeBpebbqq = new TFile(jetAlgoOutputTreeBpebbqq.c_str(), "recreate");
	  	TTree TreeBpebbqq("TreeBpebbqq","a bpebbqqimple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(11, setTrainpebbqq, setTestpebbqq, rtdCut); ///arguments: topology, set for training events, set for testing events
	  	analysis(inputFilespebbqq, 11, weightpebbqq, jetAlgoText, jetAlgo, genJetAlgo, TreeBpebbqqTrain, TreeBpebbqqTest, TreeBpebbqqMerge, TreeBpebbqq, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(11, setTrainpebbqq, setTestpebbqq, TreeBpebbqqTrain, TreeBpebbqqTest, TreeBpebbqqMerge, rtdCut);
	  	
	  	outputTreeBpebbqqTrain->cd();
	   	TreeBpebbqqTrain.Write();
	   	
	   	outputTreeBpebbqqTest->cd();
	   	TreeBpebbqqTest.Write();

		outputTreeBpebbqq->cd();
	   	TreeBpebbqq.Write();
	  	
	  	outputTreeBpebbqqTrain->Close();
		outputTreeBpebbqqTest->Close();
		outputTreeBpebbqq->Close();
	  	////////creation of File for TMVA for back pebbqq*/

		/*////////creation of File for TMVA for back peqqH 	
		vector<const char *> inputFilespeqqH; 
		//const char *inputFilepeqqH = "analysis/FilesPostDelphes/peqqHESpreadAll.root";
		//const char *inputFilepeqqH = "analysis/FilesPostDelphes/peqqHESpreadAllILDDSiDi.root";
		const char *inputFilepeqqH = "analysis/FilesPostDelphes/peqqHESpreadAllSiD2024XCC.root";
		//const char *inputFilepeqqH = "analysis/FilesPostDelphes/peqqH380AllILDDSiDi.root";
		inputFilespeqqH.push_back(inputFilepeqqH);
		sampleIndex=0;
	  	TFile *outputTreeBpeqqHTrain = new TFile(jetAlgoOutputTreeBpeqqHTrain.c_str(), "recreate");
	  	TTree TreeBpeqqHTrain("TreeBpeqqHTrain","a bpeqqHimple Tree with bpeqqHimple variables (Train)");
	  	TFile *outputTreeBpeqqHTest = new TFile(jetAlgoOutputTreeBpeqqHTest.c_str(), "recreate");
	  	TTree TreeBpeqqHTest("TreeBpeqqHTest","a bpeqqHimple Tree with bpeqqHimple variables (Test)");
	  	TTree TreeBpeqqHMerge("TreeBpeqqHMerge","a bpeqqHimple Tree with bpeqqHimple variables (merge)");
	  	set<int> setTrainpeqqH, setTestpeqqH;
		TFile *outputTreeBpeqqH = new TFile(jetAlgoOutputTreeBpeqqH.c_str(), "recreate");
	  	TTree TreeBpeqqH("TreeBpeqqH","a bpeqqHimple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(12, setTrainpeqqH, setTestpeqqH, rtdCut); ///arguments: topology, set for training events, set for testing events
	  	analysis(inputFilespeqqH, 12, weightpeqqH, jetAlgoText, jetAlgo, genJetAlgo, TreeBpeqqHTrain, TreeBpeqqHTest, TreeBpeqqHMerge, TreeBpeqqH, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(12, setTrainpeqqH, setTestpeqqH, TreeBpeqqHTrain, TreeBpeqqHTest, TreeBpeqqHMerge, rtdCut);
	  	
	  	outputTreeBpeqqHTrain->cd();
	   	TreeBpeqqHTrain.Write();
	   	
	   	outputTreeBpeqqHTest->cd();
	   	TreeBpeqqHTest.Write();

		outputTreeBpeqqH->cd();
	   	TreeBpeqqH.Write();
	  	
	  	outputTreeBpeqqHTrain->Close();
		outputTreeBpeqqHTest->Close();
		outputTreeBpeqqH->Close();
	  	////////creation of File for TMVA for back peqqH*/

		////////creation of File for TMVA for back pett 	
		vector<const char *> inputFilespett; 
		//const char *inputFilepett = "analysis/FilesPostDelphes/pettESpreadAll.root";
		//const char *inputFilepett = "analysis/FilesPostDelphes/pettESpreadAllILDDSiDi.root";
		const char *inputFilepett = "analysis/FilesPostDelphes/pettESpreadAllSiD2024XCC.root";
		//const char *inputFilepett = "analysis/FilesPostDelphes/pett380AllILDDSiDi.root";
		inputFilespett.push_back(inputFilepett);
		sampleIndex=0;
	  	TFile *outputTreeBpettTrain = new TFile(jetAlgoOutputTreeBpettTrain.c_str(), "recreate");
	  	TTree TreeBpettTrain("TreeBpettTrain","a bpettimple Tree with bpettimple variables (Train)");
	  	TFile *outputTreeBpettTest = new TFile(jetAlgoOutputTreeBpettTest.c_str(), "recreate");
	  	TTree TreeBpettTest("TreeBpettTest","a bpettimple Tree with bpettimple variables (Test)");
	  	TTree TreeBpettMerge("TreeBpettMerge","a bpettimple Tree with bpettimple variables (merge)");
	  	set<int> setTrainpett, setTestpett;
		TFile *outputTreeBpett = new TFile(jetAlgoOutputTreeBpett.c_str(), "recreate");
	  	TTree TreeBpett("TreeBpett","a bpettimple Tree with simple variables (full)");
	  	
	  	if(fileFunction=="merge") generateSetsMerge(13, setTrainpett, setTestpett, rtdCut); ///arguments: topology, set for training events, set for testing events
	  	analysis(inputFilespett, 13, weightpett, jetAlgoText, jetAlgo, genJetAlgo, TreeBpettTrain, TreeBpettTest, TreeBpettMerge, TreeBpett, fileFunction, preselection);
	  	if(fileFunction=="merge") mergeTrees(13, setTrainpett, setTestpett, TreeBpettTrain, TreeBpettTest, TreeBpettMerge, rtdCut);
	  	
	  	outputTreeBpettTrain->cd();
	   	TreeBpettTrain.Write();
	   	
	   	outputTreeBpettTest->cd();
	   	TreeBpettTest.Write();

		outputTreeBpett->cd();
	   	TreeBpett.Write();
	  	
	  	outputTreeBpettTrain->Close();
		outputTreeBpettTest->Close();
		outputTreeBpett->Close();
	  	////////creation of File for TMVA for back pett



  	}

  	
	/*//////Pt
	TCanvas *c1001 = new TCanvas();
	TCanvas *c1002 = new TCanvas();
	TCanvas *c1003 = new TCanvas();
	TCanvas *c1004 = new TCanvas();
	//////Pt*/
	
  	/////Just to plot; no file creation
  	if(fileFunction == "plot")
	{
		TTree TreeSTrain("TreeSTrain","a simple Tree with simple variables (Train)");
		TTree TreeSTest("TreeSTest","a simple Tree with simple variables (Test)");
		TTree TreeSMerge("TreeSMerge","a simple Tree with simple variables (merge)");
		TTree TreeS("TreeS","a simple Tree with simple variables (full)");
		//const char *inputFileHH = "analysis/FilesPostDelphes/GammaGammaHHESpreadAll.root";
		const char *inputFileHH = "analysis/FilesPostDelphes/GammaGammaHHESpreadAllILDDSiDi.root";
	  	//const char *inputFileHH = "analysis/FilesPostDelphes/GammaGammaHH380All.root";
		//const char *inputFileHH = "analysis/FilesPostDelphes/GammaGammaHH380AllPU.root";
		//const char *inputFileqq = "analysis/FilesPostDelphes/GammaGammabbbbqqESpreadAll.root";
		const char *inputFileqq = "analysis/FilesPostDelphes/GammaGammabbbbqqESpreadAllILDDSiDi.root";
		//const char *inputFileqq = "analysis/FilesPostDelphes/GammaGammabbbbqq380All.root";
		//const char *inputFilett = "analysis/FilesPostDelphes/GammaGammattAll.root";
		const char *inputFilett = "analysis/FilesPostDelphes/GammaGammattAllILDDSiDi.root";
	  	//const char *inputFilett = "analysis/FilesPostDelphes/GammaGammattSmall.root";
		//const char *inputFileZZ = "analysis/FilesPostDelphes/GammaGammaZZESpreadAll.root";
		const char *inputFileZZ = "analysis/FilesPostDelphes/GammaGammaZZESpreadAllILDDSiDi.root";
		//const char *inputFileZZ = "analysis/FilesPostDelphes/GammaGammaZZ380All.root";
		//const char *inputFileWW = "analysis/FilesPostDelphes/GammaGammaWWESpreadAll.root";
		const char *inputFileWW = "analysis/FilesPostDelphes/GammaGammaWWESpreadAllILDDSiDi.root";
	  	//const char *inputFileWW = "analysis/FilesPostDelphes/GammaGammaWW380All.root";
		//const char *inputFileqqX = "analysis/FilesPostDelphes/eGammaqqXAll.root";
		const char *inputFileqqX = "analysis/FilesPostDelphes/eGammaqqXAllILDDSiDi.root";
		//const char *inputFileqqqqX = "analysis/FilesPostDelphes/eGammaqqqqXAll.root";
		const char *inputFileqqqqX = "analysis/FilesPostDelphes/eGammaqqqqXAllILDDSiDi.root";
		//const char *inputFileqqHX = "analysis/FilesPostDelphes/eGammaqqHXAll.root";
		const char *inputFileqqHX = "analysis/FilesPostDelphes/eGammaqqHXAllILDDSiDi.root";
		//const char *inputFileZH = "analysis/FilesPostDelphes/GammaGammaZHESpreadAll.root";
		const char *inputFileZH = "analysis/FilesPostDelphes/GammaGammaZHESpreadAllILDDSiDi.root";
		//const char *inputFileZH = "analysis/FilesPostDelphes/GammaGammaZH380AllILDDSiDi.root";
		const char *inputFilepebb = "analysis/FilesPostDelphes/pebbESpreadAllILDDSiDi.root";
		//const char *inputFilepebb = "analysis/FilesPostDelphes/pebb380AllILDDSiDi.root";
		const char *inputFilepebbqq = "analysis/FilesPostDelphes/pebbqqESpreadAllILDDSiDi.root";
		//const char *inputFilepebbqq = "analysis/FilesPostDelphes/pebbqq380AllILDDSiDi.root";
		const char *inputFilepeqqH = "analysis/FilesPostDelphes/peqqHESpreadAllILDDSiDi.root";
		//const char *inputFilepeqqH = "analysis/FilesPostDelphes/peqqH380AllILDDSiDi.root";
		const char *inputFilepett = "analysis/FilesPostDelphes/pettESpreadAllILDDSiDi.root";
		//const char *inputFilepett = "analysis/FilesPostDelphes/pett380AllILDDSiDi.root";
		//analysis(inputFileHH, 1, weightHH, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		/*analysis(inputFileqq, 2, weightqq, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		analysis(inputFilett, 3, weightttbar, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		analysis(inputFileZZ, 4, weightZZ, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		analysis(inputFileWW, 5, weightWW, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		analysis(inputFileqqX, 6, weightqqX, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		analysis(inputFileqqqqX, 7, weightqqqqX, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		analysis(inputFileqqHX, 8, weightqqHX, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		analysis(inputFileZH, 9, weightZH, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		analysis(inputFilepebb, 10, weightpebb, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		analysis(inputFilepebbqq, 11, weightpebbqq, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);*/
		//analysis(inputFilepeqqH, 12, weightpeqqH, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);
		/*analysis(inputFilepett, 13, weightpett, jetAlgoText, jetAlgo, genJetAlgo, TreeSTrain, TreeSTest, TreeSMerge, TreeS, fileFunction, preselection);*/
		/*TTree TreeDummy("dummy","dummy");
		analysis("hh", -999, 0.001225, jetAlgoText, jetAlgo, genJetAlgo, TreeDummy, TreeDummy, TreeDummy, fileFunction, preselection); /////Just to give time for the hists. to load; delete if not plotting.*/
	}
	/////Just to plot; no file creation
	
	/*trueBPairMass(inputFileHH, 1, 0.001225);
	trueBPairMass(inputFileqq, 2, 0.3062);
	trueBPairMass(inputFilett, 3, 0.503);*/

}
