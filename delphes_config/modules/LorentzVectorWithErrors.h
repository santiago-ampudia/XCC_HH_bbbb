#ifndef LORENTZVECTORWITHERRORS
#define LORENTZVECTORWITHERRORS

#include "TLorentzVector.h"
using namespace std;


class LorentzVectorWithErrors {
public:

  LorentzVectorWithErrors() : m_LV(TLorentzVector()),
			      m_sigEnergy(0.),
			      m_sigBetaX(0.), 
			      m_sigBetaY(0.), 
			      m_sigBetaZ(0.),
			      m_chi2perNDF(1e20) {}

  LorentzVectorWithErrors(const TLorentzVector& LV, Double_t sigEnergy,
			  Double_t sigBetaX, Double_t sigBetaY, Double_t sigBetaZ, Double_t chi2perNDF=1e20, Int_t inform=-1) : m_LV(LV),
													       m_sigEnergy(sigEnergy),
													       m_sigBetaX(sigBetaX), 
													       m_sigBetaY(sigBetaY), 
													       m_sigBetaZ(sigBetaZ),
													       m_chi2perNDF(chi2perNDF),
													       m_inform(inform) {}

  TLorentzVector getLV() const {return m_LV;}
  Double_t getSigEnergy() const {return m_sigEnergy;}
  Double_t getSigBetaX() const {return m_sigBetaX;}
  Double_t getSigBetaY() const {return m_sigBetaY;}
  Double_t getSigBetaZ() const {return m_sigBetaZ;}
  Double_t getChi2perNDF() const {return m_chi2perNDF;}
  Int_t getInform() const {return m_inform;}


protected:

  
  TLorentzVector m_LV;
  Double_t m_sigEnergy;
  Double_t m_sigBetaX;
  Double_t m_sigBetaY;
  Double_t m_sigBetaZ;
  Double_t m_chi2perNDF;
  Int_t m_inform;


};

#endif
