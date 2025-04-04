#ifndef EPCONSTRAINHH
#define EPCONSTRAINHH

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "TMatrixD.h"
#include "TString.h"
#include "TObject.h"

#include "Math/Minimizer.h"
#include "LorentzVectorWithErrors.h"
#include "npsol_job.h"



using namespace std;

class epConstrainHH : public npsol_job {

public:
  //
  // Constructor and destructor
  //
  epConstrainHH(const vector<LorentzVectorWithErrors>* const mLV, double ecm, bool enableExtraTries, double nSigVar, bool majorPrint=false);
  ~epConstrainHH(); 

  //
  //
  // Minimization and chi squared
  //
  vector<LorentzVectorWithErrors>* NumericalMinimization();

  ClassDef(epConstrainHH, 1);
    

private:   

  void SetInitialFitParameters(vector<double>& variableValuesBest, double factorSig12Bx=0, double factorSig34Bx=0, double factorSig12By=0, double factorSig34By=0,
					     double factorSig12Bz=0, double factorSig34Bz=0);
  void npsolInit();

  long int objfun_user(long int* mode, const long int* n_objfun, const double* x_objfun, 
	     double* objf, double* objgrd, const long int* nstate) override;

  long int confun_user(long int* mode, const long int* ncnln, const long int* n_confun, 
		  const long int* nrowj, const long int* needc, const double* x_confun, double* c_confun, 
		  TMatrixD* cjac_confun_user, const long int* nstate) override;

  long int do_fit();

  Int_t n_do_fit;
  int n_npsol_ep;
  int nclin_ep;
  int ncnln_ep;
  TArrayD x_fit_ini;
  TArrayD bbl_ep;
  TArrayD bbu_ep;
  TMatrixD aalin_ep;

  vector<LorentzVectorWithErrors> m_mLV;
  vector<double> m_xx0;
  vector<double> m_xx;
  vector<double> m_minVal;
  vector<double> m_maxVal;
  vector<double> m_yy0;
  vector<double> m_sig0;
  vector<double> m_xe;
  vector<double> m_cov_matrix;
  vector<TString> m_name;
  int m_number_freeparam;
  double m_ecm;
  bool m_enableExtraTries;
  double m_nSigVar;
  bool m_majorPrint;
  TMatrixD m_bM;
  TMatrixD m_energyM;
  TMatrixD m_sqrtsM;
  vector<map<TString,int>*> m_vecStringInd;
  double m_chi2;
  int m_inform;

  //ClassDef(epConstrainHH, 1);
  
};

#endif
