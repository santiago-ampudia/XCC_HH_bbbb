#ifndef PXYCONSTRAINHH
#define PXYCONSTRAINHH

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "TMatrixD.h"
#include "TString.h"

#include "Math/Minimizer.h"
#include "LorentzVectorWithErrors.h"
#include "npsol_job.h"



using namespace std;

class pxyConstrainHH : public npsol_job {

public:
  //
  // Constructor and destructor
  //
  pxyConstrainHH(const vector<LorentzVectorWithErrors>* const mLV, double ecm, bool enableExtraTries, double nSigVar, bool majorPrint=false);
  ~pxyConstrainHH(); 

  //
  //
  // Minimization and chi squared
  //
  vector<LorentzVectorWithErrors>* NumericalMinimization();
    

private:   

  void SetInitialFitParameters(double* variableValuesBest, double factorSig12Bx=0, double factorSig34Bx=0, double factorSig12By=0, double factorSig34By=0,
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

  const vector<LorentzVectorWithErrors>* const m_mLV;
  double* m_xx0;
  double* m_xx;
  double* m_minVal;
  double* m_maxVal;
  TString* m_name;
  double* m_yy0;
  double* m_sig0;
  const double* m_xe;
  double* m_cov_matrix;
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
  
};

#endif
