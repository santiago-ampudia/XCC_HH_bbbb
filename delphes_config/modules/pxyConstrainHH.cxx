#include "pxyConstrainHH.h"
#include "Fit/ParameterSettings.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMatrixDSym.h"
//------------------------------------------------------------------------------------------------
pxyConstrainHH::pxyConstrainHH(const vector<LorentzVectorWithErrors>* const mLV, double ecm, bool enableExtraTries, double nSigVar, bool majorPrint)
  : npsol_job(16,0,2),n_do_fit(0),n_npsol_ep(16),nclin_ep(0),ncnln_ep(2),x_fit_ini(n_npsol_ep),bbl_ep(n_npsol_ep+nclin_ep+ncnln_ep),bbu_ep(n_npsol_ep+nclin_ep+ncnln_ep),aalin_ep(max(1,nclin_ep),n_npsol_ep),
    m_mLV(mLV),m_xx0(new double[n_npsol_ep]),m_xx(new double[n_npsol_ep]),m_minVal(new double[n_npsol_ep]),m_maxVal(new double[n_npsol_ep]),
    m_name(new TString[n_npsol_ep]),
    m_yy0(new double[n_npsol_ep]),m_sig0(new double[n_npsol_ep]),m_xe(new double[n_npsol_ep]),m_cov_matrix(0),m_number_freeparam(0),
    m_ecm(ecm),m_enableExtraTries(enableExtraTries),m_nSigVar(nSigVar),m_majorPrint(majorPrint),
    m_bM(4,4),m_energyM(4,1),m_sqrtsM(4,1),m_vecStringInd(4),m_chi2(0.),m_inform(-1)
    

{
  cout << " pxyConstrainHH::pxyConstrainHH point 000 " << " enableExtraTries= " << enableExtraTries << endl;
  for(int i=0; i<4; i++) {
    m_vecStringInd[i]= new map<TString,int>();
    m_vecStringInd[i]->insert(std::pair<TString,int>("Energy",4*i));
    m_vecStringInd[i]->insert(std::pair<TString,int>("BetaX",4*i+1));
    m_vecStringInd[i]->insert(std::pair<TString,int>("BetaY",4*i+2));
    m_vecStringInd[i]->insert(std::pair<TString,int>("BetaZ",4*i+3));
    double energy=m_mLV->at(i).getLV().E();
    m_bM(0,i)=1;
    m_bM(1,i)=m_mLV->at(i).getLV().Px()/energy;
    m_bM(2,i)=m_mLV->at(i).getLV().Py()/energy;
    m_bM(3,i)=m_mLV->at(i).getLV().Pz()/energy;
    //    cout << " i= " << i << " m_bM(*,i)= " << m_bM(0,i) << " " <<  m_bM(1,i) << " " <<  m_bM(2,i) << " " <<  m_bM(3,i) << endl;
  }
  m_sqrtsM(0,0)=m_ecm;
  m_sqrtsM(1,0)=0;
  m_sqrtsM(2,0)=0;
  m_sqrtsM(3,0)=0;
  m_energyM=TMatrixD(TMatrixD(TMatrixD::kInverted,m_bM),TMatrixD::kMult,m_sqrtsM);
  TMatrixD testConstraint(m_bM,TMatrixD::kMult,m_energyM);
  //  cout << " pxyConstrainHH::pxyConstrainHH " << " m_ecm= " << m_ecm << " testConstraint(0:3,0)= " << testConstraint(0,0) << " " << testConstraint(1,0) << " " << testConstraint(2,0) << " " << testConstraint(3,0) << endl;
  for(int i=0; i<4; i++) {
    m_xx0[m_vecStringInd[i]->at("Energy")]=m_energyM(i,0);
    m_xx0[m_vecStringInd[i]->at("BetaX")]=m_bM(1,i);
    m_xx0[m_vecStringInd[i]->at("BetaY")]=m_bM(2,i);
    m_xx0[m_vecStringInd[i]->at("BetaZ")]=m_bM(3,i);
    
    m_yy0[m_vecStringInd[i]->at("Energy")]=m_mLV->at(i).getLV().E();
    m_yy0[m_vecStringInd[i]->at("BetaX")]=m_bM(1,i);
    m_yy0[m_vecStringInd[i]->at("BetaY")]=m_bM(2,i);
    m_yy0[m_vecStringInd[i]->at("BetaZ")]=m_bM(3,i);
    
    m_sig0[m_vecStringInd[i]->at("Energy")]=m_mLV->at(i).getSigEnergy();
    m_sig0[m_vecStringInd[i]->at("BetaX")]=m_mLV->at(i).getSigBetaX();
    m_sig0[m_vecStringInd[i]->at("BetaY")]=m_mLV->at(i).getSigBetaY();
    m_sig0[m_vecStringInd[i]->at("BetaZ")]=m_mLV->at(i).getSigBetaZ();
    
    m_minVal[m_vecStringInd[i]->at("Energy")]=1.e-4;
    m_minVal[m_vecStringInd[i]->at("BetaX")]=-1.;
    m_minVal[m_vecStringInd[i]->at("BetaY")]=-1.;
    m_minVal[m_vecStringInd[i]->at("BetaZ")]=-1.;
    
    m_maxVal[m_vecStringInd[i]->at("Energy")]=0.5*m_ecm;
    m_maxVal[m_vecStringInd[i]->at("BetaX")]=1.;
    m_maxVal[m_vecStringInd[i]->at("BetaY")]=1.;
    m_maxVal[m_vecStringInd[i]->at("BetaZ")]=1.;
    
    m_name[m_vecStringInd[i]->at("Energy")]=TString("Energy")+TString::Format("%d",i+1);
    m_name[m_vecStringInd[i]->at("BetaX")]=TString("BetaX")+TString::Format("%d",i+1);
    m_name[m_vecStringInd[i]->at("BetaY")]=TString("BetaY")+TString::Format("%d",i+1);
    m_name[m_vecStringInd[i]->at("BetaZ")]=TString("BetaZ")+TString::Format("%d",i+1);
    
    //    cout << " pxyConstrainHH::pxyConstrainHH " << " i= " << i << " measured energy= " << m_mLV->at(i).getLV().E() << " initial fit energy,beta x,y,z= " << m_energyM(i,0) << " " << m_bM(1,i) << " " << m_bM(2,i) << " " << m_bM(3,i) << endl;
    //    cout << " pxyConstrainHH::pxyConstrainHH " << " i= " << i << " m_xx0[4*i+0]= " << m_xx0[4*i+0] << " m_xx0[4*i+1]= " << m_xx0[4*i+1] << " m_xx0[4*i+2]= " << m_xx0[4*i+2] << " m_xx0[4*i+3]= " << m_xx0[4*i+3] << endl;
    //    cout << " pxyConstrainHH::pxyConstrainHH " << " i= " << i << " m_sig0[4*i+0]= " << m_sig0[4*i+0] << " m_sig0[4*i+1]= " << m_sig0[4*i+1] << " m_sig0[4*i+2]= " << m_sig0[4*i+2] << " m_sig0[4*i+3]= " << m_sig0[4*i+3] << endl;
    
  }
  npsolInit();
}

//------------------------------------------------------------------------------------------------
pxyConstrainHH::~pxyConstrainHH()
{
  //  cout << " ~pxyConstrainHH point 000 " << endl;
}


//------------------------------------------------------------------------------------------------
vector<LorentzVectorWithErrors>* pxyConstrainHH::NumericalMinimization()
{
  cout << " Minimising beam parameters ...\n" << endl;

  cout << " pxyConstrainHH::NumericalMinimization " << " n_npsol_ep= " << n_npsol_ep << endl;

  //
  // Create minimizer giving a name and a name (optionally) for the specific algorithm
  //
  // disable fit for now; just use initial values for energies

 
  
  int nExtraSig12BxTries=0;
  int nExtraSig34BxTries=0;
  int nExtraSig12ByTries=0;
  int nExtraSig34ByTries=0;
  int nExtraSig12BzTries=0;
  int nExtraSig34BzTries=0;

  
  if(m_enableExtraTries) {
    nExtraSig12BxTries=4;
    nExtraSig34BxTries=4;
    nExtraSig12ByTries=0;
    nExtraSig34ByTries=0;
    nExtraSig12BzTries=0;
    nExtraSig34BzTries=0;
  }
  double maxExtraSig12Bx=2;
  double maxExtraSig34Bx=2;
  double maxExtraSig12By=2;
  double maxExtraSig34By=2;
  double maxExtraSig12Bz=2;
  double maxExtraSig34Bz=2;

  double chi2Best=1.e20;
  double* variableValuesBest = new double[n_npsol_ep];
  
  //  cout << " pxyConstrainHH::NumericalMinimization " << " m_enableExtraTries= " << m_enableExtraTries
  //     << " nExtraSig12BxTries= " << nExtraSig12BxTries << " nExtraSig34BxTries= " << nExtraSig34BxTries
  //       << " nExtraSig12ByTries= " << nExtraSig12ByTries << " nExtraSig34ByTries= " << nExtraSig34ByTries
  //       << " nExtraSig12BzTries= " << nExtraSig12BzTries << " nExtraSig34BzTries= " << nExtraSig34BzTries
  //       << " maxExtraSig12Bx= " << maxExtraSig12Bx
  //       << " op.DefaultErrorDef()= " <<  op.DefaultErrorDef() << " op.DefaultPrintLevel()= " <<  op.DefaultPrintLevel()
  //       << " op.DefaultMaxFunctionCalls()= " << op.DefaultMaxFunctionCalls() << " DefaultMaxIterations()= " << op.DefaultMaxIterations() << endl;
  for(int iSig12Bx=0; iSig12Bx < nExtraSig12BxTries+1; iSig12Bx++) {
    Double_t factorSig12Bx=(iSig12Bx == 0 ? 1 : iSig12Bx*maxExtraSig12Bx/nExtraSig12BxTries);
    for(int iSig34Bx=0; iSig34Bx < nExtraSig34BxTries+1; iSig34Bx++) {
      Double_t factorSig34Bx=(iSig34Bx == 0 ? 1 : iSig34Bx*maxExtraSig34Bx/nExtraSig34BxTries);
      for(int iSig12By=0; iSig12By < nExtraSig12ByTries+1; iSig12By++) {
	Double_t factorSig12By=(iSig12By == 0 ? factorSig34Bx : iSig12By*maxExtraSig12By/nExtraSig12ByTries);
	for(int iSig34By=0; iSig34By < nExtraSig34ByTries+1; iSig34By++) {
	  Double_t factorSig34By=(iSig34By == 0 ? factorSig12By : iSig34By*maxExtraSig34By/nExtraSig34ByTries);
	  for(int iSig12Bz=0; iSig12Bz < nExtraSig12BzTries+1; iSig12Bz++) {
	    Double_t factorSig12Bz=(iSig12Bz == 0 ? factorSig34By : iSig12Bz*maxExtraSig12Bz/nExtraSig12BzTries);
	    for(int iSig34Bz=0; iSig34Bz < nExtraSig34BzTries+1; iSig34Bz++) {
	      Double_t factorSig34Bz=(iSig34Bz == 0 ? factorSig12Bz : iSig34Bz*maxExtraSig34Bz/nExtraSig34BzTries);
	      SetInitialFitParameters(0,factorSig12Bx,factorSig34Bx,factorSig12By,factorSig34By,
				      factorSig12Bz,factorSig34Bz);
	      //	      cout << " iSig12Bx= " << iSig12Bx << " factorSig12Bx= " << factorSig12Bx <<  " iSig34Bx= " << iSig34Bx << " factorSig34Bx= " << factorSig34Bx
	      //		   << " iSig12By= " << iSig12By << " factorSig12By= " << factorSig12By <<  " iSig34By= " << iSig34By << " factorSig34By= " << factorSig34By
	      //		   << endl
	      //		   << " iSig12Bz= " << iSig12Bz << " factorSig12Bz= " << factorSig12Bz <<  " iSig34Bz= " << iSig34Bz << " factorSig34Bz= " << factorSig34Bz
	      //		   << endl;

	      //
	      // Run minimisation procedure
	      //
	      /*	      for(int i=0; i<n_npsol_ep; i++) {
		cout << " before fit i= " << i << " " << m_name[i] << " value= " << x_fit_ini[i] << " lower= " << bbl_ep[i] << " upper= " << bbu_ep[i] << endl;
		}*/
	      m_inform=do_fit();


	      int ndof = n_npsol_ep-nclin_ep-ncnln_ep;
	      m_chi2=get_objf();

  
	      cout << " pxyConstrainHH after 1st fit " << " m_inform= " << m_inform << " m_chi2= " << m_chi2 << " ndof= " << ndof << " chi2/ndf= " << m_chi2/ndof << endl;
	      if(m_chi2 < chi2Best) {
		cout << " improved m_chi2= " << m_chi2 << " chi2Best= " << chi2Best << " improved m_chi2/ndf= " << m_chi2/ndof    << endl;
		chi2Best=m_chi2;
		const TArrayD* x_npsol_output=get_x_npsol_output();
		for(int i=0; i<n_npsol_ep; i++) {
		  variableValuesBest[i]=(*x_npsol_output)[i];
		}
	      }
	    }
	  }
	}
      }
    }
  }




  //
  // Set parameters to those from best fit and minimize again
  //

  //  cout << " pxyConstrainHH::NumericalMinimization point 001 " << endl;
  if(m_enableExtraTries) {
    SetInitialFitParameters(variableValuesBest);
    for(int i=0; i<n_npsol_ep; i++) {
      cout << " before final fit i= " << i << " " << m_name[i] << " value= " << x_fit_ini[i] << " lower= " << bbl_ep[i] << " upper= " << bbu_ep[i] << endl;
    }
    m_inform=do_fit();
  }

  int ndof = n_npsol_ep-nclin_ep-ncnln_ep;
  m_chi2=get_objf();
  
  TMatrixD uCholesky(*(get_r_npsol_output()));
  TArrayD xxFinal(*(get_x_npsol_output()));
  TMatrixDSym hessianMat(TMatrixDSym(TMatrixDSym::kUnit,TMatrixDSym(n_npsol_ep)).SimilarityT(uCholesky));
  TMatrixDSym covMat(TMatrixDSym(TMatrixDSym::kInverted,hessianMat));

  cout << " after final fit " << " m_chi2= " << m_chi2 << " ndof= " << ndof << " chi2/ndf= " << m_chi2/ndof << endl;
  /*  for(int i=0; i<n_npsol_ep; i++) {
    cout << " after final fit i= " << i << " " << m_name[i] << " value= " << xxFinal[i] << " lower= " << bbl_ep[i] << " upper= " << bbu_ep[i] << endl;
    }*/
  
  vector<LorentzVectorWithErrors>* result=new vector<LorentzVectorWithErrors>(4);
  for(int i=0; i<4; i++) {
    Double_t ee=xxFinal[m_vecStringInd[i]->at("Energy")];
    Double_t px=ee*xxFinal[m_vecStringInd[i]->at("BetaX")];
    Double_t py=ee*xxFinal[m_vecStringInd[i]->at("BetaY")];
    Double_t pz=ee*xxFinal[m_vecStringInd[i]->at("BetaZ")];
    TLorentzVector tL(px,py,pz,ee);
    Double_t sigEnergy=sqrt(covMat(m_vecStringInd[i]->at("Energy"),m_vecStringInd[i]->at("Energy")));
    Double_t sigBetaX=sqrt(covMat(m_vecStringInd[i]->at("BetaX"),m_vecStringInd[i]->at("BetaX")));
    Double_t sigBetaY=sqrt(covMat(m_vecStringInd[i]->at("BetaY"),m_vecStringInd[i]->at("BetaY")));
    Double_t sigBetaZ=sqrt(covMat(m_vecStringInd[i]->at("BetaZ"),m_vecStringInd[i]->at("BetaZ")));
    result->at(i)=LorentzVectorWithErrors(tL,sigEnergy,sigBetaX,sigBetaY,sigBetaZ,m_chi2/ndof,m_inform);
  }

  
  



  return result;
  
}


void pxyConstrainHH::SetInitialFitParameters(double* variableValuesBest, double factorSig12Bx, double factorSig34Bx, double factorSig12By, double factorSig34By,
					     double factorSig12Bz, double factorSig34Bz)
{


  //  std::cout << " SetInitialFitParameters point 000 " << std::endl;

  

  if(variableValuesBest) {
    m_bM(0,0)=1;
    m_bM(1,0)=variableValuesBest[m_vecStringInd[0]->at("BetaX")];
    m_bM(2,0)=variableValuesBest[m_vecStringInd[0]->at("BetaY")];
    m_bM(3,0)=variableValuesBest[m_vecStringInd[0]->at("BetaZ")];
    m_bM(0,1)=1;
    m_bM(1,1)=variableValuesBest[m_vecStringInd[1]->at("BetaX")];
    m_bM(2,1)=variableValuesBest[m_vecStringInd[1]->at("BetaY")];
    m_bM(3,1)=variableValuesBest[m_vecStringInd[1]->at("BetaZ")];
  
    m_bM(0,2)=1;
    m_bM(1,2)=variableValuesBest[m_vecStringInd[2]->at("BetaX")];
    m_bM(2,2)=variableValuesBest[m_vecStringInd[2]->at("BetaY")];
    m_bM(3,2)=variableValuesBest[m_vecStringInd[2]->at("BetaZ")];
    m_bM(0,3)=1;
    m_bM(1,3)=variableValuesBest[m_vecStringInd[3]->at("BetaX")];
    m_bM(2,3)=variableValuesBest[m_vecStringInd[3]->at("BetaY")];
    m_bM(3,3)=variableValuesBest[m_vecStringInd[3]->at("BetaZ")];
  }
  else {
    m_bM(0,0)=1;
    m_bM(1,0)=m_xx0[m_vecStringInd[0]->at("BetaX")]+m_nSigVar*m_sig0[m_vecStringInd[0]->at("BetaX")]*(factorSig12Bx-1.);
    m_bM(2,0)=m_xx0[m_vecStringInd[0]->at("BetaY")]+m_nSigVar*m_sig0[m_vecStringInd[0]->at("BetaY")]*(factorSig12By-1.);
    m_bM(3,0)=m_xx0[m_vecStringInd[0]->at("BetaZ")]+m_nSigVar*m_sig0[m_vecStringInd[0]->at("BetaZ")]*(factorSig12Bz-1.);
    m_bM(0,1)=1;
    m_bM(1,1)=m_xx0[m_vecStringInd[1]->at("BetaX")]+m_nSigVar*m_sig0[m_vecStringInd[1]->at("BetaX")]*(factorSig12Bx-1.);
    m_bM(2,1)=m_xx0[m_vecStringInd[1]->at("BetaY")]+m_nSigVar*m_sig0[m_vecStringInd[1]->at("BetaY")]*(factorSig12By-1.);
    m_bM(3,1)=m_xx0[m_vecStringInd[1]->at("BetaZ")]+m_nSigVar*m_sig0[m_vecStringInd[1]->at("BetaZ")]*(factorSig12Bz-1.);
  
    m_bM(0,2)=1;
    m_bM(1,2)=m_xx0[m_vecStringInd[2]->at("BetaX")]+m_nSigVar*m_sig0[m_vecStringInd[2]->at("BetaX")]*(factorSig34Bx-1.);
    m_bM(2,2)=m_xx0[m_vecStringInd[2]->at("BetaY")]+m_nSigVar*m_sig0[m_vecStringInd[2]->at("BetaY")]*(factorSig34By-1.);
    m_bM(3,2)=m_xx0[m_vecStringInd[2]->at("BetaZ")]+m_nSigVar*m_sig0[m_vecStringInd[2]->at("BetaZ")]*(factorSig34Bz-1.);
    m_bM(0,3)=1;
    m_bM(1,3)=m_xx0[m_vecStringInd[3]->at("BetaX")]+m_nSigVar*m_sig0[m_vecStringInd[3]->at("BetaX")]*(factorSig34Bx-1.);
    m_bM(2,3)=m_xx0[m_vecStringInd[3]->at("BetaY")]+m_nSigVar*m_sig0[m_vecStringInd[3]->at("BetaY")]*(factorSig34By-1.);
    m_bM(3,3)=m_xx0[m_vecStringInd[3]->at("BetaZ")]+m_nSigVar*m_sig0[m_vecStringInd[3]->at("BetaZ")]*(factorSig34Bz-1.);
  }
  
    

  //  m_energyM=TMatrixD(TMatrixD(TMatrixD::kInverted,m_bM),TMatrixD::kMult,m_sqrtsM);
  for(int i=0; i<4; i++) {
    //    x_fit_ini[m_vecStringInd[i]->at("Energy"))=m_energyM(i,0);
    x_fit_ini[m_vecStringInd[i]->at("BetaX")]=m_bM(1,i);
    x_fit_ini[m_vecStringInd[i]->at("BetaY")]=m_bM(2,i);
    x_fit_ini[m_vecStringInd[i]->at("BetaZ")]=m_bM(3,i);
  }
      

}

void pxyConstrainHH::npsolInit() 
{
  //  std::cout << " pxyConstrainHH::npsolInit()  point 001 " << std::endl;
  set_option("derivative level 3");
  set_option("optimality tolerance 1.e-7");
  set_option("function precision 3.e-16");
  set_option("Hessian Yes");
  set_option("major iteration limit 5000");
  if(m_majorPrint) { 
    set_option("verify level  3");
    set_option("major print  level 90");
  }
  else {
    set_option("verify level  -1");
    set_option("major print  level -1");
  }

  bbl_ep[n_npsol_ep]=0;
  bbu_ep[n_npsol_ep]=0;

  bbl_ep[n_npsol_ep+1]=0;
  bbu_ep[n_npsol_ep+1]=0;


  //  std::cout << " pxyConstrainHH::npsolInit()  point 001 " << std::endl;
  for(Int_t i=0; i<n_npsol_ep; i++) {
    x_fit_ini[i]=m_xx0[i];
    aalin_ep(0,i)=0;
    bbl_ep[i]=m_minVal[i];
    bbu_ep[i]=m_maxVal[i];
    //    std::cout << " pxyConstrainHH::npsolInit()  " << " i= " << i << " x_fit_ini[i]= " << x_fit_ini[i] << " bbl_ep,bbu_ep= " << bbl_ep[i] << " " << bbu_ep[i] << std::endl;
  }
  //  std::cout << " pxyConstrainHH::npsolInit()  point 002 " << std::endl;

  
}
long int pxyConstrainHH::do_fit()
{
  n_do_fit++;
  if(n_do_fit <=5) {
    std::cout << " entry to pxyConstrainHH::do_fit n_do_fit= " << n_do_fit << std::endl;
  }
  if(n_do_fit == 2) {

    set_option("derivative level 3");
    set_option("optimality tolerance 1.e-7");
    set_option("function precision 3.e-16");
    set_option("Hessian Yes");
    set_option("major iteration limit 5000");
    set_option("verify level  -1");
    set_option("major print  level -1");
  }


  
  
  set_bbl(&bbl_ep);
  set_bbu(&bbu_ep);
  set_aalin(&aalin_ep);



  return run(&x_fit_ini);
}

long int pxyConstrainHH::objfun_user(long int* mode, const long int* n_objfun, const double* xx, 
				     double* objf, double* objgrd, const long int* nstate)
{
  //  std::cout << " entry to objfun_user " << " mode= " << *mode << " n_objfun= " << *n_objfun << std::endl;
  *objf=0;
  for(int i=0; i<*n_objfun; i++) {
    *objf += pow((xx[i]-m_yy0[i])/m_sig0[i],2);
    objgrd[i]=2.*(xx[i]-m_yy0[i])/pow(m_sig0[i],2);
  }
  
  

  return 0;
}



long int pxyConstrainHH::confun_user(long int* mode, const long int* ncnln, const long int* n_confun, 
				    const long int* nrowj, const long int* needc, const double* x_confun, double* c_confun, 
				    TMatrixD* cjac_confun_user, const long int* nstate)
{
  //  std::cout << " entry to confun_user " << " mode= " << *mode << " ncnln= " << *ncnln << " n_confun= " << *n_confun << std::endl;

  for(int i=0; i<*ncnln; i++) {
    c_confun[i]=0;
    for(int j=0; j<*n_confun; j++) {
      (*cjac_confun_user)(i,j)=0;
    }
  }

  for(int i=0; i<*n_confun/4; i++) {
    c_confun[0] += x_confun[m_vecStringInd[i]->at("Energy")]*x_confun[m_vecStringInd[i]->at("BetaX")];
    c_confun[1] += x_confun[m_vecStringInd[i]->at("Energy")]*x_confun[m_vecStringInd[i]->at("BetaY")];

    (*cjac_confun_user)(0,m_vecStringInd[i]->at("Energy"))=x_confun[m_vecStringInd[i]->at("BetaX")];
    (*cjac_confun_user)(0,m_vecStringInd[i]->at("BetaX"))=x_confun[m_vecStringInd[i]->at("Energy")];
	
    (*cjac_confun_user)(1,m_vecStringInd[i]->at("Energy"))=x_confun[m_vecStringInd[i]->at("BetaY")];
    (*cjac_confun_user)(1,m_vecStringInd[i]->at("BetaY"))=x_confun[m_vecStringInd[i]->at("Energy")];
	
	
  }

    
  

  return 0;
}
