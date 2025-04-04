#include <iostream>
#include "npsol_job.h"


npsol_job* g_npsol_job;

long int confun(long int* mode, const long int* ncnln, const long int* n_npsol, 
	   const long int* nrowj, const long int* needc, const double* x_confun, double* c_confun, 
	   double* cjac_confun_array_fort, const long int* nstate)
{
  TMatrixD* cjac_confun_user=g_npsol_job->cjac_confun_user;
  TMatrixD* cjac_confun_user_fort=g_npsol_job->cjac_confun_user_fort;

  cjac_confun_user->SetMatrixArray(cjac_confun_array_fort,"F");

  g_npsol_job->confun_user(mode, ncnln, n_npsol, nrowj, needc, x_confun, c_confun, 
			   cjac_confun_user, nstate);

  cjac_confun_user_fort->Transpose(*cjac_confun_user);
  memcpy(cjac_confun_array_fort,cjac_confun_user_fort->GetMatrixArray(),sizeof(double)*cjac_confun_user_fort->GetNoElements());
  
  return 0;   

}
long int objfun(long int* mode, const long int* n_objfun, const double* x_objfun, 
		double* objf, double* objgrd, const long int* nstate)
{
  g_npsol_job->objfun_user(mode, n_objfun, x_objfun, 
			   objf, objgrd, nstate);

  return 0;   

}


npsol_job::npsol_job(long int n_npsol_d, long int nclin_d, long int ncnln_d)
{
  cout << " entry to npsol_job constructor n_npsol,nclin,ncnln= " << n_npsol_d << " " << nclin_d << " " << ncnln_d << endl;
  g_npsol_job=this;
  n_npsol=n_npsol_d;
  nclin=nclin_d;
  ncnln=ncnln_d;
  bbl = new TArrayD(n_npsol+nclin+ncnln);
  bbu = new TArrayD(n_npsol+nclin+ncnln);
  //  cout << " point 003 npsol_job constructor bbl->GetSize()= " << bbl->GetSize() << " bbu->GetSize()= "  << bbu->GetSize() << endl;
  TArrayD& bbl_dr = *bbl;
  TArrayD& bbu_dr = *bbu;
  for(int i=0;i<n_npsol;i++) {
    //    cout << " point 003.10 npsol_job constructor i= " << i << endl;
    bbl_dr[i]=-1.e10;
    //    cout << " point 003.11 npsol_job constructor i= " << i << " bbl_dr[i]= " << bbl_dr[i] << endl;
    bbu_dr[i]=1.e10;
    //    cout << " point 003.12 npsol_job constructor i= " << i << " bbu_dr[i]= " << bbu_dr[i] << endl;
  } 
  nrowa=max(1L,nclin);
  nrowj=max(1L,ncnln);
  nrowr=n_npsol;
  aalin = new TMatrixD(nrowa,n_npsol);
  aalin_fort = new TMatrixD(TMatrixD::kTransposed,*aalin);
  istate = new TArrayL(n_npsol+nclin+ncnln);
  c_npsol = new TArrayD(nrowj);
  cjac_npsol_input = new TMatrixD(nrowj,n_npsol);
  cjac_npsol_fort = new TMatrixD(TMatrixD::kTransposed,*cjac_npsol_input);
  cjac_npsol_output = new TMatrixD(nrowj,n_npsol);
  clamda_input = new TArrayD(n_npsol+nclin+ncnln);
  clamda_fort = new TArrayD(n_npsol+nclin+ncnln);
  clamda_output = new TArrayD(n_npsol+nclin+ncnln);
  objgrd = new TArrayD(n_npsol);
  r_npsol_input = new TMatrixD(nrowr,n_npsol);
  r_npsol_fort = new TMatrixD(TMatrixD::kTransposed,*r_npsol_input);
  r_npsol_output = new TMatrixD(nrowr,n_npsol);
  x_npsol_input = new TArrayD(n_npsol);
  x_npsol_fort = new TArrayD(n_npsol);
  x_npsol_output = new TArrayD(n_npsol);
  leniw=3*n_npsol+nclin+2*ncnln;
  iw = new TArrayL(leniw);
  if(nclin == 0 && ncnln == 0 ) {
    lenw=20*n_npsol;
  }
  else if(ncnln == 0) {
    lenw=2*n_npsol*n_npsol+20*n_npsol+11*nclin;
  }
  else {
    lenw=2*n_npsol*n_npsol+n_npsol*nclin+2*n_npsol*ncnln+20*n_npsol+11*nclin+21*ncnln;
  }
  w_npsol = new TArrayD(lenw);
  
  //  cout << " point 006 npsol_job constructor " << " lenw= " << lenw << " nrowj= " << nrowj << " nrowr= " << nrowr << endl;
  cjac_confun_user = new TMatrixD(nrowj,n_npsol);
  cjac_confun_user_fort = new TMatrixD(TMatrixD::kTransposed,*cjac_confun_user);

}





long int npsol_job::run(const TArrayD* x_npsol_input_d)
{

  *x_npsol_input=*x_npsol_input_d;
  aalin_fort->Transpose(*aalin);
  cjac_npsol_fort->Transpose(*cjac_npsol_input);
  *clamda_fort=*clamda_input;
  r_npsol_fort->Transpose(*r_npsol_input);
  *x_npsol_fort=*x_npsol_input;
  cout << " npsol_job before call to run,   n_npsol,nclin,ncnln= " << n_npsol << " " << nclin << " " << ncnln
       <<" nrowa,nrowj,nrowr= " << nrowa << " " << nrowj << " " << nrowr
       << " leniw= " << leniw << " lenw= " << lenw << endl;

  npsol_(&n_npsol, &nclin, &ncnln, 
	 &nrowa, &nrowj, &nrowr, aalin_fort->GetMatrixArray(), bbl->GetArray(), 
	 bbu->GetArray(), &confun, &objfun, &inform, &iter, istate->GetArray(), 
	 c_npsol->GetArray(), cjac_npsol_fort->GetMatrixArray(), clamda_fort->GetArray(), &objf, 
	 objgrd->GetArray(), r_npsol_fort->GetMatrixArray(), x_npsol_fort->GetArray(), iw->GetArray(), &leniw, 
	 w_npsol->GetArray(), &lenw); 


  *x_npsol_output=*x_npsol_fort;
  r_npsol_output->Transpose(*r_npsol_fort);
  *clamda_output=*clamda_fort;
  cjac_npsol_output->Transpose(*cjac_npsol_fort);
  
  return inform;

}


void npsol_job::set_option(const char* cstr_d)
{
  size_t cstr_len=strlen(cstr_d);
  char* cstr = new char[cstr_len+1];
  strcpy(cstr,cstr_d);
  npoptn_(cstr, (ftnlen)cstr_len);
} 


long int npsol_job::objfun_user(long int* mode, const long int* n_objfun, const double* x_objfun, 
		double* objf, double* objgrd, const long int* nstate)
{
  return 0;
}

long int npsol_job::confun_user(long int* mode, const long int* ncnln, const long int* n_confun, 
		const long int* nrowj, const long int* needc, const double* x_confun, double* c_confun, 
		TMatrixD* cjac_confun_user, const long int* nstate)
{

  return 0;
}

