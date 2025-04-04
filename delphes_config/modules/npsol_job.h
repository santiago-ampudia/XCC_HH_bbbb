#ifndef Jet_Fits_npsol_job
#define Jet_Fits_npsol_job
#include <iostream>
#include "f2c.h"
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TArrayL.h"
#include <string>
using namespace std;




class npsol_job {

 private:


  long int n_npsol;
  long int nclin;
  long int ncnln;
  long int nrowa;
  long int nrowj;
  long int nrowr;
  TMatrixD* aalin;
  TMatrixD* aalin_fort;
  TArrayD* bbl;
  TArrayD* bbu;
  long int inform;
  long int iter;
  TArrayL* istate;
  TArrayD* c_npsol;
  TMatrixD* cjac_npsol_input;
  TMatrixD* cjac_npsol_fort;
  TMatrixD* cjac_npsol_output;
  TArrayD* clamda_input;
  TArrayD* clamda_fort;
  TArrayD* clamda_output;
  double objf;
  TArrayD* objgrd;
  TMatrixD* r_npsol_input;
  TMatrixD* r_npsol_fort;
  TMatrixD* r_npsol_output;
  TArrayD* x_npsol_input;
  TArrayD* x_npsol_fort;
  TArrayD* x_npsol_output;
  TArrayL* iw;
  long int leniw;
  TArrayD* w_npsol;
  long int lenw;



  

 protected: 


  long int run(const TArrayD* x_npsol_input_d);


 
 public: 
  
  TMatrixD* cjac_confun_user;
  TMatrixD* cjac_confun_user_fort;

  npsol_job(long int n_npsol_d, long int nclin_d, long int ncnln_d);

  ~npsol_job() { cout << "npsol_job destructor" << endl; }

  virtual long int objfun_user(long int* mode, const long int* n_objfun, const double* x_objfun, 
	     double* objf, double* objgrd, const long int* nstate);


  virtual long int confun_user(long int* mode, const long int* ncnln, const long int* n_confun, 
		  const long int* nrowj, const long int* needc, const double* x_confun, double* c_confun, 
		  TMatrixD* cjac_confun_user, const long int* nstate);

  void set_option(const char* cstr);

  void set_aalin(const TMatrixD* aalin_d) {*aalin = *aalin_d;}

  void set_bbl(const TArrayD* bbl_d) {*bbl = *bbl_d;}

  void set_bbu(const TArrayD* bbu_d) {*bbu = *bbu_d;}

  void set_istate(const TArrayL* istate_d) {*istate = *istate_d;}

  void set_cjac_npsol_input(const TMatrixD* cjac_npsol_input_d) {*cjac_npsol_input = *cjac_npsol_input_d;}

  void set_clamda_input(const TArrayD* clamda_input_d) {*clamda_input = *clamda_input_d;}

  void set_r_npsol_input(const TMatrixD* r_npsol_input_d) {*r_npsol_input = *r_npsol_input_d;}

  double get_objf() {return objf;}

  long int get_iter() {return iter;}

  const TArrayL* get_istate() {return istate;}

  const TArrayD* get_c_npsol() {return c_npsol;}

  const TMatrixD* get_cjac_npsol_output() {return cjac_npsol_output;}

  const TArrayD* get_clamda_output() {return clamda_output;}

  const TArrayD* get_objgrd() {return objgrd;}

  const TMatrixD* get_r_npsol_output() {return r_npsol_output;}

  const TArrayD* get_x_npsol_output() {return x_npsol_output;}

};



typedef long int (*confun_fp) (long int*, const long int*, const long int*, 
			  const long int*, const long int*, const double*, double*, 
			  double*, const long int*);
typedef long int (*objfun_fp) (long int*, const long int*, const double*, 
			  double*, double*, const long int*); 

extern "C" {
  long int npoptn_(char*, ftnlen);
  long int npsol_(long int*, long int*, long int*, 
	     long int*, long int*, long int*, double*, double*, 
	     double*, confun_fp, objfun_fp, long int*, long int*, long int*, 
	     double*, double*, double*, double*, 
	     double*, double*, double*, long int*, long int*, 
	     double*, long int*);
  long int objfun(long int*, const long int*, const double*, 
	     double*, double*, const long int*); 

  long int confun(long int*, const long int*, const long int*, 
	     const long int*, const long int*, const double*, double*, 
	     double*, const long int*);
}
#endif
