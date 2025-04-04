// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME epConstrainHH_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "epConstrainHH.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void delete_epConstrainHH(void *p);
   static void deleteArray_epConstrainHH(void *p);
   static void destruct_epConstrainHH(void *p);
   static void streamer_epConstrainHH(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::epConstrainHH*)
   {
      ::epConstrainHH *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::epConstrainHH >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("epConstrainHH", ::epConstrainHH::Class_Version(), "epConstrainHH.h", 20,
                  typeid(::epConstrainHH), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::epConstrainHH::Dictionary, isa_proxy, 16,
                  sizeof(::epConstrainHH) );
      instance.SetDelete(&delete_epConstrainHH);
      instance.SetDeleteArray(&deleteArray_epConstrainHH);
      instance.SetDestructor(&destruct_epConstrainHH);
      instance.SetStreamerFunc(&streamer_epConstrainHH);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::epConstrainHH*)
   {
      return GenerateInitInstanceLocal(static_cast<::epConstrainHH*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::epConstrainHH*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr epConstrainHH::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *epConstrainHH::Class_Name()
{
   return "epConstrainHH";
}

//______________________________________________________________________________
const char *epConstrainHH::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::epConstrainHH*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int epConstrainHH::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::epConstrainHH*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *epConstrainHH::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::epConstrainHH*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *epConstrainHH::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::epConstrainHH*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void epConstrainHH::Streamer(TBuffer &R__b)
{
   // Stream an object of class epConstrainHH.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b >> n_do_fit;
      R__b >> n_npsol_ep;
      R__b >> nclin_ep;
      R__b >> ncnln_ep;
      x_fit_ini.Streamer(R__b);
      bbl_ep.Streamer(R__b);
      bbu_ep.Streamer(R__b);
      aalin_ep.Streamer(R__b);
      {
         vector<LorentzVectorWithErrors> &R__stl =  m_mLV;
         R__stl.clear();
         TClass *R__tcl1 = TBuffer::GetClass(typeid(class LorentzVectorWithErrors));
         if (R__tcl1==0) {
            Error("m_mLV streamer","Missing the TClass object for class LorentzVectorWithErrors!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            LorentzVectorWithErrors R__t;
            R__b.StreamObject(&R__t,R__tcl1);
            R__stl.push_back(R__t);
         }
      }
      {
         vector<double> &R__stl =  m_xx0;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<double> &R__stl =  m_xx;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<double> &R__stl =  m_minVal;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<double> &R__stl =  m_maxVal;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<double> &R__stl =  m_yy0;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<double> &R__stl =  m_sig0;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<double> &R__stl =  m_xe;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<double> &R__stl =  m_cov_matrix;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            double R__t;
            R__b >> R__t;
            R__stl.push_back(R__t);
         }
      }
      {
         vector<TString> &R__stl =  m_name;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            TString R__t;
            R__t.Streamer(R__b);
            R__stl.push_back(R__t);
         }
      }
      R__b >> m_number_freeparam;
      R__b >> m_ecm;
      R__b >> m_enableExtraTries;
      R__b >> m_nSigVar;
      R__b >> m_majorPrint;
      m_bM.Streamer(R__b);
      m_energyM.Streamer(R__b);
      m_sqrtsM.Streamer(R__b);
      {
         vector<map<TString,int>*> &R__stl =  m_vecStringInd;
         R__stl.clear();
         TClass *R__tcl1 = TBuffer::GetClass(typeid(class std::map<class TString, int> *));
         if (R__tcl1==0) {
            Error("m_vecStringInd streamer","Missing the TClass object for class std::map<class TString, int> *!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            map<TString,int>* R__t;
            R__t = (map<TString,int>*)R__b.ReadObjectAny(R__tcl1);
            R__stl.push_back(R__t);
         }
      }
      R__b >> m_chi2;
      R__b >> m_inform;
      R__b.CheckByteCount(R__s, R__c, epConstrainHH::IsA());
   } else {
      R__c = R__b.WriteVersion(epConstrainHH::IsA(), kTRUE);
      R__b << n_do_fit;
      R__b << n_npsol_ep;
      R__b << nclin_ep;
      R__b << ncnln_ep;
      x_fit_ini.Streamer(R__b);
      bbl_ep.Streamer(R__b);
      bbu_ep.Streamer(R__b);
      aalin_ep.Streamer(R__b);
      {
         vector<LorentzVectorWithErrors> &R__stl =  m_mLV;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
         TClass *R__tcl1 = TBuffer::GetClass(typeid(class LorentzVectorWithErrors));
         if (R__tcl1==0) {
            Error("m_mLV streamer","Missing the TClass object for class LorentzVectorWithErrors!");
            return;
         }
            vector<LorentzVectorWithErrors>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b.StreamObject((LorentzVectorWithErrors*)&(*R__k),R__tcl1);
            }
         }
      }
      {
         vector<double> &R__stl =  m_xx0;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<double> &R__stl =  m_xx;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<double> &R__stl =  m_minVal;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<double> &R__stl =  m_maxVal;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<double> &R__stl =  m_yy0;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<double> &R__stl =  m_sig0;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<double> &R__stl =  m_xe;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<double> &R__stl =  m_cov_matrix;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<double>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << (*R__k);
            }
         }
      }
      {
         vector<TString> &R__stl =  m_name;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<TString>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            ((TString&)(*R__k)).Streamer(R__b);
            }
         }
      }
      R__b << m_number_freeparam;
      R__b << m_ecm;
      R__b << m_enableExtraTries;
      R__b << m_nSigVar;
      R__b << m_majorPrint;
      m_bM.Streamer(R__b);
      m_energyM.Streamer(R__b);
      m_sqrtsM.Streamer(R__b);
      {
         vector<map<TString,int>*> &R__stl =  m_vecStringInd;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
         TClass *R__tcl1 = TBuffer::GetClass(typeid(class std::map<class TString, int> *));
         if (R__tcl1==0) {
            Error("m_vecStringInd streamer","Missing the TClass object for class std::map<class TString, int> *!");
            return;
         }
            vector<map<TString,int>*>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b.WriteObjectAny((*R__k),R__tcl1);
            }
         }
      }
      R__b << m_chi2;
      R__b << m_inform;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_epConstrainHH(void *p) {
      delete (static_cast<::epConstrainHH*>(p));
   }
   static void deleteArray_epConstrainHH(void *p) {
      delete [] (static_cast<::epConstrainHH*>(p));
   }
   static void destruct_epConstrainHH(void *p) {
      typedef ::epConstrainHH current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_epConstrainHH(TBuffer &buf, void *obj) {
      ((::epConstrainHH*)obj)->::epConstrainHH::Streamer(buf);
   }
} // end of namespace ROOT for class ::epConstrainHH

namespace ROOT {
   static TClass *vectorlEmaplETStringcOintgRmUgR_Dictionary();
   static void vectorlEmaplETStringcOintgRmUgR_TClassManip(TClass*);
   static void *new_vectorlEmaplETStringcOintgRmUgR(void *p = nullptr);
   static void *newArray_vectorlEmaplETStringcOintgRmUgR(Long_t size, void *p);
   static void delete_vectorlEmaplETStringcOintgRmUgR(void *p);
   static void deleteArray_vectorlEmaplETStringcOintgRmUgR(void *p);
   static void destruct_vectorlEmaplETStringcOintgRmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<map<TString,int>*>*)
   {
      vector<map<TString,int>*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<map<TString,int>*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<map<TString,int>*>", -2, "vector", 383,
                  typeid(vector<map<TString,int>*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEmaplETStringcOintgRmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<map<TString,int>*>) );
      instance.SetNew(&new_vectorlEmaplETStringcOintgRmUgR);
      instance.SetNewArray(&newArray_vectorlEmaplETStringcOintgRmUgR);
      instance.SetDelete(&delete_vectorlEmaplETStringcOintgRmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEmaplETStringcOintgRmUgR);
      instance.SetDestructor(&destruct_vectorlEmaplETStringcOintgRmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<map<TString,int>*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<map<TString,int>*>","std::__1::vector<std::__1::map<TString, int, std::__1::less<TString>, std::__1::allocator<std::__1::pair<TString const, int>>>*, std::__1::allocator<std::__1::map<TString, int, std::__1::less<TString>, std::__1::allocator<std::__1::pair<TString const, int>>>*>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<map<TString,int>*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEmaplETStringcOintgRmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<map<TString,int>*>*>(nullptr))->GetClass();
      vectorlEmaplETStringcOintgRmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEmaplETStringcOintgRmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEmaplETStringcOintgRmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<map<TString,int>*> : new vector<map<TString,int>*>;
   }
   static void *newArray_vectorlEmaplETStringcOintgRmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<map<TString,int>*>[nElements] : new vector<map<TString,int>*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEmaplETStringcOintgRmUgR(void *p) {
      delete (static_cast<vector<map<TString,int>*>*>(p));
   }
   static void deleteArray_vectorlEmaplETStringcOintgRmUgR(void *p) {
      delete [] (static_cast<vector<map<TString,int>*>*>(p));
   }
   static void destruct_vectorlEmaplETStringcOintgRmUgR(void *p) {
      typedef vector<map<TString,int>*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<map<TString,int>*>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 383,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::__1::vector<double, std::__1::allocator<double>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETStringgR_Dictionary();
   static void vectorlETStringgR_TClassManip(TClass*);
   static void *new_vectorlETStringgR(void *p = nullptr);
   static void *newArray_vectorlETStringgR(Long_t size, void *p);
   static void delete_vectorlETStringgR(void *p);
   static void deleteArray_vectorlETStringgR(void *p);
   static void destruct_vectorlETStringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TString>*)
   {
      vector<TString> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TString>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TString>", -2, "vector", 383,
                  typeid(vector<TString>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETStringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TString>) );
      instance.SetNew(&new_vectorlETStringgR);
      instance.SetNewArray(&newArray_vectorlETStringgR);
      instance.SetDelete(&delete_vectorlETStringgR);
      instance.SetDeleteArray(&deleteArray_vectorlETStringgR);
      instance.SetDestructor(&destruct_vectorlETStringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TString> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TString>","std::__1::vector<TString, std::__1::allocator<TString>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TString>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETStringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TString>*>(nullptr))->GetClass();
      vectorlETStringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETStringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETStringgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TString> : new vector<TString>;
   }
   static void *newArray_vectorlETStringgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TString>[nElements] : new vector<TString>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETStringgR(void *p) {
      delete (static_cast<vector<TString>*>(p));
   }
   static void deleteArray_vectorlETStringgR(void *p) {
      delete [] (static_cast<vector<TString>*>(p));
   }
   static void destruct_vectorlETStringgR(void *p) {
      typedef vector<TString> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TString>

namespace ROOT {
   static TClass *vectorlELorentzVectorWithErrorsgR_Dictionary();
   static void vectorlELorentzVectorWithErrorsgR_TClassManip(TClass*);
   static void *new_vectorlELorentzVectorWithErrorsgR(void *p = nullptr);
   static void *newArray_vectorlELorentzVectorWithErrorsgR(Long_t size, void *p);
   static void delete_vectorlELorentzVectorWithErrorsgR(void *p);
   static void deleteArray_vectorlELorentzVectorWithErrorsgR(void *p);
   static void destruct_vectorlELorentzVectorWithErrorsgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<LorentzVectorWithErrors>*)
   {
      vector<LorentzVectorWithErrors> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<LorentzVectorWithErrors>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<LorentzVectorWithErrors>", -2, "vector", 383,
                  typeid(vector<LorentzVectorWithErrors>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlELorentzVectorWithErrorsgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<LorentzVectorWithErrors>) );
      instance.SetNew(&new_vectorlELorentzVectorWithErrorsgR);
      instance.SetNewArray(&newArray_vectorlELorentzVectorWithErrorsgR);
      instance.SetDelete(&delete_vectorlELorentzVectorWithErrorsgR);
      instance.SetDeleteArray(&deleteArray_vectorlELorentzVectorWithErrorsgR);
      instance.SetDestructor(&destruct_vectorlELorentzVectorWithErrorsgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<LorentzVectorWithErrors> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<LorentzVectorWithErrors>","std::__1::vector<LorentzVectorWithErrors, std::__1::allocator<LorentzVectorWithErrors>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<LorentzVectorWithErrors>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlELorentzVectorWithErrorsgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<LorentzVectorWithErrors>*>(nullptr))->GetClass();
      vectorlELorentzVectorWithErrorsgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlELorentzVectorWithErrorsgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlELorentzVectorWithErrorsgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<LorentzVectorWithErrors> : new vector<LorentzVectorWithErrors>;
   }
   static void *newArray_vectorlELorentzVectorWithErrorsgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<LorentzVectorWithErrors>[nElements] : new vector<LorentzVectorWithErrors>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlELorentzVectorWithErrorsgR(void *p) {
      delete (static_cast<vector<LorentzVectorWithErrors>*>(p));
   }
   static void deleteArray_vectorlELorentzVectorWithErrorsgR(void *p) {
      delete [] (static_cast<vector<LorentzVectorWithErrors>*>(p));
   }
   static void destruct_vectorlELorentzVectorWithErrorsgR(void *p) {
      typedef vector<LorentzVectorWithErrors> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<LorentzVectorWithErrors>

namespace {
  void TriggerDictionaryInitialization_epConstrainHH_dict_Impl() {
    static const char* headers[] = {
"epConstrainHH.h",
nullptr
    };
    static const char* includePaths[] = {
"/opt/homebrew/Cellar/root/6.32.04/include/root",
"/Users/ampudia/DelphesTutorial/delphes/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "epConstrainHH_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$epConstrainHH.h")))  epConstrainHH;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "epConstrainHH_dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "epConstrainHH.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"epConstrainHH", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("epConstrainHH_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_epConstrainHH_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_epConstrainHH_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_epConstrainHH_dict() {
  TriggerDictionaryInitialization_epConstrainHH_dict_Impl();
}
