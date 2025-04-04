// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME pxyConstrainHH_dict
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
#include "pxyConstrainHH.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *pxyConstrainHH_Dictionary();
   static void pxyConstrainHH_TClassManip(TClass*);
   static void delete_pxyConstrainHH(void *p);
   static void deleteArray_pxyConstrainHH(void *p);
   static void destruct_pxyConstrainHH(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::pxyConstrainHH*)
   {
      ::pxyConstrainHH *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::pxyConstrainHH));
      static ::ROOT::TGenericClassInfo 
         instance("pxyConstrainHH", "pxyConstrainHH.h", 19,
                  typeid(::pxyConstrainHH), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &pxyConstrainHH_Dictionary, isa_proxy, 0,
                  sizeof(::pxyConstrainHH) );
      instance.SetDelete(&delete_pxyConstrainHH);
      instance.SetDeleteArray(&deleteArray_pxyConstrainHH);
      instance.SetDestructor(&destruct_pxyConstrainHH);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::pxyConstrainHH*)
   {
      return GenerateInitInstanceLocal(static_cast<::pxyConstrainHH*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::pxyConstrainHH*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *pxyConstrainHH_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::pxyConstrainHH*>(nullptr))->GetClass();
      pxyConstrainHH_TClassManip(theClass);
   return theClass;
   }

   static void pxyConstrainHH_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_pxyConstrainHH(void *p) {
      delete (static_cast<::pxyConstrainHH*>(p));
   }
   static void deleteArray_pxyConstrainHH(void *p) {
      delete [] (static_cast<::pxyConstrainHH*>(p));
   }
   static void destruct_pxyConstrainHH(void *p) {
      typedef ::pxyConstrainHH current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::pxyConstrainHH

namespace {
  void TriggerDictionaryInitialization_pxyConstrainHH_dict_Impl() {
    static const char* headers[] = {
"pxyConstrainHH.h",
nullptr
    };
    static const char* includePaths[] = {
"/opt/homebrew/Cellar/root/6.32.04/include/root",
"/Users/ampudia/DelphesTutorial/delphes/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "pxyConstrainHH_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$pxyConstrainHH.h")))  pxyConstrainHH;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "pxyConstrainHH_dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "pxyConstrainHH.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"pxyConstrainHH", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("pxyConstrainHH_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_pxyConstrainHH_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_pxyConstrainHH_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_pxyConstrainHH_dict() {
  TriggerDictionaryInitialization_pxyConstrainHH_dict_Impl();
}
