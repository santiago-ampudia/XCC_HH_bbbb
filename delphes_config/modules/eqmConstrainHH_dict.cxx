// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME eqmConstrainHH_dict
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
#include "eqmConstrainHH.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *eqmConstrainHH_Dictionary();
   static void eqmConstrainHH_TClassManip(TClass*);
   static void delete_eqmConstrainHH(void *p);
   static void deleteArray_eqmConstrainHH(void *p);
   static void destruct_eqmConstrainHH(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::eqmConstrainHH*)
   {
      ::eqmConstrainHH *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::eqmConstrainHH));
      static ::ROOT::TGenericClassInfo 
         instance("eqmConstrainHH", "eqmConstrainHH.h", 19,
                  typeid(::eqmConstrainHH), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &eqmConstrainHH_Dictionary, isa_proxy, 0,
                  sizeof(::eqmConstrainHH) );
      instance.SetDelete(&delete_eqmConstrainHH);
      instance.SetDeleteArray(&deleteArray_eqmConstrainHH);
      instance.SetDestructor(&destruct_eqmConstrainHH);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::eqmConstrainHH*)
   {
      return GenerateInitInstanceLocal(static_cast<::eqmConstrainHH*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::eqmConstrainHH*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *eqmConstrainHH_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::eqmConstrainHH*>(nullptr))->GetClass();
      eqmConstrainHH_TClassManip(theClass);
   return theClass;
   }

   static void eqmConstrainHH_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_eqmConstrainHH(void *p) {
      delete (static_cast<::eqmConstrainHH*>(p));
   }
   static void deleteArray_eqmConstrainHH(void *p) {
      delete [] (static_cast<::eqmConstrainHH*>(p));
   }
   static void destruct_eqmConstrainHH(void *p) {
      typedef ::eqmConstrainHH current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::eqmConstrainHH

namespace {
  void TriggerDictionaryInitialization_eqmConstrainHH_dict_Impl() {
    static const char* headers[] = {
"eqmConstrainHH.h",
nullptr
    };
    static const char* includePaths[] = {
"/opt/homebrew/Cellar/root/6.32.04/include/root",
"/Users/ampudia/DelphesTutorial/delphes/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "eqmConstrainHH_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$eqmConstrainHH.h")))  eqmConstrainHH;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "eqmConstrainHH_dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "eqmConstrainHH.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"eqmConstrainHH", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("eqmConstrainHH_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_eqmConstrainHH_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_eqmConstrainHH_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_eqmConstrainHH_dict() {
  TriggerDictionaryInitialization_eqmConstrainHH_dict_Impl();
}
