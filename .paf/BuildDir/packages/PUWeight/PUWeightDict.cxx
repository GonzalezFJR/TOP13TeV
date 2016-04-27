// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PUWeightDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "PUWeight.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *PUWeight_Dictionary();
   static void PUWeight_TClassManip(TClass*);
   static void delete_PUWeight(void *p);
   static void deleteArray_PUWeight(void *p);
   static void destruct_PUWeight(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PUWeight*)
   {
      ::PUWeight *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::PUWeight));
      static ::ROOT::TGenericClassInfo 
         instance("PUWeight", "PUWeight.h", 55,
                  typeid(::PUWeight), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &PUWeight_Dictionary, isa_proxy, 0,
                  sizeof(::PUWeight) );
      instance.SetDelete(&delete_PUWeight);
      instance.SetDeleteArray(&deleteArray_PUWeight);
      instance.SetDestructor(&destruct_PUWeight);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PUWeight*)
   {
      return GenerateInitInstanceLocal((::PUWeight*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::PUWeight*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *PUWeight_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::PUWeight*)0x0)->GetClass();
      PUWeight_TClassManip(theClass);
   return theClass;
   }

   static void PUWeight_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_PUWeight(void *p) {
      delete ((::PUWeight*)p);
   }
   static void deleteArray_PUWeight(void *p) {
      delete [] ((::PUWeight*)p);
   }
   static void destruct_PUWeight(void *p) {
      typedef ::PUWeight current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::PUWeight

namespace {
  void TriggerDictionaryInitialization_PUWeightDict_Impl() {
    static const char* headers[] = {
"PUWeight.h",
0
    };
    static const char* includePaths[] = {
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/PUWeight/../PUWeight/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/PUWeight/../mt2/",
"/nfs/fanae/PAF_releases/head/include",
"/mnt_pool/fanae105/root_releases/root-6.06.02.slc6/include",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/PUWeight/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PUWeightDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$PUWeight.h")))  PUWeight;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PUWeightDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "PUWeight.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"PUWeight", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PUWeightDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PUWeightDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PUWeightDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PUWeightDict() {
  TriggerDictionaryInitialization_PUWeightDict_Impl();
}
