// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LeptonSFDict

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
#include "LeptonSF.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *LeptonSF_Dictionary();
   static void LeptonSF_TClassManip(TClass*);
   static void *new_LeptonSF(void *p = 0);
   static void *newArray_LeptonSF(Long_t size, void *p);
   static void delete_LeptonSF(void *p);
   static void deleteArray_LeptonSF(void *p);
   static void destruct_LeptonSF(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LeptonSF*)
   {
      ::LeptonSF *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::LeptonSF));
      static ::ROOT::TGenericClassInfo 
         instance("LeptonSF", "LeptonSF.h", 19,
                  typeid(::LeptonSF), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &LeptonSF_Dictionary, isa_proxy, 0,
                  sizeof(::LeptonSF) );
      instance.SetNew(&new_LeptonSF);
      instance.SetNewArray(&newArray_LeptonSF);
      instance.SetDelete(&delete_LeptonSF);
      instance.SetDeleteArray(&deleteArray_LeptonSF);
      instance.SetDestructor(&destruct_LeptonSF);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LeptonSF*)
   {
      return GenerateInitInstanceLocal((::LeptonSF*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LeptonSF*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *LeptonSF_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::LeptonSF*)0x0)->GetClass();
      LeptonSF_TClassManip(theClass);
   return theClass;
   }

   static void LeptonSF_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_LeptonSF(void *p) {
      return  p ? new(p) ::LeptonSF : new ::LeptonSF;
   }
   static void *newArray_LeptonSF(Long_t nElements, void *p) {
      return p ? new(p) ::LeptonSF[nElements] : new ::LeptonSF[nElements];
   }
   // Wrapper around operator delete
   static void delete_LeptonSF(void *p) {
      delete ((::LeptonSF*)p);
   }
   static void deleteArray_LeptonSF(void *p) {
      delete [] ((::LeptonSF*)p);
   }
   static void destruct_LeptonSF(void *p) {
      typedef ::LeptonSF current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LeptonSF

namespace {
  void TriggerDictionaryInitialization_LeptonSFDict_Impl() {
    static const char* headers[] = {
"LeptonSF.h",
0
    };
    static const char* includePaths[] = {
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/LeptonSF/../TreeAnalysisTop/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/LeptonSF/../PUWeight/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/LeptonSF/../mt2/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/LeptonSF/../LeptonSF/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/LeptonSF/../BTagSFUtil/",
"/nfs/fanae/PAF_releases/head/include",
"/mnt_pool/fanae105/root_releases/root-6.06.00.slc6/include",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/LeptonSF/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LeptonSFDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$LeptonSF.h")))  LeptonSF;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LeptonSFDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "LeptonSF.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"LeptonSF", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LeptonSFDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LeptonSFDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LeptonSFDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LeptonSFDict() {
  TriggerDictionaryInitialization_LeptonSFDict_Impl();
}
