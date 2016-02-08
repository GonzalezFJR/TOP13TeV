// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME mt2Dict

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
#include "mt2.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *mt2_Dictionary();
   static void mt2_TClassManip(TClass*);
   static void *new_mt2(void *p = 0);
   static void *newArray_mt2(Long_t size, void *p);
   static void delete_mt2(void *p);
   static void deleteArray_mt2(void *p);
   static void destruct_mt2(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::mt2*)
   {
      ::mt2 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::mt2));
      static ::ROOT::TGenericClassInfo 
         instance("mt2", "mt2.h", 15,
                  typeid(::mt2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &mt2_Dictionary, isa_proxy, 0,
                  sizeof(::mt2) );
      instance.SetNew(&new_mt2);
      instance.SetNewArray(&newArray_mt2);
      instance.SetDelete(&delete_mt2);
      instance.SetDeleteArray(&deleteArray_mt2);
      instance.SetDestructor(&destruct_mt2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::mt2*)
   {
      return GenerateInitInstanceLocal((::mt2*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::mt2*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *mt2_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::mt2*)0x0)->GetClass();
      mt2_TClassManip(theClass);
   return theClass;
   }

   static void mt2_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_mt2(void *p) {
      return  p ? new(p) ::mt2 : new ::mt2;
   }
   static void *newArray_mt2(Long_t nElements, void *p) {
      return p ? new(p) ::mt2[nElements] : new ::mt2[nElements];
   }
   // Wrapper around operator delete
   static void delete_mt2(void *p) {
      delete ((::mt2*)p);
   }
   static void deleteArray_mt2(void *p) {
      delete [] ((::mt2*)p);
   }
   static void destruct_mt2(void *p) {
      typedef ::mt2 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::mt2

namespace {
  void TriggerDictionaryInitialization_mt2Dict_Impl() {
    static const char* headers[] = {
"mt2.h",
0
    };
    static const char* includePaths[] = {
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/mt2/../TreeAnalysisTop/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/mt2/../PUWeight/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/mt2/../mt2/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/mt2/../LeptonSF/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/mt2/../BTagSFUtil/",
"/nfs/fanae/PAF_releases/head/include",
"/mnt_pool/fanae105/root_releases/root-6.06.00.slc6/include",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/mt2/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "mt2Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$mt2.h")))  mt2;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "mt2Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "mt2.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"mt2", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("mt2Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_mt2Dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_mt2Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_mt2Dict() {
  TriggerDictionaryInitialization_mt2Dict_Impl();
}
