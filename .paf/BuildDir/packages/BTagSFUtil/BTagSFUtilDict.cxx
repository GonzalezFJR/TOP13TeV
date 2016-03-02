// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME BTagSFUtilDict

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
#include "BTagSFUtil.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *BTagSFUtil_Dictionary();
   static void BTagSFUtil_TClassManip(TClass*);
   static void delete_BTagSFUtil(void *p);
   static void deleteArray_BTagSFUtil(void *p);
   static void destruct_BTagSFUtil(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::BTagSFUtil*)
   {
      ::BTagSFUtil *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::BTagSFUtil));
      static ::ROOT::TGenericClassInfo 
         instance("BTagSFUtil", "BTagSFUtil.h", 10,
                  typeid(::BTagSFUtil), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &BTagSFUtil_Dictionary, isa_proxy, 0,
                  sizeof(::BTagSFUtil) );
      instance.SetDelete(&delete_BTagSFUtil);
      instance.SetDeleteArray(&deleteArray_BTagSFUtil);
      instance.SetDestructor(&destruct_BTagSFUtil);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::BTagSFUtil*)
   {
      return GenerateInitInstanceLocal((::BTagSFUtil*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::BTagSFUtil*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *BTagSFUtil_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::BTagSFUtil*)0x0)->GetClass();
      BTagSFUtil_TClassManip(theClass);
   return theClass;
   }

   static void BTagSFUtil_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrapper around operator delete
   static void delete_BTagSFUtil(void *p) {
      delete ((::BTagSFUtil*)p);
   }
   static void deleteArray_BTagSFUtil(void *p) {
      delete [] ((::BTagSFUtil*)p);
   }
   static void destruct_BTagSFUtil(void *p) {
      typedef ::BTagSFUtil current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::BTagSFUtil

namespace {
  void TriggerDictionaryInitialization_BTagSFUtilDict_Impl() {
    static const char* headers[] = {
"BTagSFUtil.h",
0
    };
    static const char* includePaths[] = {
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/BTagSFUtil/../TreeAnalysisTop/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/BTagSFUtil/../PUWeight/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/BTagSFUtil/../mt2/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/BTagSFUtil/../LeptonSF/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/BTagSFUtil/../BTagSFUtil/",
"/nfs/fanae/PAF_releases/head/include",
"/mnt_pool/fanae105/root_releases/root-6.06.00.slc6/include",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/BTagSFUtil/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "BTagSFUtilDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$BTagSFUtil.h")))  BTagSFUtil;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "BTagSFUtilDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "BTagSFUtil.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"BTagSFUtil", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("BTagSFUtilDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_BTagSFUtilDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_BTagSFUtilDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_BTagSFUtilDict() {
  TriggerDictionaryInitialization_BTagSFUtilDict_Impl();
}
