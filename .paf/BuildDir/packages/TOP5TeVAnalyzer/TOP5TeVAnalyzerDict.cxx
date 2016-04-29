// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TOP5TeVAnalyzerDict

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
#include "TOP5TeVAnalyzer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TOP5TeVAnalyzer(void *p = 0);
   static void *newArray_TOP5TeVAnalyzer(Long_t size, void *p);
   static void delete_TOP5TeVAnalyzer(void *p);
   static void deleteArray_TOP5TeVAnalyzer(void *p);
   static void destruct_TOP5TeVAnalyzer(void *p);
   static void streamer_TOP5TeVAnalyzer(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TOP5TeVAnalyzer*)
   {
      ::TOP5TeVAnalyzer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TOP5TeVAnalyzer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TOP5TeVAnalyzer", ::TOP5TeVAnalyzer::Class_Version(), "TOP5TeVAnalyzer.h", 187,
                  typeid(::TOP5TeVAnalyzer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TOP5TeVAnalyzer::Dictionary, isa_proxy, 16,
                  sizeof(::TOP5TeVAnalyzer) );
      instance.SetNew(&new_TOP5TeVAnalyzer);
      instance.SetNewArray(&newArray_TOP5TeVAnalyzer);
      instance.SetDelete(&delete_TOP5TeVAnalyzer);
      instance.SetDeleteArray(&deleteArray_TOP5TeVAnalyzer);
      instance.SetDestructor(&destruct_TOP5TeVAnalyzer);
      instance.SetStreamerFunc(&streamer_TOP5TeVAnalyzer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TOP5TeVAnalyzer*)
   {
      return GenerateInitInstanceLocal((::TOP5TeVAnalyzer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TOP5TeVAnalyzer*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TOP5TeVAnalyzer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TOP5TeVAnalyzer::Class_Name()
{
   return "TOP5TeVAnalyzer";
}

//______________________________________________________________________________
const char *TOP5TeVAnalyzer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TOP5TeVAnalyzer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TOP5TeVAnalyzer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TOP5TeVAnalyzer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TOP5TeVAnalyzer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TOP5TeVAnalyzer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TOP5TeVAnalyzer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TOP5TeVAnalyzer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TOP5TeVAnalyzer::Streamer(TBuffer &R__b)
{
   // Stream an object of class TOP5TeVAnalyzer.

   PAFChainItemSelector::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TOP5TeVAnalyzer(void *p) {
      return  p ? new(p) ::TOP5TeVAnalyzer : new ::TOP5TeVAnalyzer;
   }
   static void *newArray_TOP5TeVAnalyzer(Long_t nElements, void *p) {
      return p ? new(p) ::TOP5TeVAnalyzer[nElements] : new ::TOP5TeVAnalyzer[nElements];
   }
   // Wrapper around operator delete
   static void delete_TOP5TeVAnalyzer(void *p) {
      delete ((::TOP5TeVAnalyzer*)p);
   }
   static void deleteArray_TOP5TeVAnalyzer(void *p) {
      delete [] ((::TOP5TeVAnalyzer*)p);
   }
   static void destruct_TOP5TeVAnalyzer(void *p) {
      typedef ::TOP5TeVAnalyzer current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TOP5TeVAnalyzer(TBuffer &buf, void *obj) {
      ((::TOP5TeVAnalyzer*)obj)->::TOP5TeVAnalyzer::Streamer(buf);
   }
} // end of namespace ROOT for class ::TOP5TeVAnalyzer

namespace {
  void TriggerDictionaryInitialization_TOP5TeVAnalyzerDict_Impl() {
    static const char* headers[] = {
"TOP5TeVAnalyzer.h",
0
    };
    static const char* includePaths[] = {
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TOP5TeVAnalyzer/../TreeAnalysisTop/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TOP5TeVAnalyzer/../TOP5TeVAnalyzer/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TOP5TeVAnalyzer/../PUWeight/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TOP5TeVAnalyzer/../mt2/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TOP5TeVAnalyzer/../LeptonSF/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TOP5TeVAnalyzer/../BTagSFUtil/",
"/nfs/fanae/PAF_releases/head/include",
"/mnt_pool/fanae105/root_releases/root-6.06.02.slc6/include",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TOP5TeVAnalyzer/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TOP5TeVAnalyzerDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TOP5TeVAnalyzer.h")))  TOP5TeVAnalyzer;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TOP5TeVAnalyzerDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TOP5TeVAnalyzer.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TOP5TeVAnalyzer", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TOP5TeVAnalyzerDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TOP5TeVAnalyzerDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TOP5TeVAnalyzerDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TOP5TeVAnalyzerDict() {
  TriggerDictionaryInitialization_TOP5TeVAnalyzerDict_Impl();
}
