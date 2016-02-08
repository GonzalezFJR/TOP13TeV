// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TreeAnalysisTopDict

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
#include "TreeAnalysisTop.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TreeAnalysisTop(void *p = 0);
   static void *newArray_TreeAnalysisTop(Long_t size, void *p);
   static void delete_TreeAnalysisTop(void *p);
   static void deleteArray_TreeAnalysisTop(void *p);
   static void destruct_TreeAnalysisTop(void *p);
   static void streamer_TreeAnalysisTop(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TreeAnalysisTop*)
   {
      ::TreeAnalysisTop *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TreeAnalysisTop >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TreeAnalysisTop", ::TreeAnalysisTop::Class_Version(), "TreeAnalysisTop.h", 187,
                  typeid(::TreeAnalysisTop), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TreeAnalysisTop::Dictionary, isa_proxy, 16,
                  sizeof(::TreeAnalysisTop) );
      instance.SetNew(&new_TreeAnalysisTop);
      instance.SetNewArray(&newArray_TreeAnalysisTop);
      instance.SetDelete(&delete_TreeAnalysisTop);
      instance.SetDeleteArray(&deleteArray_TreeAnalysisTop);
      instance.SetDestructor(&destruct_TreeAnalysisTop);
      instance.SetStreamerFunc(&streamer_TreeAnalysisTop);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TreeAnalysisTop*)
   {
      return GenerateInitInstanceLocal((::TreeAnalysisTop*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TreeAnalysisTop*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TreeAnalysisTop::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TreeAnalysisTop::Class_Name()
{
   return "TreeAnalysisTop";
}

//______________________________________________________________________________
const char *TreeAnalysisTop::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TreeAnalysisTop*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TreeAnalysisTop::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TreeAnalysisTop*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TreeAnalysisTop::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TreeAnalysisTop*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TreeAnalysisTop::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TreeAnalysisTop*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TreeAnalysisTop::Streamer(TBuffer &R__b)
{
   // Stream an object of class TreeAnalysisTop.

   PAFChainItemSelector::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TreeAnalysisTop(void *p) {
      return  p ? new(p) ::TreeAnalysisTop : new ::TreeAnalysisTop;
   }
   static void *newArray_TreeAnalysisTop(Long_t nElements, void *p) {
      return p ? new(p) ::TreeAnalysisTop[nElements] : new ::TreeAnalysisTop[nElements];
   }
   // Wrapper around operator delete
   static void delete_TreeAnalysisTop(void *p) {
      delete ((::TreeAnalysisTop*)p);
   }
   static void deleteArray_TreeAnalysisTop(void *p) {
      delete [] ((::TreeAnalysisTop*)p);
   }
   static void destruct_TreeAnalysisTop(void *p) {
      typedef ::TreeAnalysisTop current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TreeAnalysisTop(TBuffer &buf, void *obj) {
      ((::TreeAnalysisTop*)obj)->::TreeAnalysisTop::Streamer(buf);
   }
} // end of namespace ROOT for class ::TreeAnalysisTop

namespace {
  void TriggerDictionaryInitialization_TreeAnalysisTopDict_Impl() {
    static const char* headers[] = {
"TreeAnalysisTop.h",
0
    };
    static const char* includePaths[] = {
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TreeAnalysisTop/../TreeAnalysisTop/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TreeAnalysisTop/../PUWeight/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TreeAnalysisTop/../mt2/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TreeAnalysisTop/../LeptonSF/",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TreeAnalysisTop/../BTagSFUtil/",
"/nfs/fanae/PAF_releases/head/include",
"/mnt_pool/fanae105/root_releases/root-6.06.00.slc6/include",
"/mnt_pool/fanae105/user/juanr/TOP13TeV/.paf/BuildDir/packages/TreeAnalysisTop/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TreeAnalysisTopDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TreeAnalysisTop.h")))  TreeAnalysisTop;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TreeAnalysisTopDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TreeAnalysisTop.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TreeAnalysisTop", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TreeAnalysisTopDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TreeAnalysisTopDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TreeAnalysisTopDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TreeAnalysisTopDict() {
  TriggerDictionaryInitialization_TreeAnalysisTopDict_Impl();
}
