// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
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
#include "linear_least_squares.h"
#include "Utilities.h"
#include "Utilities.h"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_Dict_Impl() {
    static const char* headers[] = {
"linear_least_squares.h",
"Utilities.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/include/",
"/mnt/c/Users/marratia/Linux/buildroot/include",
"/mnt/c/Users/marratia/Linux/kalmanalignmentsvn/tags/release-1.0/util/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "linear_least_squares.h"
#include "Utilities.h"
#include "Utilities.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ function setXYZ;

#pragma link C++ function GetSubMatrix;
#pragma link C++ function SetSubMatrix;
#pragma link C++ function GetSubVector;
#pragma link C++ function SetSubVector;

#pragma link C++ function CLHEPtoROOT;

#pragma link C++ function get_object;
#pragma link C++ function get_tree;
#pragma link C++ function init_align_tree;

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"CLHEPtoROOT", payloadCode, "@",
"GetSubMatrix", payloadCode, "@",
"GetSubVector", payloadCode, "@",
"SetSubMatrix", payloadCode, "@",
"SetSubVector", payloadCode, "@",
"get_object", payloadCode, "@",
"get_tree", payloadCode, "@",
"init_align_tree", payloadCode, "@",
"setXYZ", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Dict() {
  TriggerDictionaryInitialization_Dict_Impl();
}
