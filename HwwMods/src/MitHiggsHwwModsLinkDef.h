// $Id: MitHiggsHwwModsLinkDef.h,v 1.3 2009/08/19 10:55:42 ceballos Exp $

#ifndef MITHIGGS_HWWMODS_LINKDEF_H
#define MITHIGGS_HWWMODS_LINKDEF_H
#include "MitHiggs/HwwMods/interface/HwwAnaMod.h"
#include "MitHiggs/HwwMods/interface/HwwEvtPreSelMod.h"
#include "MitHiggs/HwwMods/interface/HwwMCEvtSelMod.h"
#include "MitHiggs/HwwMods/interface/HwwMakeNtupleMod.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::HwwAnaMod+;
#pragma link C++ class mithep::HwwEvtPreSelMod+;
#pragma link C++ class mithep::HwwMCEvtSelMod+;
#pragma link C++ class mithep::HwwMakeNtupleMod+;
#endif
