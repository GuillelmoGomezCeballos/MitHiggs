// $Id: MitHiggsHwwwModsLinkDef.h,v 1.8 2009/04/29 15:10:30 loizides Exp $

#ifndef MITHIGGS_HWWWMODS_LINKDEF_H
#define MITHIGGS_HWWWMODS_LINKDEF_H
#include "MitHiggs/HwwwMods/interface/H2lAnaMod.h"
#include "MitHiggs/HwwwMods/interface/H2lNtupMod.h"
#include "MitHiggs/HwwwMods/interface/HwwwAnaMod.h"
#include "MitHiggs/HwwwMods/interface/HwwwNtupMod.h"
#include "MitHiggs/HwwwMods/interface/HwwwSignalMCMod.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::H2lAnaMod+;
#pragma link C++ class mithep::H2lNtupMod+;
#pragma link C++ class mithep::HwwwAnaMod+;
#pragma link C++ class mithep::HwwwNtupMod+;
#pragma link C++ class mithep::HwwwSignalMCMod+;
#endif
