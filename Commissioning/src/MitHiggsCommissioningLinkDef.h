// $Id $

#ifndef COMMISSIONING_LINKDEF_H
#define COMMISSIONING_LINKDEF_H
#include "MitHiggs/Commissioning/interface/ElectronCommissioning.h"
#include "MitHiggs/Commissioning/interface/MuonCommissioning.h"
#include "MitHiggs/Commissioning/interface/JetCommissioning.h"

#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::ElectronCommissioning+;
#pragma link C++ class mithep::MuonCommissioning+;
#pragma link C++ class mithep::JetCommissioning+;

#endif
