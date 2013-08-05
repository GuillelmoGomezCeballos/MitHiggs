// $Id: MitHiggsUtilsLinkDef.h,v 1.1 2012/02/07 16:00:06 sixie Exp $

#ifndef MITHIGGS_UTILS_LINKDEF_H
#define MITHIGGS_UTILS_LINKDEF_H

#include "MitHiggs/Utils/interface/MitStyle.h"
#include "MitHiggs/Utils/interface/PlotUtils.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::PlotUtils;
#pragma link C++ class mithep::EfficiencyUtils;
#pragma link C++ class mithep::MitStyle;
#endif
