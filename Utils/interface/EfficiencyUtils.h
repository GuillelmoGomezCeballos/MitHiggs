//--------------------------------------------------------------------------------------------------
// $Id: EfficiencyUtils.h,v 1.1 2012/02/07 16:00:06 sixie Exp $
//
// EfficiencyUtils
//
// Math utility functions.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_UTILS_EFFICIENCYUTILS_H
#define MITHIGGS_UTILS_EFFICIENCYUTILS_H

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <utility>
#include <TGraphAsymmErrors.h>
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitCommon/DataFormats/interface/TH3DAsymErr.h"

namespace mithep
{
  class EfficiencyUtils {
    public:
      static TGraphAsymmErrors* createEfficiencyGraph(TH1F* numerator, TH1F* denominator,
                                               std::string histname, std::vector<Double_t> bins, 
                                               Int_t errorType,
                                               Double_t xlow, Double_t xhigh, 
                                               Double_t ylow, Double_t yhigh);
      static TH2DAsymErr* createEfficiencyHist2D(TH2F* numerator, TH2F* denominator,
                                          std::string histname, 
                                          std::vector<double> xbins, std::vector<double> ybins, 
                                                 Int_t errorType, Bool_t printDebug );
      static void createEfficiencyHist2D(TH2F* numerator, TH2F* denominator, 
                                         std::string histname, 
                                         std::vector<double> xbins, std::vector<double> ybins, 
                                         Int_t errorType, TFile *file );


       static TH3DAsymErr* createEfficiencyHist3D(TH3F* numerator, TH3F* denominator,
                                          std::string histname, 
                                          std::vector<double> xbins, std::vector<double> ybins, 
                                          std::vector<double> zbins, 
                                          Int_t errorType );
 
    ClassDef(EfficiencyUtils, 0) // utitily functions
  };
}


#endif
