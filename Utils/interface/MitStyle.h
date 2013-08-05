#ifndef MITHIGGS_UTILS_MITSTYLE_H
#define MITHIGGS_UTILS_MITSTYLE_H

#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>

namespace mithep
{
  class MitStyle {
    public:
      static TCanvas* MakeCanvas(const char* name, const char *title, int dX = 500, int dY = 500);
      static void     InitSubPad(TPad* pad, int i);
      static void     InitHist  (TH1 *hist, const char *xtit, 
                                 const char *ytit  = "Number of Entries",
                                 EColor color = kBlack);
      static void     SetStyle  ();
       
      ClassDef(MitStyle, 0) // Mit Plotting Style        
  };
}

#endif
