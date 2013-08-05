//--------------------------------------------------------------------------------------------------
// $Id: HwwAnaMod.h,v 1.2 2009/06/15 15:00:19 loizides Exp $
//
// HwwAnaMod
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITANA_PHYSICSMOD_HWWEVTSELMOD_H
#define MITANA_PHYSICSMOD_HWWEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

class TH1D;

namespace mithep 
{
  class HwwAnaMod : public BaseMod
  {
    public:
      HwwAnaMod(const char *name="HwwAnaMod", 
                const char *title="H->WW analysis module");

      void         SetCleanJetsName(const char *n)     { fCleanJetsName = n;    }
      void         SetCleanLeptonsName(const char *n)  { fCleanLeptonsName = n; }
      void         SetCleanMetName(const char *n)      { fCleanMetName = n;     }

    protected:
      void         Process();
      void         SlaveBegin();

      TString      fCleanJetsName;        //clean jets name (input)
      TString      fCleanLeptonsName;     //clean leptons name (input)
      TString      fCleanMetName;         //clean missing et name (input)
      TH1D        *fNAccCounters;         //!history of cuts
      TH1D        *fNAllJets;             //!number of all jets
      TH1D        *fNCenJets;             //!number of central
      TH1D        *fDelPhiLepLep;         //!delta phi between two leptons
      TH1D        *fDelPhiDilMet;         //!delta phi between dilepton and met
      TH1D        *fMtH;                  //!transverse higgs mass

    ClassDef(HwwAnaMod,1) // HWW analysis modules
  };
}
#endif
