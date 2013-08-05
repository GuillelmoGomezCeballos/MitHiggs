//--------------------------------------------------------------------------------------------------
// $Id: H2lNtupMod.h,v 1.2 2009/06/15 15:00:20 loizides Exp $
//
// H2lNtupMod
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_HWWWMODS_H2LNTUPMOD_H
#define MITHIGGS_HWWWMODS_H2LNTUPMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

class TNtuple;

namespace mithep 
{
  class H2lNtupMod : public BaseMod
  {
    public:
      H2lNtupMod(const char *name="H2lNtupMod", 
                  const char *title="H->WWW(2l) ntuple maker module");

      void         SetCleanJetsName(const char *n)     { fCleanJetsName    = n; }
      void         SetCleanLeptonsName(const char *n)  { fCleanLeptonsName = n; }
      void         SetCleanMetName(const char *n)      { fCleanMetName     = n; }
      void         SetMinJetPt(Double_t pt)            { fMinJetPt   = pt;      }
      void         SetMinMaxPt(Double_t pt)            { fMinMaxPt   = pt;      }
      void         SetMinPt(Double_t pt)               { fMinPt      = pt;      }
      void         SetMissEt(Double_t met)             { fMinMissEt  = met;     }  
      void         SetWeight(Double_t w)               { fWeight = w;           }

    protected:
      void         Process();
      void         SlaveBegin();

      TString      fCleanJetsName;        //clean jets name (input)
      TString      fCleanLeptonsName;     //clean leptons name (input)
      TString      fCleanMetName;         //clean missing et name (input)
      Double_t     fMinPt;                //minimum pt for leptons
      Double_t     fMinMaxPt;             //minimum pt for maximum lepton
      Double_t     fMinMissEt;            //minimum missing Et
      Double_t     fMinJetPt;             //minimum jet pt
      Double_t     fWeight;               //sample weight (1 if not specified)
      TNtuple     *fNtuple;               //!the ntuple

    ClassDef(H2lNtupMod,1) // HWWW(2l) ntuple maker module
  };
}
#endif
