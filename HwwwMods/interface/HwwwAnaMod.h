//--------------------------------------------------------------------------------------------------
// $Id: HwwwAnaMod.h,v 1.9 2009/06/15 15:00:20 loizides Exp $
//
// HwwwAnaMod
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_HWWWMODS_HWWWANAMOD_H
#define MITHIGGS_HWWWMODS_HWWWANAMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

class TH1D;

namespace mithep 
{
  class Particle;
  class Met;
  class HwwwAnaMod : public BaseMod
  {
    public:
      HwwwAnaMod(const char *name="HwwwAnaMod", 
                 const char *title="H->WWW analysis module");

      void         SetCleanJetsName(const char *n)     { fCleanJetsName = n;    }
      void         SetCleanLeptonsName(const char *n)  { fCleanLeptonsName = n; }
      void         SetCleanMetName(const char *n)      { fCleanMetName = n;     }
      void         SetDileptonMinMass(Double_t m)      { fDilMinMass = m;       }
      void         SetMinMaxPt(Double_t pt)            { fMinMaxPt   = pt;      }
      void         SetMinPt(Double_t pt)               { fMinPt      = pt;      }
      void         SetMissEt(Double_t met)             { fMinMissEt  = met;     }  

    protected:
      void         ComputeMt(const Particle *p0, const Particle *p1, const Met *met);
      void         Process();
      void         SlaveBegin();

      TString      fCleanJetsName;        //clean jets name (input)
      TString      fCleanLeptonsName;     //clean leptons name (input)
      TString      fCleanMetName;         //clean missing et name (input)
      Double_t     fMinPt;                //minimum pt for leptons
      Double_t     fMinMaxPt;             //minimum pt for maximum lepton
      Double_t     fMinMissEt;            //minimum missing Et
      Double_t     fDilMinMass;           //minimum dilepton mass
      TH1D        *fNAccCounters;         //!history of cuts
      TH1D        *fNAllJets;             //!number of all jets
      TH1D        *fNCenJets;             //!number of central
      TH1D        *fDelPhiLepLep;         //!delta phi between two leptons
      TH1D        *fDelPhiDilMet;         //!delta phi between dilepton and met
      TH1D        *fMtH;                  //!transverse higgs mass
      TH1D        *fMissEtBefore;         //!missing et distribution (before cuts)
      TH1D        *fMissEtAfter;          //!missing et distribution (after cuts)
      TH1D        *fDiLeptonM;            //!dilepton mass distribution
      TH1D        *fL1L2Mass;             //!mass distribution for mindelphi leptons
      TH1D        *fL1L2DelPhi;           //!delphi distribution for mindelphi leptons
      TH1D        *fL1L2DelMet;           //!delphi distribution between dileptons and met

    ClassDef(HwwwAnaMod,1) // HWWW analysis module
  };
}
#endif
