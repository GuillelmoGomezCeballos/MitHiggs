//--------------------------------------------------------------------------------------------------
// $Id: H2lAnaMod.h,v 1.5 2009/06/15 15:00:20 loizides Exp $
//
// H2lAnaMod
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_HWWWMODS_H2LANAMOD_H
#define MITHIGGS_HWWWMODS_H2LANAMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

class TH1D;

namespace mithep 
{
  class H2lAnaMod : public BaseMod
  {
    public:
      H2lAnaMod(const char *name="H2lAnaMod", 
                const char *title="H->WWW (2l) analysis module");
      ~H2lAnaMod() {}

      void         SetCleanJetsName(const char *n)     { fCleanJetsName = n;    }
      void         SetCleanLeptonsName(const char *n)  { fCleanLeptonsName = n; }
      void         SetCleanMetName(const char *n)      { fCleanMetName = n;     }
      void         SetMaxZMass(Double_t m)             { fMaxZMass = m;         }
      void         SetMinMaxPt(Double_t pt)            { fMinMaxPt   = pt;      }
      void         SetMinDilMass(Double_t m)           { fDilMinMass = m;       }
      void         SetMinJetEt(Double_t jpt)           { fMinJetPt  = jpt;      }  
      void         SetMinPt(Double_t pt)               { fMinPt      = pt;      }
      void         SetMinZMass(Double_t m)             { fMinZMass = m;         }
      void         SetMissEt(Double_t met)             { fMinMissEt  = met;     }  

    protected:
      void         Process();
      void         SlaveBegin();

      TString      fCleanJetsName;        //clean jets name (input)
      TString      fCleanLeptonsName;     //clean leptons name (input)
      TString      fCleanMetName;         //clean missing et name (input)
      Double_t     fMinPt;                //minimum pt for leptons
      Double_t     fMinMaxPt;             //minimum pt for maximum lepton
      Double_t     fDilMinMass;           //minimum dilepton mass
      Double_t     fMinMissEt;            //minimum missing Et
      Double_t     fMinJetPt;             //minimum jet Pt to be accepted as a jet
      Double_t     fMinZMass;             //minimum Z mass
      Double_t     fMaxZMass;             //maximum Z mass
      TH1D        *fNAccCounters;         //!history of cuts
      TH1D        *fNLepCounts;           //!same sign lepton counts
      TH1D        *fNAllJets;             //!number of all jets
      TH1D        *fNCenJets;             //!number of central
      TH1D        *fAllDiLepMass;         //!dilepton mass for all dilepton pairs
      TH1D        *fElElMass;             //!electron-electron mass for all pairs
      TH1D        *fElMuMass;             //!electron-muon mass for all pairs
      TH1D        *fMuMuMass;             //!muon-muon mass for all pairs
      TH1D        *fMissEtBefore;         //!missing et distribution (before cuts)
      TH1D        *fMissEtAfter;          //!missing et distribution (after cuts)
      
    ClassDef(H2lAnaMod,1) // H2l analysis module
  };
}
#endif
