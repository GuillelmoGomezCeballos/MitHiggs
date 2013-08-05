//--------------------------------------------------------------------------------------------------
// $Id: HwwwSignalMCMod.h,v 1.5 2009/06/15 15:00:20 loizides Exp $
//
// HwwwSignalMCMod
//
// This module collects interesting generator information and publishes collections
// for subsequent modules.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_HWWWMODS_SIGNALMCMOD_H
#define MITHIGGS_HWWWMODS_SIGNALMCMOD_H

#include "MitPhysics/Mods/interface/GeneratorMod.h" 

class TH1D;
class TH2D;

namespace mithep 
{
  class Particle;
  class HwwwSignalMCMod : public GeneratorMod
  {
    public:
      HwwwSignalMCMod(const char *name="HwwwSignalMCMod", 
                      const char *title="HwwwSignalMC module");
      ~HwwwSignalMCMod() {}

      void          SetDileptonMinMass(Double_t m) { fDilMinMass = m;   }
      void          SetDoSelect(Bool_t b)          { fDoSelect   = b;   }
      void          SetDoSameSign(Bool_t b)        { fDoSS       = b;   }
      void          SetMinMaxPt(Double_t pt)       { fMinMaxPt   = pt;  }
      void          SetMinPt(Double_t pt)          { fMinPt      = pt;  }
      void          SetMissEt(Double_t met)        { fMinMissEt  = met; }  

    protected:
      void          ComputeMt(const Particle *l0, const Particle *l1, const Particle *met);
      void          Process();
      void          SlaveBegin();

      Bool_t        fDoSelect;         //=true then select events
      Bool_t        fDoSS;             //=true then accept 2 SS lepton events
      Double_t      fMinPt;            //minimum pt for leptons
      Double_t      fMinMaxPt;         //minimum pt for maximum lepton
      Double_t      fMinMissEt;        //minimum missing Et
      Double_t      fDilMinMass;       //minimum dilepton mass
      TH1D         *fNAccCounters;     //!cut counters
      TH1D         *fMcHMt;            //!higgs      transverse mass
      TH1D         *fMcWWMt;           //!higgs->WW  transverse mass
      TH1D         *fMcWWMt2;          //!higgs->WW  transverse mass check
      TH1D         *fMcWWWMt;          //!higgs->WWW transverse mass
      TH1D         *fMcMissEt;         //!missing Et (for higgs->WW)
      TH1D         *fMcMissEtAll;      //!missing Et (for all neutrinos)
      TH1D         *fMcLeptonPt;       //!lepton pt distribution (3l)
      TH1D         *fMc2lDiLeptonM;    //!di-lepton mass distribution (2l)
      TH1D         *fMc3lDiLeptonM;    //!di-lepton mass distribution (3l)
      TH1D         *fMcHists[1024];    //!histograms potentially used
      
      ClassDef(HwwwSignalMCMod, 1) // Signal MC analysis module
  };
}
#endif
