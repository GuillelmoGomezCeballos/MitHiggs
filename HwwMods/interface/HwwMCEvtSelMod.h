//--------------------------------------------------------------------------------------------------
// $Id: HwwMCEvtSelMod.h,v 1.4 2009/06/15 15:00:19 loizides Exp $
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef HWWMODS_HWWMCEVTSELMOD_H
#define HWWMODS_HWWMCEVTSELMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCParticleCol.h"

class TH1D;

namespace mithep 
{
  class HwwMCEvtSelMod : public BaseMod 
  {
    public:
      HwwMCEvtSelMod(const char *name="HwwMCEvtSelMod", 
                     const char *title="HWW MC based event selection module");
      ~HwwMCEvtSelMod() {}

      const char              *GetMCPartName()              const { return fMCPartName; }
      void                     SetMCPartName(const char *n)       { fMCPartName=n; }

      enum EProdType {
        kUnknown = -1,
        kGluGlu,
        kVBF,
        kAll
      };

      enum EDecayType {
        kUnknownDecay = -1,
        kElEl,
        kElMu,
        kElTauHad,
        kElTauEl,
        kElTauMu,
        kMuMu,
        kMuTauHad,
        kMuTauEl,
        kMuTauMu,
        kTauElTauEl,
        kTauElTauMu,
        kTauMuTauMu,
        kTauElTauHad,
        kTauMuTauHad,
        kTauHadTauHad,
        kAllDecays
      };

    protected:
      MCParticle::EPartType    ClassifyTau(const MCParticle *tau);
      void                     Process();
      void                     SlaveBegin();
      void                     Terminate();

      TString                  fMCPartName;       //branch name of MCParticle collection
      const MCParticleCol     *fMCParticles;      //!generated particle branch

      const MCParticle        *fHiggsPart;         //!
      const MCParticle        *fHiggsMo1;          //!
      const MCParticle        *fHiggsMo2;          //!
      const MCParticle        *fHiggsDo1;          //!
      const MCParticle        *fHiggsDo2;          //!
      const MCParticle        *fLepton1;           //!
      const MCParticle        *fLepton2;           //!
      const MCParticle        *fNeutrino1;         //!
      const MCParticle        *fNeutrino2;         //!

      EProdType                fProdType;          //!
      EDecayType               fDecayType;         //!
      FourVector               fMissingEnergy;     //!

      TH1D                    *fHiggsPtHist;       //!pt higgs histogram
      TH1D                    *fHiggsEtaHist;      //!eta higgs histogram
      TH1D                    *fHWWPtHist;         //!pt w-boson from higgs histogram
      TH1D                    *fHWWEtaHist;        //!eta w-boson from higgs histogram
      TH1D                    *fWWPtHist;          //!pt w-boson bkg histogram
      TH1D                    *fWWEtaHist;         //!eta w-boson bkg histogram
      TH1D                    *fHiggsMassHist;     //!mass higgs histogram
      TH1D                    *fHiggsTransMHist;   //!transverse mass higgs histogram

      ClassDef(HwwMCEvtSelMod,1) // HWW MC based event selection module
  };
}
#endif
