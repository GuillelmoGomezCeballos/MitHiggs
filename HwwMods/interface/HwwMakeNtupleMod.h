//------------------------------------------------------------------------------
// $Id: HwwMakeNtupleMod.h,v 1.40 2013/07/13 08:06:26 ceballos Exp $
//
// HwwMakeNtupleMod
//
// A Module to produce a flat ntuple (2l+X analyses)
//
//
// Authors: ceballos
//------------------------------------------------------------------------------

#ifndef HWWMODS_HwwMakeNtupleMod_H
#define HWWMODS_HwwMakeNtupleMod_H
#include <iomanip>
#include <iostream>
#include <fstream>
#include "TRandom.h"

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include <TFile.h>
#include "TH2D.h"
#include "Smurf/Core/SmurfTree.h"
#include "MitPhysics/Utils/interface/PUReweighting.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/MuonIDMVA.h"
#include "MitPhysics/Utils/interface/MVAMet.h"
#include "MitAna/DataTree/interface/LheWeightCol.h"

class TH1D;
class TH2D;
class TTree;

namespace mithep 
{
  class HwwMakeNtupleMod : public BaseMod
  {
    public:
    HwwMakeNtupleMod(const char *name="HwwMakeNtupleMod", 
		 const char *title="Example analysis module with all branches");
      ~HwwMakeNtupleMod();

      const char *GetJetBranchName()                      { return fCaloJetName0;                 }
      const char *GetPFJetBranchName()                    { return fPFJetName0;                   }
      const char *GetPFMetName()                          { return fPFMetName;                    }
      Double_t    GetProcessID()                          { return fDecay;                        }
      Double_t    GetJetScaleSyst()                       { return fJetScaleSyst;                 } 
      Double_t    GetPtJetCut()                           { return fPtJetCut;                     }
      Double_t    GetEtaJetCut()                          { return fEtaJetCut;                    }

      void   SetFakeRatePredictionType(Int_t d)           { fFakeRatePredictionType       = d;    }
      void   SetFillNtupleType(Int_t d)                   { fFillNtupleType               = d;    }
      void   SetJetBranchName(const char *name)           { fCaloJetName0                 = name; }
      void   SetPFJetBranchName(const char *name)         { fPFJetName0                   = name; }
      void   SetPFTauName(TString s)                      { fPFTauName                    = s;    } 
      void   SetPFMetName(TString s)                      { fPFMetName                    = s;    } 
      void   SetProcessID(Double_t x)                     { fDecay                        = x;    }
      void   SetJetScaleSyst(Double_t x)                  { fJetScaleSyst                 = x;    }
      void   SetPtJetCut(Double_t x)                      { fPtJetCut                     = x;    }
      void   SetEtaJetCut(Double_t x)                     { fEtaJetCut                    = x;    }
      void   SetCleanJetsNoPtCutName (TString s)          { fCleanJetsNoPtCutName         = s;    }
      void   SetIsData(Bool_t b)                          { fIsData                       = b;    }	 
      void   SetFillPhotonTemplate(Bool_t b)              { fFillPhotonTemplate           = b;    }	 
      void   SetDoPileupReweighting(Bool_t b)             { fDoPileupReweighting          = b;    }	 
      void   SetOutputName(const char *f)                 { fOutputName                   = f;    }
      void   SetMuonFakeName(const char *name)  	  { fMuonFakeName    	          = name; }
      void   SetElectronFakeName(const char *name)  	  { fElectronFakeName	          = name; }
      void   SetLeptonFakeName(const char *name)  	  { fLeptonFakeName  	          = name; }
      void   SetIntRadius(Double_t dr)                    { fIntRadius                    = dr;   }
      void   SetIs42x(Bool_t b)                           { fIs42x                        = b;    }
      void   SetCorrectedJetsName(TString s)              { fCorrectedJetsName            = s;    }   
      void   SetMVAElVersion(int i)                       { fMVAElVersion                 = i;    }
      void   SetMVAMuVersion(int i)                       { fMVAMuVersion                 = i;    }
      void   SetRhoType(RhoUtilities::RhoType type)       { fTheRhoType                   = type; };
      void   SetAddLheWeights(bool b)                     { fAddLheWeights                = b; };

    protected:
      Int_t             fFakeRatePredictionType;       //Which kind of Fake Rate prediction
                                                       // 0 == leptons, 1 == leptons+fake, 2 == fakes
      Int_t             fFillNtupleType;               //Which kind of fill ntuple we want
                                                       // 0-2 == full, 3-5 == smaller
      Double_t          fPtJetCut;
      Double_t          fEtaJetCut;
      TString           fPFTauName;
      TString           fPFMetName;
      const PFMetCol   *fPFMetStd;
      TString           fMuonName;	               // name of muon collection
      TString           fElectronName;	               // name of electron collection
      TString           fTrackName;	               // name of track collection
      TString           fPFCandidatesName;
      const PFCandidateCol *fPFCandidates;
      TString           fCleanJetsNoPtCutName;         // name of clean central jets collection with no pt requirement
      TString           fMCqqHsName;                   // name of qq jets at gen level
      TString           fEvtHdrName;	               // name of event header branch
      TString           fVertexName;	               // name of vertex collection
      TString           fMCEvInfoName;                 //event info branch name
      const             MCEventInfo *fMCEventInfo;     //!event info branch pointer
      TString           fCaloJetName0;
      TString           fPFJetName0;
      const CaloJetCol  *fCaloJet0;
      const PFJetCol    *fPFJet0;
      const MuonCol     *fMuons;                       // Muon branch
      const ElectronCol *fElectrons;                   // Electron branch
      const TrackCol    *fTracks;                      // Track branch     
      const EventHeader *fEventHeader;                 // event header for current event
      const VertexCol   *fVertices;                    // Vertices branches
      const PileupInfoCol *fPileupInfos;               // PileupInfo branch
      const MCParticleCol *fParticles;                 // MC particle collection handle
      Double_t          fDecay;                        // code for MC process ID
      Double_t          fJetScaleSyst;
      Bool_t            fIsData;                       //=true then it does nothing (def=0)
      Bool_t            fFillPhotonTemplate;           //=true if filling photon met templates
      Bool_t            fDoPileupReweighting;          //=true if filling photon met templates
      TString           fPileupEnergyDensityName;
      const PileupEnergyDensityCol *fPileupEnergyDensity;
      TFile            *fOutputFile;	   	       //output file handle
      TString           fOutputName;	   	       //output file name
      TString           fMuonFakeName;
      TString           fElectronFakeName;
      TString           fLeptonFakeName;
      Double_t          fIntRadius;
      Bool_t            fIs42x;
      ElectronIDMVA    *fElectronIDMVA;
      MuonTools        *fMuonTools;	      // interface to tools for muon ID
      MuonIDMVA        *fMuonIDMVA;	      // helper class for MuonMVA
      TString           fCorrectedJetsName;
      Int_t             fMVAElVersion;
      Int_t             fMVAMuVersion;
      RhoUtilities::RhoType fTheRhoType;
      MVAMet           *fMVAMet;
      Bool_t            fAddLheWeights;
      TString           fLheWeightsName;
      const LheWeightCol *fLheWeights;
      Int_t             fNEventsSelected;              //selected events

      PUReweighting    *fPUReweighting;  
      SmurfTree         fSmurfTree;
      vector<TH2D*>     fDYNNLOKFactorHists;           //vector of hist for Drell-Yan NNLO Kfactor 

      TH2D             *LoadHisto(const char *fname, TFile *file) const;

      void      Begin();
      void      Process();
      void      SlaveBegin();
      void      SlaveTerminate();
      void      Terminate();	  

      ClassDef(HwwMakeNtupleMod,1) // TAM example analysis module
  };
}
#endif
