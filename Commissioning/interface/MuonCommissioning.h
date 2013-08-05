//--------------------------------------------------------------------------------------------------
// $Id $
//
// MuonCommissioning
//
// 
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_COMMISSIONING_MUONCOMMISSIONING_H
#define MITHIGGS_COMMISSIONING_MUONCOMMISSIONING_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"

class TH1D;
class TH2D;
class TH3D;
class TTree;

namespace mithep 
{
  class MuonCommissioning : public BaseMod
  {
    public:
      MuonCommissioning(const char *name="MuonCommissioning", 
                       const char *title="Electron Commissioning Module");

      const char     *GetMuonBranchName()          const { return fMuonBranchName;             }
      const char     *GetMetName()                 const { return fMetName;                    }
      const char     *GetCleanJetsName()           const { return fCleanJetsName;              }
      const Double_t  GetSampleType()              const { return fSampleType;                 }

      void   SetIsData (Bool_t b)                        { fIsData                       = b;  }
      void   SetMuonBranchName (TString s)               { fMuonBranchName               = s;  }
      void   SetMetName (TString s)                      { fMetName                      = s;  }
      void   SetCleanJetsName (TString s)                { fCleanJetsName                = s;  }
      void   SetSampleType (Double_t n)                  { fSampleType                   = n;  }

    protected:
      void                     Process();
      void                     SlaveBegin();
      void                     SlaveTerminate();

      Bool_t                   fIsData;
      Double_t                 fSampleType;
      TString                  fVertexName;	        //name of vertex collection
      TString                  fConversionBranchName;   //name of electron collection (input)
      TString                  fMuonBranchName;
      TString                  fMetName;                    
      TString                  fCleanJetsName;              
      const VertexCol         *fVertices;                   //!vertices branches
      const DecayParticleCol  *fConversions;                //!conversion collection
      const MuonCol           *fMuons;
      const MetCol            *fMet;                        //!Missing Et
      const TrackCol          *fTracks;                    
      const PhotonCol         *fPhotons;

      Int_t                    fMuonCount;
      Int_t                    fCleanMuonCount;
      Int_t                    fRealMuonCount;
      Float_t                  fMuonTreeVariables[45];  //array holds Muon tree vars
      TTree                   *fMuonTree;               //Muon tree
      


      TH1D *fNMuons  ;
      TH1D *fNGlobalMuons  ;
      TH1D *fNTrackerMuons  ;
      TH1D *fNStandaloneMuons  ;
      TH1D *fNCaloMuons  ;

      //Kinematics
      TH1D *fGlobalMuonPt  ;
      TH1D *fGlobalMuonEta  ;
      TH1D *fGlobalMuonPhi  ;
      TH1D *fTrackerMuonPt  ;
      TH1D *fTrackerMuonEta  ;
      TH1D *fTrackerMuonPhi  ;
      TH1D *fStandaloneMuonPt  ;
      TH1D *fStandaloneMuonEta  ;
      TH1D *fStandaloneMuonPhi  ;
      TH1D *fCaloMuonPt  ;
      TH1D *fCaloMuonEta  ;
      TH1D *fCaloMuonPhi  ;

      //GlobalMuon Critical Plots
      TH1D *fGlobalMuonNHits  ;
      TH1D *fGlobalMuonD0  ;
      TH1D *fGlobalMuonEmEnergy  ;
      TH1D *fGlobalMuonEmS9Energy  ;
      TH1D *fGlobalMuonHadEnergy  ;
      TH1D *fGlobalMuonHadS9Energy  ;
      TH1D *fGlobalMuonHoEnergy  ;
      TH1D *fGlobalMuonHoS9Energy  ;
      TH1D *fGlobalMuonIsoR03SumPt  ;
      TH1D *fGlobalMuonIsoR03EmEt  ;
      TH1D *fGlobalMuonIsoR03HadEt  ;
      TH1D *fGlobalMuonIsoR03HoEt  ;
      TH1D *fGlobalMuonIsoR05SumPt  ;
      TH1D *fGlobalMuonIsoR05EmEt  ;
      TH1D *fGlobalMuonIsoR05HadEt  ;
      TH1D *fGlobalMuonIsoR05HoEt  ;
      TH1D *fGlobalMuonNChambers  ;
      TH1D *fGlobalMuonNSegments  ;
      TH1D *fGlobalMuonTrackChi2OverNdof  ;


      //TrackerMuon Critical Plots
      TH1D *fTrackerMuonNHits  ;
      TH1D *fTrackerMuonD0  ;
      TH1D *fTrackerMuonEmEnergy  ;
      TH1D *fTrackerMuonEmS9Energy  ;
      TH1D *fTrackerMuonHadEnergy  ;
      TH1D *fTrackerMuonHadS9Energy  ;
      TH1D *fTrackerMuonHoEnergy  ;
      TH1D *fTrackerMuonHoS9Energy  ;
      TH1D *fTrackerMuonIsoR03SumPt  ;
      TH1D *fTrackerMuonIsoR03EmEt  ;
      TH1D *fTrackerMuonIsoR03HadEt  ;
      TH1D *fTrackerMuonIsoR03HoEt  ;
      TH1D *fTrackerMuonIsoR05SumPt  ;
      TH1D *fTrackerMuonIsoR05EmEt  ;
      TH1D *fTrackerMuonIsoR05HadEt  ;
      TH1D *fTrackerMuonIsoR05HoEt  ;
      TH1D *fTrackerMuonNChambers  ;
      TH1D *fTrackerMuonNSegments  ;
      TH1D *fTrackerMuonTrackChi2OverNdof  ;

   
      TH1D *fCaloMet  ;
      TH1D *fTCMet  ;
      TH1D *fPFMet  ;
      
      TH1D                    *fPrimaryVertexBeamSpotX;
      TH1D                    *fPrimaryVertexBeamSpotY;
      TH1D                    *fPrimaryVertexBeamSpotZ;
      TH1D                    *fTrackPt;
      TH1D                    *fTrackEta;
      TH1D                    *fTrackPhi;
      TH1D                    *fTrackNHits;

      ClassDef(MuonCommissioning,1);
  };
}
#endif
