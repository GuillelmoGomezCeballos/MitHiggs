//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronCommissioning
//
// 
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_COMMISSIONING_ELECTRONCOMMISSIONING_H
#define MITHIGGS_COMMISSIONING_ELECTRONCOMMISSIONING_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitAna/DataTree/interface/BasicClusterCol.h"
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
  class ElectronCommissioning : public BaseMod
  {
    public:
      ElectronCommissioning(const char *name="ElectronCommissioning", 
                       const char *title="Electron Commissioning Module");

      const char     *GetElectronBranchName()      const { return fElectronBranchName;         }
      const char     *GetMetName()                 const { return fMetName;                    }
      const char     *GetCleanJetsName()           const { return fCleanJetsName;              }
      const Double_t  GetSampleType()              const { return fSampleType;                 }

      void   SetIsData (Bool_t b)                        { fIsData                       = b;  }
      void   SetElectronBranchName (TString s)           { fElectronBranchName           = s;  }
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
      TString                  fElectronBranchName;
      TString                  fMetName;                    
      TString                  fCleanJetsName;              
      const VertexCol         *fVertices;                   //!vertices branches
      const DecayParticleCol  *fConversions;                //!conversion collection
      const ElectronCol       *fElectrons;
      const MuonCol           *fMuons;
      const PhotonCol         *fPhotons;
      const MetCol            *fMet;                        //!Missing Et
      const TrackCol          *fTracks;                    
      const TrackCol          *fGsfTracks;                    
      const SuperClusterCol   *fBarrelSuperClusters;
      const SuperClusterCol   *fEndcapSuperClusters;
      const BasicClusterCol   *fBarrelBasicClusters;
      const BasicClusterCol   *fEndcapBasicClusters;

      Int_t                    fElectronsCount;
      Int_t                    fCleanElectronsCount;
      Int_t                    fRealElectronsCount;
      Float_t                  fElectronTreeVariables[54];  //array holds electron tree vars
      TTree                   *fElectronTree;               //electron tree
      


      TH1D *fNBarrelSuperClusters ;
      TH1D *fBarrelSuperClusterEt ;
      TH1D *fBarrelSuperClusterEta  ;
      TH1D *fBarrelSuperClusterPhi  ;
      TH1D *fBarrelSuperClusterESeedOverE  ;
      TH1D *fBarrelSuperClusterEtaWidth ;
      TH1D *fBarrelSuperClusterPhiWidth ;
      TH1D *fBarrelSuperClusterPreshowerEnergy ;
      TH1D *fBarrelSuperClusterSize ;

      TH1D *fNEndcapSuperClusters ;
      TH1D *fEndcapSuperClusterEt ;
      TH1D *fEndcapSuperClusterEta  ;
      TH1D *fEndcapSuperClusterPhi  ;
      TH1D *fEndcapSuperClusterESeedOverE  ;
      TH1D *fEndcapSuperClusterEtaWidth ;
      TH1D *fEndcapSuperClusterPhiWidth ;
      TH1D *fEndcapSuperClusterPreshowerEnergy ;
      TH1D *fEndcapSuperClusterSize ;


      TH1D *fNElectrons  ;

      //Kinematics
      TH1D *fElectronECALDrivenPt  ;
      TH1D *fElectronECALDrivenEta  ;
      TH1D *fElectronECALDrivenPhi  ;
      TH1D *fElectronTrackerDrivenPt  ;
      TH1D *fElectronTrackerDrivenEta  ;
      TH1D *fElectronTrackerDrivenPhi  ;


      //EGamma Critical Plots
      TH1D *fElectronMva ;
      TH1D *fElectronESeed  ;
      TH1D *fElectronEBrem  ;
      TH1D *fElectronETot ;
      TH1D *fElectronGsfPOut ;
      TH1D *fElectronGsfPIn ;
      TH1D *fElectronGsfPInMinusPOut ;
      TH1D *fElectronESuperClusterOverP ;
      TH1D *fElectronDeltaEtaSuperClusterTrackAtVtx ;
      TH1D *fElectronDeltaPhiSuperClusterTrackAtVtx ;
      TH1D *fElectronESeedClusterOverPout ;
      TH1D *fElectronHadronicOverEm ;
      TH1D *fElectronD0;
      TH1D *fElectronTrackIsolationDr04 ;
      TH1D *fElectronEcalRecHitIsoDr04 ;
      TH1D *fElectronHcalTowerSumEtDr04 ;
      TH1D *fElectronEcalRecHitIsoDr03 ;
      TH1D *fElectronHcalTowerSumEtDr03 ;
      TH1D *fElectronTrackIsolationDr03 ;
      TH1D *fElectronCoviEtaiEta ;

      //Additional EGamma GsfTrack Related Plots
      TH1D *fElectronFBrem  ;
      TH1D *fElectronGsfPtOverKFPt  ;

      //Additional Electron Related Plots
      TH1D *fElectronClassification ;
      TH1D *fElectronCovEtaEta ;
      TH1D *fElectronDeltaEtaSeedClusterTrackAtCalo ;
      TH1D *fElectronDeltaPhiSeedClusterTrackAtCalo ;
      TH1D *fElectronE15 ;
      TH1D *fElectronE25Max ;
      TH1D *fElectronE55 ;
      TH1D *fElectronESeedClusterOverPIn ;
      TH1D *fElectronFracSharedHits ;
      TH1D *fElectronHcalDepth1OverEcal ;
      TH1D *fElectronHcalDepth2OverEcal ;
      TH1D *fElectronNumberOfClusters ;
      TH1D *fElectronHcalDepth1TowerSumEtDr04 ;
      TH1D *fElectronHcalDepth2TowerSumEtDr04 ;
      TH1D *fElectronHcalDepth1TowerSumEtDr03 ;
      TH1D *fElectronHcalDepth2TowerSumEtDr03 ;
      TH1D *fElectronIDLikelihood ;

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

      ClassDef(ElectronCommissioning,1) // TAM example analysis module
  };
}
#endif
