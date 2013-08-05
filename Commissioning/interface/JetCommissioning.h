//--------------------------------------------------------------------------------------------------
// $Id $
//
// JetCommissioning
//
// 
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITHIGGS_COMMISSIONING_JETCOMMISSIONING_H
#define MITHIGGS_COMMISSIONING_JETCOMMISSIONING_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/TrackJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
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
  class JetCommissioning : public BaseMod
  {
    public:
      JetCommissioning(const char *name="JetCommissioning", 
                       const char *title="Jet Commissioning Module");

      const char     *GetJetBranchName()      const { return fJetBranchName;         }
      const char     *GetMetName()                 const { return fMetName;                    }
      const char     *GetCleanJetsName()           const { return fCleanJetsName;              }
      const Double_t  GetSampleType()              const { return fSampleType;                 }

      void   SetIsData (Bool_t b)                        { fIsData                       = b;  }
      void   SetJetBranchName (TString s)           { fJetBranchName           = s;  }
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
      TString                  fJetBranchName;
      TString                  fMetName;                    
      TString                  fCleanJetsName;              
      const VertexCol         *fVertices;                   //!vertices branches
      const DecayParticleCol  *fConversions;                //!conversion collection
      const CaloJetCol        *fCaloJets;
      const PFJetCol          *fPFJets;
      const TrackJetCol       *fTrackJets;
      const MuonCol           *fMuons;
      const MetCol            *fMet;                        //!Missing Et
      const TrackCol          *fTracks;                    

      Int_t                    fJetsCount;
      Int_t                    fCleanJetsCount;
      Int_t                    fRealJetsCount;
      Float_t                  fJetTreeVariables[48];  //array holds electron tree vars
      TTree                   *fJetTree;               //electron tree
      

      //*********************************************************************************************
      //CaloJet Variables
      //*********************************************************************************************  
      TH1D *fNCaloJets;
      TH1D *fJetPt  ;
      TH1D *fJetEta  ;
      TH1D *fJetPhi  ;
      TH1D *fJetEMF;
      TH1D *fJetN90Hits;
      TH1D *fJetFHPD;
      TH1D *fJetFRBX;
      TH1D *fJetEnergyFractionH;
      TH1D *fJetEnergyFractionEm;
      TH1D *fJetSigmaEta;
      TH1D *fJetSigmaPhi;
      TH1D *fJetNConstituents;
      TH1D *fJetNTracks;
      TH1D *fJetChargedFraction;

      //*********************************************************************************************
      //PFJet Variables
      //*********************************************************************************************
      TH1D *fNPFJets;
      TH1D *fPFJetPt  ;
      TH1D *fPFJetEta  ;
      TH1D *fPFJetPhi  ;
      TH1D *fPFJetSigmaEta;
      TH1D *fPFJetSigmaPhi;
      TH1D *fPFJetChargedFraction;
      TH1D *fPFJetChargedEMFraction;
      TH1D *fPFJetChargedHadronFraction;
      TH1D *fPFJetNeutralHadronFraction;
      TH1D *fPFJetNeutralEMFraction;
      TH1D *fPFJetChargedMultiplicity;
  

      //*********************************************************************************************
      //TrackJet Variables
      //*********************************************************************************************
      TH1D *fNTrackJets;
      TH1D *fTrackJetPt  ;
      TH1D *fTrackJetEta  ;
      TH1D *fTrackJetPhi  ;
      TH1D *fTrackJetSigmaEta;
      TH1D *fTrackJetSigmaPhi;
      TH1D *fTrackJetNTracks;
      TH1D *fTrackJetMaxTrackPtOverPt;


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

      ClassDef(JetCommissioning,1) // TAM example analysis module
  };
}
#endif
