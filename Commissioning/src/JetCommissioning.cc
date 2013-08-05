// $Id: JetCommissioning.cc,v 1.1 2010/04/02 14:08:15 sixie Exp $

#include "MitHiggs/Commissioning/interface/JetCommissioning.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include <map>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TTree.h>

using namespace mithep;

ClassImp(mithep::JetCommissioning)

//--------------------------------------------------------------------------------------------------
  JetCommissioning::JetCommissioning(const char *name, const char *title) : 
  BaseMod(name,title),
  fVertexName("PrimaryVertexes"),
  fMetName("NotSet"),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fVertices(0),
  fConversions(0),
  fMet(0)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
void JetCommissioning::Process()
{
  // Process entries of the tree. 
  LoadEventObject(fVertexName,     fVertices);
//   LoadBranch(Names::gkMuonBrn);
  LoadBranch(Names::gkTrackBrn);
  LoadBranch(Names::gkCaloJetBrn);
  LoadBranch(Names::gkPFJetBrn);
  LoadBranch(Names::gkTrackJetBrn);

  bool printdebug = false;

  ObjArray<Electron> *CleanElectrons = dynamic_cast<ObjArray<Electron>* >
    (FindObjThisEvt(ModNames::gkCleanElectronsName));
  if (CleanElectrons) {
  } else {
    cout << "Error: Clean Electron Collection " << ModNames::gkCleanElectronsName << " could not be loaded.\n";
  }

  const MetCol *caloMetCol = GetObjThisEvt<MetCol>(ModNames::gkCleanCaloMetName);
  const Met *caloMet = 0;
  if (caloMetCol) {
    caloMet = caloMetCol->At(0);
  } else {
    cout << "Error: Met Collection " << ModNames::gkCleanCaloMetName << " could not be loaded.\n";
  }

  const MetCol *TCMetCol = GetObjThisEvt<MetCol>("pubTCMet");
  const Met *TCMet = 0;
  if (TCMetCol) {
    TCMet = TCMetCol->At(0);
  } else {
    cout << "Error: Met Collection " << "pubTCMet" << " could not be loaded.\n";
  }

  const MetCol *PFMetCol = GetObjThisEvt<MetCol>("pubPFMet");
  const Met *PFMet = 0;
  if (PFMetCol) {
    PFMet = PFMetCol->At(0);
  } else {
    cout << "Error: Met Collection " << "pubPFMet" << " could not be loaded.\n";
  }



  //Event Selection Variables
  Bool_t passDiJetSelection = kFALSE;
  

  //***********************************************************************************************
  //Met
  //***********************************************************************************************
  //cout << "MET : " << caloMet << " " << TCMet << " " << PFMet << endl;

  fCaloMet->Fill(min(caloMet->Pt(),99.9));
  fTCMet->Fill(min(TCMet->Pt(),99.9));
  fPFMet->Fill(min(PFMet->Pt(),99.9));


  //***********************************************************************************************
  //Vertex
  //***********************************************************************************************
  const Vertex *primaryVertex = 0;

  for (UInt_t i=0; i < fVertices->GetEntries(); ++i) {
    if (!primaryVertex || (primaryVertex && primaryVertex->NTracksFit() < fVertices->At(i)->NTracksFit())) {
      primaryVertex = fVertices->At(i);
    }
  }
  if (fVertices->GetEntries() == 1) {
    fPrimaryVertexBeamSpotX->Fill(primaryVertex->X()); 
    fPrimaryVertexBeamSpotY->Fill(primaryVertex->Y()); 
    fPrimaryVertexBeamSpotZ->Fill(primaryVertex->Z()); 
  }
  if (!primaryVertex) return;

  //***********************************************************************************************
  //All Tracks
  //***********************************************************************************************
  for (UInt_t i=0; i< fTracks->GetEntries(); ++i) {
    fTrackNHits->Fill(fTracks->At(i)->NHits());
    if (fTracks->At(i)->NHits() > 8) {
      fTrackPt->Fill(fTracks->At(i)->Pt());
      fTrackEta->Fill(fTracks->At(i)->Eta());
      fTrackPhi->Fill(fTracks->At(i)->Phi());
    }
  }


  //***********************************************************************************************
  //Event Selection Cuts
  //***********************************************************************************************
  //Dijet Selection
  if (0 == 0) {
    passDiJetSelection = kTRUE;
  }
  




  //***********************************************************************************************
  //CaloJets
  //***********************************************************************************************
  Int_t NCaloJets = 0;
  for (UInt_t i=0; i< fCaloJets->GetEntries(); ++i) {


    //kinematics
    fJetPt->Fill(fCaloJets->At(i)->Pt());

    //10GeV cut
    if (fCaloJets->At(i)->Pt() < 10.0) continue;

    NCaloJets++;
    fJetEta->Fill(fCaloJets->At(i)->Eta());
    fJetPhi->Fill(fCaloJets->At(i)->Phi());
    
    fJetEMF->Fill(fCaloJets->At(i)->RestrictedEMF());
    fJetN90Hits->Fill(fCaloJets->At(i)->HitsInN90());
    fJetFHPD->Fill(fCaloJets->At(i)->FHPD());
    fJetFRBX->Fill(fCaloJets->At(i)->FRBX());
    fJetEnergyFractionH->Fill(fCaloJets->At(i)->EnergyFractionH());
    fJetEnergyFractionEm->Fill(fCaloJets->At(i)->EnergyFractionEm());
    fJetSigmaEta->Fill(fCaloJets->At(i)->SigmaEta());
    fJetSigmaPhi->Fill(fCaloJets->At(i)->SigmaPhi());
    
    fJetNConstituents->Fill(fCaloJets->At(i)->NConstituents());

    //get info on charged component
    Double_t ChargedPt = 0;
    Int_t NCharged = 0;
    for (UInt_t t=0; t< fTracks->GetEntries(); ++t) {
      if (fTracks->At(t)->Quality().Quality(TrackQuality::highPurity) && fTracks->At(t)->Pt() > 0.3) {
        ChargedPt += fTracks->At(t)->Pt();
        NCharged++;
      }
    }
    fJetNTracks->Fill(NCharged);
    fJetChargedFraction->Fill(ChargedPt/fCaloJets->At(i)->Pt());

  }
  fNCaloJets->Fill(NCaloJets);





  //***********************************************************************************************
  //PFJets
  //***********************************************************************************************
  Int_t NPFJets = 0;
  for (UInt_t i=0; i< fPFJets->GetEntries(); ++i) {

    //5GeV cut
    if (fPFJets->At(i)->Pt() < 5.0) continue;

    NPFJets++;
    //kinematics
    fPFJetPt->Fill(fPFJets->At(i)->Pt());
    fPFJetEta->Fill(fPFJets->At(i)->Eta());
    fPFJetPhi->Fill(fPFJets->At(i)->Phi());    
    fPFJetSigmaEta->Fill(fPFJets->At(i)->SigmaEta());
    fPFJetSigmaPhi->Fill(fPFJets->At(i)->SigmaPhi());

    fPFJetChargedFraction->Fill((fPFJets->At(i)->ChargedEmEnergy() + fPFJets->At(i)->ChargedHadronEnergy()) / fPFJets->At(i)->E());
    fPFJetChargedEMFraction->Fill(fPFJets->At(i)->ChargedEmEnergy()/fPFJets->At(i)->E());
    fPFJetChargedHadronFraction->Fill(fPFJets->At(i)->ChargedHadronEnergy()/fPFJets->At(i)->E());
    fPFJetNeutralHadronFraction->Fill(fPFJets->At(i)->NeutralHadronEnergy()/fPFJets->At(i)->E());
    fPFJetNeutralEMFraction->Fill(fPFJets->At(i)->NeutralEmEnergy()/fPFJets->At(i)->E());
    fPFJetChargedMultiplicity->Fill(fPFJets->At(i)->ChargedMultiplicity());

  }

  fNPFJets->Fill(NPFJets);


  //***********************************************************************************************
  //TrackJets
  //***********************************************************************************************
  Int_t NTrackJets = 0;
  for (UInt_t i=0; i< fTrackJets->GetEntries(); ++i) {

    if (fTrackJets->At(i)->Pt() < 6.0) continue;

    NTrackJets++;
    //kinematics
    fTrackJetPt->Fill(fTrackJets->At(i)->Pt());
    fTrackJetEta->Fill(fTrackJets->At(i)->Eta());
    fTrackJetPhi->Fill(fTrackJets->At(i)->Phi());    
    fTrackJetSigmaEta->Fill(fTrackJets->At(i)->SigmaEta());
    fTrackJetSigmaPhi->Fill(fTrackJets->At(i)->SigmaPhi());

    fTrackJetNTracks->Fill(fTrackJets->At(i)->NTracks());
    Double_t maxTrackPt = 0;
    for (UInt_t t=0; t< fTrackJets->At(i)->NTracks(); ++t) {
      if (maxTrackPt < fTrackJets->At(i)->Trk(t)->Pt()) {
        maxTrackPt = fTrackJets->At(i)->Trk(t)->Pt();
      }
    }
    fTrackJetMaxTrackPtOverPt->Fill(maxTrackPt/fTrackJets->At(i)->Pt());

  }
  fNTrackJets->Fill(NTrackJets);



}

//--------------------------------------------------------------------------------------------------
void JetCommissioning::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.
  ReqEventObject(fVertexName,                fVertices, kTRUE);
  ReqBranch(Names::gkMuonBrn,                fMuons);
  ReqBranch(Names::gkTrackBrn,               fTracks);
  ReqBranch(Names::gkCaloJetBrn,             fCaloJets);
  ReqBranch(Names::gkPFJetBrn,               fPFJets);
  ReqBranch(Names::gkTrackJetBrn,            fTrackJets);


  //*********************************************************************************************
  //CaloJet Variables
  //*********************************************************************************************  
  AddTH1(fNCaloJets,"hNCaloJets",";Number of CaloJets",20,-0.5,19.5);

  AddTH1(fJetPt  ,"hJetPt", ";p_{T} [GeV/c];Number of Events",200,0,50);
  AddTH1(fJetEta  ,"hJetEta", ";#eta ;Number of Events",200,-3,3);
  AddTH1(fJetPhi  ,"hJetPhi", ";#phi;Number of Events",200,-3.2,3.2);

  //see PAS JME 09-008
  //JetCuts
  //EMF > 0.01, n90hits > 4, fHPD < 0.98, fRBX < 0.98, sigmaEta < 0.01, sigmaPhi < 0.01

  AddTH1(fJetEMF,"hJetEMF",";Jet EM Fraction",500,0,1);
  AddTH1(fJetN90Hits,"hJetN90Hits",";Jet N90Hits",500,0,1);
  AddTH1(fJetFHPD,"hJetFHPD",";Jet HPD Fraction",500,0,1);
  AddTH1(fJetFRBX,"hJetFRBX",";Jet RBX Fraction",500,0,1);
  AddTH1(fJetEnergyFractionH,"hJetEnergyFractionH",";Had Energy Fraction",500,0,1);
  AddTH1(fJetEnergyFractionEm,"hJetEnergyFractionEm",";Em Energy Fraction",500,0,1);
  AddTH1(fJetSigmaEta,"hJetSigmaEta",";SigmaEta",500,0,1);
  AddTH1(fJetSigmaPhi,"hJetSigmaPhi",";SigmaPhi",500,0,1);

  AddTH1(fJetNConstituents,"hJetNConstituents",";Jet NConstituents",500,0,1);
  AddTH1(fJetNTracks,"hJetNTracks",";Jet NTracks",500,0,1);
  AddTH1(fJetChargedFraction,"hJetChargedFraction",";Jet Charged Fraction",500,0,1);
  //energy density: tower energy sum / number of towers



  //*********************************************************************************************
  //PFJet Variables
  //*********************************************************************************************
  AddTH1(fNPFJets,"hNPFJets",";Number of PFJets",20,-0.5,19.5);
  AddTH1(fPFJetPt  ,"hPFJetPt", ";p_{T} [GeV/c];Number of Events",200,0,50);
  AddTH1(fPFJetEta  ,"hPFJetEta", ";#eta ;Number of Events",200,-3,3);
  AddTH1(fPFJetPhi  ,"hPFJetPhi", ";#phi;Number of Events",200,-3.2,3.2);

  //see AN-2010/003 for PF Jet ID
  AddTH1(fPFJetSigmaEta,"hPFJetSigmaEta",";SigmaEta",500,0,1);
  AddTH1(fPFJetSigmaPhi,"hPFJetSigmaPhi",";SigmaPhi",500,0,1);
  AddTH1(fPFJetChargedFraction,"hPFJetChargedFraction",";Charged Fraction",500,0,1);
  AddTH1(fPFJetChargedEMFraction,"hPFJetChargedEMFraction",";Charged EM Fraction",500,0,1);
  AddTH1(fPFJetChargedHadronFraction,"hPFJetChargedHadronFraction",";Charged Hadron Fraction",500,0,1);
  AddTH1(fPFJetNeutralHadronFraction,"hPFJetNeutralHadronFraction",";Neutral Hadron Fraction",500,0,1);
  AddTH1(fPFJetNeutralEMFraction,"hPFJetNeutralEMFraction",";Neutral EM Fraction",500,0,1);
  AddTH1(fPFJetChargedMultiplicity,"hPFJetChargedMultiplicity",";Charged Multiplicity",500,0,1);
  

  //*********************************************************************************************
  //TrackJet Variables
  //*********************************************************************************************
  AddTH1(fNTrackJets,"hNTrackJets",";Number of TrackJets",20,-0.5,19.5);

  AddTH1(fTrackJetPt  ,"hTrackJetPt", ";p_{T} [GeV/c];Number of Events",200,0,50);
  AddTH1(fTrackJetEta  ,"hTrackJetEta", ";#eta ;Number of Events",200,-3,3);
  AddTH1(fTrackJetPhi  ,"hTrackJetPhi", ";#phi;Number of Events",200,-3.2,3.2);
  AddTH1(fTrackJetSigmaEta,"hTrackJetSigmaEta",";SigmaEta",500,0,1);
  AddTH1(fTrackJetSigmaPhi,"hTrackJetSigmaPhi",";SigmaPhi",500,0,1);
  AddTH1(fTrackJetNTracks,"hTrackJetNTracks",";Neutral EM Fraction",100,-0.5,99.5);
  AddTH1(fTrackJetMaxTrackPtOverPt,"hTrackJetMaxTrackPtOverPt",";MaxTrackPtOverPt",200,0,1);





  AddTH1(fCaloMet  ,"hCaloMet", ";Met [GeV/c];Number of Events",500,0,100);
  AddTH1(fTCMet  ,"hTCMet", ";Met [GeV/c];Number of Events",500,0,100);
  AddTH1(fPFMet  ,"hPFMet", ";Met [GeV/c];Number of Events",500,0,100);

  AddTH1(fPrimaryVertexBeamSpotX  ,"hPrimaryVertexBeamSpotX", ";x [cm];Number of Events",200,0,0.5);
  AddTH1(fPrimaryVertexBeamSpotY  ,"hPrimaryVertexBeamSpotY", ";y [cm];Number of Events",200,0,0.5);
  AddTH1(fPrimaryVertexBeamSpotZ  ,"hPrimaryVertexBeamSpotZ", ";z [cm];Number of Events",200,-15,15);

  AddTH1(fTrackPt  ,"hTrackPt", ";p_{T} [GeV/c];Number of Events",200,0,50);
  AddTH1(fTrackEta  ,"hTrackEta", ";#eta ;Number of Events",200,-3,3);
  AddTH1(fTrackPhi  ,"hTrackPhi", ";#phi;Number of Events",200,-3.2,3.2);
  AddTH1(fTrackNHits  ,"hTrackNHits", ";NHits;Number of Events",30,-0.5,29.5);

}

//--------------------------------------------------------------------------------------------------
void JetCommissioning::SlaveTerminate()
{

}
