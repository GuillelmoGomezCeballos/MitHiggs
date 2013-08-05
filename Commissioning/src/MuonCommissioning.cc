// $Id: MuonCommissioning.cc,v 1.1 2010/04/02 14:08:15 sixie Exp $

#include "MitHiggs/Commissioning/interface/MuonCommissioning.h"
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

ClassImp(mithep::MuonCommissioning)

//--------------------------------------------------------------------------------------------------
  MuonCommissioning::MuonCommissioning(const char *name, const char *title) : 
  BaseMod(name,title),
  fVertexName("PrimaryVertexes"),
  fMuonBranchName(Names::gkMuonBrn),
  fMetName("NotSet"),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fVertices(0),
  fConversions(0),
  fMuons(0),
  fMet(0)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
void MuonCommissioning::Process()
{
  // Process entries of the tree. 
  LoadEventObject(fVertexName,     fVertices);
  LoadBranch(fMuonBranchName);
  LoadBranch(Names::gkMuonBrn);
  LoadBranch(Names::gkTrackBrn);

  bool printdebug = false;

  ObjArray<Muon> *CleanMuons = dynamic_cast<ObjArray<Muon>* >
    (FindObjThisEvt(ModNames::gkCleanMuonsName));
  if (CleanMuons) {
  } else {
    cout << "Error: Clean Electron Collection " << ModNames::gkCleanMuonsName << " could not be loaded.\n";
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
  Bool_t passWSelection = kFALSE;
  Bool_t passZSelection = kFALSE;
  Bool_t passJPsiSelection = kFALSE;
  

  //***********************************************************************************************
  //Met
  //***********************************************************************************************
  //cout << "MET : " << caloMet << " " << TCMet << " " << PFMet << endl;

  fCaloMet->Fill(caloMet->Pt());
  fTCMet->Fill(TCMet->Pt());
  fPFMet->Fill(PFMet->Pt());


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
  //WSelection
  if (caloMet->Pt() > 35.0 && CleanMuons->GetEntries() == 1 && CleanMuons->At(0)->Pt() > 20.0) {
    passWSelection = kTRUE;
  }
  

  //Z Selection
  for (UInt_t i=0; i< CleanMuons->GetEntries(); ++i) {
    if (!(CleanMuons->At(i)->Pt() > 20.0 && CleanMuons->At(i)->AbsEta() < 2.5))
      continue;
      
    for (UInt_t j=i+1; j< CleanMuons->GetEntries(); ++j) {
      if (!(CleanMuons->At(j)->Pt() > 20.0 && CleanMuons->At(j)->AbsEta() < 2.5))
        continue;

      CompositeParticle *ZBoson = new CompositeParticle;
      ZBoson->AddDaughter(CleanMuons->At(i));
      ZBoson->AddDaughter(CleanMuons->At(j));
      if (ZBoson->Mass() > 70 && ZBoson->Mass()<110) {
        passZSelection = kTRUE;
      }
      delete ZBoson;
    }
  }


  //***********************************************************************************************
  //Muons
  //***********************************************************************************************

  fNMuons->Fill(fMuons->GetEntries());

//   if (fMuons->GetEntries() > 0 || fMuons->GetEntries() > 0 || fPhotons->GetEntries() > 0) {
//     cout << "LumiSec: " << GetEventHeader()->LumiSec() << " RunNum: " << GetEventHeader()->RunNum()
//          << " EventNum: " << GetEventHeader()->EvtNum() << endl;
//   }
  
  Int_t NGlobalMuons = 0;
  Int_t NTrackerMuons = 0;
  Int_t NStandaloneMuons = 0;
  Int_t NCaloMuons = 0;

  for (UInt_t i=0; i< fMuons->GetEntries(); ++i) {
    fMuonCount++;
    
//     cout << "Muon " << fMuons->At(i)->Pt() << " " << fMuons->At(i)->Eta() << " " << fMuons->At(i)->Phi() 
//          << fMuons->At(i)->FBrem() << " " << fMuons->At(i)->CoviEtaiEta() << " " 
//          << endl;  
    
    

    //Kinematics
    if (fMuons->At(i)->IsGlobalMuon()) {
      NGlobalMuons++;
      fGlobalMuonPt->Fill(fMuons->At(i)->Pt());
      fGlobalMuonEta->Fill(fMuons->At(i)->Eta());
      fGlobalMuonPhi->Fill(fMuons->At(i)->Phi());
    }
    else if (fMuons->At(i)->IsTrackerMuon()) {
      NTrackerMuons++;
      fTrackerMuonPt->Fill(fMuons->At(i)->Pt());
      fTrackerMuonEta->Fill(fMuons->At(i)->Eta());
      fTrackerMuonPhi->Fill(fMuons->At(i)->Phi());
    }
    else if (fMuons->At(i)->IsStandaloneMuon()) {
      NStandaloneMuons++;
      fStandaloneMuonPt->Fill(fMuons->At(i)->Pt());
      fStandaloneMuonEta->Fill(fMuons->At(i)->Eta());
      fStandaloneMuonPhi->Fill(fMuons->At(i)->Phi());
    } else if (fMuons->At(i)->IsCaloMuon()) {
      NCaloMuons++;
      fCaloMuonPt->Fill(fMuons->At(i)->Pt());
      fCaloMuonEta->Fill(fMuons->At(i)->Eta());
      fCaloMuonPhi->Fill(fMuons->At(i)->Phi());
    }

    //Muon Plots
    if (fMuons->At(i)->IsGlobalMuon()) {
      fGlobalMuonNHits->Fill(fMuons->At(i)->GlobalTrk()->NHits());
      fGlobalMuonD0->Fill(fMuons->At(i)->GlobalTrk()->D0Corrected(*primaryVertex));
      fGlobalMuonEmEnergy->Fill(fMuons->At(i)->EmEnergy());
      fGlobalMuonEmS9Energy->Fill(fMuons->At(i)->EmS9Energy());
      fGlobalMuonHadEnergy->Fill(fMuons->At(i)->HadEnergy());
      fGlobalMuonHadS9Energy->Fill(fMuons->At(i)->HadS9Energy());
      fGlobalMuonHoEnergy->Fill(fMuons->At(i)->HoEnergy());
      fGlobalMuonHoS9Energy->Fill(fMuons->At(i)->HoS9Energy());
      fGlobalMuonIsoR03SumPt->Fill(fMuons->At(i)->IsoR03SumPt());
      fGlobalMuonIsoR03EmEt->Fill(fMuons->At(i)->IsoR03EmEt());
      fGlobalMuonIsoR03HadEt->Fill(fMuons->At(i)->IsoR03HadEt());
      fGlobalMuonIsoR03HoEt->Fill(fMuons->At(i)->IsoR03HoEt());
      fGlobalMuonIsoR05SumPt->Fill(fMuons->At(i)->IsoR05SumPt());
      fGlobalMuonIsoR05EmEt->Fill(fMuons->At(i)->IsoR05EmEt());
      fGlobalMuonIsoR05HadEt->Fill(fMuons->At(i)->IsoR05HadEt());
      fGlobalMuonIsoR05HoEt->Fill(fMuons->At(i)->IsoR05HoEt());
      fGlobalMuonNChambers->Fill(fMuons->At(i)->NChambers());
      fGlobalMuonNSegments->Fill(fMuons->At(i)->NSegments());
      fGlobalMuonTrackChi2OverNdof->Fill(fMuons->At(i)->GlobalTrk()->Chi2()/fMuons->At(i)->GlobalTrk()->Ndof());
    }

    //Muon Plots
    if (fMuons->At(i)->IsTrackerMuon()) {
      fTrackerMuonNHits->Fill(fMuons->At(i)->TrackerTrk()->NHits());
      fTrackerMuonD0->Fill(fMuons->At(i)->TrackerTrk()->D0Corrected(*primaryVertex));
      fTrackerMuonEmEnergy->Fill(fMuons->At(i)->EmEnergy());
      fTrackerMuonEmS9Energy->Fill(fMuons->At(i)->EmS9Energy());
      fTrackerMuonHadEnergy->Fill(fMuons->At(i)->HadEnergy());
      fTrackerMuonHadS9Energy->Fill(fMuons->At(i)->HadS9Energy());
      fTrackerMuonHoEnergy->Fill(fMuons->At(i)->HoEnergy());
      fTrackerMuonHoS9Energy->Fill(fMuons->At(i)->HoS9Energy());
      fTrackerMuonIsoR03SumPt->Fill(fMuons->At(i)->IsoR03SumPt());
      fTrackerMuonIsoR03EmEt->Fill(fMuons->At(i)->IsoR03EmEt());
      fTrackerMuonIsoR03HadEt->Fill(fMuons->At(i)->IsoR03HadEt());
      fTrackerMuonIsoR03HoEt->Fill(fMuons->At(i)->IsoR03HoEt());
      fTrackerMuonIsoR05SumPt->Fill(fMuons->At(i)->IsoR05SumPt());
      fTrackerMuonIsoR05EmEt->Fill(fMuons->At(i)->IsoR05EmEt());
      fTrackerMuonIsoR05HadEt->Fill(fMuons->At(i)->IsoR05HadEt());
      fTrackerMuonIsoR05HoEt->Fill(fMuons->At(i)->IsoR05HoEt());
      fTrackerMuonNChambers->Fill(fMuons->At(i)->NChambers());
      fTrackerMuonNSegments->Fill(fMuons->At(i)->NSegments());
      fTrackerMuonTrackChi2OverNdof->Fill(fMuons->At(i)->TrackerTrk()->Chi2()/fMuons->At(i)->TrackerTrk()->Ndof());
    }
    
    //initilize 
    for (int k=0; k < 45; ++k) {
      fMuonTreeVariables[k] = 0.0;
    }

    //Kinematics
    fMuonTreeVariables[0] = 1.0;
    fMuonTreeVariables[1] = fSampleType;

    //Kinematics
    fMuonTreeVariables[2] = fMuons->At(i)->Pt();
    fMuonTreeVariables[3] = fMuons->At(i)->Eta();
    fMuonTreeVariables[4] = fMuons->At(i)->Phi();

    //Type
    Int_t type = 0;
    if (fMuons->At(i)->IsGlobalMuon()) type += 1;
    if (fMuons->At(i)->IsTrackerMuon()) type += 10;
    if (fMuons->At(i)->IsStandaloneMuon()) type += 100;
    if (fMuons->At(i)->IsCaloMuon()) type += 1000;
    fMuonTreeVariables[5] = type;
      
    //EGamma Muon Critical Variables
    
    fMuonTreeVariables[6] = fMuons->At(i)->BestTrk()->NHits();
    fMuonTreeVariables[7] = fMuons->At(i)->BestTrk()->D0Corrected(*primaryVertex);
    fMuonTreeVariables[8] = fMuons->At(i)->EmEnergy();
    fMuonTreeVariables[9] = fMuons->At(i)->EmS9Energy();
    fMuonTreeVariables[10] = fMuons->At(i)->HadEnergy();
    fMuonTreeVariables[11] = fMuons->At(i)->HadS9Energy();
    fMuonTreeVariables[12] = fMuons->At(i)->HoEnergy();
    fMuonTreeVariables[13] = fMuons->At(i)->HoS9Energy();
    fMuonTreeVariables[14] = fMuons->At(i)->IsoR03SumPt();
    fMuonTreeVariables[15] = fMuons->At(i)->IsoR03EmEt();
    fMuonTreeVariables[16] = fMuons->At(i)->IsoR03HadEt();
    fMuonTreeVariables[17] = fMuons->At(i)->IsoR03HoEt();
    fMuonTreeVariables[18] = fMuons->At(i)->IsoR05SumPt();
    fMuonTreeVariables[19] = fMuons->At(i)->IsoR05EmEt();
    fMuonTreeVariables[20] = fMuons->At(i)->IsoR05HadEt();
    fMuonTreeVariables[21] = fMuons->At(i)->IsoR05HoEt();
    fMuonTreeVariables[22] = fMuons->At(i)->NChambers();
    fMuonTreeVariables[23] = fMuons->At(i)->NSegments();
    fMuonTreeVariables[24] = fMuons->At(i)->BestTrk()->Chi2()/fMuons->At(i)->BestTrk()->Ndof();

    //Muon ID
    fMuonTreeVariables[25] = (fMuons->At(i)->Quality().Quality(MuonQuality::TrackerMuonArbitrated) ? 1 : 0);
    fMuonTreeVariables[26] = (fMuons->At(i)->Quality().Quality(MuonQuality::AllArbitrated) ? 1 : 0);
    fMuonTreeVariables[27] = (fMuons->At(i)->Quality().Quality(MuonQuality::GlobalMuonPromptTight) ? 1 : 0);
    fMuonTreeVariables[28] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMLastStationLoose) ? 1 : 0);
    fMuonTreeVariables[29] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMLastStationTight) ? 1 : 0);
    fMuonTreeVariables[30] = (fMuons->At(i)->Quality().Quality(MuonQuality::TM2DCompatibilityLoose) ? 1 : 0);
    fMuonTreeVariables[31] = (fMuons->At(i)->Quality().Quality(MuonQuality::TM2DCompatibilityTight) ? 1 : 0);
    fMuonTreeVariables[32] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMOneStationLoose) ? 1 : 0);
    fMuonTreeVariables[33] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMOneStationTight) ? 1 : 0);
    fMuonTreeVariables[34] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMLastStationOptimizedLowPtLoose) ? 1 : 0);
    fMuonTreeVariables[35] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMLastStationOptimizedLowPtTight) ? 1 : 0);
    fMuonTreeVariables[36] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMLastStationAngLoose) ? 1 : 0);
    fMuonTreeVariables[37] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMLastStationAngTight) ? 1 : 0);
    fMuonTreeVariables[38] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMOneStationAngLoose) ? 1 : 0);
    fMuonTreeVariables[39] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMOneStationAngTight) ? 1 : 0);
    fMuonTreeVariables[40] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMLastStationOptimizedBarrelLowPtLoose) ? 1 : 0);
    fMuonTreeVariables[41] = (fMuons->At(i)->Quality().Quality(MuonQuality::TMLastStationOptimizedBarrelLowPtTight) ? 1 : 0);
    fMuonTreeVariables[42] = (passWSelection ? 1 : 0);
    fMuonTreeVariables[43] = (passZSelection ?  1 : 0);
    fMuonTreeVariables[44] = (passJPsiSelection ? 1 : 0);
    fMuonTree->Fill();

  }

  fNGlobalMuons->Fill( NGlobalMuons);
  fNTrackerMuons->Fill( NTrackerMuons);
  fNStandaloneMuons->Fill( NStandaloneMuons);
  fNCaloMuons ->Fill(NCaloMuons);


}

//--------------------------------------------------------------------------------------------------
void MuonCommissioning::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.
//   ReqEventObject(fConversionBranchName,      fConversions, kTRUE);
  ReqEventObject(fVertexName,                fVertices, kTRUE);
  ReqBranch(fMuonBranchName,                 fMuons);
  ReqBranch(Names::gkMuonBrn,                fMuons);
  ReqBranch(Names::gkPhotonBrn,              fPhotons);
  ReqBranch(Names::gkTrackBrn,               fTracks);


  //*********************************************************************************************
  //Muon Ntuple
  //*********************************************************************************************
  fMuonTree = new TTree("MuonCommissioningTree", "MuonCommissioningTree");
  char* TreeFormat;
  TreeFormat = "weight/F:sampleType/F:pt/F:eta/F:phi/F:muonType/F:nhits/F:d0/F:EmEnergy/F:EmS9Energy/F:HadEnergy/F:HadS9Energy/F:HoEnergy/F:HoS9Energy/F:IsoR03SumPt/F:IsoR03EmEt/F:IsoR03HadEt/F:IsoR03HoEt/F:IsoR05SumPt/F:IsoR05EmEt/F:IsoR05HadEt/F:IsoR05HoEt/F:NChambers/F:NSegments/F:Chi2PerDoF/F:TrackerMuonArbitrated/F:AllArbitrated /F:GlobalMuonPromptTight/F:TMLastStationLoose/F:TMLastStationTight/F:TM2DCompatibilityLoose/F:TM2DCompatibilityTight/F:TMOneStationLoose/F:TMOneStationTight/F:TMLastStationOptimizedLowPtLoose/F:TMLastStationOptimizedLowPtTight/F:TMLastStationAngLoose/F:TMLastStationAngTight/F:TMOneStationAngLoose/F:TMOneStationAngTight/F:TMLastStationOptimizedBarrelLowPtLoose/F:TMLastStationOptimizedBarrelLowPtTight/F:passWSelection/F:passZSelection/F:passJPsiSelection/F";
  fMuonTree->Branch("muon_branch", &fMuonTreeVariables,TreeFormat);
  AddOutput(fMuonTree);


                                                    


  //*********************************************************************************************
  //Muon Histograms
  //*********************************************************************************************
   AddTH1(fNMuons  ,"hNMuons", ";# of Muons;Number of Events",10,-0.5,9.5);
   AddTH1(fNGlobalMuons  ,"hNGlobalMuons", ";# of Global Muons;Number of Events",10,-0.5,9.5);
   AddTH1(fNTrackerMuons  ,"hNTrackerMuons", ";# of Tracker Muons;Number of Events",10,-0.5,9.5);
   AddTH1(fNStandaloneMuons  ,"hNStandaloneMuons", ";# of Standalone Muons;Number of Events",10,-0.5,9.5);
   AddTH1(fNCaloMuons  ,"hNCaloMuons", ";# of Calo Muons;Number of Events",10,-0.5,9.5);

   //Kinematics
   AddTH1(fGlobalMuonPt  ,"hGlobalMuonPt", ";Muon Pt [GeV/c];Number of Events",200,0,200);
   AddTH1(fGlobalMuonEta  ,"hGlobalMuonEta", ";Muon #eta;Number of Events",200,-5,5);
   AddTH1(fGlobalMuonPhi  ,"hGlobalMuonPhi", ";Muon Pt [GeV/c];Number of Events",200,0,200);
   AddTH1(fTrackerMuonPt  ,"hTrackerMuonPt", ";Muon Pt [GeV/c];Number of Events",200,0,200);
   AddTH1(fTrackerMuonEta  ,"hTrackerMuonEta", ";Muon #eta;Number of Events",200,-5,5);
   AddTH1(fTrackerMuonPhi  ,"hTrackerMuonPhi", ";Muon Pt [GeV/c];Number of Events",200,0,200);
   AddTH1(fStandaloneMuonPt  ,"hStandaloneMuonPt", ";Muon Pt [GeV/c];Number of Events",200,0,200);
   AddTH1(fStandaloneMuonEta  ,"hStandaloneMuonEta", ";Muon #eta;Number of Events",200,-5,5);
   AddTH1(fStandaloneMuonPhi  ,"hStandaloneMuonPhi", ";Muon Pt [GeV/c];Number of Events",200,0,200);
   AddTH1(fCaloMuonPt  ,"hCaloMuonPt", ";Muon Pt [GeV/c];Number of Events",200,0,200);
   AddTH1(fCaloMuonEta  ,"hCaloMuonEta", ";Muon #eta;Number of Events",200,-5,5);
   AddTH1(fCaloMuonPhi  ,"hCaloMuonPhi", ";Muon Pt [GeV/c];Number of Events",200,0,200);

   //GlobalMuon Critical Plots
   AddTH1(fGlobalMuonNHits  ,"hGlobalMuonNHits", ";# of Hits;Number of GlobalMuons",100,-0.5,99.5);
   AddTH1(fGlobalMuonD0  ,"hGlobalMuonD0", "; D0 [cm];Number of GlobalMuons",500,0,5);
   AddTH1(fGlobalMuonEmEnergy  , "hGlobalMuonEmEnergy", "", 500, 0, 10); 
   AddTH1(fGlobalMuonEmS9Energy  , "hGlobalMuonEmS9Energy", "", 500, 0, 10); 
   AddTH1(fGlobalMuonHadEnergy  , "hGlobalMuonHadEnergy", "", 500, 0, 10);
   AddTH1(fGlobalMuonHadS9Energy  , "hGlobalMuonHadS9Energy", "", 500, 0, 10);
   AddTH1(fGlobalMuonHoEnergy  , "hGlobalMuonHoEnergy", "", 500, 0, 10); 
   AddTH1(fGlobalMuonHoS9Energy  , "hGlobalMuonHoS9Energy", "", 500, 0, 10);
   AddTH1(fGlobalMuonIsoR03SumPt  , "hGlobalMuonIsoR03SumPt", "", 500, 0, 10);
   AddTH1(fGlobalMuonIsoR03EmEt  , "hGlobalMuonIsoR03EmEt", "", 500, 0, 10);
   AddTH1(fGlobalMuonIsoR03HadEt  , "hGlobalMuonIsoR03HadEt", "", 500, 0, 10);
   AddTH1(fGlobalMuonIsoR03HoEt  , "hGlobalMuonIsoR03HoEt", "", 500, 0, 10);
   AddTH1(fGlobalMuonIsoR05SumPt  , "hGlobalMuonIsoR05SumPt", "", 500, 0, 10);
   AddTH1(fGlobalMuonIsoR05EmEt  , "hGlobalMuonIsoR05EmEt", "", 500, 0, 10);
   AddTH1(fGlobalMuonIsoR05HadEt  , "hGlobalMuonIsoR05HadEt", "", 500, 0, 10);
   AddTH1(fGlobalMuonIsoR05HoEt  , "hGlobalMuonIsoR05HoEt", "", 500, 0, 10);
   AddTH1(fGlobalMuonNChambers  , "hGlobalMuonNChambers", "", 20, -0.5, 19.5);
   AddTH1(fGlobalMuonNSegments  , "hGlobalMuonNSegments", "", 20, -0.5, 19.5);
   AddTH1(fGlobalMuonTrackChi2OverNdof  , "hGlobalMuonTrackChi2OverNdof", "", 500, 0, 20);


   //TrackerMuon Critical Plots
   AddTH1(fTrackerMuonNHits  ,"hTrackerMuonNHits", ";# of Hits;Number of TrackerMuons",100,-0.5,99.5);
   AddTH1(fTrackerMuonD0  ,"hTrackerMuonD0", "; D0 [cm];Number of TrackerMuons",500,0,5);
   AddTH1(fTrackerMuonEmEnergy  , "hTrackerMuonEmEnergy", "", 500, 0, 10); 
   AddTH1(fTrackerMuonEmS9Energy  , "hTrackerMuonEmS9Energy", "", 500, 0, 10); 
   AddTH1(fTrackerMuonHadEnergy  , "hTrackerMuonHadEnergy", "", 500, 0, 10);
   AddTH1(fTrackerMuonHadS9Energy  , "hTrackerMuonHadS9Energy", "", 500, 0, 10);
   AddTH1(fTrackerMuonHoEnergy  , "hTrackerMuonHoEnergy", "", 500, 0, 10); 
   AddTH1(fTrackerMuonHoS9Energy  , "hTrackerMuonHoS9Energy", "", 500, 0, 10);
   AddTH1(fTrackerMuonIsoR03SumPt  , "hTrackerMuonIsoR03SumPt", "", 500, 0, 10);
   AddTH1(fTrackerMuonIsoR03EmEt  , "hTrackerMuonIsoR03EmEt", "", 500, 0, 10);
   AddTH1(fTrackerMuonIsoR03HadEt  , "hTrackerMuonIsoR03HadEt", "", 500, 0, 10);
   AddTH1(fTrackerMuonIsoR03HoEt  , "hTrackerMuonIsoR03HoEt", "", 500, 0, 10);
   AddTH1(fTrackerMuonIsoR05SumPt  , "hTrackerMuonIsoR05SumPt", "", 500, 0, 10);
   AddTH1(fTrackerMuonIsoR05EmEt  , "hTrackerMuonIsoR05EmEt", "", 500, 0, 10);
   AddTH1(fTrackerMuonIsoR05HadEt  , "hTrackerMuonIsoR05HadEt", "", 500, 0, 10);
   AddTH1(fTrackerMuonIsoR05HoEt  , "hTrackerMuonIsoR05HoEt", "", 500, 0, 10);
   AddTH1(fTrackerMuonNChambers  , "hTrackerMuonNChambers", "", 20, -0.5, 19.5);
   AddTH1(fTrackerMuonNSegments  , "hTrackerMuonNSegments", "", 20, -0.5, 19.5);
   AddTH1(fTrackerMuonTrackChi2OverNdof  , "hTrackerMuonTrackChi2OverNdof", "", 500, 0, 20);

   
   AddTH1(fCaloMet  ,"hCaloMet", ";Met [GeV/c];Number of Events",200,0,40);
   AddTH1(fTCMet  ,"hTCMet", ";Met [GeV/c];Number of Events",200,0,40);
   AddTH1(fPFMet  ,"hPFMet", ";Met [GeV/c];Number of Events",200,0,40);
   
   AddTH1(fPrimaryVertexBeamSpotX  ,"hPrimaryVertexBeamSpotX", ";x [cm];Number of Events",200,0,0.5);
   AddTH1(fPrimaryVertexBeamSpotY  ,"hPrimaryVertexBeamSpotY", ";y [cm];Number of Events",200,0,0.5);
   AddTH1(fPrimaryVertexBeamSpotZ  ,"hPrimaryVertexBeamSpotZ", ";z [cm];Number of Events",200,-15,15);
   
   AddTH1(fTrackPt  ,"hTrackPt", ";p_{T} [GeV/c];Number of Events",200,0,50);
   AddTH1(fTrackEta  ,"hTrackEta", ";#eta ;Number of Events",200,-3,3);
   AddTH1(fTrackPhi  ,"hTrackPhi", ";#phi;Number of Events",200,-3.2,3.2);
   AddTH1(fTrackNHits  ,"hTrackNHits", ";NHits;Number of Events",30,-0.5,29.5);
   
}

//--------------------------------------------------------------------------------------------------
void MuonCommissioning::SlaveTerminate()
{
  cout << "Total Muons Processed : " << fMuonCount << endl;
  cout << "Total Real Muons Processed : " << fRealMuonCount << endl;
  cout << "Total Clean Muons Processed : " << fCleanMuonCount << endl;

}
