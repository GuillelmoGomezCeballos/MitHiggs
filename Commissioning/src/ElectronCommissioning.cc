// $Id: ElectronCommissioning.cc,v 1.1 2010/04/02 14:08:15 sixie Exp $

#include "MitHiggs/Commissioning/interface/ElectronCommissioning.h"
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

ClassImp(mithep::ElectronCommissioning)

//--------------------------------------------------------------------------------------------------
  ElectronCommissioning::ElectronCommissioning(const char *name, const char *title) : 
  BaseMod(name,title),
  fVertexName("PrimaryVertexes"),
  fElectronBranchName(Names::gkElectronBrn),
  fMetName("NotSet"),
  fCleanJetsName(ModNames::gkCleanJetsName),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fVertices(0),
  fConversions(0),
  fElectrons(0),
  fMet(0)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
void ElectronCommissioning::Process()
{
  // Process entries of the tree. 
  LoadEventObject(fVertexName,     fVertices);
  LoadBranch(fElectronBranchName);
  LoadBranch(Names::gkMuonBrn);
  LoadBranch(Names::gkTrackBrn);
  LoadBranch(Names::gkGsfTrackBrn);
  LoadBranch(Names::gkPhotonBrn);
//   LoadBranch(fConversionBranchName);
  LoadBranch(Names::gkBarrelBasicClusterBrn);
  LoadBranch(Names::gkEndcapBasicClusterBrn);
  LoadBranch(Names::gkBarrelSuperClusterBrn);
  LoadBranch(Names::gkEndcapSuperClusterBrn);

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
  if (caloMet->Pt() > 35.0 && CleanElectrons->GetEntries() == 1 && CleanElectrons->At(0)->Pt() > 20.0) {
    passWSelection = kTRUE;
  }
  

  //Z Selection
  for (UInt_t i=0; i< CleanElectrons->GetEntries(); ++i) {
    if (!(CleanElectrons->At(i)->Pt() > 20.0 && CleanElectrons->At(i)->AbsEta() < 2.5))
      continue;
      
    for (UInt_t j=i+1; j< CleanElectrons->GetEntries(); ++j) {
      if (!(CleanElectrons->At(j)->Pt() > 20.0 && CleanElectrons->At(j)->AbsEta() < 2.5))
        continue;

      CompositeParticle *ZBoson = new CompositeParticle;
      ZBoson->AddDaughter(CleanElectrons->At(i));
      ZBoson->AddDaughter(CleanElectrons->At(j));
      if (ZBoson->Mass() > 70 && ZBoson->Mass()<110) {
        passZSelection = kTRUE;
      }
      delete ZBoson;
    }
  }



  //***********************************************************************************************
  //Super Clusters
  //***********************************************************************************************
  fNBarrelSuperClusters->Fill(fBarrelSuperClusters->GetEntries());
  for (UInt_t i=0; i< fBarrelSuperClusters->GetEntries(); ++i) {

    //Additional EGamma Supercluster Related Plots
    fBarrelSuperClusterEta->Fill(fBarrelSuperClusters->At(i)->Eta());
    fBarrelSuperClusterPhi->Fill(fBarrelSuperClusters->At(i)->Phi());
    fBarrelSuperClusterEt->Fill(fBarrelSuperClusters->At(i)->Et());
    fBarrelSuperClusterESeedOverE->Fill(fBarrelSuperClusters->At(i)->Seed()->Energy()/fBarrelSuperClusters->At(i)->Energy());

    //Additional Variables
    fBarrelSuperClusterEtaWidth->Fill(fBarrelSuperClusters->At(i)->EtaWidth());
    fBarrelSuperClusterPhiWidth->Fill(fBarrelSuperClusters->At(i)->PhiWidth());
    fBarrelSuperClusterPreshowerEnergy->Fill(fBarrelSuperClusters->At(i)->PreshowerEnergy());
    fBarrelSuperClusterSize->Fill(fBarrelSuperClusters->At(i)->ClusterSize());

    //H/E variables
  }

  fNEndcapSuperClusters->Fill(fEndcapSuperClusters->GetEntries());
  for (UInt_t i=0; i< fEndcapSuperClusters->GetEntries(); ++i) {

    //Additional EGamma Supercluster Related Plots
    fEndcapSuperClusterEta->Fill(fEndcapSuperClusters->At(i)->Eta());
    fEndcapSuperClusterPhi->Fill(fEndcapSuperClusters->At(i)->Phi());
    fEndcapSuperClusterEt->Fill(fEndcapSuperClusters->At(i)->Et());
    fEndcapSuperClusterESeedOverE->Fill(fEndcapSuperClusters->At(i)->Seed()->Energy()/fEndcapSuperClusters->At(i)->Energy());

    //Additional Variables
    fEndcapSuperClusterEtaWidth->Fill(fEndcapSuperClusters->At(i)->EtaWidth());
    fEndcapSuperClusterPhiWidth->Fill(fEndcapSuperClusters->At(i)->PhiWidth());
    fEndcapSuperClusterPreshowerEnergy->Fill(fEndcapSuperClusters->At(i)->PreshowerEnergy());
    fEndcapSuperClusterSize->Fill(fEndcapSuperClusters->At(i)->ClusterSize());
 
    //H/E variables
 }


  //***********************************************************************************************
  //Electrons
  //***********************************************************************************************

  fNElectrons->Fill(fElectrons->GetEntries());

//   if (fElectrons->GetEntries() > 0 || fMuons->GetEntries() > 0 || fPhotons->GetEntries() > 0) {
//     cout << "LumiSec: " << GetEventHeader()->LumiSec() << " RunNum: " << GetEventHeader()->RunNum()
//          << " EventNum: " << GetEventHeader()->EvtNum() << endl;
//   }
  
  for (UInt_t i=0; i< fElectrons->GetEntries(); ++i) {
    fElectronsCount++;
    
//     cout << "Electron " << fElectrons->At(i)->Pt() << " " << fElectrons->At(i)->Eta() << " " << fElectrons->At(i)->Phi() 
//          << fElectrons->At(i)->FBrem() << " " << fElectrons->At(i)->CoviEtaiEta() << " " 
//          << endl;  
    
    

    //Kinematics
    if (fElectrons->At(i)->IsEcalDriven()) {
      fElectronECALDrivenPt->Fill(fElectrons->At(i)->Pt());
      fElectronECALDrivenEta->Fill(fElectrons->At(i)->Eta());
      fElectronECALDrivenPhi->Fill(fElectrons->At(i)->Phi());
    }
    if (fElectrons->At(i)->IsTrackerDriven()) {
      fElectronTrackerDrivenPt->Fill(fElectrons->At(i)->Pt());
      fElectronTrackerDrivenEta->Fill(fElectrons->At(i)->Eta());
      fElectronTrackerDrivenPhi->Fill(fElectrons->At(i)->Phi());
    }

    //EGamma Critical Plots
    fElectronMva->Fill(fElectrons->At(i)->Mva());
    fElectronESeed->Fill(fElectrons->At(i)->SCluster()->Seed()->Energy());
    fElectronEBrem->Fill(fElectrons->At(i)->SCluster()->Energy()-fElectrons->At(i)->SCluster()->Seed()->Energy());
    fElectronETot->Fill(fElectrons->At(i)->SCluster()->Energy());
    fElectronGsfPOut->Fill(fElectrons->At(i)->POut());
    fElectronGsfPIn->Fill(fElectrons->At(i)->PIn());
    fElectronGsfPInMinusPOut->Fill(fElectrons->At(i)->POut()-fElectrons->At(i)->PIn());
    fElectronESuperClusterOverP->Fill(fElectrons->At(i)->ESuperClusterOverP());
    fElectronDeltaEtaSuperClusterTrackAtVtx->Fill(fElectrons->At(i)->DeltaEtaSuperClusterTrackAtVtx());
    fElectronDeltaPhiSuperClusterTrackAtVtx->Fill(fElectrons->At(i)->DeltaPhiSuperClusterTrackAtVtx());
    fElectronESeedClusterOverPout->Fill(fElectrons->At(i)->ESeedClusterOverPout());
    fElectronHadronicOverEm->Fill(fElectrons->At(i)->HadronicOverEm());
    fElectronD0->Fill(fElectrons->At(i)->GsfTrk()->D0Corrected(*primaryVertex));
    fElectronTrackIsolationDr04->Fill(fElectrons->At(i)->TrackIsolationDr04());
    fElectronEcalRecHitIsoDr04->Fill(fElectrons->At(i)->EcalRecHitIsoDr04());
    fElectronHcalTowerSumEtDr04->Fill(fElectrons->At(i)->HcalTowerSumEtDr04());
    fElectronTrackIsolationDr03->Fill(fElectrons->At(i)->TrackIsolationDr03());
    fElectronEcalRecHitIsoDr03->Fill(fElectrons->At(i)->EcalRecHitIsoDr03());
    fElectronHcalTowerSumEtDr03->Fill(fElectrons->At(i)->HcalTowerSumEtDr03());
    fElectronCoviEtaiEta->Fill(fElectrons->At(i)->CoviEtaiEta());



    //Additional EGamma GsfTrack Related Plots
    fElectronFBrem->Fill(fElectrons->At(i)->FBrem());
    if (fElectrons->At(i)->HasTrackerTrk()) {
      fElectronGsfPtOverKFPt->Fill(fElectrons->At(i)->GsfTrk()->Pt() / fElectrons->At(i)->TrackerTrk()->Pt());
    } else {
      fElectronGsfPtOverKFPt->Fill(0.0);
    }

    //Additional Electron Related Plots
    fElectronClassification->Fill(fElectrons->At(i)->Classification());
    fElectronCovEtaEta->Fill(fElectrons->At(i)->CovEtaEta());
    fElectronDeltaEtaSeedClusterTrackAtCalo->Fill(fElectrons->At(i)->DeltaEtaSeedClusterTrackAtCalo());
    fElectronDeltaPhiSeedClusterTrackAtCalo->Fill(fElectrons->At(i)->DeltaPhiSeedClusterTrackAtCalo());
    fElectronE15->Fill(fElectrons->At(i)->E15());
    fElectronE25Max->Fill(fElectrons->At(i)->E25Max());
    fElectronE55->Fill(fElectrons->At(i)->E55());
    fElectronESeedClusterOverPIn->Fill(fElectrons->At(i)->ESeedClusterOverPIn());
    fElectronFracSharedHits->Fill(fElectrons->At(i)->FracSharedHits());
    fElectronHcalDepth1OverEcal->Fill(fElectrons->At(i)->HcalDepth1OverEcal());
    fElectronHcalDepth2OverEcal->Fill(fElectrons->At(i)->HcalDepth2OverEcal());
    fElectronNumberOfClusters->Fill(fElectrons->At(i)->NumberOfClusters());
    fElectronHcalDepth1TowerSumEtDr04->Fill(fElectrons->At(i)->HcalDepth1TowerSumEtDr04());
    fElectronHcalDepth2TowerSumEtDr04->Fill(fElectrons->At(i)->HcalDepth2TowerSumEtDr04());
    fElectronHcalDepth1TowerSumEtDr03->Fill(fElectrons->At(i)->HcalDepth1TowerSumEtDr03());
    fElectronHcalDepth2TowerSumEtDr03->Fill(fElectrons->At(i)->HcalDepth2TowerSumEtDr03());
    fElectronIDLikelihood->Fill(fElectrons->At(i)->IDLikelihood());


    //initilize 
    for (int k=0; k < 53; ++k) {
      fElectronTreeVariables[k] = 0.0;
    }

    //Kinematics
    fElectronTreeVariables[0] = 1.0;
    fElectronTreeVariables[1] = fSampleType;

    //Kinematics
    fElectronTreeVariables[2] = fElectrons->At(i)->Pt();
    fElectronTreeVariables[3] = fElectrons->At(i)->Eta();
    fElectronTreeVariables[4] = fElectrons->At(i)->Phi();

    //Type
    Int_t type = 0;
    if (fElectrons->At(i)->IsEcalDriven()) type += 1;
    if (fElectrons->At(i)->IsTrackerDriven()) type += 10;
    fElectronTreeVariables[5] = type;
    
    //EB/EE
    Int_t location = 0;
    if (fElectrons->At(i)->IsEB()) location = 1;
    else if (fElectrons->At(i)->IsEE()) location = 2;
    else if (fElectrons->At(i)->IsEBEEGap()) location = 3;
    else if (fElectrons->At(i)->IsEBEtaGap()) location = 4;
    else if (fElectrons->At(i)->IsEBPhiGap()) location = 5;
    else if (fElectrons->At(i)->IsEEDeeGap()) location = 6;
    else if (fElectrons->At(i)->IsEERingGap()) location = 7;
    fElectronTreeVariables[6] = location;
  
    //EGamma Electron Critical Variables
    fElectronTreeVariables[7] = fElectrons->At(i)->SCluster()->Seed()->Energy(); 
    fElectronTreeVariables[8] = fElectrons->At(i)->SCluster()->Energy();
    fElectronTreeVariables[9] = fElectrons->At(i)->POut();
    fElectronTreeVariables[10] = fElectrons->At(i)->PIn();
    fElectronTreeVariables[11] = fElectrons->At(i)->ESuperClusterOverP();
    fElectronTreeVariables[12] = fElectrons->At(i)->DeltaEtaSuperClusterTrackAtVtx();
    fElectronTreeVariables[13] = fElectrons->At(i)->DeltaPhiSuperClusterTrackAtVtx();
    fElectronTreeVariables[14] = fElectrons->At(i)->ESeedClusterOverPout();
    fElectronTreeVariables[15] = fElectrons->At(i)->HadronicOverEm();
    fElectronTreeVariables[16] = fElectrons->At(i)->GsfTrk()->D0Corrected(*primaryVertex);
    fElectronTreeVariables[17] = fElectrons->At(i)->TrackIsolationDr04();
    fElectronTreeVariables[18] = fElectrons->At(i)->EcalRecHitIsoDr04();
    fElectronTreeVariables[19] = fElectrons->At(i)->HcalTowerSumEtDr04();
    fElectronTreeVariables[20] = fElectrons->At(i)->TrackIsolationDr03();
    fElectronTreeVariables[21] = fElectrons->At(i)->EcalRecHitIsoDr03();
    fElectronTreeVariables[22] = fElectrons->At(i)->HcalTowerSumEtDr03();
    fElectronTreeVariables[23] = fElectrons->At(i)->CoviEtaiEta();
    fElectronTreeVariables[24] = fElectrons->At(i)->Mva();

    //Additional EGamma GsfTrack Related Plots
    
    fElectronTreeVariables[25] = fElectrons->At(i)->FBrem();
    fElectronTreeVariables[26] = fElectrons->At(i)->GsfTrk()->Pt();
    if (fElectrons->At(i)->HasTrackerTrk()) {
      fElectronTreeVariables[27] = fElectrons->At(i)->TrackerTrk()->Pt();
    } else {
      fElectronTreeVariables[27] = -99;
    }

    //Additional Electron Related Plots
    fElectronTreeVariables[28] = fElectrons->At(i)->Classification();
    fElectronTreeVariables[29] = fElectrons->At(i)->CovEtaEta();
    fElectronTreeVariables[30] = fElectrons->At(i)->DeltaEtaSeedClusterTrackAtCalo();
    fElectronTreeVariables[31] = fElectrons->At(i)->DeltaPhiSeedClusterTrackAtCalo();
    fElectronTreeVariables[32] = fElectrons->At(i)->E15();
    fElectronTreeVariables[33] = fElectrons->At(i)->E25Max();
    fElectronTreeVariables[34] = fElectrons->At(i)->E55();
    fElectronTreeVariables[35] = fElectrons->At(i)->ESeedClusterOverPIn();
    fElectronTreeVariables[36] = fElectrons->At(i)->FracSharedHits();
    fElectronTreeVariables[37] = fElectrons->At(i)->HcalDepth1OverEcal();
    fElectronTreeVariables[38] = fElectrons->At(i)->HcalDepth2OverEcal();
    fElectronTreeVariables[39] = fElectrons->At(i)->NumberOfClusters();
    fElectronTreeVariables[40] = fElectrons->At(i)->HcalDepth1TowerSumEtDr04();
    fElectronTreeVariables[41] = fElectrons->At(i)->HcalDepth2TowerSumEtDr04();
    fElectronTreeVariables[42] = fElectrons->At(i)->HcalDepth1TowerSumEtDr03();
    fElectronTreeVariables[43] = fElectrons->At(i)->HcalDepth2TowerSumEtDr03();
    fElectronTreeVariables[44] = fElectrons->At(i)->IDLikelihood();


    //Seed Cluster variables
    fElectronTreeVariables[45] = fElectrons->At(i)->SCluster()->Seed()->E1x3();
    fElectronTreeVariables[46] = fElectrons->At(i)->SCluster()->Seed()->E3x1();
    fElectronTreeVariables[47] = fElectrons->At(i)->SCluster()->Seed()->E2x2();    
    fElectronTreeVariables[48] = fElectrons->At(i)->SCluster()->Seed()->E3x3();
    fElectronTreeVariables[49] = fElectrons->At(i)->SCluster()->Seed()->EMax();
    fElectronTreeVariables[50] = fElectrons->At(i)->SCluster()->Seed()->NHits();
    
    //Event Selection Cuts
    fElectronTreeVariables[51] = (passWSelection ? 1 : 0);
    fElectronTreeVariables[52] = (passZSelection ?  1 : 0);
    fElectronTreeVariables[53] = (passJPsiSelection ? 1 : 0);
    fElectronTree->Fill();

  }




}

//--------------------------------------------------------------------------------------------------
void ElectronCommissioning::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.
//   ReqEventObject(fConversionBranchName,      fConversions, kTRUE);
  ReqEventObject(fVertexName,                fVertices, kTRUE);
  ReqBranch(fElectronBranchName,             fElectrons);
  ReqBranch(Names::gkMuonBrn,                fMuons);
  ReqBranch(Names::gkPhotonBrn,              fPhotons);
  ReqBranch(Names::gkTrackBrn,               fTracks);
  ReqBranch(Names::gkGsfTrackBrn,            fGsfTracks);
  ReqBranch(Names::gkBarrelSuperClusterBrn,  fBarrelSuperClusters);
  ReqBranch(Names::gkEndcapSuperClusterBrn,  fEndcapSuperClusters);
  ReqBranch(Names::gkBarrelBasicClusterBrn,  fBarrelBasicClusters);
  ReqBranch(Names::gkEndcapBasicClusterBrn,  fEndcapBasicClusters);


  //*********************************************************************************************
  //Electron Ntuple
  //*********************************************************************************************
  fElectronTree = new TTree("ElectronCommissioningTree", "ElectronCommissioningTree");
  char* TreeFormat;
  TreeFormat = "weight/F:sampleType/F:pt/F:eta/F:phi/F:electronType/F:electronLocation/F:seedEnergy/F:superclusterEnergy/F:pout/F:pin/F:ESuperClusterOverP/F:DeltaEtaSuperClusterTrackAtVtx/F:DeltaPhiSuperClusterTrackAtVtx/F:ESeedClusterOverPout/F:HadronicOverEm/F:D0/F:TrackIsolationDr04/F:EcalRecHitIsoDr04/F:HcalTowerSumEtDr04/F:TrackIsolationDr03/F:EcalRecHitIsoDr03/F:HcalTowerSumEtDr03/F:CoviEtaiEta/F:Mva/F:FBrem/F:GsfTrackPt/F:TrackerTrackPt/F:Classification/F:CovEtaEta/F:DeltaEtaSeedClusterTrackAtCalo/F:DeltaPhiSeedClusterTrackAtCalo/F:E15/F:E25Max/F:E55/F:ESeedClusterOverPIn/F:FracSharedHits/F:HcalDepth1OverEcal/F:HcalDepth2OverEcal/F:NumberOfClusters/F:HcalDepth1TowerSumEtDr04/F:HcalDepth2TowerSumEtDr04/F:HcalDepth1TowerSumEtDr03/F:HcalDepth2TowerSumEtDr03/F:IDLikelihood/F:E1x3/F:E3x1/F:E2x2/F:E3x3/F:EMax/F:NHits/F:passWSelection/F:passZSelection/F:passJPsiSelection/F";
  fElectronTree->Branch("electron_branch", &fElectronTreeVariables,TreeFormat);
  AddOutput(fElectronTree);


  //*********************************************************************************************
  //Super Cluster Histograms
  //*********************************************************************************************
  AddTH1(fNBarrelSuperClusters ,"hNBarrelSuperClusters", ";# of Super Clusters;Number of Events",50,-0.5,49.5);
  AddTH1(fBarrelSuperClusterEt ,"hBarrelSuperClusterEt", ";Super Cluster Et [GeV/c];Number of Superclusters",200,0,200);
  AddTH1(fBarrelSuperClusterEta  ,"hBarrelSuperClusterEta", ";Super Cluster #eta;Number of Superclusters",200,-5,5);
  AddTH1(fBarrelSuperClusterPhi  ,"hBarrelSuperClusterPhi", ";Super Cluster #phi;Number of Superclusters",200,-3.5,3.5);
  AddTH1(fBarrelSuperClusterESeedOverE  ,"hBarrelSuperClusterESeedOverE", ";ESeed/ESuperCluster;Number of Superclusters",200,-3.5,3.5);  
  AddTH1(fBarrelSuperClusterEtaWidth , "hBarrelSuperClusterEtaWidth", ";CoviEtaiEta; # of Superclusters" ,500,0,0.1);
  AddTH1(fBarrelSuperClusterPhiWidth , "hBarrelSuperClusterPhiWidth", ";CoviPhiiPhi; # of Superclusters" ,500,0,0.1);
  AddTH1(fBarrelSuperClusterPreshowerEnergy ,"hBarrelSuperClusterPreshowerEnergy", ";Preshower Energy [GeV];Number of Superclsuters",200,0,200);
  AddTH1(fBarrelSuperClusterSize ,"hBarrelSuperClusterSize", ";# of Constituent Clusters;Number of Superclsuters",20,-0.5,19.5);

  AddTH1(fNEndcapSuperClusters ,"hNEndcapSuperClusters", ";# of Super Clusters;Number of Events",50,-0.5,49.5);
  AddTH1(fEndcapSuperClusterEt ,"hEndcapSuperClusterEt", ";Super Cluster Et [GeV/c];Number of Superclusters",200,0,200);
  AddTH1(fEndcapSuperClusterEta  ,"hEndcapSuperClusterEta", ";Super Cluster #eta;Number of Superclusters",200,-5,5);
  AddTH1(fEndcapSuperClusterPhi  ,"hEndcapSuperClusterPhi", ";Super Cluster #phi;Number of Superclusters",200,-3.5,3.5);
  AddTH1(fEndcapSuperClusterESeedOverE  ,"hEndcapSuperClusterESeedOverE", ";ESeed/ESuperCluster;Number of Superclusters",200,-3.5,3.5);  
  AddTH1(fEndcapSuperClusterEtaWidth , "hEndcapSuperClusterEtaWidth", ";CoviEtaiEta; # of Superclusters" ,500,0,0.1);
  AddTH1(fEndcapSuperClusterPhiWidth , "hEndcapSuperClusterPhiWidth", ";CoviPhiiPhi; # of Superclusters" ,500,0,0.1);
  AddTH1(fEndcapSuperClusterPreshowerEnergy ,"hEndcapSuperClusterPreshowerEnergy", ";Preshower Energy [GeV];Number of Superclsuters",200,0,200);
  AddTH1(fEndcapSuperClusterSize ,"hEndcapSuperClusterSize", ";# of Constituent Clusters;Number of Superclsuters",20,-0.5,19.5);

  //*********************************************************************************************
  //Electron Histograms
  //*********************************************************************************************
   AddTH1(fNElectrons  ,"hNElectrons", ";# of Electrons;Number of Events",10,-0.5,9.5);

   //Kinematics
   AddTH1(fElectronECALDrivenPt  ,"hElectronECALDrivenPt", ";Electron Pt [GeV/c];Number of Events",200,0,200);
   AddTH1(fElectronECALDrivenEta  ,"hElectronECALDrivenEta", ";Electron #eta;Number of Events",200,-5,5);
   AddTH1(fElectronECALDrivenPhi  ,"hElectronECALDrivenPhi", ";Electron Pt [GeV/c];Number of Events",200,0,200);
   AddTH1(fElectronTrackerDrivenPt  ,"hElectronTrackerDrivenPt", ";Electron #phi;Number of Events",200,-3.5,3.5);
   AddTH1(fElectronTrackerDrivenEta  ,"hElectronTrackerDrivenEta", ";Electron #eta;Number of Events",200,-5,5);
   AddTH1(fElectronTrackerDrivenPhi  ,"hElectronTrackerDrivenPhi", ";Electron #phi;Number of Events",200,-3.5,3.5);

   //EGamma Critical Plots
   AddTH1(fElectronMva , "hElectronMva", ";Mva; # of Electrons" ,500,0,1);
   AddTH1(fElectronESeed  ,"hElectronESeed", "; Seed Cluster Energy [GeV];Number of Events",500,0,200);
   AddTH1(fElectronEBrem  ,"hElectronEBrem", "; Bremmed Energy [GeV];Number of Events",500,0,200);
   AddTH1(fElectronETot ,"hElectronETot", "; Energy [GeV];Number of Events",500,0,200);
   AddTH1(fElectronGsfPOut ,"hElectronGsfPOut", "; GsfTrack P Out [GeV];Number of Events",500,0,200);
   AddTH1(fElectronGsfPIn ,"hElectronGsfPIn", "; GsfTrack P In [GeV];Number of Events",500,0,200);
   AddTH1(fElectronGsfPInMinusPOut ,"hElectronGsfPInMinusPOut", "; GsfTrack PIn - POut [GeV];Number of Events",500,0,200);
   AddTH1(fElectronESuperClusterOverP , "hElectronESuperClusterOverP", ";ESuperClusterOverP; # of Electrons" ,500,0,20);
   AddTH1(fElectronDeltaEtaSuperClusterTrackAtVtx , "hElectronDeltaEtaSuperClusterTrackAtVtx", ";DeltaEtaSuperClusterTrackAtVtx; # of Electrons" ,500,-0.2,0.2);
   AddTH1(fElectronDeltaPhiSuperClusterTrackAtVtx , "hElectronDeltaPhiSuperClusterTrackAtVtx", ";DeltaPhiSuperClusterTrackAtVtx; # of Electrons" ,500,-0.5,0.5);
   AddTH1(fElectronESeedClusterOverPout , "hElectronESeedClusterOverPout", ";ESeedClusterOverPout; # of Electrons" ,500,0,100);
   AddTH1(fElectronHadronicOverEm , "hElectronHadronicOverEm", ";HadronicOverEm; # of Electrons" ,500,0,5);
   AddTH1(fElectronD0, "hElectronD0", ";Electron Track D0 [cm]; # of Electrons" ,500,0,5);
   AddTH1(fElectronTrackIsolationDr04 , "hElectronTrackIsolationDr04", ";TrackIsolationDr04; # of Electrons" ,500,0,100);
   AddTH1(fElectronEcalRecHitIsoDr04 , "hElectronEcalRecHitIsoDr04", ";EcalRecHitIsoDr04; # of Electrons" ,500,0,10);
   AddTH1(fElectronHcalTowerSumEtDr04 , "hElectronHcalTowerSumEtDr04", ";HcalTowerSumEtDr04; # of Electrons" ,500,0,20);
   AddTH1(fElectronEcalRecHitIsoDr03 , "hElectronEcalRecHitIsoDr03", ";EcalRecHitIsoDr03; # of Electrons" ,500,0,10);
   AddTH1(fElectronHcalTowerSumEtDr03 , "hElectronHcalTowerSumEtDr03", ";HcalTowerSumEtDr03; # of Electrons" ,500,0,10);
   AddTH1(fElectronTrackIsolationDr03 , "hElectronTrackIsolationDr03", ";TrackIsolationDr03; # of Electrons" ,500,0,10);
   AddTH1(fElectronCoviEtaiEta , "hElectronCoviEtaiEta", ";CoviEtaiEta; # of Electrons" ,500,0,0.1);

   //Additional EGamma GsfTrack Related Plots
   AddTH1(fElectronFBrem  ,"hElectronFBrem", "; Brem Fraction;Number of Events",500,0,2.0);
   AddTH1(fElectronGsfPtOverKFPt  ,"hElectronGsfPtOverKFPt", "; GsfTrack Pt / Kalman Filter Track Pt;Number of Events",500,0,5.0);

   //Additional Electron Related Plots
   AddTH1(fElectronClassification , "hElectronClassification", ";Classification; # of Electrons" ,500,0,10);
   AddTH1(fElectronCovEtaEta , "hElectronCovEtaEta", ";CovEtaEta; # of Electrons" ,500,0,0.1);
   AddTH1(fElectronDeltaEtaSeedClusterTrackAtCalo , "hElectronDeltaEtaSeedClusterTrackAtCalo", 
          ";DeltaEtaSeedClusterTrackAtCalo; # of Electrons" ,500,-0.2,0.2);
   AddTH1(fElectronDeltaPhiSeedClusterTrackAtCalo , "hElectronDeltaPhiSeedClusterTrackAtCalo", 
          ";DeltaPhiSeedClusterTrackAtCalo; # of Electrons" ,500,-0.5,0.5);
   AddTH1(fElectronE15 , "hElectronE15", ";E15; # of Electrons" ,1000,0,100);
   AddTH1(fElectronE25Max , "hElectronE25Max", ";E25Max; # of Electrons" ,1000,0,1000);
   AddTH1(fElectronE55 , "hElectronE55", ";E55; # of Electrons" ,1000,0,1000);
   AddTH1(fElectronESeedClusterOverPIn , "hElectronESeedClusterOverPIn", 
 	 ";ESeedClusterOverPIn; # of Electrons" ,500,0,20);
   AddTH1(fElectronFracSharedHits , "hElectronFracSharedHits", 
 	 ";FracSharedHits; # of Electrons" ,500,0,1);
   AddTH1(fElectronHcalDepth1OverEcal , "hElectronHcalDepth1OverEcal", 
 	 ";HcalDepth1OverEcal; # of Electrons" ,500,0,5);
   AddTH1(fElectronHcalDepth2OverEcal , "hElectronHcalDepth2OverEcal", 
 	 ";HcalDepth2OverEcal; # of Electrons" ,500,0,2);
   AddTH1(fElectronNumberOfClusters , "hElectronNumberOfClusters", 
 	 ";NumberOfClusters; # of Electrons" ,500,-0.5,499.5);
   AddTH1(fElectronHcalDepth1TowerSumEtDr04 , "hElectronHcalDepth1TowerSumEtDr04", 
 	 ";HcalDepth1TowerSumEtDr04; # of Electrons" ,500,0,20);
   AddTH1(fElectronHcalDepth2TowerSumEtDr04 , "hElectronHcalDepth2TowerSumEtDr04", 
 	 ";HcalDepth2TowerSumEtDr04; # of Electrons" ,500,0,10);
   AddTH1(fElectronHcalDepth1TowerSumEtDr03 , "hElectronHcalDepth1TowerSumEtDr03", 
 	 ";HcalDepth1TowerSumEtDr03; # of Electrons" ,500,0,10);
   AddTH1(fElectronHcalDepth2TowerSumEtDr03 , "hElectronHcalDepth2TowerSumEtDr03", 
 	 ";HcalDepth2TowerSumEtDr03; # of Electrons" ,500,0,10);
   AddTH1(fElectronIDLikelihood , "hElectronIDLikelihood", 
          ";IDLikelihood; # of Electrons" ,500,0,1);


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
void ElectronCommissioning::SlaveTerminate()
{
  cout << "Total Electrons Processed : " << fElectronsCount << endl;
  cout << "Total Real Electrons Processed : " << fRealElectronsCount << endl;
  cout << "Total Clean Electrons Processed : " << fCleanElectronsCount << endl;

}
