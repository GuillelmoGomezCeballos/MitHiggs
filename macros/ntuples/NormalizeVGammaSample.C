#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <utility>
#include "TFile.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Smurf/Core/SmurfTree.h"
#include <vector>
#include <algorithm>
#include "MitAna/Utils/interface/SimpleTable.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include "Smurf/Analysis/HWWlvlv/factors.h"
#include <TROOT.h>
#include <TChain.h>

void NormalizeVGammaSample(TString  InputFilename0 = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/wgamma.root",
                           TString  InputFilename1 = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/wgammafo.root",
			   bool applyCorrection = true
  ) {

  Double_t eventCount0 = 0;

  TH1D *hDRatioPhotonElectron = new TH1D("hDRatioPhotonElectron","hDRatioPhotonElectron",10,0.0,2.5); hDRatioPhotonElectron->Sumw2();

  TH1D *hDLepSel[10];
  hDLepSel[0] = new TH1D("hDLepSel_0","hDLepSel_0",20,0.0,100.0);
  hDLepSel[1] = new TH1D("hDLepSel_1","hDLepSel_1",20,0.0,100.0);
  hDLepSel[2] = new TH1D("hDLepSel_2","hDLepSel_2",10,0.0,2.5);
  hDLepSel[3] = new TH1D("hDLepSel_3","hDLepSel_3",10,0.0,2.5);
  hDLepSel[4] = new TH1D("hDLepSel_4","hDLepSel_4",20,0.0,100.0);
  hDLepSel[5] = new TH1D("hDLepSel_5","hDLepSel_5",20,0.0,100.0);
  hDLepSel[6] = new TH1D("hDLepSel_6","hDLepSel_6",10,0.0,2.5);
  hDLepSel[7] = new TH1D("hDLepSel_7","hDLepSel_7",10,0.0,2.5);
  for(int i=0; i<8; i++) hDLepSel[i]->Sumw2();

  TH1D *hDpSel[20];
  hDpSel[0]  = new TH1D("hDpSel_0" ,"hDpSel_0" ,20,0.0,100.0);
  hDpSel[1]  = new TH1D("hDpSel_1" ,"hDpSel_1" ,20,0.0,100.0);
  hDpSel[2]  = new TH1D("hDpSel_2" ,"hDpSel_2" ,10,0.0,2.5);
  hDpSel[3]  = new TH1D("hDpSel_3" ,"hDpSel_3" ,10,0.0,2.5);
  hDpSel[4]  = new TH1D("hDpSel_4" ,"hDpSel_4" ,40,0.0,200.0);
  hDpSel[5]  = new TH1D("hDpSel_5" ,"hDpSel_5" ,40,0.0,200.0);
  hDpSel[6]  = new TH1D("hDpSel_6" ,"hDpSel_6" ,40,0.0,200.0);
  hDpSel[7]  = new TH1D("hDpSel_7" ,"hDpSel_7" ,36,0.0,160.0);
  hDpSel[8]  = new TH1D("hDpSel_8" ,"hDpSel_8" ,40,0.0,200.0);
  hDpSel[9]  = new TH1D("hDpSel_9" ,"hDpSel_9" ,40,0.0,200.0);
  hDpSel[10] = new TH1D("hDpSel_10","hDpSel_10",20,0.0,100.0);
  hDpSel[11] = new TH1D("hDpSel_11","hDpSel_11",20,0.0,100.0);
  hDpSel[12] = new TH1D("hDpSel_12","hDpSel_12",10,0.0,2.5);
  hDpSel[13] = new TH1D("hDpSel_13","hDpSel_13",10,0.0,2.5);
  hDpSel[14] = new TH1D("hDpSel_14","hDpSel_14",40,0.0,200.0);
  hDpSel[15] = new TH1D("hDpSel_15","hDpSel_15",40,0.0,200.0);
  hDpSel[16] = new TH1D("hDpSel_16","hDpSel_16",40,0.0,200.0);
  hDpSel[17] = new TH1D("hDpSel_17","hDpSel_17",36,0.0,160.0);
  hDpSel[18] = new TH1D("hDpSel_18","hDpSel_18",40,0.0,200.0);
  hDpSel[19] = new TH1D("hDpSel_19","hDpSel_19",40,0.0,200.0);
  for(int i=0; i<20; i++) hDpSel[i]->Sumw2();

  SmurfTree vgammaEvent;
  vgammaEvent.LoadTree(InputFilename0.Data());
  vgammaEvent.InitTree(0);
  vgammaEvent.tree_->SetName("tree");

  TFile *fRatioPhotonElectron = TFile::Open("/data/smurf/data/Run2012_Summer12_SmurfV9_53X/auxiliar/ratio_photon_electron.root");
  TH1D *fhDRatioPhotonElectron = (TH1D*)(fRatioPhotonElectron->Get("hDRatioPhotonElectron"));
  assert(fhDRatioPhotonElectron);
  fhDRatioPhotonElectron->SetDirectory(0);
  fRatioPhotonElectron->Close();
  delete fRatioPhotonElectron;

  for (int n=0;n<vgammaEvent.tree_->GetEntries();n++) { 
    vgammaEvent.tree_->GetEntry(n);

    bool lid = (vgammaEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection 
           &&  (vgammaEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection;
    if (!lid) continue;

    if( vgammaEvent.lep1_.pt() <= 20	    	     ) continue; // cut on leading lepton pt
    if( vgammaEvent.lep2_.pt() <= 10	    	     ) continue; // cut on trailing lepton pt
    if( vgammaEvent.dilep_.M() <= 12                 ) continue; // cut on minimum dilepton mass

    Double_t add = vgammaEvent.sfWeightPU_*vgammaEvent.sfWeightEff_*vgammaEvent.sfWeightTrig_;
    Double_t myWeight = vgammaEvent.scale1fb_*add;    

    if(!(TMath::Abs(vgammaEvent.lep1McId_) == 11 || TMath::Abs(vgammaEvent.lep1McId_) == 13)) hDLepSel[0]->Fill(TMath::Min(vgammaEvent.lep1_.pt(),99.999),myWeight);
    if(!(TMath::Abs(vgammaEvent.lep2McId_) == 11 || TMath::Abs(vgammaEvent.lep2McId_) == 13)) hDLepSel[0]->Fill(TMath::Min(vgammaEvent.lep2_.pt(),99.999),myWeight);
    if(!(TMath::Abs(vgammaEvent.lep1McId_) == 11 || TMath::Abs(vgammaEvent.lep1McId_) == 13)) hDLepSel[2]->Fill(TMath::Min(TMath::Abs(vgammaEvent.lep1_.eta()),2.499),myWeight);
    if(!(TMath::Abs(vgammaEvent.lep2McId_) == 11 || TMath::Abs(vgammaEvent.lep2McId_) == 13)) hDLepSel[2]->Fill(TMath::Min(TMath::Abs(vgammaEvent.lep2_.eta()),2.499),myWeight);
    if((TMath::Abs(vgammaEvent.lep1McId_) == 11 || TMath::Abs(vgammaEvent.lep1McId_) == 13)) hDLepSel[4]->Fill(TMath::Min(vgammaEvent.lep1_.pt(),99.999),myWeight);
    if((TMath::Abs(vgammaEvent.lep2McId_) == 11 || TMath::Abs(vgammaEvent.lep2McId_) == 13)) hDLepSel[4]->Fill(TMath::Min(vgammaEvent.lep2_.pt(),99.999),myWeight);
    if((TMath::Abs(vgammaEvent.lep1McId_) == 11 || TMath::Abs(vgammaEvent.lep1McId_) == 13)) hDLepSel[6]->Fill(TMath::Min(TMath::Abs(vgammaEvent.lep1_.eta()),2.499),myWeight);
    if((TMath::Abs(vgammaEvent.lep2McId_) == 11 || TMath::Abs(vgammaEvent.lep2McId_) == 13)) hDLepSel[6]->Fill(TMath::Min(TMath::Abs(vgammaEvent.lep2_.eta()),2.499),myWeight);

    eventCount0 += myWeight;

  }

  cout << "eventCount0 Event Count 1fb^-1 : " << eventCount0 << endl;
 
  Double_t eventCount1 = 0;

  SmurfTree vgammafoEvent;
  vgammafoEvent.LoadTree(InputFilename1.Data());
  vgammafoEvent.InitTree(0);
  vgammafoEvent.tree_->SetName("tree");

  for (int n=0;n<vgammafoEvent.tree_->GetEntries();n++) { 
    vgammafoEvent.tree_->GetEntry(n);

    bool lid = ((vgammafoEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection 
            &&  (vgammafoEvent.cuts_ & SmurfTree::Lep2LooseEleV2) == SmurfTree::Lep2LooseEleV2) ||
	       ((vgammafoEvent.cuts_ & SmurfTree::Lep1LooseEleV2) == SmurfTree::Lep1LooseEleV2 
           &&   (vgammafoEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection);
    if (!lid) continue;

    if( vgammafoEvent.lep1_.pt() <= 20	    	     ) continue; // cut on leading lepton pt
    if( vgammafoEvent.lep2_.pt() <= 10	    	     ) continue; // cut on trailing lepton pt
    if( TMath::Abs(vgammafoEvent.lep1_.eta()) >= 2.5 ) continue; // cut on leading lepton pt
    if( TMath::Abs(vgammafoEvent.lep2_.eta()) >= 2.5 ) continue; // cut on trailing lepton pt
    if( vgammafoEvent.dilep_.M() <= 12               ) continue; // cut on minimum dilepton mass

    Double_t add = vgammafoEvent.sfWeightPU_*vgammafoEvent.sfWeightEff_*vgammafoEvent.sfWeightTrig_;
    Double_t myWeight = vgammafoEvent.scale1fb_*add;    

    hDpSel[0]->Fill(TMath::Min(vgammafoEvent.lep1_.pt(),99.999),myWeight);
    hDpSel[1]->Fill(TMath::Min(vgammafoEvent.lep2_.pt(),99.999),myWeight);
    hDpSel[2]->Fill(TMath::Min(TMath::Abs(vgammafoEvent.lep1_.eta()),2.499),myWeight);
    hDpSel[3]->Fill(TMath::Min(TMath::Abs(vgammafoEvent.lep2_.eta()),2.499),myWeight);
    hDpSel[4]->Fill(TMath::Min((double)vgammafoEvent.mt_,199.999),myWeight);
    hDpSel[5]->Fill(TMath::Min((double)vgammafoEvent.dilep_.M(),199.999),myWeight);
    hDpSel[6]->Fill(TMath::Min((double)vgammafoEvent.met_,199.999),myWeight);
    hDpSel[7]->Fill(TMath::Min((double)vgammafoEvent.dPhi_*180/TMath::Pi(),179.999),myWeight);
    hDpSel[8]->Fill(TMath::Min((double)vgammafoEvent.mt1_,199.999),myWeight);
    hDpSel[9]->Fill(TMath::Min((double)vgammafoEvent.mt2_,199.999),myWeight);
    if(!(TMath::Abs(vgammafoEvent.lep1McId_) == 11 || TMath::Abs(vgammafoEvent.lep1McId_) == 13) && applyCorrection == true) myWeight = myWeight * ratioPhotonElectron(fhDRatioPhotonElectron,TMath::Abs(vgammafoEvent.lep1_.eta()));
    if(!(TMath::Abs(vgammafoEvent.lep2McId_) == 11 || TMath::Abs(vgammafoEvent.lep2McId_) == 13) && applyCorrection == true) myWeight = myWeight * ratioPhotonElectron(fhDRatioPhotonElectron,TMath::Abs(vgammafoEvent.lep2_.eta()));
    hDpSel[10]->Fill(TMath::Min(vgammafoEvent.lep1_.pt(),99.999),myWeight);
    hDpSel[11]->Fill(TMath::Min(vgammafoEvent.lep2_.pt(),99.999),myWeight);
    hDpSel[12]->Fill(TMath::Min(TMath::Abs(vgammafoEvent.lep1_.eta()),2.499),myWeight);
    hDpSel[13]->Fill(TMath::Min(TMath::Abs(vgammafoEvent.lep2_.eta()),2.499),myWeight);
    hDpSel[14]->Fill(TMath::Min((double)vgammafoEvent.mt_,199.999),myWeight);
    hDpSel[15]->Fill(TMath::Min((double)vgammafoEvent.dilep_.M(),199.999),myWeight);
    hDpSel[16]->Fill(TMath::Min((double)vgammafoEvent.met_,199.999),myWeight);
    hDpSel[17]->Fill(TMath::Min((double)vgammafoEvent.dPhi_*180/TMath::Pi(),179.999),myWeight);
    hDpSel[18]->Fill(TMath::Min((double)vgammafoEvent.mt1_,199.999),myWeight);
    hDpSel[19]->Fill(TMath::Min((double)vgammafoEvent.mt2_,199.999),myWeight);

    if(!(TMath::Abs(vgammafoEvent.lep1McId_) == 11 || TMath::Abs(vgammafoEvent.lep1McId_) == 13)) hDLepSel[1]->Fill(TMath::Min(vgammafoEvent.lep1_.pt(),99.999),myWeight);
    if(!(TMath::Abs(vgammafoEvent.lep2McId_) == 11 || TMath::Abs(vgammafoEvent.lep2McId_) == 13)) hDLepSel[1]->Fill(TMath::Min(vgammafoEvent.lep2_.pt(),99.999),myWeight);
    if(!(TMath::Abs(vgammafoEvent.lep1McId_) == 11 || TMath::Abs(vgammafoEvent.lep1McId_) == 13)) hDLepSel[3]->Fill(TMath::Min(TMath::Abs(vgammafoEvent.lep1_.eta()),2.499),myWeight);
    if(!(TMath::Abs(vgammafoEvent.lep2McId_) == 11 || TMath::Abs(vgammafoEvent.lep2McId_) == 13)) hDLepSel[3]->Fill(TMath::Min(TMath::Abs(vgammafoEvent.lep2_.eta()),2.499),myWeight);
    if((TMath::Abs(vgammafoEvent.lep1McId_) == 11 || TMath::Abs(vgammafoEvent.lep1McId_) == 13)) hDLepSel[5]->Fill(TMath::Min(vgammafoEvent.lep1_.pt(),99.999),myWeight);
    if((TMath::Abs(vgammafoEvent.lep2McId_) == 11 || TMath::Abs(vgammafoEvent.lep2McId_) == 13)) hDLepSel[5]->Fill(TMath::Min(vgammafoEvent.lep2_.pt(),99.999),myWeight);
    if((TMath::Abs(vgammafoEvent.lep1McId_) == 11 || TMath::Abs(vgammafoEvent.lep1McId_) == 13)) hDLepSel[7]->Fill(TMath::Min(TMath::Abs(vgammafoEvent.lep1_.eta()),2.499),myWeight);
    if((TMath::Abs(vgammafoEvent.lep2McId_) == 11 || TMath::Abs(vgammafoEvent.lep2McId_) == 13)) hDLepSel[7]->Fill(TMath::Min(TMath::Abs(vgammafoEvent.lep2_.eta()),2.499),myWeight);

    eventCount1 += myWeight;

  }

  cout << "eventCount1 Event Count 1fb^-1 : " << eventCount1 << endl;
 
  double NormalizationWeight = eventCount0/eventCount1;
  cout << "Normalized  Event Count 1fb^-1 : " << NormalizationWeight << endl;

  hDLepSel[0]->Divide(hDLepSel[1]);
  hDLepSel[2]->Divide(hDLepSel[3]);
  hDLepSel[4]->Divide(hDLepSel[5]);
  hDLepSel[6]->Divide(hDLepSel[7]);

  hDRatioPhotonElectron->Add(hDLepSel[2]);
  TFile *weightPE = new TFile("ratio_photon_electron.root", "RECREATE");
  weightPE->cd();
  hDRatioPhotonElectron->Write();
  weightPE->Close();

  //*************************************************************************************************
  // Create new normalized tree
  //*************************************************************************************************
  SmurfTree newEvent;
  newEvent.LoadTree(InputFilename1.Data());
  newEvent.InitTree(0);
  newEvent.tree_->SetName("tree");

  TString OutputFilename = InputFilename1;
  OutputFilename.ReplaceAll(".root","_newweight.root");
  cout << "creating new root file: " << OutputFilename << endl;
  TFile *outputFile = new TFile(OutputFilename.Data(), "RECREATE");
  outputFile->cd();

  for(int i=0; i<8; i++) hDLepSel[i]->Write();
  for(int i=0; i<20; i++) hDpSel[i]->Scale(1./hDpSel[i]->GetSumOfWeights());
  for(int i=0; i<20; i++) hDpSel[i]->Write();

  TTree *normalizedTree = newEvent.tree_->CloneTree(0);  
  for (int n=0;n<newEvent.tree_->GetEntries();n++) { 
    newEvent.tree_->GetEntry(n);
    
    newEvent.scale1fb_ = newEvent.scale1fb_*NormalizationWeight;
    if((vgammafoEvent.cuts_ & SmurfTree::Lep1LooseEleV2) == SmurfTree::Lep1LooseEleV2) newEvent.cuts_ |= SmurfTree::Lep1FullSelection;
    if((vgammafoEvent.cuts_ & SmurfTree::Lep2LooseEleV2) == SmurfTree::Lep2LooseEleV2) newEvent.cuts_ |= SmurfTree::Lep2FullSelection;
    if((vgammafoEvent.cuts_ & SmurfTree::Lep3LooseEleV2) == SmurfTree::Lep3LooseEleV2) newEvent.cuts_ |= SmurfTree::Lep3FullSelection;
    normalizedTree->Fill(); 
  }

  normalizedTree->Write();
  outputFile->Close();

}
