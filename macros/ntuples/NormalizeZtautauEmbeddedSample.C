#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <utility>
#include "TFile.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1F.h"
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

//*************************************************************************************************
//Main part of the macro
//*************************************************************************************************
void NormalizeZtautauEmbeddedSample(const string dyttMCFile    = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/dyll.root",
                                    const string InputFilename = "/data/smurf/data/Run2012_Summer12_SmurfV9_53X/mitf-alljets/data_ztt.root"
  ) {

  const int NPoints = 11;
  TH1D* hDOpt0[NPoints];
  hDOpt0[0]  = new TH1D("hDOpt0_0", "electron pt", 100, 0, 100); 	hDOpt0[0] ->Sumw2();
  hDOpt0[1]  = new TH1D("hDOpt0_1", "muon pt", 100, 0, 100); 	        hDOpt0[1] ->Sumw2();
  hDOpt0[2]  = new TH1D("hDOpt0_2", "projected MET", 100, 0, 100); 	hDOpt0[2] ->Sumw2();
  hDOpt0[3]  = new TH1D("hDOpt0_3", "projected trackMET", 100, 0, 100); hDOpt0[3] ->Sumw2();
  hDOpt0[4]  = new TH1D("hDOpt0_4", "min(pMET,pTrackMET)", 100, 0, 100);hDOpt0[4] ->Sumw2();
  hDOpt0[5]  = new TH1D("hDOpt0_5", "MT", 100, 0, 100); 	        hDOpt0[5] ->Sumw2();
  hDOpt0[6]  = new TH1D("hDOpt0_6", "MLL", 100, 0, 100);        	hDOpt0[6] ->Sumw2();
  hDOpt0[7]  = new TH1D("hDOpt0_7", "PTLL", 100, 0, 100); 	        hDOpt0[7] ->Sumw2();
  hDOpt0[8]  = new TH1D("hDOpt0_8", "Njets",  10,-0.5, 9.5);            hDOpt0[8] ->Sumw2();
  hDOpt0[9]  = new TH1D("hDOpt0_9", "TopVeto",   2,-0.5, 1.5);          hDOpt0[9] ->Sumw2();
  hDOpt0[10] = new TH1D("hDOpt0_10","Nvtx", 60, 0,  60); 	        hDOpt0[10]->Sumw2();
  TH1D* hDOpt1[NPoints];
  hDOpt1[0]  = new TH1D("hDOpt1_0", "electron pt", 100, 0, 100);	hDOpt1[0] ->Sumw2();
  hDOpt1[1]  = new TH1D("hDOpt1_1", "muon pt", 100, 0, 100);		hDOpt1[1] ->Sumw2();
  hDOpt1[2]  = new TH1D("hDOpt1_2", "projected MET", 100, 0, 100);	hDOpt1[2] ->Sumw2();
  hDOpt1[3]  = new TH1D("hDOpt1_3", "projected trackMET", 100, 0, 100); hDOpt1[3] ->Sumw2();
  hDOpt1[4]  = new TH1D("hDOpt1_4", "min(pMET,pTrackMET)", 100, 0, 100);hDOpt1[4] ->Sumw2();
  hDOpt1[5]  = new TH1D("hDOpt1_5", "MT", 100, 0, 100); 		hDOpt1[5] ->Sumw2();
  hDOpt1[6]  = new TH1D("hDOpt1_6", "MLL", 100, 0, 100);		hDOpt1[6] ->Sumw2();
  hDOpt1[7]  = new TH1D("hDOpt1_7", "PTLL", 100, 0, 100);		hDOpt1[7] ->Sumw2();
  hDOpt1[8]  = new TH1D("hDOpt1_8", "Njets",  10,-0.5, 9.5);		hDOpt1[8] ->Sumw2();
  hDOpt1[9]  = new TH1D("hDOpt1_9", "TopVeto",   2,-0.5, 1.5);  	hDOpt1[9] ->Sumw2();
  hDOpt1[10] = new TH1D("hDOpt1_10","Nvtx", 60, 0,  60);		hDOpt1[10]->Sumw2();

  //*************************************************************************************************
  //Count dytt MC Normalization
  //*************************************************************************************************
  Double_t DYMCEventcount_emu = 0;

  SmurfTree dyttEvent;
  dyttEvent.LoadTree(dyttMCFile.c_str());
  dyttEvent.InitTree(0);
  dyttEvent.tree_->SetName("tree");

  for (int n=0;n<dyttEvent.tree_->GetEntries();n++) { 
    dyttEvent.tree_->GetEntry(n);

    if (!(((dyttEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
       && ((dyttEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
    if( dyttEvent.lid3_ != 0              ) continue; // cut on third lepton
    if( dyttEvent.lq1_*dyttEvent.lq2_ > 0 ) continue; // cut on opposite-dytt leptons
    if( dyttEvent.dilep_.mass() <= 12.0   ) continue; // cut on low dilepton mass
    if( dyttEvent.lep1_.pt() <= 20	  ) continue; // cut on leading lepton pt
    if( dyttEvent.lep2_.pt() <= 10	  ) continue; // cut on trailing lepton pt

    //MC -> Data Corrections
    
    Double_t add = dyttEvent.sfWeightPU_*dyttEvent.sfWeightEff_*dyttEvent.sfWeightTrig_;
    Double_t myWeight = dyttEvent.scale1fb_*(1.0)*add;    

    if (dyttEvent.type_ == SmurfTree::em || dyttEvent.type_ == SmurfTree::me) {
      DYMCEventcount_emu += myWeight;
      if(dyttEvent.type_ == SmurfTree::em){
        hDOpt0[0]->Fill(TMath::Min((double)dyttEvent.lep1_.Pt(),99.999),myWeight);
        hDOpt0[1]->Fill(TMath::Min((double)dyttEvent.lep2_.Pt(),99.999),myWeight);
      } else {
        hDOpt0[0]->Fill(TMath::Min((double)dyttEvent.lep2_.Pt(),99.999),myWeight);
        hDOpt0[1]->Fill(TMath::Min((double)dyttEvent.lep1_.Pt(),99.999),myWeight);
      }
      hDOpt0[2]->Fill(TMath::Min((double)dyttEvent.pmet_,99.999),myWeight);
      hDOpt0[3]->Fill(TMath::Min((double)dyttEvent.pTrackMet_,99.999),myWeight);
      hDOpt0[4]->Fill(TMath::Min((double)dyttEvent.pmet_,(double)dyttEvent.pTrackMet_),myWeight);
      hDOpt0[5]->Fill(TMath::Min((double)dyttEvent.mt_,99.999),myWeight);
      hDOpt0[6]->Fill(TMath::Min((double)dyttEvent.dilep_.mass(),99.999),myWeight);
      hDOpt0[7]->Fill(TMath::Min((double)dyttEvent.dilep_.Pt(),99.999),myWeight);
      hDOpt0[8]->Fill(TMath::Min((double)dyttEvent.njets_,9.499),myWeight);
      hDOpt0[9]->Fill((double)(dyttEvent.cuts_ & SmurfTree::TopVeto) ==  SmurfTree::TopVeto,myWeight);
      hDOpt0[10]->Fill(TMath::Min((double)dyttEvent.nvtx_,59.999),myWeight);
    }

  }

  cout << "DYtt Event Count (e-mu) 1fb^-1 : " << DYMCEventcount_emu << endl;

  //*************************************************************************************************
  //Count embedded Normalization
  //*************************************************************************************************
  Double_t EmbeddedEventcount_emu = 0;

  SmurfTree embeddedEvent;
  embeddedEvent.LoadTree(InputFilename.c_str());
  embeddedEvent.InitTree(0);
  embeddedEvent.tree_->SetName("tree");

  for (int n=0;n<embeddedEvent.tree_->GetEntries();n++) { 
    embeddedEvent.tree_->GetEntry(n);

    if (!(((embeddedEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
       && ((embeddedEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
    if( embeddedEvent.lid3_ != 0                  ) continue; // cut on third lepton
    if( embeddedEvent.lq1_*embeddedEvent.lq2_ > 0 ) continue; // cut on opposite-embedded leptons
    if( embeddedEvent.dilep_.mass() <= 12.0       ) continue; // cut on low dilepton mass
    if( embeddedEvent.lep1_.pt() <= 20	    	  ) continue; // cut on leading lepton pt
    if( embeddedEvent.lep2_.pt() <= 10	    	  ) continue; // cut on trailing lepton pt

    if (embeddedEvent.type_ == SmurfTree::em || embeddedEvent.type_ == SmurfTree::me) {
      
      Double_t myWeight = embeddedEvent.scale1fb_*embeddedEvent.sfWeightEff_*embeddedEvent.sfWeightTrig_;    
      EmbeddedEventcount_emu = EmbeddedEventcount_emu + myWeight;
      if(embeddedEvent.type_ == SmurfTree::em){
        hDOpt1[0]->Fill(TMath::Min((double)embeddedEvent.lep1_.Pt(),99.999),myWeight);
        hDOpt1[1]->Fill(TMath::Min((double)embeddedEvent.lep2_.Pt(),99.999),myWeight);
      } else {
        hDOpt1[0]->Fill(TMath::Min((double)embeddedEvent.lep2_.Pt(),99.999),myWeight);
        hDOpt1[1]->Fill(TMath::Min((double)embeddedEvent.lep1_.Pt(),99.999),myWeight);
      }
      hDOpt1[2]->Fill(TMath::Min((double)embeddedEvent.pmet_,99.999),myWeight);
      hDOpt1[3]->Fill(TMath::Min((double)embeddedEvent.pTrackMet_,99.999),myWeight);
      hDOpt1[4]->Fill(TMath::Min((double)embeddedEvent.pmet_,(double)embeddedEvent.pTrackMet_),myWeight);
      hDOpt1[5]->Fill(TMath::Min((double)embeddedEvent.mt_,99.999),myWeight);
      hDOpt1[6]->Fill(TMath::Min((double)embeddedEvent.dilep_.mass(),99.999),myWeight);
      hDOpt1[7]->Fill(TMath::Min((double)embeddedEvent.dilep_.Pt(),99.999),myWeight);
      hDOpt1[8]->Fill(TMath::Min((double)embeddedEvent.njets_,9.499),myWeight);
      hDOpt1[9]->Fill((double)(embeddedEvent.cuts_ & SmurfTree::TopVeto) ==  SmurfTree::TopVeto,myWeight);
      hDOpt1[10]->Fill(TMath::Min((double)embeddedEvent.nvtx_,59.999),myWeight);
    }

  }

  cout << "Embedded Event Count (e-mu) : " << EmbeddedEventcount_emu << endl;
  Double_t NormalizationWeight = DYMCEventcount_emu / EmbeddedEventcount_emu ;
  cout << "Normalization Weight = " << NormalizationWeight << endl;
  for(int i=0; i<NPoints; i++) hDOpt1[i]->Scale(NormalizationWeight);

  TString OutputFilename = "test.root";
  TFile *outputFile = new TFile(OutputFilename.Data(), "RECREATE");
  outputFile->cd();
  for(int i=0; i<NPoints; i++) hDOpt0[i]->Write();
  for(int i=0; i<NPoints; i++) hDOpt1[i]->Write();
  outputFile->Close();
}
