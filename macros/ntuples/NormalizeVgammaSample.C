//
//root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/NormalizeHwwNtuple.C+\(\"/home/sixie/CMSSW_3_1_2/src/HwwNtuple_s8-h160ww2l-id11_noskim_0000.root\",\"s8-h160ww2l-id11\",\"test.root\"\) 


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
void NormalizeVgammaSample(int nsel = 0) {
  string Input1Filename = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets/wgamma.root";
  string Input2Filename = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets/wgammafo.root";
  if(nsel == 1 || nsel == 2){
    Input1Filename = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets/zgamma.root";
    Input2Filename = "/data/smurf/data/Run2012_Summer12_SmurfV9_52X/mitf-alljets/zgammafo.root";
  }
  //*************************************************************************************************
  //Count Wgamma MC Normalization
  //*************************************************************************************************
  Double_t WGMCEventcount  = 0;
  Double_t WGMCEventcountE = 0;

  SmurfTree vgammaEvent;
  vgammaEvent.LoadTree(Input1Filename.c_str());
  vgammaEvent.InitTree(0);
  vgammaEvent.tree_->SetName("tree");

  for (int n=0;n<vgammaEvent.tree_->GetEntries();n++) { 
    vgammaEvent.tree_->GetEntry(n);

    bool cuts = true;
    if      (nsel == 0){
      cuts = cuts && vgammaEvent.lep1_.pt() > 20 && vgammaEvent.lep2_.pt() > 10 && vgammaEvent.lid3_ == 0 &&
      (((vgammaEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) &&
       ((vgammaEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection));
    }
    else if(nsel == 1){
      cuts = cuts && vgammaEvent.lep1_.pt() > 20 && vgammaEvent.lep2_.pt() > 10 && vgammaEvent.lep3_.pt() > 10 && vgammaEvent.lid3_ != 0 &&
       (((vgammaEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) &&
        ((vgammaEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) &&
        ((vgammaEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection));
    }
    else if(nsel == 2){
      cuts = cuts && vgammaEvent.lep1_.pt() > 20 && vgammaEvent.lep2_.pt() > 10 && vgammaEvent.lep3_.pt() > 10 && vgammaEvent.lid3_ != 0 &&
       (((vgammaEvent.cuts_ & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) &&
        ((vgammaEvent.cuts_ & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) &&
        ((vgammaEvent.cuts_ & SmurfTree::Lep3FullSelection) == SmurfTree::Lep3FullSelection)) &&
	TMath::Abs((vgammaEvent.lep1_+vgammaEvent.lep2_+vgammaEvent.lep3_).M()-91.1976) < 10.0;
    }

    if(cuts == false) continue;
    
    Double_t add = vgammaEvent.sfWeightPU_*vgammaEvent.sfWeightEff_*vgammaEvent.sfWeightTrig_;
    Double_t myWeight = vgammaEvent.scale1fb_*add;    
    
    WGMCEventcount  += myWeight;
    WGMCEventcountE += myWeight*myWeight;
  }

  printf("WG   Event Count 1fb^-1 : %f +/- %f\n",WGMCEventcount,sqrt(WGMCEventcountE));
 
  //*************************************************************************************************
  //Count embedded Normalization
  //*************************************************************************************************
  Double_t WGFOMCEventcount  = 0;
  Double_t WGFOMCEventcountE = 0;

  SmurfTree vgammafoEvent;
  vgammafoEvent.LoadTree(Input2Filename.c_str());
  vgammafoEvent.InitTree(0);
  vgammafoEvent.tree_->SetName("tree");

  for (int n=0;n<vgammafoEvent.tree_->GetEntries();n++) { 
    vgammafoEvent.tree_->GetEntry(n);

    bool cuts = true;
    if     (nsel == 0){
      cuts = cuts && vgammafoEvent.lep1_.pt() > 20 && vgammafoEvent.lep2_.pt() > 10 && vgammafoEvent.lid3_ == 0 &&
      (((vgammafoEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) ||
       ((vgammafoEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection));
    } 
    else if(nsel == 1){
      cuts = cuts && vgammafoEvent.lep1_.pt() > 20 && vgammafoEvent.lep2_.pt() > 10 && vgammafoEvent.lep3_.pt() > 10 && vgammafoEvent.lid3_ != 0 &&
       (((vgammafoEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) ||
        ((vgammafoEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
        ((vgammafoEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection));
    }
    else if(nsel == 2){
      cuts = cuts && vgammafoEvent.lep1_.pt() > 20 && vgammafoEvent.lep2_.pt() > 10 && vgammafoEvent.lep3_.pt() > 10 && vgammafoEvent.lid3_ != 0 &&
       (((vgammafoEvent.cuts_ & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) ||
        ((vgammafoEvent.cuts_ & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
        ((vgammafoEvent.cuts_ & SmurfTree::Lep3FullSelection) != SmurfTree::Lep3FullSelection)) &&
	TMath::Abs((vgammafoEvent.lep1_+vgammafoEvent.lep2_+vgammafoEvent.lep3_).M()-91.1976) < 10.0;
    }

    if(cuts == false) continue;
    
    Double_t add = vgammafoEvent.sfWeightPU_*vgammafoEvent.sfWeightEff_*vgammafoEvent.sfWeightTrig_;
    Double_t myWeight = vgammafoEvent.scale1fb_*add;    

    WGFOMCEventcount  += myWeight;
    WGFOMCEventcountE += myWeight*myWeight;
  }

  printf("WGFO Event Count 1fb^-1 : %f +/- %f\n",WGFOMCEventcount,sqrt(WGFOMCEventcountE));

  printf("Ratio : %f --> %0.5f +/- %0.5f\n",WGFOMCEventcount/WGMCEventcount,WGMCEventcount/WGFOMCEventcount,sqrt(WGMCEventcountE)/WGFOMCEventcount);
}
