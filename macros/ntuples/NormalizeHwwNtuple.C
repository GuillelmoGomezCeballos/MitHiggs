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
#include "/home/ceballos/releases/CMSSW_5_3_14/src/Smurf/Core/SmurfTree.h"
#include <vector>
#include <algorithm>
#include "MitAna/Utils/interface/SimpleTable.h"
#include <TROOT.h>
#include <TChain.h>


//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, int nsel)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  //TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwMakeNtupleMod");
  //assert(dir);

  TTree* t;
  
  if     (nsel==0) t = (TTree*) gROOT->FindObject("HwwTree0");
  else if(nsel==1) t = (TTree*) gROOT->FindObject("HwwTree1");
  else if(nsel==2) t = (TTree*) gROOT->FindObject("HwwTree2");
  else if(nsel==3) t = (TTree*) gROOT->FindObject("HwwTree3");
  else if(nsel==4) t = (TTree*) gROOT->FindObject("HwwTree4");
  else             assert(1);
  t->SetName("tree");
  assert(t);

  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}

//--------------------------------------------------------------------------------------------------
// Get Total Number of Events in the sample
//--------------------------------------------------------------------------------------------------
Double_t getNormalizationWeight(TString filename, TString datasetName) {
  // Get Normalization Weight Factor

  //Get Number of Events in the Sample
  TFile *file = new TFile(filename.Data(),"READ");
  if (!file) {
    cout << "Could not open file " << filename << endl;
    return 0;
  }

  //TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  //if (!dir) {
  //  cout << "Could not find directory AnaFwkMod"
  //       << " in file " << filename << endl;
  //  delete file;
  //  return 0;
  //}

  TH1D *hist = (TH1D*) gROOT->FindObject("hDAllEvents");
  if (!hist) {
    cout << "Could not find histogram hDAllEvents in directory AnaFwkMod"
         << " in file " << filename << endl;
    //delete dir;
    delete file;
    return 0;
  }
  Double_t NEvents = hist->Integral();
  cout << "Original events in the sample: " << NEvents << endl;

  //Get CrossSection
  mithep::SimpleTable xstab("$CMSSW_BASE/src/MitPhysics/data/xs.dat");
  Double_t CrossSection = xstab.Get(datasetName.Data());
  Double_t Weight = CrossSection / NEvents;
  // weight for data is always 1 (-1 to make a trick for fakes)
  if(CrossSection < 0) Weight = -1.0;

  // Only gg->H and qqH samples
  if(datasetName.Contains("-gf-") ||
     datasetName.Contains("-vbf-")){
    // Only high mass samples
    if(datasetName.Contains("h250") || datasetName.Contains("h300") || datasetName.Contains("h350") || datasetName.Contains("h400") || 
       datasetName.Contains("h450") || datasetName.Contains("h500") || datasetName.Contains("h550") || datasetName.Contains("h600") ||
       datasetName.Contains("h700") || datasetName.Contains("h800") || datasetName.Contains("h900") || datasetName.Contains("h1000")){
      TH1D *histBeforeWeight = (TH1D*) gROOT->FindObject("hDHKFactor_1");
      TH1D *histAfterWeight  = (TH1D*) gROOT->FindObject("hDHKFactor_2");
      cout << "Nevents before/after reweighting: " << histBeforeWeight->GetSumOfWeights() << " " << 
                                                      histAfterWeight->GetSumOfWeights() << endl;
      Weight = Weight * histBeforeWeight->GetSumOfWeights() / histAfterWeight->GetSumOfWeights();
    }
  }

  //delete dir;
  delete file;
  return Weight;

}

//*************************************************************************************************
//Main part of the macro
//*************************************************************************************************
void NormalizeHwwNtuple(const string InputFilename, const string datasetName, 
                        const string OutputFilename, const int nsel) {

  //TTree* HwwTree = getTreeFromFile(InputFilename.c_str(),nsel);
  //assert(HwwTree);

  SmurfTree sigEvent;
  sigEvent.LoadTree(InputFilename.c_str(),nsel);
  sigEvent.InitTree(0);
  sigEvent.tree_->SetName("tree");

  Double_t normalizationWeight = getNormalizationWeight(InputFilename, datasetName);

  //*************************************************************************************************
  //Create new normalized tree
  //*************************************************************************************************
  TFile *outputFile = new TFile(OutputFilename.c_str(), "RECREATE");
  outputFile->cd();

  TTree *normalizedTree = sigEvent.tree_->CloneTree(0);
  
  cout << "Events in the ntuple: " << sigEvent.tree_->GetEntries() << endl;
  
  for (int n=0;n<sigEvent.tree_->GetEntries();n++) { 
    sigEvent.tree_->GetEntry(n);
    if     (normalizationWeight < 0){
      sigEvent.scale1fb_ = sigEvent.scale1fb_ * 1.0;
    }
    else if(nsel == 0){
      sigEvent.scale1fb_ = sigEvent.scale1fb_ * normalizationWeight *  1.0;
    }
    else if(nsel == 1){
      sigEvent.scale1fb_ = sigEvent.scale1fb_ * normalizationWeight *  1.0;
    }
    else {
      sigEvent.scale1fb_ = sigEvent.scale1fb_ * normalizationWeight *  1.0;
    }
    normalizedTree->Fill(); 
  }

  normalizedTree->Write();
  outputFile->Close();
}
