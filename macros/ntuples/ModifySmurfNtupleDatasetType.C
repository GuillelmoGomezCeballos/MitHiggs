//
//root -l -b -q $CMSSW_BASE/src/Smurf/ProcessingAndSkimming/FixProjectedMet.C+\(\"/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP/mitf-alljets/data_2l.root\",\"/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP_Fixed/mitf-alljets/data_2l.root\",-1\) 


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
#include <TROOT.h>
#include <TChain.h>
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"


//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname)
{
  bool verbose = true;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  //TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwMakeNtupleMod");
  //assert(dir);

  TTree* t;
  t = (TTree*) inf->Get("HwwTree0");
  if(!t) t = (TTree*) inf->Get("HwwTree1");
  if(!t) t = (TTree*) inf->Get("HwwTree2");
  if(!t) t = (TTree*) inf->Get("tree");
  assert(t);
  t->SetName("tree");

  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}


//*************************************************************************************************
//Main part of the macro
//*************************************************************************************************
void ModifySmurfNtupleDatasetType(const string InputFilename, const string OutputFilename, 
                                  Int_t type, string TargetDatasetType ) {

  TTree* HwwTree = getTreeFromFile(InputFilename.c_str());
  assert(HwwTree);
  SmurfTree sigEvent;
  sigEvent.LoadTree(InputFilename.c_str(),type);
  sigEvent.InitTree(0);
  sigEvent.tree_->SetName("tree");

  //*************************************************************************************************
  //Create new normalized tree
  //*************************************************************************************************
  cout << "Output File : " << OutputFilename << endl;
  TFile *outputFile = new TFile(OutputFilename.c_str(), "RECREATE");
  outputFile->cd();
  TTree *normalizedTree = sigEvent.tree_->CloneTree(0);
  
  cout << "Events in the ntuple: " << sigEvent.tree_->GetEntries() << endl;
  
  for (int n=0;n<sigEvent.tree_->GetEntries();n++) { 
    sigEvent.tree_->GetEntry(n);
    if (n%100000==0) cout << "Processed Event " << n << endl;    

    if (TargetDatasetType == "qqww_MCAtNLO_ScaleUp") {
      sigEvent.dstype_ = SmurfTree::ggww;
    } 
    if (TargetDatasetType == "qqww_MCAtNLO_ScaleDown") {
      sigEvent.dstype_ = SmurfTree::ggzz;
    } 

    normalizedTree->Fill(); 
  }
  cout << "Events in output ntuple: " << normalizedTree->GetEntries() << endl;
  normalizedTree->Write();
  outputFile->Close();
}
