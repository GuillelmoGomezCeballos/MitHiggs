//
//root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/RandomizeHwwNtuple.C+\(\"/home/sixie/CMSSW_3_1_2/src/HwwNtuple_s8-h160ww2l-id11_noskim_0000.root\"\) 


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
#include "MitNtupleEvent.C"
#include <vector>
#include <algorithm>


//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get("HwwTree");
  assert(t);

  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}

//*************************************************************************************************
//Main part of the macro
//*************************************************************************************************
void RandomizeHwwNtuple(const string InputFilename, const string OutputFilename) {

  TTree* HwwTree = getTreeFromFile(InputFilename.c_str());
  assert(HwwTree);

  MitNtupleEvent event(HwwTree);

  //*************************************************************************************************
  //Create randomized list of indices
  //*************************************************************************************************
  vector<Int_t> indices;
  for (Int_t i=0; i<HwwTree->GetEntries(); ++i) {
    indices.push_back(i);
  }
  random_shuffle(indices.begin(),indices.end());


  //*************************************************************************************************
  //Create new randomized tree
  //*************************************************************************************************
  TFile *outputFile = new TFile(OutputFilename.c_str(), "RECREATE");
  outputFile->cd();

  TTree *randomizedTree = HwwTree->CloneTree(0);
  for (int n=0;n<HwwTree->GetEntries();n++) { 
    event.GetEntry(indices[n]);
    randomizedTree->Fill();
  }
  randomizedTree->Write();
  cout << "Randomized Tree Size: " << randomizedTree->GetEntries() << endl;
  outputFile->Close();
}

