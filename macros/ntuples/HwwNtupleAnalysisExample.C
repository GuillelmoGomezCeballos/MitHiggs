//
//root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/HwwNtupleAnalysisExample.C+\(\"/home/sixie/ntuples/HwwNtuple/filler/009/HwwNtuple_all_noskim.root\"\) 


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
void HwwNtupleAnalysisExample(const string InputFilename) {

  TTree* HwwTree = getTreeFromFile(InputFilename.c_str());
  assert(HwwTree);

  MitNtupleEvent event(HwwTree);

  //*************************************************************************************************
  //Fill new randomized tree
  //*************************************************************************************************
  for (int n=0;n<HwwTree->GetEntries();n++) { 
    event.GetEntry(n);

    float eventweight = event.H_weight;
    int finalstatetype = (int)event.H_ltype;
    int decay = (int)event.H_decay;

    cout << eventweight << " " << finalstatetype << " " << decay << endl;
  }
}

