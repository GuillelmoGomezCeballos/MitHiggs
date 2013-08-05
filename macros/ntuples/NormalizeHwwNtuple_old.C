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
#include "MitNtupleEvent.C"
#include <vector>
#include <algorithm>
#include "MitAna/Utils/interface/SimpleTable.h"


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
  else            assert(1);
  t->SetName("HwwTree");
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
Double_t getNormalizationWeight(string filename, string datasetName) {
  // Get Normalization Weight Factor

  //Get Number of Events in the Sample
  TFile *file = new TFile(filename.c_str(),"READ");
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
  Double_t CrossSection = xstab.Get(datasetName.c_str());
  Double_t Weight = CrossSection / NEvents;
  // weight for data is always 1
  if(CrossSection < 0) Weight = 1.0;

  //delete dir;
  delete file;
  return Weight;

}


//*************************************************************************************************
//Main part of the macro
//*************************************************************************************************
void NormalizeHwwNtuple(const string InputFilename, const string datasetName, 
                        const string OutputFilename, const int nsel) {

  TTree* HwwTree = getTreeFromFile(InputFilename.c_str(),nsel);
  assert(HwwTree);

  MitNtupleEvent event(HwwTree);

  Double_t normalizationWeight = getNormalizationWeight(InputFilename, datasetName);

  //*************************************************************************************************
  //Create new normalized tree
  //*************************************************************************************************
  TFile *outputFile = new TFile(OutputFilename.c_str(), "RECREATE");
  outputFile->cd();

  TTree *normalizedTree = HwwTree->CloneTree(0);
  
  cout << "Events in the ntuple: " << HwwTree->GetEntries() << endl;
  
  for (int n=0;n<HwwTree->GetEntries();n++) { 
    event.GetEntry(n);
    event.H_weight = event.H_weight * normalizationWeight;
    normalizedTree->Fill(); 
  }

  normalizedTree->Write();
  outputFile->Close();
}

