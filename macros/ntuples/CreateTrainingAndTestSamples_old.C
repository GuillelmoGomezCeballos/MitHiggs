//
//root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/CreateTrainingAndTestSamples.C+\(\"/home/sixie/CMSSW_3_1_2/src/test.root\",0.5\) 


#include <iostream>
//#include <stdlib.h>
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
void CreateTrainingAndTestSamples(const string InputFilename, Double_t TestSampleFraction) {

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


  if ( TestSampleFraction > 1.0 || TestSampleFraction < 0.0 ) {
    cerr << "TestSampleFraction = " << TestSampleFraction << " is not in the range [0,1]. " 
         << endl;
    assert( TestSampleFraction >= 1.0 && TestSampleFraction <= 0.0 );
  }
  cout << "Input Tree Size : " << HwwTree->GetEntries() << endl;

  //Remove the '.root' from end of filename
  size_t p;
  p = InputFilename.find(".root");  
  string tmpInputFilename = InputFilename.substr(0,p);
  
  //change 'unrandomized' to 'randomized' in the file name.
  p = tmpInputFilename.find("_unrandomized");
  if (p != string::npos) {
    tmpInputFilename = tmpInputFilename.substr(0,p) + "_randomized" ;
  }

  //*************************************************************************************************
  //Create Randomized Sample Tree
  //*************************************************************************************************
  TFile *allSampleFile = new TFile((tmpInputFilename+".all.root").c_str(), "RECREATE");
  allSampleFile->cd();
  TTree *allSampleTree = HwwTree->CloneTree(0);

  for (Int_t n=0;n<HwwTree->GetEntries();n++) { 
    event.GetEntry(indices[n]);
    allSampleTree->Fill(); 
  }  
  allSampleTree->Write();
  cout << "All Tree Size: " << allSampleTree->GetEntries() << endl;
  allSampleFile->Close();

 

  //*************************************************************************************************
  //Create Test Sample Tree
  //*************************************************************************************************
  //For some reason I need to make another TTree, otherwise when I try to clone it, root crashes.
  TTree* HwwTree2 = getTreeFromFile(InputFilename.c_str());
  assert(HwwTree2);
  assert(HwwTree2->GetEntries() == 
         HwwTree->GetEntries());
  MitNtupleEvent event2(HwwTree2);

  TFile *testSampleFile = new TFile((tmpInputFilename+".test.root").c_str(), "RECREATE");
  testSampleFile->cd();

  Int_t testSampleSize = Int_t(TestSampleFraction * double(HwwTree->GetEntries()));
  TTree *testSampleTree = HwwTree2->CloneTree(0);

  for (Int_t n=0;n<testSampleSize;n++) { 
    event2.GetEntry(indices[n]);
    event2.H_weight = event2.H_weight / (TestSampleFraction);
    testSampleTree->Fill(); 
  }  
  testSampleTree->Write();
  cout << "Test Tree Size: " << testSampleTree->GetEntries() << endl;
  testSampleFile->Close();



  //*************************************************************************************************
  //Create Training Sample Tree
  //*************************************************************************************************
  //For some reason I need to make another TTree, otherwise when I try to clone it, root crashes.
  TTree* HwwTree3 = getTreeFromFile(InputFilename.c_str());
  assert(HwwTree3);
  assert(HwwTree3->GetEntries() ==
         HwwTree->GetEntries());
  MitNtupleEvent event3(HwwTree3);

  TFile *trainingSampleFile = new TFile((tmpInputFilename+".training.root").c_str(), "RECREATE");
  trainingSampleFile->cd();
  TTree *trainingSampleTree = HwwTree3->CloneTree(0);

  for (Int_t n=testSampleSize;n<HwwTree->GetEntries();n++) { 
    event3.GetEntry(indices[n]);
    event3.H_weight = event3.H_weight / (1 - TestSampleFraction);
    trainingSampleTree->Fill(); 
  }
  trainingSampleTree->Write();
  cout << "Training Tree Size: " << trainingSampleTree->GetEntries() << endl;
  trainingSampleFile->Close();

}

