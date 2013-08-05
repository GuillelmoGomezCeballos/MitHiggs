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
#include "/home/ceballos/releases/CMSSW_5_2_3_patch3/src/Smurf/Core/SmurfTree.h"
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include <TChain.h>


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

  TTree* t = (TTree*)inf->Get("tree");
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

  SmurfTree sigEvent1;
  sigEvent1.LoadTree(InputFilename.c_str(),-1);
  sigEvent1.InitTree(0);
    
  //*************************************************************************************************
  //Create randomized list of indices
  //*************************************************************************************************
  vector<Int_t> indices;
  for (Int_t i=0; i<sigEvent1.tree_->GetEntries(); ++i) {
    indices.push_back(i);
  }
  random_shuffle(indices.begin(),indices.end());


  if ( TestSampleFraction > 1.0 || TestSampleFraction < 0.0 ) {
    cerr << "TestSampleFraction = " << TestSampleFraction << " is not in the range [0,1]. " 
         << endl;
    assert( TestSampleFraction >= 1.0 && TestSampleFraction <= 0.0 );
  }
  cout << "Input Tree Size : " << sigEvent1.tree_->GetEntries() << endl;

  //Remove the '.root' from end of filename
  size_t p;
  p = InputFilename.find(".root");
  string tmpInputFilename = InputFilename.substr(0,p);
  
  //change 'unrandomized' to 'randomized' in the file name.
  p = tmpInputFilename.find("_unrandomized");
  if (p != string::npos) {
    tmpInputFilename = tmpInputFilename.substr(0,p) + "_randomized" ;
  }
  
  cout << "Finishing here for now!" << endl;
  return;
  //*************************************************************************************************
  //Create Randomized Sample Tree
  //*************************************************************************************************
  TFile *allSampleFile = new TFile((tmpInputFilename+".all.root").c_str(), "RECREATE");
  allSampleFile->cd();
  TTree *allSampleTree = sigEvent1.tree_->CloneTree(0);

  for (Int_t n=0;n<sigEvent1.tree_->GetEntries();n++) { 
    if(n%10000 == 0) cout << n << endl;
    sigEvent1.tree_->GetEntry(indices[n]);
    allSampleTree->Fill(); 
  }  
  allSampleTree->Write();
  cout << "All Tree Size: " << allSampleTree->GetEntries() << endl;
  allSampleFile->Close();

  if (TestSampleFraction > 0.0){ // Only if we want to define a test sample
    //*************************************************************************************************
    //Create Test Sample Tree
    //*************************************************************************************************
    //For some reason I need to make another TTree, otherwise when I try to clone it, root crashes.
    SmurfTree sigEvent2;
    sigEvent2.LoadTree(InputFilename.c_str(),-1);
    sigEvent2.InitTree(0);

    assert(sigEvent1.tree_->GetEntries() == 
           sigEvent2.tree_->GetEntries());

    TFile *testSampleFile = new TFile((tmpInputFilename+".test.root").c_str(), "RECREATE");
    testSampleFile->cd();

    Int_t testSampleSize = Int_t(TestSampleFraction * double(sigEvent1.tree_->GetEntries()));
    TTree *testSampleTree = sigEvent2.tree_->CloneTree(0);

    for (Int_t n=0;n<testSampleSize;n++) { 
      if(n%10000 == 0) cout << n << endl;
      sigEvent2.tree_->GetEntry(indices[n]);
      sigEvent2.scale1fb_ = sigEvent2.scale1fb_ / (TestSampleFraction);
      testSampleTree->Fill(); 
    }  
    testSampleTree->Write();
    cout << "Test Tree Size: " << testSampleTree->GetEntries() << endl;
    testSampleFile->Close();

    //*************************************************************************************************
    //Create Training Sample Tree
    //*************************************************************************************************
    //For some reason I need to make another TTree, otherwise when I try to clone it, root crashes.
    SmurfTree sigEvent3;
    sigEvent3.LoadTree(InputFilename.c_str(),-1);
    sigEvent3.InitTree(0);

    assert(sigEvent1.tree_->GetEntries() == 
           sigEvent3.tree_->GetEntries());

    TFile *trainingSampleFile = new TFile((tmpInputFilename+".training.root").c_str(), "RECREATE");
    trainingSampleFile->cd();
    TTree *trainingSampleTree = sigEvent3.tree_->CloneTree(0);

    for (Int_t n=testSampleSize;n<sigEvent1.tree_->GetEntries();n++) { 
      if(n%10000 == 0) cout << n << endl;
      sigEvent3.tree_->GetEntry(indices[n]);
      sigEvent3.scale1fb_ = sigEvent3.scale1fb_ / (1 - TestSampleFraction);
      trainingSampleTree->Fill(); 
    }
    trainingSampleTree->Write();
    cout << "Training Tree Size: " << trainingSampleTree->GetEntries() << endl;
    trainingSampleFile->Close();
  }
}
