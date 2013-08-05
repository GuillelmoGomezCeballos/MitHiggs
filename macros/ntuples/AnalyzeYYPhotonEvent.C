//root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/AnalyzeYYPhotonEvent.C+\(\"/home/sixie/hist/YYTree/filler/009/YYTree_s8-zmm-id11_noskim.root\"\)



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
#include "YYPhotonEvent.C"
#include <vector>
#include <algorithm>
#include "MitAna/Utils/interface/SimpleTable.h"


//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, Bool_t loadFromDir = kTRUE)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = 0;
  if (loadFromDir) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("PhotonValidationMod");
    assert(dir);
    

    t = (TTree*)dir->Get("YangYongPhotonTree");
    assert(t);
  } else {
    t = (TTree*)inf->Get("YangYongPhotonTree");  
  }

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

  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  if (!dir) {
    cout << "Could not find directory AnaFwkMod"
         << " in file " << filename << endl;
    delete file;
    return 0;
  }

  TH1F *hist = (TH1F*)dir->Get("hDAllEvents");
  if (!hist) {
    cout << "Could not find histogram hDAllEvents in directory AnaFwkMod"
         << " in file " << filename << endl;
    delete dir;
    delete file;
    return 0;
  }  
  Double_t NEvents = hist->Integral();

  //Get CrossSection
  mithep::SimpleTable xstab("$CMSSW_BASE/src/MitPhysics/data/xs.dat");
  Double_t CrossSection = xstab.Get(datasetName.c_str());
  Double_t Weight = CrossSection / NEvents;

  delete dir;
  delete file;
  return Weight;

}


//*************************************************************************************************
//Main part of the macro
//*************************************************************************************************
void AnalyzeYYPhotonEvent(const string InputFilename) {

  
  TTree* YYTree = getTreeFromFile(InputFilename.c_str());
  assert(YYTree);    
  YYPhotonEvent photon(YYTree);
      
  TH1D *plot = new TH1D("hPlot", ";XAxis; Number of Events ", 6, -3,3);

  for (int n=0;n<YYTree->GetEntries();n++) { 
    photon.GetEntry(n);

    Double_t weight = photon.photon_branch_weight;    
    //cout << weight << " " << photon.photon_branch_NPhotons << endl;
    if (photon.photon_branch_NPhotons > 0 ) {
      plot->Fill(photon.photon_branch_Photon1Eta);
    }

  }

  TCanvas *c = new TCanvas("canvas" , "canvas" , 800, 600);
  plot->DrawCopy();
  
}

