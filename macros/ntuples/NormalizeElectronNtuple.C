//root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/NormalizeElectronNtuple.C+\(\"/home/sixie/hist/AllAnalysis/filler/011//AllAnalysis_s09-wwll10-mc3_noskim.root\",\"s09-wwll10-mc3\",\"/home/sixie/ntuples/ElectronOptimization//ElectronOptimization_s09-wwll10-mc3_noskim_normalized.root\",0,\"\"\)

//root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/NormalizeElectronNtuple.C+\(\"/home/sixie/ntuples/ElectronOptimization/ElectronOptimization_bkg_normalized.root\",\"bkg\",\"/home/sixie/ntuples/ElectronOptimization/ElectronOptimization_bkg_normalized_reweighted.root\",1,\"/home/sixie/Tools/root/ElectronIDOptimization/ElectronOptimizationNormalization.root\"\)



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
#include "MitNtupleElectron.C"
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
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("ElectronValidationMod");
    assert(dir);
    

    t = (TTree*)dir->Get("ElectronIDOptimizationTree");
    assert(t);
  } else {
    t = (TTree*)inf->Get("ElectronIDOptimizationTree");  
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
void NormalizeElectronNtuple(const string InputFilename, const string datasetName, 
                             const string OutputFilename, Int_t sampleType, 
                             const string normalizationFile = "") {


  Double_t normalizationWeight = 0;
  Bool_t useReweightFactor = kFALSE;
  TH2F *PtEtaReweightFactor = 0;
  Double_t overallWJetsNormalizationFactor = 0;

  //For Normalizing each sample individually
  if (sampleType == 0 || sampleType >= 10) {
    TTree* electronTree = getTreeFromFile(InputFilename.c_str());
    assert(electronTree);    
    MitNtupleElectron ele(electronTree);
    
    normalizationWeight = getNormalizationWeight(InputFilename, datasetName);
    //*************************************************************************************************
    //Create new normalized tree
    //*************************************************************************************************
    TFile *outputFile = new TFile(OutputFilename.c_str(), "RECREATE");
    outputFile->cd();    
    TTree *normalizedTree = electronTree->CloneTree(0);
    
    for (int n=0;n<electronTree->GetEntries();n++) { 
      ele.GetEntry(n);

      if (sampleType == 10) {
        if (ele.electron_branch_passedSelectionCut == 1) {
          if (datasetName == "s09-wwll10-mc3" || datasetName == "s09-ttbar-mc3") {
            //for mu-e there should be only 1 electron candidate
            ele.electron_branch_weight = normalizationWeight; 
            normalizedTree->Fill(); 
          } else if (datasetName == "s09-we-mc3" || datasetName == "s09-wm-mc3" || datasetName == "s09-wt-mc3") {
            //For the W->enu sample, fill only the fake electron candidates.
            ele.electron_branch_weight = normalizationWeight; 
            if (ele.electron_branch_electronType < 100) {
              normalizedTree->Fill(); 
            }
          } else {
            cout << "Warning: The specified dataset " << datasetName 
                 << " is not recognized to be one of the selection cut samples.\n";
          }
        }        
      } else {
        //For regular samples just fill all electrons
        ele.electron_branch_weight = normalizationWeight; 
        normalizedTree->Fill(); 
      }      
    }
    normalizedTree->Write();
    cout << "Original Tree Entries: " << electronTree->GetEntries() << "  Normalized Tree Entries: " << normalizedTree->GetEntries() << endl;
    outputFile->Close();     
  }
  //For Normalization of Background sample
  else if (sampleType == 1) {
    TTree* electronTree = getTreeFromFile(InputFilename.c_str(), kFALSE);
    assert(electronTree);     
    MitNtupleElectron ele(electronTree);

    //*************************************************************************************************
    //Create new reweighted tree
    //*************************************************************************************************
    TFile *outputFile = new TFile(OutputFilename.c_str(), "RECREATE");
    TTree *reweightedTree = electronTree->CloneTree(0);    

    if (normalizationFile == "") {
      //normalize according to total number of electrons expected in W+Jets/W+gamma bkg
      Double_t totalWeight = 0;
      for (int n=0;n<electronTree->GetEntries();n++) { 
        ele.GetEntry(n);
        totalWeight += ele.electron_branch_weight;
      }

      cout << "Input Bkg Ntuple Total Weight = " << totalWeight << endl;
      normalizationWeight = 1126.5 / totalWeight;
    } else {
      //do pt/eta Reweighting

      TFile *reweightInputFile = new TFile(normalizationFile.c_str(), "READ");
      assert(reweightInputFile);

      TH2F *WJetsFakeElectronPtEta = (TH2F*)reweightInputFile->Get("hWJetsFakeElectronPtEta");
      assert(WJetsFakeElectronPtEta);      
      WJetsFakeElectronPtEta = (TH2F*)(WJetsFakeElectronPtEta->Clone());
      WJetsFakeElectronPtEta->SetDirectory(0);
      reweightInputFile->Close();

      //Create histogram for reweighting factor
      PtEtaReweightFactor = (TH2F*)WJetsFakeElectronPtEta->Clone();
      PtEtaReweightFactor->SetName("PtEtaReweightFactor");
      PtEtaReweightFactor->SetDirectory(0);

      TH2F *BkgFakeElectronPtEta = new TH2F("BkgFakeElectronPtEta", ";Pt [GeV/c];#eta;" , WJetsFakeElectronPtEta->GetXaxis()->GetNbins(), WJetsFakeElectronPtEta->GetXaxis()->GetBinLowEdge(1), WJetsFakeElectronPtEta->GetXaxis()->GetBinUpEdge(WJetsFakeElectronPtEta->GetXaxis()->GetNbins()), WJetsFakeElectronPtEta->GetYaxis()->GetNbins(), WJetsFakeElectronPtEta->GetYaxis()->GetBinLowEdge(1), WJetsFakeElectronPtEta->GetYaxis()->GetBinUpEdge(WJetsFakeElectronPtEta->GetYaxis()->GetNbins()));
      BkgFakeElectronPtEta->SetDirectory(0);
      for (int n=0;n<electronTree->GetEntries();n++) { 
        ele.GetEntry(n);
        BkgFakeElectronPtEta->Fill(ele.electron_branch_pt, ele.electron_branch_eta, ele.electron_branch_weight);        
      }
      PtEtaReweightFactor->Divide(WJetsFakeElectronPtEta, BkgFakeElectronPtEta, 1.0,1.0,"B");

      useReweightFactor = kTRUE;    
    }

    cout << "Reweighting Background Ntuple\n";
    outputFile->cd();    

    //check Reweighting
    TH2F *ReweightedBkgFakeElectronPtEta = new TH2F("BkgFakeElectronPtEta", ";Pt [GeV/c];#eta;" , 200, 0, 200, 200, -3.0, 3.0);    
    ReweightedBkgFakeElectronPtEta->SetDirectory(0);
    for (int n=0;n<electronTree->GetEntries();n++) { 
      if (n%250000 == 0) cout << "Entry " << n << endl;
      ele.GetEntry(n);
      if (useReweightFactor) {
        Double_t reweightFactor = PtEtaReweightFactor->GetBinContent(PtEtaReweightFactor->GetXaxis()->FindFixBin(ele.electron_branch_pt),PtEtaReweightFactor->GetYaxis()->FindFixBin(ele.electron_branch_eta));
        ele.electron_branch_weight = ele.electron_branch_weight*reweightFactor;//*overallWJetsNormalizationFactor;
      } else {
        ele.electron_branch_weight = normalizationWeight;
      }
      reweightedTree->Fill(); 

      ReweightedBkgFakeElectronPtEta->Fill(ele.electron_branch_pt, ele.electron_branch_eta, ele.electron_branch_weight);
    }
    reweightedTree->Write();
    cout << "Original Tree Entries: " << electronTree->GetEntries() << "  Normalized Tree Entries: " << reweightedTree->GetEntries() << endl;
    outputFile->Close();
  
    TCanvas *cv = new TCanvas("BkgReweightedFakeElectronPtEta", "BkgReweightedFakeElectronPtEta", 0,0,800,600);
    ReweightedBkgFakeElectronPtEta->ProjectionX()->DrawCopy("E1");
    cv->SaveAs("BkgReweightedFakeElectronPt.gif");
    ReweightedBkgFakeElectronPtEta->ProjectionY()->DrawCopy("E1");
    cv->SaveAs("BkgReweightedFakeElectronEta.gif");
  } 
  
  //For Normalization of Signal sample
  else if (sampleType == 2) {
    TTree* electronTree = getTreeFromFile(InputFilename.c_str(), kFALSE);
    assert(electronTree);     
    MitNtupleElectron ele(electronTree);

    //*************************************************************************************************
    //Create new reweighted tree
    //*************************************************************************************************
    TFile *outputFile = new TFile(OutputFilename.c_str(), "RECREATE");
    TTree *reweightedTree = electronTree->CloneTree(0);    

    if (normalizationFile == "") {
      //normalize according to total number of electrons expected in WW -> ee nunu signal
      Double_t totalWeight = 0;
      for (int n=0;n<electronTree->GetEntries();n++) { 
        ele.GetEntry(n);
        totalWeight += ele.electron_branch_weight;
      }

      cout << "Input Bkg Ntuple Total Weight = " << totalWeight << endl;
      normalizationWeight = 1126.5 / totalWeight;
    } else {
      //do pt/eta Reweighting

      TFile *reweightInputFile = new TFile(normalizationFile.c_str(), "READ");
      assert(reweightInputFile);
      TH2F *WWSigElectronPtEta = (TH2F*)reweightInputFile->Get("hWWRealElectronPtEta");
      assert(WWSigElectronPtEta);      
      WWSigElectronPtEta = (TH2F*)(WWSigElectronPtEta->Clone());
      WWSigElectronPtEta->SetDirectory(0);
      reweightInputFile->Close();

      //Create histogram for reweighting factor
      PtEtaReweightFactor = (TH2F*)WWSigElectronPtEta->Clone();
      PtEtaReweightFactor->SetName("PtEtaReweightFactor");
      PtEtaReweightFactor->SetDirectory(0);

      TH2F *SigSampleElectronPtEta = new TH2F("SigSampleElectronPtEta", ";Pt [GeV/c];#eta;" , WWSigElectronPtEta->GetXaxis()->GetNbins(), WWSigElectronPtEta->GetXaxis()->GetBinLowEdge(1), WWSigElectronPtEta->GetXaxis()->GetBinUpEdge(WWSigElectronPtEta->GetXaxis()->GetNbins()), WWSigElectronPtEta->GetYaxis()->GetNbins(), WWSigElectronPtEta->GetYaxis()->GetBinLowEdge(1), WWSigElectronPtEta->GetYaxis()->GetBinUpEdge(WWSigElectronPtEta->GetYaxis()->GetNbins()));
      SigSampleElectronPtEta->SetDirectory(0);
      for (int n=0;n<electronTree->GetEntries();n++) { 
        ele.GetEntry(n);
        SigSampleElectronPtEta->Fill(ele.electron_branch_pt, ele.electron_branch_eta, ele.electron_branch_weight);        
      }
      PtEtaReweightFactor->Divide(WWSigElectronPtEta, SigSampleElectronPtEta, 1.0,1.0,"B");

      useReweightFactor = kTRUE;    
    }

    cout << "Reweighting Signal Ntuple\n";
    outputFile->cd();    

    //check Reweighting
    TH2F *ReweightedSigRealElectronPtEta = new TH2F("ReweightedSigRealElectronPtEta", ";Pt [GeV/c];#eta;" , 200, 0, 200, 200, -3.0, 3.0);    
    ReweightedSigRealElectronPtEta->SetDirectory(0);
    for (int n=0;n<electronTree->GetEntries();n++) { 
      if (n%250000 == 0) cout << "Entry " << n << endl;
      ele.GetEntry(n);
      if (useReweightFactor) {
        Double_t reweightFactor = PtEtaReweightFactor->GetBinContent(PtEtaReweightFactor->GetXaxis()->FindFixBin(ele.electron_branch_pt),PtEtaReweightFactor->GetYaxis()->FindFixBin(ele.electron_branch_eta));
        ele.electron_branch_weight = ele.electron_branch_weight*reweightFactor;//*overallWJetsNormalizationFactor;
      } else {
        ele.electron_branch_weight = normalizationWeight;
      }
      reweightedTree->Fill(); 
      ReweightedSigRealElectronPtEta->Fill(ele.electron_branch_pt, ele.electron_branch_eta, ele.electron_branch_weight);
    }
    reweightedTree->Write();
    cout << "Original Tree Entries: " << electronTree->GetEntries() << "  Normalized Tree Entries: " << reweightedTree->GetEntries() << endl;
    outputFile->Close();
  
    TCanvas *cv = new TCanvas("BkgReweightedFakeElectronPtEta", "BkgReweightedFakeElectronPtEta", 0,0,800,600);
    ReweightedSigRealElectronPtEta->ProjectionX()->DrawCopy("E1");
    cv->SaveAs("SigReweightedRealElectronPt.gif");
    ReweightedSigRealElectronPtEta->ProjectionY()->DrawCopy("E1");
    cv->SaveAs("SigReweightedRealElectronEta.gif");
  } 
  else {
    cout << "Warning: Specified sampleType " << sampleType << " is not recognized.\n";
  }



}

