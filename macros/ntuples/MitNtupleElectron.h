//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 14 11:34:27 2009 by ROOT version 5.22/00a
// from TTree ElectronIDOptimizationTree/ElectronIDOptimizationTree
// found on file: /home/sixie/hist/AllAnalysis/filler/011/AllAnalysis_s09-zee-mc3_noskim.root
//////////////////////////////////////////////////////////

#ifndef MitNtupleElectron_h
#define MitNtupleElectron_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class MitNtupleElectron {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         electron_branch_weight;
   Float_t         electron_branch_sampleType;
   Float_t         electron_branch_pt;
   Float_t         electron_branch_eta;
   Float_t         electron_branch_phi;
   Float_t         electron_branch_electronType;
   Float_t         electron_branch_CaloIsolation;
   Float_t         electron_branch_CovEtaEta;
   Float_t         electron_branch_CoviEtaiEta;
   Float_t         electron_branch_DeltaEtaSuperClusterTrackAtVtx;
   Float_t         electron_branch_DeltaEtaSeedClusterTrackAtCalo;
   Float_t         electron_branch_DeltaPhiSuperClusterTrackAtVtx;
   Float_t         electron_branch_DeltaPhiSeedClusterTrackAtCalo;
   Float_t         electron_branch_ESuperClusterOverP;
   Float_t         electron_branch_ESeedClusterOverPout;
   Float_t         electron_branch_ESeedClusterOverPIn;
   Float_t         electron_branch_FBrem;
   Float_t         electron_branch_HadronicOverEm;
   Float_t         electron_branch_EcalRecHitIsoDr04;
   Float_t         electron_branch_TrackIsolationDr04;
   Float_t         electron_branch_HcalTowerSumEtDr04;
   Float_t         electron_branch_HcalDepth1TowerSumEtDr04;
   Float_t         electron_branch_HcalDepth2TowerSumEtDr04;
   Float_t         electron_branch_E15;
   Float_t         electron_branch_E25Max;
   Float_t         electron_branch_E55;
   Float_t         electron_branch_FracSharedHits;
   Float_t         electron_branch_IDLikelihood;
   Float_t         electron_branch_Mva;
   Float_t         electron_branch_PIn;
   Float_t         electron_branch_POut;
   Float_t         electron_branch_NumberOfClusters;
   Float_t         electron_branch_Classification;
   Float_t         electron_branch_d0;
   Float_t         electron_branch_isConversion;
   Float_t         electron_branch_passedChargeFilter;
   Float_t         electron_branch_IsBarrel;
   Float_t         electron_branch_passedSelectionCut;
   Float_t         electron_branch_passedZeeSignalCut;
   Float_t         electron_branch_passedWJetsBkgCut;

   // List of branches
   TBranch        *b_electron_branch;   //!

   MitNtupleElectron(TTree *tree=0);
   virtual ~MitNtupleElectron();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MitNtupleElectron_cxx
MitNtupleElectron::MitNtupleElectron(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/sixie/hist/AllAnalysis/filler/011/AllAnalysis_s09-zee-mc3_noskim.root");
      if (!f) {
         f = new TFile("/home/sixie/hist/AllAnalysis/filler/011/AllAnalysis_s09-zee-mc3_noskim.root");
      }
      tree = (TTree*)gDirectory->Get("ElectronIDOptimizationTree");

   }
   Init(tree);
}

MitNtupleElectron::~MitNtupleElectron()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MitNtupleElectron::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MitNtupleElectron::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MitNtupleElectron::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("electron_branch", &electron_branch_weight, &b_electron_branch);
   Notify();
}

Bool_t MitNtupleElectron::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MitNtupleElectron::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MitNtupleElectron::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MitNtupleElectron_cxx
