//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep  7 14:10:46 2009 by ROOT version 5.22/00a
// from TTree YangYongPhotonTree/YangYongPhotonTree
// found on file: /home/sixie/hist/YYTree/filler/009/YYTree_s8-zmm-id11_noskim.root
//////////////////////////////////////////////////////////

#ifndef YYPhotonEvent_h
#define YYPhotonEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class YYPhotonEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         photon_branch_weight;
   Float_t         photon_branch_sampleType;
   Float_t         photon_branch_NPhotons;
   Float_t         photon_branch_Photon1Flag;
   Float_t         photon_branch_Photon1Pt;
   Float_t         photon_branch_Photon1Eta;
   Float_t         photon_branch_Photon1Phi;
   Float_t         photon_branch_Photon1EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon1EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon1HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon1HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon1HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon1HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon1HcalRecHitIso;
   Float_t         photon_branch_Photon1HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon1HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon1HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon1HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon1SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon1SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon1HadOverEm;
   Float_t         photon_branch_Photon2Flag;
   Float_t         photon_branch_Photon2Pt;
   Float_t         photon_branch_Photon2Eta;
   Float_t         photon_branch_Photon2Phi;
   Float_t         photon_branch_Photon2EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon2EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon2HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon2HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon2HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon2HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon2HcalRecHitIso;
   Float_t         photon_branch_Photon2HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon2HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon2HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon2HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon2SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon2SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon2HadOverEm;
   Float_t         photon_branch_Photon3Flag;
   Float_t         photon_branch_Photon3Pt;
   Float_t         photon_branch_Photon3Eta;
   Float_t         photon_branch_Photon3Phi;
   Float_t         photon_branch_Photon3EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon3EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon3HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon3HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon3HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon3HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon3HcalRecHitIso;
   Float_t         photon_branch_Photon3HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon3HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon3HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon3HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon3SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon3SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon3HadOverEm;
   Float_t         photon_branch_Photon4Flag;
   Float_t         photon_branch_Photon4Pt;
   Float_t         photon_branch_Photon4Eta;
   Float_t         photon_branch_Photon4Phi;
   Float_t         photon_branch_Photon4EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon4EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon4HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon4HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon4HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon4HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon4HcalRecHitIso;
   Float_t         photon_branch_Photon4HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon4HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon4HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon4HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon4SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon4SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon4HadOverEm;
   Float_t         photon_branch_Photon5Flag;
   Float_t         photon_branch_Photon5Pt;
   Float_t         photon_branch_Photon5Eta;
   Float_t         photon_branch_Photon5Phi;
   Float_t         photon_branch_Photon5EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon5EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon5HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon5HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon5HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon5HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon5HcalRecHitIso;
   Float_t         photon_branch_Photon5HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon5HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon5HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon5HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon5SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon5SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon5HadOverEm;
   Float_t         photon_branch_Photon6Flag;
   Float_t         photon_branch_Photon6Pt;
   Float_t         photon_branch_Photon6Eta;
   Float_t         photon_branch_Photon6Phi;
   Float_t         photon_branch_Photon6EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon6EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon6HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon6HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon6HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon6HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon6HcalRecHitIso;
   Float_t         photon_branch_Photon6HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon6HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon6HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon6HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon6SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon6SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon6HadOverEm;
   Float_t         photon_branch_Photon7Flag;
   Float_t         photon_branch_Photon7Pt;
   Float_t         photon_branch_Photon7Eta;
   Float_t         photon_branch_Photon7Phi;
   Float_t         photon_branch_Photon7EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon7EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon7HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon7HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon7HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon7HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon7HcalRecHitIso;
   Float_t         photon_branch_Photon7HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon7HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon7HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon7HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon7SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon7SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon7HadOverEm;
   Float_t         photon_branch_Photon8Flag;
   Float_t         photon_branch_Photon8Pt;
   Float_t         photon_branch_Photon8Eta;
   Float_t         photon_branch_Photon8Phi;
   Float_t         photon_branch_Photon8EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon8EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon8HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon8HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon8HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon8HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon8HcalRecHitIso;
   Float_t         photon_branch_Photon8HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon8HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon8HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon8HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon8SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon8SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon8HadOverEm;
   Float_t         photon_branch_Photon9Flag;
   Float_t         photon_branch_Photon9Pt;
   Float_t         photon_branch_Photon9Eta;
   Float_t         photon_branch_Photon9Phi;
   Float_t         photon_branch_Photon9EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon9EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon9HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon9HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon9HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon9HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon9HcalRecHitIso;
   Float_t         photon_branch_Photon9HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon9HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon9HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon9HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon9SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon9SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon9HadOverEm;
   Float_t         photon_branch_Photon10Flag;
   Float_t         photon_branch_Photon10Pt;
   Float_t         photon_branch_Photon10Eta;
   Float_t         photon_branch_Photon10Phi;
   Float_t         photon_branch_Photon10EcalRecHitIsoDr03;
   Float_t         photon_branch_Photon10EcalRecHitIsoDr04;
   Float_t         photon_branch_Photon10HcalDepth1TowerSumEtDr03;
   Float_t         photon_branch_Photon10HcalDepth1TowerSumEtDr04;
   Float_t         photon_branch_Photon10HcalDepth2TowerSumEtDr03;
   Float_t         photon_branch_Photon10HcalDepth2TowerSumEtDr04;
   Float_t         photon_branch_Photon10HcalRecHitIso;
   Float_t         photon_branch_Photon10HcalTowerSumEtDr03;
   Float_t         photon_branch_Photon10HcalTowerSumEtDr04;
   Float_t         photon_branch_Photon10HollowConeTrkIsoDr03;
   Float_t         photon_branch_Photon10HollowConeTrkIsoDr04;
   Float_t         photon_branch_Photon10SolidConeTrkIsoDr03;
   Float_t         photon_branch_Photon10SolidConeTrkIsoDr04;
   Float_t         photon_branch_Photon10HadOverEm;
   Float_t         photon_branch_NMuons;
   Float_t         photon_branch_Muon1Pt;
   Float_t         photon_branch_Muon1Eta;
   Float_t         photon_branch_Muon1Phi;
   Float_t         photon_branch_Muon2Pt;
   Float_t         photon_branch_Muon2Eta;
   Float_t         photon_branch_Muon2Phi;

   // List of branches
   TBranch        *b_photon_branch;   //!

   YYPhotonEvent(TTree *tree=0);
   virtual ~YYPhotonEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef YYPhotonEvent_cxx
YYPhotonEvent::YYPhotonEvent(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/sixie/hist/YYTree/filler/009/YYTree_s8-zmm-id11_noskim.root");
      if (!f) {
         f = new TFile("/home/sixie/hist/YYTree/filler/009/YYTree_s8-zmm-id11_noskim.root");
      }
      tree = (TTree*)gDirectory->Get("YangYongPhotonTree");

   }
   Init(tree);
}

YYPhotonEvent::~YYPhotonEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t YYPhotonEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t YYPhotonEvent::LoadTree(Long64_t entry)
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

void YYPhotonEvent::Init(TTree *tree)
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

   fChain->SetBranchAddress("photon_branch", &photon_branch_weight, &b_photon_branch);
   Notify();
}

Bool_t YYPhotonEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void YYPhotonEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t YYPhotonEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef YYPhotonEvent_cxx
