// $Id: runHwwwExample.C,v 1.17 2009/04/29 15:10:30 loizides Exp $

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Utils/interface/SimpleTable.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PlotKineMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/SelMods/interface/DilepSelMod.h"
#include "MitPhysics/Mods/interface/EffMod.h"
#include "MitHiggs/HwwwMods/interface/H2lAnaMod.h"
#include "MitHiggs/HwwwMods/interface/H2lNtupMod.h"
#include "MitHiggs/HwwMods/interface/HwwAnaMod.h"
#include "MitHiggs/HwwwMods/interface/HwwwAnaMod.h"
#include "MitHiggs/HwwwMods/interface/HwwwNtupMod.h"
#include "MitHiggs/HwwwMods/interface/HwwwSignalMCMod.h"
#endif

//--------------------------------------------------------------------------------------------------
void runHwwwExample(Bool_t      mcomp = 1,
                    const char *files = 0,
                    UInt_t nev        = 0)
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;
  gErrorIgnoreLevel = kInfo;

  const Double_t ptMin    = 10;
  const Double_t ptMaxMin = 20;
  const Double_t minEt    = 30;
  const Double_t minDi    = 12;
  const Double_t minJMet  = 15; 
  const char *muInput  = Names::gkMuonBrn;
  const char *elInput  = Names::gkElectronBrn;
  const char *jetInput = Names::gkCaloJetBrn;
  const char *metInput = Names::gkCaloMetBrn;

  // setup analysis object
  Analysis *ana = new Analysis;
  ana->SetUseHLT(0);
  if (nev>0)
    ana->SetProcessNEvents(nev);
  TString ofname(gSystem->Getenv("MIT_OUTPUTFILE"));
  if (ofname.IsNull()) 
    ana->SetOutputName("hwww-hists.root");
  else 
    ana->SetOutputName(ofname);

  if (files)
    ana->AddFile(files);

  // setup modules
  GeneratorMod *genMod = new GeneratorMod;
  genMod->SetFillHist(1);
  ana->AddSuperModule(genMod);

  if (mcomp) {
    HwwwSignalMCMod *hwwwMcMod = new HwwwSignalMCMod;
    hwwwMcMod->SetFillHist(1);
    genMod->Add(hwwwMcMod);

    HwwwSignalMCMod *hwwwSSMcMod = new HwwwSignalMCMod("SameSignMCMod");
    hwwwSSMcMod->SetFillHist(1);
    hwwwSSMcMod->SetDoSameSign(1);
    genMod->Add(hwwwSSMcMod);
  }

  MuonIDMod *muCl = new MuonIDMod;  
  muCl->SetInputName(muInput);
  muCl->SetPtMin(ptMin);
  ana->AddSuperModule(muCl);

  ElectronIDMod *elId = new ElectronIDMod;
  elId->SetInputName(elInput);
  elId->SetPtMin(ptMin);
  ana->AddSuperModule(elId);

  ElectronCleaningMod *elCl = new ElectronCleaningMod;
  elCl->SetGoodElectronsName(elId->GetOutputName());
  elCl->SetCleanMuonsName(muCl->GetOutputName());
  elId->Add(elCl);

  PublisherMod<CaloJet,Jet> *pubJet = new PublisherMod<CaloJet,Jet>("JetPub");
  pubJet->SetInputName(jetInput);
  pubJet->SetOutputName(Form("Pub%s",jetInput));
  ana->AddSuperModule(pubJet);

  JetIDMod *jetId = new JetIDMod;            
  jetId->SetInputName(pubJet->GetOutputName());
  pubJet->Add(jetId);

  JetCleaningMod *jetCl = new JetCleaningMod;
  jetCl->SetGoodJetsName(jetId->GetOutputName());
  jetCl->SetCleanPhotonsName("");
  jetCl->SetCleanElectronsName(elCl->GetOutputName());
  jetId->Add(jetCl);

  PublisherMod<CaloMet,Met> *pubMet = new PublisherMod<CaloMet,Met>("MetPub");
  pubMet->SetInputName(metInput);
  pubMet->SetOutputName(ModNames::gkCleanCaloMetName);
  ana->AddSuperModule(pubMet);

  MergeLeptonsMod *merger = new MergeLeptonsMod;
  merger->SetMuonsName(muCl->GetOutputName());
  merger->SetElectronsName(elCl->GetOutputName());
  ana->AddSuperModule(merger);

  EffMod *effMod = new EffMod;
  effMod->SetCol1Name(ModNames::gkMCLeptonsName);
  effMod->SetCol2Name(merger->GetOutputName());
  effMod->SetRadius(0.15);
  merger->Add(effMod);

  EffMod *effMod1 = new EffMod("MuonEff");
  effMod1->SetCol1Name(ModNames::gkMCLeptonsName);
  effMod1->SetCol2Name(muCl->GetOutputName());
  effMod1->SetType(MCParticle::kMu);
  effMod1->SetRadius(0.15);
  merger->Add(effMod1);

  EffMod *effMod2 = new EffMod("ElectronEff");
  effMod2->SetCol1Name(ModNames::gkMCLeptonsName);
  effMod2->SetCol2Name(elCl->GetOutputName());
  effMod2->SetType(MCParticle::kEl);
  effMod2->SetRadius(0.15);
  merger->Add(effMod2);

  Double_t w = 1.;
  TString dataset(gSystem->Getenv("MIT_DATASET"));
  if (!dataset.IsNull()) {
    mithep::SimpleTable xstab("$CMSSW_BASE/src/MitPhysics/data/xs.dat");
    if (xstab.Has(dataset))
      w = xstab.Get(dataset) * 1000; //pb into fb
  }

  H2lNtupMod *www2lNtupMod = new H2lNtupMod;
  www2lNtupMod->SetCleanJetsName(jetCl->GetOutputName());
  www2lNtupMod->SetCleanLeptonsName(merger->GetOutputName());
  www2lNtupMod->SetCleanMetName(pubMet->GetOutputName());
  www2lNtupMod->SetMinJetPt(minJMet);
  www2lNtupMod->SetMinMaxPt(ptMaxMin);
  www2lNtupMod->SetMinPt(ptMin);
  www2lNtupMod->SetMissEt(minJMet);
  www2lNtupMod->SetWeight(w);
  merger->Add(www2lNtupMod);

  HwwwNtupMod *www3lNtupMod = new HwwwNtupMod;
  www3lNtupMod->SetCleanJetsName(jetCl->GetOutputName());
  www3lNtupMod->SetCleanLeptonsName(merger->GetOutputName());
  www3lNtupMod->SetCleanMetName(pubMet->GetOutputName());
  www3lNtupMod->SetMinJetPt(minJMet);
  www3lNtupMod->SetMinMaxPt(ptMaxMin);
  www3lNtupMod->SetMinPt(ptMin);
  www3lNtupMod->SetMissEt(minJMet);
  www3lNtupMod->SetWeight(w);
  merger->Add(www3lNtupMod);

  DilepSelMod *diLepSelMod = new DilepSelMod;
  diLepSelMod->SetCleanLeptonsName(merger->GetOutputName());
  diLepSelMod->SetMinPt(ptMin);
  diLepSelMod->SetMinDilMass(minDi);
  diLepSelMod->SetMinZMass(80);
  diLepSelMod->SetMaxZMass(100);
  merger->Add(diLepSelMod);

  GenericSelMod<Particle> *llSelMod = new GenericSelMod<Particle>("TwoLeptonSel");
  llSelMod->SetInputName(merger->GetOutputName());
  llSelMod->SetPtMin(ptMin);
  llSelMod->SetMinMaxPt(ptMaxMin);
  llSelMod->SetMinCounts(2);
  diLepSelMod->Add(llSelMod);

  HwwAnaMod *h2lMod = new HwwAnaMod;
  h2lMod->SetCleanJetsName(jetCl->GetOutputName());
  h2lMod->SetCleanLeptonsName(merger->GetOutputName());
  h2lMod->SetCleanMetName(pubMet->GetOutputName());
  llSelMod->Add(h2lMod);

  H2lAnaMod *hCount2lMod = new H2lAnaMod("SameSignAnaMod");
  hCount2lMod->SetCleanJetsName(jetCl->GetOutputName());
  hCount2lMod->SetCleanLeptonsName(merger->GetOutputName());
  hCount2lMod->SetCleanMetName(pubMet->GetOutputName());
  hCount2lMod->SetMinZMass(80);
  hCount2lMod->SetMaxZMass(100);
  llSelMod->Add(hCount2lMod);

  HwwwSignalMCMod *hwwwSSMcMod2 = new HwwwSignalMCMod("SameSignMCMod2");
  hwwwSSMcMod2->SetFillHist(1);
  hwwwSSMcMod2->SetDoSameSign(1);
  if (mcomp)
    llSelMod->Add(hwwwSSMcMod2);

  GenericSelMod<Particle> *lllSelMod = new GenericSelMod<Particle>("ThreeLeptonSel");
  lllSelMod->SetInputName(merger->GetOutputName());
  lllSelMod->SetPtMin(ptMin);
  lllSelMod->SetMinMaxPt(ptMaxMin);
  lllSelMod->SetMinCounts(3);
  diLepSelMod->Add(lllSelMod);

  HwwwAnaMod *h3lMod = new HwwwAnaMod;
  h3lMod->SetCleanJetsName(jetCl->GetOutputName());
  h3lMod->SetCleanLeptonsName(merger->GetOutputName());
  h3lMod->SetCleanMetName(pubMet->GetOutputName());
  h3lMod->SetMissEt(minEt);
  h3lMod->SetDileptonMinMass(minDi);
  lllSelMod->Add(h3lMod);

  HwwwSignalMCMod *hwwwMcMod2 = new HwwwSignalMCMod("HwwwSignalMCMod2");
  hwwwMcMod2->SetFillHist(1);
  hwwwMcMod2->SetMissEt(minEt);
  hwwwMcMod2->SetDileptonMinMass(minDi);
  if (mcomp)
    lllSelMod->Add(hwwwMcMod2);

  HwwwSignalMCMod *hwwwMcMod3 = new HwwwSignalMCMod("HwwwSignalMCMod3");
  hwwwMcMod3->SetFillHist(1);
  hwwwMcMod3->SetDoSelect(1);
  hwwwMcMod3->SetMissEt(minEt);
  hwwwMcMod3->SetDileptonMinMass(minDi);
  if (mcomp)
    merger->Add(hwwwMcMod3);

  HwwwSignalMCMod *hwwwSSMcMod3 = new HwwwSignalMCMod("SameSignMCMod3");
  hwwwSSMcMod3->SetFillHist(1);
  hwwwSSMcMod3->SetDoSelect(1);
  hwwwSSMcMod3->SetDoSameSign(1);
  if (mcomp)
    merger->Add(hwwwSSMcMod3);

  DilepSelMod *diLep3SelForMC = new DilepSelMod("DiLep3SelForMC");
  diLep3SelForMC->SetCleanLeptonsName(merger->GetOutputName());
  diLep3SelForMC->SetMinPt(ptMin);
  diLep3SelForMC->SetMinDilMass(minDi);
  diLep3SelForMC->SetMinZMass(80);
  diLep3SelForMC->SetMaxZMass(100);
  hwwwMcMod3->Add(diLep3SelForMC);

  GenericSelMod<Particle> *lllSelMod2 = new GenericSelMod<Particle>("ThreeLeptonSelForMc");
  lllSelMod2->SetInputName(merger->GetOutputName());
  lllSelMod2->SetPtMin(ptMin);
  lllSelMod2->SetMinMaxPt(ptMaxMin);
  lllSelMod2->SetMinCounts(3);
  diLep3SelForMC->Add(lllSelMod2);

  HwwwAnaMod *h3lMod2 = new HwwwAnaMod("HwwwAnaMod2");
  h3lMod2->SetCleanJetsName(jetCl->GetOutputName());
  h3lMod2->SetCleanLeptonsName(merger->GetOutputName());
  h3lMod2->SetCleanMetName(pubMet->GetOutputName());
  h3lMod2->SetMissEt(minEt);
  h3lMod2->SetDileptonMinMass(minDi);
  lllSelMod2->Add(h3lMod2);

  DilepSelMod *diLep2SelForMC = new DilepSelMod("DiLep2SelForMC");
  diLep2SelForMC->SetCleanLeptonsName(merger->GetOutputName());
  diLep2SelForMC->SetMinPt(ptMin);
  diLep2SelForMC->SetMinDilMass(minDi);
  diLep2SelForMC->SetMinZMass(80);
  diLep2SelForMC->SetMaxZMass(100);
  hwwwSSMcMod3->Add(diLep2SelForMC);

  GenericSelMod<Particle> *llSelMod2 = new GenericSelMod<Particle>("TwoLeptonSelForMc");
  llSelMod2->SetInputName(merger->GetOutputName());
  llSelMod2->SetPtMin(ptMin);
  llSelMod2->SetMinMaxPt(ptMaxMin);
  llSelMod2->SetMinCounts(2);
  diLep2SelForMC->Add(llSelMod2);

  H2lAnaMod *hCount2lMod2 = new H2lAnaMod("SameSignAnaMod2");
  hCount2lMod2->SetCleanJetsName(jetCl->GetOutputName());
  hCount2lMod2->SetCleanLeptonsName(merger->GetOutputName());
  hCount2lMod2->SetCleanMetName(pubMet->GetOutputName());
  llSelMod2->Add(hCount2lMod2);

  if (1) {
    PlotKineMod<Muon> *plotAllMus = new PlotKineMod<Muon>("AllMuons");
    plotAllMus->SetInputName(muCl->GetInputName());
    plotAllMus->SetPtMax(360);
    muCl->Add(plotAllMus);

    PlotKineMod<Muon> *plotClMus = new PlotKineMod<Muon>("CleanMuons");
    plotClMus->SetLoadBranch(0);
    plotClMus->SetInputName(muCl->GetOutputName());
    plotClMus->SetPtMax(360);
    muCl->Add(plotClMus);

    PlotKineMod<Electron> *plotAllEls = new PlotKineMod<Electron>("AllElectrons");
    plotAllEls->SetInputName(elId->GetInputName());
    plotAllEls->SetPtMax(360);
    elCl->Add(plotAllEls);

    PlotKineMod<Electron> *plotClEls = new PlotKineMod<Electron>("CleanElectrons");
    plotClEls->SetLoadBranch(0);
    plotClEls->SetInputName(elCl->GetOutputName());
    plotClEls->SetPtMax(360);
    elCl->Add(plotClEls);

    PlotKineMod<Jet> *plotAllJets = new PlotKineMod<Jet>("AllJets");
    plotAllJets->SetLoadBranch(0);
    plotAllJets->SetInputName(jetId->GetInputName());
    plotAllJets->SetPtMax(360);
    jetCl->Add(plotAllJets);

    PlotKineMod<Jet> *plotClJets = new PlotKineMod<Jet>("CleanJets");
    plotClJets->SetLoadBranch(0);
    plotClJets->SetInputName(jetId->GetOutputName());
    plotClJets->SetPtMax(500);
    jetCl->Add(plotClJets);

    PlotKineMod<CaloMet> *plotAllMet = new PlotKineMod<CaloMet>("AllMet");
    plotAllMet->SetInputName(pubMet->GetInputName());
    plotAllMet->SetPtMax(500);
    pubMet->Add(plotAllMet);

    PlotKineMod<CaloMet> *plotMll = new PlotKineMod<CaloMet>("Metll");
    plotMll->SetInputName(pubMet->GetInputName());
    plotMll->SetPtMax(500);
    llSelMod->Add(plotMll);

    PlotKineMod<CaloMet> *plotMlll = new PlotKineMod<CaloMet>("Metlll");
    plotMlll->SetInputName(pubMet->GetInputName());
    plotMlll->SetPtMax(500);
    lllSelMod->Add(plotMlll);
  }

  // run the analysis after successful initialisation
  TStopwatch timer;
  ana->Run(!gROOT->IsBatch());
  timer.Stop();
  timer.Print();
}
