#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "EWKAna/Ntupler/interface/HwwNtuplerMod.hh"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/PDFProducerMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/CaloMetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/PFTauIDMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/PartonFlavorHistoryMod.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h" 
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#include "MitHiggs/EwkMods/interface/SignalFilterForJettinessMod.h"
#include "MitHiggs/HwwMods/interface/HwwMakeNtupleMod.h"

#endif

using namespace mithep;


//==================================================================================================
/*
 * Run on a BAMBU fileset
 *
 * Example usage:
 * root -l -q -b $CMSSW_BASE/src/MitHiggs/macros/runMacros/runAllNtupler.C+\(\"0000\",\"noskim\",\"r11a-dmu-m10-v1\",\"cern/filefi/022\",\"/home/mitprod/catalog\",\"SmurfNtuple\",10000,-1\) 
 *
 * Output file name has standard format: <dataset>_<fileset>_ntuple.root
 *
 */
void runHwwNtupler(
  const char *fileset    = "0000",
  const char *skim       = "noskim",
  const char *dataset    = "r11a-dmu-pr-v4",
  const char *book       = "cern/filefi/025",
  const char *catalogDir = "/home/mitprod/catalog",
  const char *outputName = "SmurfNtuple",
  int   nEvents          = 1000,
  int   sampleID         = -1
  )
{
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //******************************************************************
  //Set up Options
  //******************************************************************
  bool doHWWNtupler       = true;
  bool doSmurfNtupler     = true;
  bool usePDFProducer     = false;
  string pdfSetName = "";

  bool isData             = false;
  bool isDataMuonElectron = false;
  bool isDataDMuon        = false;
  bool isDataSMuon        = false;
  bool isDataDElectron    = false;
  bool isDataSElectron    = false;
  bool isDataSPhoton      = false;
  bool applyMllGenCut     = false;
  bool applyZVetoGenCut   = false;

  double ptJetCut = 30.0;
  double etaJetCut = 5.0;

  int processid = -999999999;
  TString fInputFilenameKF = "/home/ceballos/releases/CMSSW_4_2_2/src/MitPhysics/data/HWW_KFactors_PowhegToNNLL_160_7TeV.dat";
  TString MCType = "kMCTypeUndef";
  Bool_t applyPartonFlavorFilter = kFALSE;
  Bool_t applyISRFilter          = kFALSE;
  Bool_t applyVVFilter           = kFALSE;
  Int_t fakeRatePredictionType = 0;
  Bool_t useSelectGenLeptons   = kTRUE;
  Bool_t isPhotonControlSample = kFALSE;

  int fDecay = sampleID; 
  if (sampleID < 0) fDecay = (-1) * (abs(sampleID) % 10);
  if (sampleID == 10333) fDecay = 3;
  
  Int_t runSkim = 0; if (sampleID <= -1 && sampleID >= -5) runSkim = 1; 
  if (sampleID <= -11 && sampleID >= -15) runSkim = 2; 
  if (sampleID <= -21 && sampleID >= -25) runSkim = 3; 

  if (fDecay < 0) {
    isData = true;
    if (fDecay == -1) isDataDMuon = true;
    if (fDecay == -2) isDataDElectron = true;
    if (fDecay == -3) isDataMuonElectron = true;
    if (fDecay == -4) isDataSMuon = true;
    if (fDecay == -5) isDataSElectron = true;
    if (fDecay == -6) isDataSPhoton = true;
  }

  if (fDecay == 14) applyVVFilter = kTRUE;
  if (fDecay > 20000) {
    fakeRatePredictionType = 1;
    useSelectGenLeptons    = kFALSE;
  }
  if (isDataSPhoton || fDecay == 10666 || fDecay == 11666) {
    isPhotonControlSample = kTRUE;
    doHWWNtupler = kFALSE;
  }

  cout << "Summarize Run Options: " << "fDecay == " << fDecay << " "
       << "runSkim == " << runSkim << " "
       << endl;

  //******************************************************************
  //Modules
  //******************************************************************

  // Generator info
  GeneratorMod *GeneratorMod1 = new GeneratorMod;
  GeneratorMod1->SetPrintDebug(kFALSE);
  GeneratorMod1->SetPtLeptonMin(0.0);
  GeneratorMod1->SetEtaLeptonMax(2.7);
  GeneratorMod1->SetPtPhotonMin(15.0);
  GeneratorMod1->SetEtaPhotonMax(2.7);
  GeneratorMod1->SetPtRadPhotonMin(10.0);
  GeneratorMod1->SetEtaRadPhotonMax(2.7);
  GeneratorMod1->SetIsData(isData);
  GeneratorMod1->SetFillHist(!isData);
  if(applyMllGenCut == kTRUE){
    GeneratorMod1->SetPdgIdCut(23);
    GeneratorMod1->SetMassMaxCut(50.);
  }
  else if(applyZVetoGenCut == kTRUE){
    GeneratorMod1->SetPdgIdCut(23);
    GeneratorMod1->SetMassMinCut(20000.);
    GeneratorMod1->SetMassMaxCut(20000.);
  }
  GeneratorMod1->SetApplyISRFilter(applyISRFilter);
  GeneratorMod1->SetApplyVVFilter(applyVVFilter);
  GeneratorMod1->SetAllowWWEvents(kTRUE);
  GeneratorMod1->SetAllowWZEvents(kFALSE);
  GeneratorMod1->SetAllowZZEvents(kFALSE);

  PartonFlavorHistoryMod *PartonFlavorHistoryMod1 = new PartonFlavorHistoryMod;
  PartonFlavorHistoryMod1->SetMCSampleType(MCType);
  PartonFlavorHistoryMod1->SetApplyPartonFlavorFilter(applyPartonFlavorFilter);

  HLTMod *hltmod = new HLTMod;
  if(isData == true && isDataMuonElectron == true) {
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",150000,161176);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",150000,161176);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",161179,163261);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",161179,163261);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v3",163262,164237);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",163262,164237);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v4",165085,165888);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v4",165085,165888);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v5",165900,166967);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v5",165900,166967);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v6",166968,170053);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v6",166968,170053);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v8",170054,173198);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3",170054,173198);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4",173199,178380);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4",173199,178380);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7",178381,179889);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7",178381,179889);
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8",179890,999999);
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8",179890,999999);
  }
  else if(isData == true && isDataDMuon == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&HLT_DoubleMu7_v1",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&HLT_DoubleMu7_v1",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&HLT_DoubleMu7_v2",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&HLT_Mu13_Mu8_v2" ,165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&HLT_Mu13_Mu8_v2" ,165900,167043);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&HLT_Mu13_Mu8_v4" ,167044,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&HLT_Mu13_Mu8_v6" ,170054,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&HLT_Mu13_Mu8_v7" ,173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&HLT_Mu17_Mu8_v10" ,178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&HLT_Mu17_TkMu8_v3" ,178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&HLT_Mu17_Mu8_v11" ,179890,999999);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&HLT_Mu17_TkMu8_v4" ,179890,999999);
  }
  else if(isData == true && isDataSMuon == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&!HLT_DoubleMu7_v1&HLT_Mu15_v2",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&!HLT_DoubleMu7_v1&HLT_Mu15_v2",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&HLT_Mu24_v2",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&HLT_IsoMu17_v6",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&!HLT_Mu13_Mu8_v2&HLT_Mu30_v3",165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&!HLT_Mu13_Mu8_v2&HLT_IsoMu17_v8",165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&HLT_Mu30_v3",165900,167043);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&HLT_IsoMu17_v9",165900,167043);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v4&HLT_Mu30_v5",167044,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v4&HLT_IsoMu17_eta2p1_v1",167044,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&HLT_Mu40_v5",170054,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&HLT_IsoMu24_v8",170054,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&HLT_Mu40_eta2p1_v1",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&HLT_IsoMu30_eta2p1_v3",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&HLT_IsoMu24_eta2p1_v3",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&HLT_Mu40_eta2p1_v4",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&HLT_IsoMu30_eta2p1_v6",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&HLT_IsoMu24_eta2p1_v6",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&HLT_Mu40_eta2p1_v5",179890,999999);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&HLT_IsoMu30_eta2p1_v7",179890,999999);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&HLT_IsoMu24_eta2p1_v7",179890,999999);
  }
  else if(isData == true && isDataDElectron == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&!HLT_Mu24_v2&!HLT_IsoMu17_v6&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&!HLT_Mu13_Mu8_v2&!HLT_Mu30_v3&!HLT_IsoMu17_v8&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&!HLT_Mu13_Mu8_v4&!HLT_Mu30_v3&!HLT_Mu30_v5&!HLT_IsoMu17_v9&!HLT_IsoMu17_eta2p1_v1&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",165900,166967);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&!HLT_Mu13_Mu8_v4&!HLT_Mu30_v3&!HLT_Mu30_v5&!HLT_IsoMu17_v9&!HLT_IsoMu17_eta2p1_v1&HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",166968,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&!HLT_Mu40_v5&!HLT_IsoMu24_v8&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",170054,170759);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&!HLT_Mu40_v5&!HLT_IsoMu24_v8&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",170760,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&!HLT_Mu40_eta2p1_v1&!HLT_IsoMu30_eta2p1_v3&!HLT_IsoMu24_eta2p1_v3&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&!HLT_Mu40_eta2p1_v4&!HLT_IsoMu30_eta2p1_v6&!HLT_IsoMu24_eta2p1_v6&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&!HLT_Mu40_eta2p1_v5&!HLT_IsoMu30_eta2p1_v7&!HLT_IsoMu24_eta2p1_v7&HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",179890,999999);
  }
  else if(isData == true && isDataSElectron == true) {
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v1&!HLT_Mu17_Ele8_CaloIdL_v1&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1&HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",150000,161176);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v2&!HLT_Mu17_Ele8_CaloIdL_v2&!HLT_DoubleMu7_v1&!HLT_Mu15_v2&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2&HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",161179,163261);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v3&!HLT_Mu17_Ele8_CaloIdL_v3&!HLT_DoubleMu7_v2&!HLT_Mu24_v2&!HLT_IsoMu17_v6&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3&HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",163262,164237);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v4&!HLT_Mu17_Ele8_CaloIdL_v4&!HLT_Mu13_Mu8_v2&!HLT_Mu30_v3&!HLT_IsoMu17_v8&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4&HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",165085,165888);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&!HLT_Mu13_Mu8_v4&!HLT_Mu30_v3&!HLT_Mu30_v5&!HLT_IsoMu17_v9&!HLT_IsoMu17_eta2p1_v1&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5&HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4",165900,166967);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdL_v5&!HLT_Mu17_Ele8_CaloIdL_v5&!HLT_Mu8_Ele17_CaloIdL_v6&!HLT_Mu17_Ele8_CaloIdL_v6&!HLT_Mu13_Mu8_v2&!HLT_Mu13_Mu8_v4&!HLT_Mu30_v3&!HLT_Mu30_v5&!HLT_IsoMu17_v9&!HLT_IsoMu17_eta2p1_v1&!HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6&HLT_Ele52_CaloIdVT_TrkIdT_v3",166968,170053);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3&!HLT_Mu17_Ele8_CaloIdL_v8&!HLT_Mu13_Mu8_v6&!HLT_Mu40_v5&!HLT_IsoMu24_v8&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7&HLT_Ele65_CaloIdVT_TrkIdT_v3",170054,173198);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4&!HLT_Mu13_Mu8_v7&!HLT_Mu40_eta2p1_v1&!HLT_IsoMu30_eta2p1_v3&!HLT_IsoMu24_eta2p1_v3&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8&HLT_Ele65_CaloIdVT_TrkIdT_v4",173199,178380);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7&!HLT_Mu17_Mu8_v10&!HLT_Mu17_TkMu8_v3&!HLT_Mu40_eta2p1_v4&!HLT_IsoMu30_eta2p1_v6&!HLT_IsoMu24_eta2p1_v6&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9&HLT_Ele80_CaloIdVT_TrkIdT_v2",178381,179889);
    hltmod->AddTrigger("!HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8&!HLT_Mu17_Mu8_v11&!HLT_Mu17_TkMu8_v4&!HLT_Mu40_eta2p1_v5&!HLT_IsoMu30_eta2p1_v7&!HLT_IsoMu24_eta2p1_v7&!HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10&HLT_Ele80_CaloIdVT_TrkIdT_v3",179890,999999);
  }
  else if(isData == true && isDataSPhoton == true) {
    hltmod->AddTrigger("HLT_Photon20_CaloIdVL_IsoL_v1");
    hltmod->AddTrigger("HLT_Photon20_CaloIdVL_IsoL_v2");
    hltmod->AddTrigger("HLT_Photon20_CaloIdVL_IsoL_v3");
    hltmod->AddTrigger("HLT_Photon20_CaloIdVL_IsoL_v4");
    hltmod->AddTrigger("HLT_Photon20_CaloIdVL_IsoL_v5");
    hltmod->AddTrigger("HLT_Photon20_CaloIdVL_IsoL_v6");
    hltmod->AddTrigger("HLT_Photon20_CaloIdVL_IsoL_v7");
    hltmod->AddTrigger("HLT_Photon20_CaloIdVL_IsoL_v8");
    hltmod->AddTrigger("HLT_Photon20_CaloIdVL_IsoL_v9");

    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_v1");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_v2");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_v3");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_v4");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_v5");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_v6");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_v7");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_v8");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v1");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v2");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v3");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v4");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v5");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v6");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v7");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v8");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v9");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v10");
    hltmod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v11");

    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v1");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v2");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v3");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v4");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v5");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v6");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v7");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v8");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v9");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_v1");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_v2");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_v3");
    hltmod->AddTrigger("HLT_Photon50_CaloIdVL_v4");

    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_v1");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_v2");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_v3");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_v4");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_v5");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_v6");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_v7");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v1");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v2");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v3");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v4");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v5");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v6");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v7");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v8");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v9");
    hltmod->AddTrigger("HLT_Photon75_CaloIdVL_IsoL_v10");

    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_IsoL_v1");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_IsoL_v2");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_IsoL_v3");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_IsoL_v4");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_IsoL_v5");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_IsoL_v6");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_IsoL_v7");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_v1");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_v2");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_v3");
    hltmod->AddTrigger("HLT_Photon90_CaloIdVL_v4");
    hltmod->AddTrigger("HLT_Photon125_v1");
    hltmod->AddTrigger("HLT_Photon125_v2");
    hltmod->AddTrigger("HLT_Photon135_v1");
    hltmod->AddTrigger("HLT_Photon135_v2");
    hltmod->AddTrigger("HLT_Photon200_NoHE_v1");
    hltmod->AddTrigger("HLT_Photon200_NoHE_v2");
    hltmod->AddTrigger("HLT_Photon200_NoHE_v3");
    hltmod->AddTrigger("HLT_Photon200_NoHE_v4");
    hltmod->AddTrigger("HLT_Photon400_v1");
    hltmod->AddTrigger("HLT_Photon400_v2");
  } 
  else if(sampleID == 10333 || sampleID == 11666 || sampleID == 10666 ) {
    hltmod->AddTrigger("HLT_Mu15_v1");
    hltmod->AddTrigger("!HLT_Mu15_v1");
  }
  else {
    hltmod->AddTrigger("HLT_Mu15_v2");
    hltmod->AddTrigger("HLT_IsoMu17_v5");
    hltmod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2");
    hltmod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2");
    hltmod->AddTrigger("HLT_DoubleMu7_v1");
    hltmod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2");
    hltmod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2"); 
    hltmod->AddTrigger("!HLT_Mu15_v2");
    hltmod->SetAbortIfNotAccepted(kFALSE);
  }
  hltmod->SetTrigObjsName("myhltobjs");


  //------------------------------------------------------------------------------------------------
  // Dilepton Skim
  //------------------------------------------------------------------------------------------------
  MuonIDMod *muTightId = new MuonIDMod;  
  muTightId->SetClassType ("GlobalTracker");
  muTightId->SetIDType    ("WWMuIdV3");
  muTightId->SetIsoType   ("PFIso");
  muTightId->SetApplyD0Cut(kTRUE);
  muTightId->SetApplyDZCut(kTRUE);
  muTightId->SetWhichVertex(0);
  muTightId->SetOutputName("SkimMuons");

  ElectronIDMod       *electronTightId       = new ElectronIDMod;
  electronTightId->SetIDType(TString("VBTFWorkingPointLowPtId"));
  electronTightId->SetIsoType(TString("PFIso"));
  electronTightId->SetApplyConversionFilterType1(kTRUE);
  electronTightId->SetApplyConversionFilterType2(kFALSE);
  electronTightId->SetChargeFilter(kFALSE);
  electronTightId->SetApplyD0Cut(kTRUE);
  electronTightId->SetApplyDZCut(kTRUE);
  electronTightId->SetWhichVertex(0);
  electronTightId->SetNExpectedHitsInnerCut(0);
  electronTightId->SetOutputName("SkimElectrons");


  MuonIDMod *muDenominator = new MuonIDMod;  
  muDenominator->SetClassType ("GlobalTracker");
  muDenominator->SetIDType    ("WWMuIdV3");
  muDenominator->SetIsoType   ("PFIsoEffectiveAreaCorrected");
  muDenominator->SetApplyD0Cut(kTRUE);
  muDenominator->SetApplyDZCut(kTRUE);
  muDenominator->SetD0Cut(0.20);
  muDenominator->SetDZCut(0.10);
  muDenominator->SetPFIsoCut(0.40); 
  muDenominator->SetOutputName("SkimDenominatorMuons");
  muDenominator->SetWhichVertex(0);

  ElectronIDMod *electronDenominator = new ElectronIDMod;
  electronDenominator->SetIDType("VBTFWorkingPointFakeableId");
  electronDenominator->SetIsoType("TrackJura");
  electronDenominator->SetTrackIsoCut(0.2);
  electronDenominator->SetEcalJurIsoCut(0.2);
  electronDenominator->SetHcalIsoCut(0.2);
  electronDenominator->SetApplyConversionFilterType1(kTRUE);
  electronDenominator->SetApplyConversionFilterType2(kFALSE);
  electronDenominator->SetChargeFilter              (kFALSE);
  electronDenominator->SetApplyD0Cut                (kTRUE);
  electronDenominator->SetApplyDZCut                (kTRUE);
  electronDenominator->SetNExpectedHitsInnerCut     (0);
  electronDenominator->SetD0Cut(0.02);
  electronDenominator->SetDZCut(0.10);
  electronDenominator->SetOutputName("SkimDenominatorElectrons");

  MergeLeptonsMod *mergedTight = new MergeLeptonsMod;
  mergedTight->SetMuonsName    (muTightId->GetOutputName());
  mergedTight->SetElectronsName(electronTightId->GetOutputName());
  mergedTight->SetOutputName("mergedTightLeptons");

  MergeLeptonsMod *mergedLoose = new MergeLeptonsMod;
  mergedLoose->SetMuonsName    (muDenominator->GetOutputName());
  mergedLoose->SetElectronsName(electronDenominator->GetOutputName());
  mergedLoose->SetOutputName("mergedLooseLeptons");

  GenericSelMod<mithep::Particle> *selModTight = new GenericSelMod<mithep::Particle>;
  selModTight->SetPtMin(0.0);
  selModTight->SetMinCounts(1);
  selModTight->SetColName(mergedTight->GetOutputName());

  GenericSelMod<mithep::Particle> *selModLoose = new GenericSelMod<mithep::Particle>;
  selModLoose->SetPtMin(0.0);
  selModLoose->SetMinCounts(1);
  selModLoose->SetColName(mergedLoose->GetOutputName());

  GenericSelMod<mithep::Particle> *selModDoubleLoose = new GenericSelMod<mithep::Particle>;
  selModDoubleLoose->SetPtMin(0.0);
  selModDoubleLoose->SetMinCounts(2);
  selModDoubleLoose->SetColName(mergedLoose->GetOutputName());


  //------------------------------------------------------------------------------------------------
  // Run RunLumiSelectionMod
  //------------------------------------------------------------------------------------------------
  RunLumiSelectionMod *runLumiSelection = new RunLumiSelectionMod;      
  runLumiSelection->SetAcceptMC(!isData);
  runLumiSelection->SetAcceptAll(kTRUE);
  runLumiSelection->AddJSONFile("certifiedUCSD.json");
  runLumiSelection->SetAbortIfNotAccepted(kFALSE);

  //------------------------------------------------------------------------------------------------
  // KFactor Producer Mainly for MC@NLO
  //------------------------------------------------------------------------------------------------
  HKFactorProducer *HKFactorProducer1 = new HKFactorProducer;
  HKFactorProducer1->SetProcessID(0); //0 means don't use weights. 998 uses weights
  HKFactorProducer1->SetFillHist(!isData);

  //------------------------------------------------------------------------------------------------
  // PV filter selection
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);
  if (sampleID == 10333 || sampleID > 10000) { goodPVFilterMod->SetVertexesName("DAPrimaryVertexes"); cout << "use DAPrimaryVertexes\n";}
  else goodPVFilterMod->SetVertexesName("PrimaryVertexes");

  Bool_t isFastSim = kFALSE;

  PDFProducerMod *PDFProducerMod1 = new PDFProducerMod;
  PDFProducerMod1->SetPrintDebug(kFALSE);
  PDFProducerMod1->SetFillHist(kTRUE);
  PDFProducerMod1->SetRunPDF(usePDFProducer);
  PDFProducerMod1->SetPDFName("cteq65.LHgrid");
  //PDFProducerMod1->SetPDFName("cteq6ll.LHpdf");
  PDFProducerMod1->SetIsData(isData);

  //***********************************************************************************8
  //Lepton Selection
  //***********************************************************************************8

  // Object ID and Cleaning Sequence
  MuonIDMod *muonID1 = new MuonIDMod;
  muonID1->SetClassType("GlobalTracker");
  muonID1->SetIDType("WWMuIdV3");
  muonID1->SetIsoType("PFIso");
  muonID1->SetApplyD0Cut(kTRUE);
  muonID1->SetApplyDZCut(kTRUE);
  muonID1->SetWhichVertex(0);

  ElectronIDMod *electronID1 = new ElectronIDMod;
  electronID1->SetIDType("VBTFWorkingPointLowPtId");
  electronID1->SetIsoType("PFIso");
  electronID1->SetApplyConversionFilterType1(kTRUE);
  electronID1->SetApplyConversionFilterType2(kFALSE);
  electronID1->SetChargeFilter(kFALSE);
  electronID1->SetApplyD0Cut(kTRUE);
  electronID1->SetApplyDZCut(kTRUE);
  electronID1->SetWhichVertex(0);
  electronID1->SetNExpectedHitsInnerCut(0);

  PhotonIDMod *photonIDMod1 = new PhotonIDMod;
  photonIDMod1->SetIsoType("MITPUCorrected");
  photonIDMod1->SetHadOverEmMax(0.05);
  photonIDMod1->SetApplyPixelSeed(kFALSE);
  photonIDMod1->SetApplyElectronVetoConvRecovery(kTRUE);
  photonIDMod1->SetApplyConversionId(kTRUE);

  PFTauIDMod *pftauIDMod1 = new PFTauIDMod;
  pftauIDMod1->SetPFTausName("HPSTaus");
  pftauIDMod1->SetIsHPSSel(kFALSE);


  //***********************************************************************************8
  //Jet Selection
  //***********************************************************************************8

  const char *jetInput1 = "AKt5PFJets";
  PublisherMod<PFJet,Jet> *pubJet1 = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet1->SetInputName(jetInput1);
  pubJet1->SetOutputName(Form("Pub%s",jetInput1));

  JetCorrectionMod *jetCorr1_ntuple = new JetCorrectionMod;
  jetCorr1_ntuple->AddCorrectionFromFile("MitPhysics/data/START42_V12_AK5PF_L1FastJet.txt"); 
  jetCorr1_ntuple->AddCorrectionFromFile("MitPhysics/data/START42_V12_AK5PF_L2Relative.txt"); 
  jetCorr1_ntuple->AddCorrectionFromFile("MitPhysics/data/START42_V12_AK5PF_L3Absolute.txt");
  if(isData == true){
    jetCorr1_ntuple->AddCorrectionFromFile("MitPhysics/data/START42_V12_AK5PF_L2L3Residual.txt");  
  }
  jetCorr1_ntuple->SetInputName(pubJet1->GetOutputName());
  jetCorr1_ntuple->SetCorrectedName("CorrectedJets_ntuple");

  JetIDMod *theJetID1_ntuple = new JetIDMod;
  theJetID1_ntuple->SetInputName(jetCorr1_ntuple->GetOutputName());
  theJetID1_ntuple->SetPtCut(20.0);
  theJetID1_ntuple->SetEtaMaxCut(etaJetCut);
  theJetID1_ntuple->SetJetEEMFractionMinCut(0.00);
  theJetID1_ntuple->SetOutputName("GoodJets_ntuple");
  theJetID1_ntuple->SetApplyBetaCut(kFALSE);



  //***********************************************************************************8
  //Cleaning
  //***********************************************************************************8
  ElectronCleaningMod *electronCleaning1 = new ElectronCleaningMod;
  PhotonCleaningMod *photonCleaningMod1 = new PhotonCleaningMod;
  PFTauCleaningMod *pftauCleaningMod1 = new PFTauCleaningMod;

  JetCleaningMod *theJetCleaning1_ntuple = new JetCleaningMod;
  theJetCleaning1_ntuple->SetGoodJetsName("GoodJets_ntuple");
  theJetCleaning1_ntuple->SetCleanJetsName("CleanJets_ntuple");
  if (isPhotonControlSample) {
    theJetCleaning1_ntuple->SetApplyPhotonRemoval(kTRUE);
  }

  MergeLeptonsMod *merger1 = new MergeLeptonsMod;
  merger1->SetMuonsName(muonID1->GetOutputName());
  merger1->SetElectronsName(electronCleaning1->GetOutputName());

  JetIDMod *theJetID2_ntuple = new JetIDMod;
  theJetID2_ntuple->SetInputName(jetCorr1_ntuple->GetOutputName());
  theJetID2_ntuple->SetPtCut(0.0);
  theJetID2_ntuple->SetEtaMaxCut(etaJetCut);
  theJetID2_ntuple->SetJetEEMFractionMinCut(0.00);
  theJetID2_ntuple->SetOutputName("GoodJetsNoPtCut_ntuple");
  theJetID2_ntuple->SetApplyBetaCut(kFALSE);

  JetCleaningMod *theJetCleaning2_ntuple = new JetCleaningMod;
  theJetCleaning2_ntuple->SetGoodJetsName("GoodJetsNoPtCut_ntuple");
  theJetCleaning2_ntuple->SetCleanJetsName("CleanJetsNoPtCut_ntuple");
  if (isPhotonControlSample) {
    theJetCleaning2_ntuple->SetApplyPhotonRemoval(kTRUE);
  }

  //***********************************************************************************8
  //MET
  //***********************************************************************************8
  const char *metInput = "TCMet";
  PublisherMod<Met,Met> *pubMet = new PublisherMod<Met,Met>("MetPub");
  pubMet->SetInputName(metInput);
  pubMet->SetOutputName(Form("Pub%s",metInput));

  PublisherMod<CaloMet> *pubCaloMet = new PublisherMod<CaloMet>;
  pubCaloMet->SetName("CaloMetPub");
  pubCaloMet->SetInputName("CorMuonMet");
  pubCaloMet->SetOutputName("pubCaloMet");

  CaloMetCorrectionMod *metCaloCorr = new CaloMetCorrectionMod;
  metCaloCorr->SetInputName(pubCaloMet->GetOutputName());
  metCaloCorr->SetCorrectedJetsName(jetCorr1_ntuple->GetOutputName());
  metCaloCorr->SetOutputName("pubCaloCorrectedMet");

  const char *metPFInput = "PFMet";
  PublisherMod<PFMet,Met> *pubPFMet = new PublisherMod<PFMet,Met>("MetPFPub");
  pubPFMet->SetInputName(metPFInput);
  pubPFMet->SetOutputName(Form("Pub%s",metPFInput));

  //***********************************************************************************8
  //Fakeable Objects Definition
  //***********************************************************************************8

  // Lepton ID with loose requirements
  MuonIDMod *muonIDFakeable = new MuonIDMod;
  muonIDFakeable->SetClassType("GlobalTracker");
  muonIDFakeable->SetIDType("WWMuIdV3");
  muonIDFakeable->SetIsoType("PFIsoEffectiveAreaCorrected");
  muonIDFakeable->SetApplyD0Cut(kTRUE);
  muonIDFakeable->SetApplyDZCut(kTRUE);
  muonIDFakeable->SetD0Cut(0.20);
  muonIDFakeable->SetPFIsoCut(0.40);
  muonIDFakeable->SetCleanMuonsName("CleanMuonsFakeable");
  muonIDFakeable->SetWhichVertex(0);

  ElectronIDMod *electronIDFakeable = new ElectronIDMod;
  electronIDFakeable->SetIDType("VBTFWorkingPointFakeableId");
  electronIDFakeable->SetIsoType("TrackJura");
  electronIDFakeable->SetTrackIsoCut(0.2);
  electronIDFakeable->SetEcalJurIsoCut(0.2);
  electronIDFakeable->SetHcalIsoCut(0.2);
  electronIDFakeable->SetApplyConversionFilterType1(kTRUE);
  electronIDFakeable->SetApplyConversionFilterType2(kFALSE);
  electronIDFakeable->SetChargeFilter(kFALSE);
  electronIDFakeable->SetApplyD0Cut(kTRUE);
  electronIDFakeable->SetApplyDZCut(kTRUE);
  electronIDFakeable->SetNExpectedHitsInnerCut(0);
  electronIDFakeable->SetD0Cut(0.02);
  electronIDFakeable->SetGoodElectronsName("GoodElectronsFakeable");
  electronIDFakeable->SetWhichVertex(0);

  ElectronCleaningMod *electronCleaningFakeable = new ElectronCleaningMod;
  electronCleaningFakeable->SetCleanMuonsName(muonIDFakeable->GetOutputName());
  electronCleaningFakeable->SetGoodElectronsName(electronIDFakeable->GetOutputName());
  electronCleaningFakeable->SetCleanElectronsName("CleanElectronsFakeable");

  MergeLeptonsMod *mergerFakeable = new MergeLeptonsMod;
  mergerFakeable->SetMuonsName(muonIDFakeable->GetOutputName());
  mergerFakeable->SetElectronsName(electronCleaningFakeable->GetOutputName());
  mergerFakeable->SetMergedName("MergedLeptonsFakeable");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile =  TString("test/") + TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  printf("\nRoot output: %s\n\n",rootFile.Data());  


  //***********************************************************************************8
  //Smurf Ntupler
  //***********************************************************************************8
  TString rootFileHwwMake0 = TString("test/");
  rootFileHwwMake0 += TString(outputName);
  rootFileHwwMake0 += TString("_smurf0_") + TString(dataset) + TString("_") + TString(skim); 
  if (TString(fileset) != TString(""))
    rootFileHwwMake0 += TString("_") + TString(fileset);
  rootFileHwwMake0 += TString(".root");
  HwwMakeNtupleMod *HwwMakeNtupleMod0 = new HwwMakeNtupleMod;
  HwwMakeNtupleMod0->SetCleanJetsName("CleanJets_ntuple");
  HwwMakeNtupleMod0->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  HwwMakeNtupleMod0->SetMetName(pubPFMet->GetOutputName());
  HwwMakeNtupleMod0->SetJetScaleSyst(0.0);
  HwwMakeNtupleMod0->SetPtJetCut(ptJetCut);
  HwwMakeNtupleMod0->SetEtaJetCut(etaJetCut);
  HwwMakeNtupleMod0->SetProcessID(fDecay);
  HwwMakeNtupleMod0->SetIsData(isData);
  HwwMakeNtupleMod0->SetFakeRatePredictionType(0);
  HwwMakeNtupleMod0->SetMuonFakeName("CleanMuonsFakeable");
  HwwMakeNtupleMod0->SetElectronFakeName("CleanElectronsFakeable");
  HwwMakeNtupleMod0->SetLeptonFakeName("MergedLeptonsFakeable");
  HwwMakeNtupleMod0->SetFillNtupleType(0);
  HwwMakeNtupleMod0->SetOutputName(rootFileHwwMake0);


  TString rootFileHwwMake1 = TString("test/");
  rootFileHwwMake1 += TString(outputName);
  rootFileHwwMake1 += TString("_smurf1_") + TString(dataset) + TString("_") + TString(skim); 
  if (TString(fileset) != TString(""))
    rootFileHwwMake1 += TString("_") + TString(fileset);
  rootFileHwwMake1 += TString(".root");
  HwwMakeNtupleMod *HwwMakeNtupleMod1 = new HwwMakeNtupleMod;
  HwwMakeNtupleMod1->SetCleanJetsName("CleanJets_ntuple");
  HwwMakeNtupleMod1->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  HwwMakeNtupleMod1->SetMetName(pubPFMet->GetOutputName());
  HwwMakeNtupleMod1->SetJetScaleSyst(0.0);
  HwwMakeNtupleMod1->SetPtJetCut(ptJetCut);
  HwwMakeNtupleMod1->SetEtaJetCut(etaJetCut);
  HwwMakeNtupleMod1->SetProcessID(fDecay);
  HwwMakeNtupleMod1->SetIsData(isData);
  HwwMakeNtupleMod1->SetFakeRatePredictionType(1);
  HwwMakeNtupleMod1->SetMuonFakeName("CleanMuonsFakeable");
  HwwMakeNtupleMod1->SetElectronFakeName("CleanElectronsFakeable");
  HwwMakeNtupleMod1->SetLeptonFakeName("MergedLeptonsFakeable");
  HwwMakeNtupleMod1->SetFillNtupleType(1);
  HwwMakeNtupleMod1->SetOutputName(rootFileHwwMake1);


  TString rootFileHwwMake2 = TString("test/");
  rootFileHwwMake2 += TString(outputName);
  rootFileHwwMake2 += TString("_smurf2_") + TString(dataset) + TString("_") + TString(skim); 
  if (TString(fileset) != TString(""))
    rootFileHwwMake2 += TString("_") + TString(fileset);
  rootFileHwwMake2 += TString(".root");
  HwwMakeNtupleMod *HwwMakeNtupleMod2 = new HwwMakeNtupleMod;
  HwwMakeNtupleMod2->SetCleanJetsName("CleanJets_ntuple");
  HwwMakeNtupleMod2->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  HwwMakeNtupleMod2->SetMetName(pubPFMet->GetOutputName());
  HwwMakeNtupleMod2->SetJetScaleSyst(0.0);
  HwwMakeNtupleMod2->SetPtJetCut(ptJetCut);
  HwwMakeNtupleMod2->SetEtaJetCut(etaJetCut);
  HwwMakeNtupleMod2->SetProcessID(fDecay);
  HwwMakeNtupleMod2->SetIsData(isData);
  HwwMakeNtupleMod2->SetFakeRatePredictionType(2);
  HwwMakeNtupleMod2->SetMuonFakeName("CleanMuonsFakeable");
  HwwMakeNtupleMod2->SetElectronFakeName("CleanElectronsFakeable");
  HwwMakeNtupleMod2->SetLeptonFakeName("MergedLeptonsFakeable");
  HwwMakeNtupleMod2->SetFillNtupleType(2);
  HwwMakeNtupleMod2->SetOutputName(rootFileHwwMake2);


  TString rootFileHwwMakePhoton = TString("test/");
  rootFileHwwMakePhoton += TString(outputName);
  rootFileHwwMakePhoton += TString("_smurfPhoton_") + TString(dataset) + TString("_") + TString(skim); 
  if (TString(fileset) != TString(""))
    rootFileHwwMakePhoton += TString("_") + TString(fileset);
  rootFileHwwMakePhoton += TString(".root");
  HwwMakeNtupleMod *HwwMakeNtupleModPhoton = new HwwMakeNtupleMod;
  HwwMakeNtupleModPhoton->SetCleanJetsName("CleanJets_ntuple");
  HwwMakeNtupleModPhoton->SetCleanJetsNoPtCutName("CleanJetsNoPtCut_ntuple");
  HwwMakeNtupleModPhoton->SetMetName(pubPFMet->GetOutputName());
  HwwMakeNtupleModPhoton->SetJetScaleSyst(0.0);
  HwwMakeNtupleModPhoton->SetPtJetCut(ptJetCut);
  HwwMakeNtupleModPhoton->SetEtaJetCut(etaJetCut);
  HwwMakeNtupleModPhoton->SetProcessID(fDecay);
  HwwMakeNtupleModPhoton->SetIsData(isData);
  HwwMakeNtupleModPhoton->SetFakeRatePredictionType(0);
  HwwMakeNtupleModPhoton->SetFillPhotonTemplate(kTRUE);
  HwwMakeNtupleModPhoton->SetMuonFakeName("CleanMuonsFakeable");
  HwwMakeNtupleModPhoton->SetElectronFakeName("CleanElectronsFakeable");
  HwwMakeNtupleModPhoton->SetLeptonFakeName("MergedLeptonsFakeable");
  HwwMakeNtupleModPhoton->SetFillNtupleType(0);
  HwwMakeNtupleModPhoton->SetOutputName(rootFileHwwMakePhoton);

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // Chain modules together
  if(isData == false){
    GeneratorMod1->Add(PartonFlavorHistoryMod1);
    PartonFlavorHistoryMod1->Add(runLumiSelection);
    runLumiSelection->Add(HKFactorProducer1);
    HKFactorProducer1->Add(goodPVFilterMod);
  }
  else {
    GeneratorMod1->Add(runLumiSelection);
    runLumiSelection->Add(goodPVFilterMod);
  }

  //------------------------------------------------------------------------------------------------
  // Run Lepton + Denominator Skim
  //------------------------------------------------------------------------------------------------

  //loose+loose
  if (runSkim == 1) {
    goodPVFilterMod->Add(muTightId);
    muTightId->Add(electronTightId);
    electronTightId->Add(muDenominator);
    muDenominator->Add(electronDenominator);
    electronDenominator->Add(mergedTight);
    mergedTight->Add(mergedLoose);
    mergedLoose->Add(selModDoubleLoose);
    selModDoubleLoose->Add(muonID1); 
  } 
  
  //loose
  else if (runSkim == 2) {
    goodPVFilterMod->Add(muTightId);
    muTightId->Add(electronTightId);
    electronTightId->Add(muDenominator);
    muDenominator->Add(electronDenominator);
    electronDenominator->Add(mergedTight);
    mergedTight->Add(mergedLoose);
    mergedLoose->Add(selModLoose);
    selModLoose->Add(muonID1); 
  }
  //tight+loose
  else if (runSkim == 3) {
    goodPVFilterMod->Add(muTightId);
    muTightId->Add(electronTightId);
    electronTightId->Add(muDenominator);
    muDenominator->Add(electronDenominator);
    electronDenominator->Add(mergedTight);
    mergedTight->Add(mergedLoose);
    mergedLoose->Add(selModTight);
    selModTight->Add(selModDoubleLoose); 
    selModDoubleLoose->Add(muonID1);
  } else {
   goodPVFilterMod->Add(muonID1);
  }


  //Standard Sequence
  muonID1->Add(electronID1);
  electronID1->Add(photonIDMod1);
  photonIDMod1->Add(pftauIDMod1);
  pftauIDMod1->Add(pubJet1);
  pubJet1->Add(jetCorr1_ntuple);
  jetCorr1_ntuple->Add(theJetID1_ntuple);
  theJetID1_ntuple->Add(electronCleaning1);
  electronCleaning1->Add(photonCleaningMod1);
  photonCleaningMod1->Add(pftauCleaningMod1);
  pftauCleaningMod1->Add(theJetCleaning1_ntuple);
  theJetCleaning1_ntuple->Add(theJetID2_ntuple);
  theJetID2_ntuple->Add(theJetCleaning2_ntuple);
  theJetCleaning2_ntuple->Add(merger1);
  merger1->Add(pubMet);
  pubMet->Add(pubCaloMet);
  pubCaloMet->Add(metCaloCorr);
  metCaloCorr->Add(pubPFMet);
  pubPFMet->Add(muonIDFakeable);
  muonIDFakeable->Add(electronIDFakeable);
  electronIDFakeable->Add(electronCleaningFakeable);
  electronCleaningFakeable->Add(mergerFakeable);

  //SmurfNtupler
  if (doSmurfNtupler) {
    mergerFakeable->Add(hltmod);
    
    if (isPhotonControlSample) {
      hltmod->Add(HwwMakeNtupleModPhoton);
    } else {
      hltmod->Add(HwwMakeNtupleMod0);
      HwwMakeNtupleMod0->Add(HwwMakeNtupleMod1);
      HwwMakeNtupleMod1->Add(HwwMakeNtupleMod2);
    }
  }
  
  //------------------------------------------------------------------------------------------------
  //
  // setup analysis object
  //
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kFALSE);
  if(nEvents >= 0) 
    ana->SetProcessNEvents(nEvents);

  ana->AddSuperModule(GeneratorMod1);
  ana->SetPrintScale(100);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n\n",book,dataset,skim,fileset);
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset);
  else 
    d = c->FindDataset(book,skimdataset.Data(),fileset);
  ana->AddDataset(d);

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(0);

  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}

