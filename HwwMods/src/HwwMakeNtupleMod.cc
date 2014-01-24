// $Id: HwwMakeNtupleMod.cc,v 1.142 2013/11/19 17:32:25 ceballos Exp $

#include "MitHiggs/HwwMods/interface/HwwMakeNtupleMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TTree.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/CompositeParticleCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/GeneratorTools.h"
#include "MitAna/DataTree/interface/GenericParticleCol.h"

using namespace mithep;
ClassImp(mithep::HwwMakeNtupleMod)

//--------------------------------------------------------------------------------------------------
HwwMakeNtupleMod::HwwMakeNtupleMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fFakeRatePredictionType(0),
  fFillNtupleType(0),
  fPtJetCut(30.0),
  fEtaJetCut(3.0),
  fPFTauName("random"),
  fPFMetName("PFMet"),
  fPFMetStd(0),
  fMuonName(Names::gkMuonBrn),
  fElectronName(Names::gkElectronBrn),
  fTrackName(Names::gkTrackBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFCandidates(0),
  fCleanJetsNoPtCutName("random0"),
  fMCqqHsName(ModNames::gkMCqqHsName),
  fEvtHdrName(Names::gkEvtHeaderBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fMCEvInfoName(Names::gkMCEvtInfoBrn),
  fMCEventInfo(0),
  fCaloJetName0("AKt5Jets"),
  fPFJetName0("AKt5PFJets"),
  fCaloJet0(0),
  fPFJet0(0),
  fMuons(0),
  fElectrons(0),
  fEventHeader(0),
  fVertices(0),
  fPileupInfos(0),
  fParticles(0),
  fDecay(0),
  fJetScaleSyst(0.0),
  fIsData(kFALSE),
  fFillPhotonTemplate(kFALSE),
  fDoPileupReweighting(kFALSE),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPileupEnergyDensity(0),
  fOutputFile(0),
  fOutputName("ntuple.root"),
  fMuonFakeName("random"),
  fElectronFakeName("random"),
  fLeptonFakeName("random"),
  fIntRadius(0.0),
  fIs42x(kFALSE),
  fElectronIDMVA(0),
  fMuonTools(0),
  fMuonIDMVA(0),
  fCorrectedJetsName("dummy"),
  fMVAElVersion(0),
  fMVAMuVersion(0),
  fTheRhoType(RhoUtilities::CMS_RHO_RHOKT6PFJETS),
  fAddLheWeights(kFALSE),
  fLheWeightsName("LheWeights"),
  fLheWeights(0),
  fNEventsSelected(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
HwwMakeNtupleMod::~HwwMakeNtupleMod() 
{
  // Destructor.
}

//--------------------------------------------------------------------------------------------------
void HwwMakeNtupleMod::Begin()
{
  // Run startup code on the client machine. For this module, we don't do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void HwwMakeNtupleMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  LoadBranch(fEvtHdrName);
  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fCaloJetName0);
  LoadBranch(fPFJetName0);
  LoadBranch(fTrackName);
  LoadBranch(fPFCandidatesName);
  if (!fIsData) LoadBranch(Names::gkPileupInfoBrn);
  if (!fIsData) LoadBranch(Names::gkMCPartBrn);
  LoadEventObject(fPileupEnergyDensityName, fPileupEnergyDensity);
  LoadBranch(fTrackName);
  if(fAddLheWeights) LoadBranch(fLheWeightsName);

  //************************************************************************************************
  //Get NNLO Weight
  //************************************************************************************************
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight",0);

  //************************************************************************************************
  //Obtain all the good objects from the event cleaning module
  //************************************************************************************************
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  ElectronOArr *OriginalCleanElectrons  = GetObjThisEvt<ElectronOArr>(ModNames::gkCleanElectronsName);
  MuonOArr  *OriginalCleanMuons         = GetObjThisEvt<MuonOArr>(ModNames::gkCleanMuonsName);
  PhotonOArr  *OriginalCleanPhotons     = GetObjThisEvt<PhotonOArr>(ModNames::gkCleanPhotonsName);
  PFTauOArr  *CleanTaus                 = GetObjThisEvt<PFTauOArr>(fPFTauName);

  JetOArr *OriginalCleanJetsNoPtCut     = GetObjThisEvt<JetOArr>(fCleanJetsNoPtCutName);
  ParticleOArr *OriginalLeptons         = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
  LoadBranch(fPFMetName);
  const PFMet *pfMetStatic              = fPFMetStd->At(0);
  MetOArr *GenMet                       = GetObjThisEvt<MetOArr>(ModNames::gkMCMETName);

  //************************************************************************************************
  //Obtain MC Particle Collections
  //************************************************************************************************
  MCParticleOArr *GenLeptons   = 0;
  MCParticleOArr *GenTaus      = 0;
  MCParticleOArr *GenPhotons   = 0;
  ObjArray<MCParticle> *GenLeptonsAndTaus = new ObjArray<MCParticle>;
  MCParticleCol *GenBosons   = 0;    
  if(fIsData == kFALSE){
    GenLeptons   = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCLeptonsName);
    GenTaus      = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCTausName);
    GenPhotons   = GetObjThisEvt<MCParticleOArr>(ModNames::gkMCPhotonsName);
    GenLeptonsAndTaus = new ObjArray<MCParticle>;
    for (UInt_t i=0; i<GenLeptons->GetEntries(); i++)
      GenLeptonsAndTaus->Add(GenLeptons->At(i));
    for (UInt_t i=0; i<GenTaus->GetEntries(); i++)
      GenLeptonsAndTaus->Add(GenTaus->At(i));
    GenBosons = GetObjThisEvt<MCParticleCol>(ModNames::gkMCBosonsName);    
  }

  //***********************************************************************************************
  //Import Fakeable object Collections
  //***********************************************************************************************
  MuonOArr  *CleanMuonsFakeable        = GetObjThisEvt<MuonOArr>(fMuonFakeName);
  ElectronOArr *CleanElectronsFakeable = GetObjThisEvt<ElectronOArr>(fElectronFakeName);
  ParticleOArr *leptonsFakeable        = GetObjThisEvt<ParticleOArr>(fLeptonFakeName);

  ParticleOArr *leptonsOnlyFake = new ObjArray<Particle>;
  MuonOArr *CleanMuonsOnlyFake = new ObjArray<Muon>;
  ElectronOArr *CleanElectronsOnlyFake = new ObjArray<Electron>;
  if(fFakeRatePredictionType == 1 || fFakeRatePredictionType == 2){
    for(UInt_t i=0; i<leptonsFakeable->GetEntries(); i++){
      Bool_t isOnlyFake = kTRUE;
      for(UInt_t j=0; j<OriginalLeptons->GetEntries(); j++) {
	if(OriginalLeptons->At(j) == leptonsFakeable->At(i)) {
  	  isOnlyFake = kFALSE;
          break;
	}
      }
      if(isOnlyFake == kTRUE) leptonsOnlyFake->Add(leptonsFakeable->At(i));
    }
    leptonsOnlyFake->Sort();

    for(UInt_t i=0; i<CleanMuonsFakeable->GetEntries(); i++){
      Bool_t isOnlyFake = kTRUE;
      for(UInt_t j=0; j<OriginalCleanMuons->GetEntries(); j++) {
	if(OriginalCleanMuons->At(j) == CleanMuonsFakeable->At(i)) {
  	  isOnlyFake = kFALSE;
          break;
	}
      }
      if(isOnlyFake == kTRUE) CleanMuonsOnlyFake->Add(CleanMuonsFakeable->At(i));
    }
    CleanMuonsOnlyFake->Sort();

    for(UInt_t i=0; i<CleanElectronsFakeable->GetEntries(); i++){
      Bool_t isOnlyFake = kTRUE;
      for(UInt_t j=0; j<OriginalCleanElectrons->GetEntries(); j++) {
	if(OriginalCleanElectrons->At(j) == CleanElectronsFakeable->At(i)) {
  	  isOnlyFake = kFALSE;
          break;
	}
      }
      if(isOnlyFake == kTRUE) CleanElectronsOnlyFake->Add(CleanElectronsFakeable->At(i));
    }
    CleanElectronsOnlyFake->Sort();
  }

  //***********************************************************************************************
  //-----------------------------------------------------------------------------------------------
  //MAIN LOOP OVER FAKE EVENT HEADERS
  //-----------------------------------------------------------------------------------------------
  UInt_t nLoop = 1;
  if(fFakeRatePredictionType == 1) nLoop = leptonsOnlyFake->GetEntries();
  //***********************************************************************************************
  for (UInt_t in=0; in<nLoop; in++) {
 
    //********************************************************************************************
    //make lepton collections
    //********************************************************************************************

    ObjArray<Particle> *leptons = new ObjArray<Particle>;
    MuonOArr *CleanMuons = new ObjArray<Muon>;
    ElectronOArr *CleanElectrons = new ObjArray<Electron>;
    if(fFakeRatePredictionType == 0 || fFakeRatePredictionType == 1 || fFakeRatePredictionType == 3){
      for (UInt_t j=0;j<OriginalLeptons->GetEntries() ; j++) {
        leptons->Add(OriginalLeptons->At(j));
        leptons->At(leptons->GetEntries()-1)->SetIsFakeable(kFALSE);
      }

      for (UInt_t j=0;j<OriginalCleanMuons->GetEntries() ; j++) {
        CleanMuons->Add(OriginalCleanMuons->At(j));
      }

      for (UInt_t j=0;j<OriginalCleanElectrons->GetEntries() ; j++) {
        CleanElectrons->Add(OriginalCleanElectrons->At(j));
      }
    }

    if (fFakeRatePredictionType == 1) {
      leptons->Add(leptonsOnlyFake->At(in));
      leptons->At(leptons->GetEntries()-1)->SetIsFakeable(kTRUE);
      if (leptonsOnlyFake->At(in)->ObjType() == kMuon) {
        for (UInt_t j=0;j<CleanMuonsOnlyFake->GetEntries() ; j++) {
          if (MathUtils::DeltaR(leptonsOnlyFake->At(in)->Mom(), CleanMuonsOnlyFake->At(j)->Mom()) < 0.05) {
            CleanMuons->Add(CleanMuonsOnlyFake->At(j));
            break;
          }
        }
      }
      if (leptonsOnlyFake->At(in)->ObjType() == kElectron) {
        for (UInt_t j=0;j<CleanElectronsOnlyFake->GetEntries() ; j++) {
          if (MathUtils::DeltaR(leptonsOnlyFake->At(in)->Mom(), CleanElectronsOnlyFake->At(j)->Mom()) < 0.05) {
            CleanElectrons->Add(CleanElectronsOnlyFake->At(j));
            break;
          }
        }
      }
    }

    if (fFakeRatePredictionType == 2) {
      for (UInt_t j=0;j<leptonsOnlyFake->GetEntries() ; j++) {
        leptons->Add(leptonsOnlyFake->At(j));
        leptons->At(leptons->GetEntries()-1)->SetIsFakeable(kTRUE);
      	if (leptonsOnlyFake->At(j)->ObjType() == kMuon) {
      	  for (UInt_t k=0;k<CleanMuonsOnlyFake->GetEntries() ; k++) {
      	    if (MathUtils::DeltaR(leptonsOnlyFake->At(j)->Mom(), CleanMuonsOnlyFake->At(k)->Mom()) < 0.05) {
      	      CleanMuons->Add(CleanMuonsOnlyFake->At(k));
      	      break;
      	    }
      	  }
      	}
      	if (leptonsOnlyFake->At(j)->ObjType() == kElectron) {
      	  for (UInt_t k=0;k<CleanElectronsOnlyFake->GetEntries() ; k++) {
      	    if (MathUtils::DeltaR(leptonsOnlyFake->At(j)->Mom(), CleanElectronsOnlyFake->At(k)->Mom()) < 0.05) {
      	      CleanElectrons->Add(CleanElectronsOnlyFake->At(k));
      	      break;
      	    }
      	  }
      	}
      }
    }

    if(fFakeRatePredictionType == 3 && OriginalCleanPhotons->GetEntries() >= 1){
      leptons->Add(OriginalCleanPhotons->At(0));
      leptons->At(leptons->GetEntries()-1)->SetIsFakeable(kTRUE);
    }

    leptons->Sort();
    CleanMuons->Sort();
    CleanElectrons->Sort();

    //********************************************************************************************
    //make jets collection
    //********************************************************************************************
    ObjArray<Jet> *CleanJetsNoPtCut = new ObjArray<Jet>;
    if (fFakeRatePredictionType == 0 ) {
      for (UInt_t j=0;j<OriginalCleanJetsNoPtCut->GetEntries() ; j++) {
        if(fFillPhotonTemplate == kFALSE) CleanJetsNoPtCut->Add(OriginalCleanJetsNoPtCut->At(j));
        else if(OriginalCleanPhotons->GetEntries() >= 1) {
	  if(MathUtils::DeltaR(OriginalCleanPhotons->At(0)->Mom(), OriginalCleanJetsNoPtCut->At(j)->Mom()) > 0.3) {
	    CleanJetsNoPtCut->Add(OriginalCleanJetsNoPtCut->At(j));
          }
        }
      }
    }
    else if (fFakeRatePredictionType == 1) {
      for (UInt_t j=0;j<OriginalCleanJetsNoPtCut->GetEntries() ; j++) {
        if(MathUtils::DeltaR(leptonsOnlyFake->At(in)->Mom(), OriginalCleanJetsNoPtCut->At(j)->Mom()) > 0.3) {
	  CleanJetsNoPtCut->Add(OriginalCleanJetsNoPtCut->At(j));
        }
      }
    }
    else if (fFakeRatePredictionType == 2) {
      for (UInt_t k=0;k<leptons->GetEntries() ; k++) {
      	for (UInt_t j=0;j<OriginalCleanJetsNoPtCut->GetEntries() ; j++) {
	  Bool_t isLepton = kFALSE;
      	  if(MathUtils::DeltaR(leptons->At(k)->Mom(), OriginalCleanJetsNoPtCut->At(j)->Mom()) < 0.3) {
	    isLepton = kTRUE;
	    break;
      	  }
	  if(isLepton == kFALSE) CleanJetsNoPtCut->Add(OriginalCleanJetsNoPtCut->At(j));
      	}
      }
    }
    else if (fFakeRatePredictionType == 3 && OriginalCleanPhotons->GetEntries() >= 1) {
      for (UInt_t j=0;j<OriginalCleanJetsNoPtCut->GetEntries() ; j++) {
        if(MathUtils::DeltaR(OriginalCleanPhotons->At(0)->Mom(), OriginalCleanJetsNoPtCut->At(j)->Mom()) > 0.3) {
	  CleanJetsNoPtCut->Add(OriginalCleanJetsNoPtCut->At(j));
        }
      }
    }
    CleanJetsNoPtCut->Sort(); 

    if(leptons->GetEntries() != (CleanMuons->GetEntries() + CleanElectrons->GetEntries()) && fFakeRatePredictionType != 3) assert(1);

    //*********************************************************************************************
    //Initialize
    //*********************************************************************************************
    double deltaPhiLLJet = 999.;
    std::vector<double> leptonsDz;
    std::vector<double> leptonsMVA;
    std::vector<double> lepDetEta;
    std::vector<int> leptonsMVAPass;
    std::vector<bool> leptonsAgreeQ;
    double zDiffMax = 0.0;
    bool hasZCand = kFALSE;

    ObjArray<Muon> *DirtyMuons = new ObjArray<Muon>;
    for (UInt_t i=0; i<fMuons->GetEntries(); i++) {
      const Muon *mu = fMuons->At(i);
      for (UInt_t nj=i+1; nj<fMuons->GetEntries(); nj++) {
        const Muon *mu2 = fMuons->At(nj);
        if(mu->Charge() != mu2->Charge()){
	  CompositeParticle looseDilepton;
	  looseDilepton.AddDaughter(mu);
	  looseDilepton.AddDaughter(mu2);
	  double zDiff = TMath::Abs(mu ->BestTrk()->DzCorrected(*fVertices->At(0))-
	                            mu2->BestTrk()->DzCorrected(*fVertices->At(0)));
	  if(mu->Pt() > 10 && mu2->Pt() > 10 &&  zDiff < 0.1 &&
	     TMath::Abs(looseDilepton.Mass()-91.1876) < 15.0) hasZCand = kTRUE;
	}
      }
      bool isCleanMuon = kFALSE;
      for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
        if(fMuons->At(i) == CleanMuons->At(j) &&
           CleanMuons->At(j)->Pt() > 10) isCleanMuon = kTRUE;
      }
      if(isCleanMuon == kTRUE) continue;

      if(MuonTools::PassSoftMuonCut(mu, fVertices, 0.2)) DirtyMuons->Add(mu);
    }
    DirtyMuons->Sort();

    for (UInt_t i=0; i<fElectrons->GetEntries(); i++) {
      const Electron *el = fElectrons->At(i);
      for (UInt_t nj=i+1; nj<fElectrons->GetEntries(); nj++) {
        const Electron *el2 = fElectrons->At(nj);
        if(el->Charge() != el2->Charge()){
	  CompositeParticle looseDilepton;
	  looseDilepton.AddDaughter(el);
	  looseDilepton.AddDaughter(el2);
	  double zDiff = TMath::Abs(el ->GsfTrk()->DzCorrected(*fVertices->At(0))-
	                            el2->GsfTrk()->DzCorrected(*fVertices->At(0)));
	  if(el->Pt() > 10 && el2->Pt() > 10 && zDiff < 0.1 &&
	     TMath::Abs(looseDilepton.Mass()-91.1876) < 15.0) hasZCand = kTRUE;
	}
      }
    }
    // Make lepton vector from muons and electrons
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
      double pDz = CleanMuons->At(j)->BestTrk()->DzCorrected(*fVertices->At(0));
      leptonsDz.push_back(pDz);
    }

    for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
      double pDz = CleanElectrons->At(j)->GsfTrk()->DzCorrected(*fVertices->At(0));
      leptonsDz.push_back(pDz);
    }

    // New MET
    PFMet* pfMet = new PFMet();

    // Preselection requirements
    if ( (!fFillPhotonTemplate && fFakeRatePredictionType != 3 && leptons->GetEntries() >= 2 && leptons->At(0)->Pt() > 20 &&
	   (
	    (leptons->GetEntries() == 2)
	    ||
	    (leptons->GetEntries() == 3)
	    )
	 )
	   || 
	 (!fFillPhotonTemplate && fFakeRatePredictionType == 3 && leptons->GetEntries() >= 2 && leptons->At(0)->Pt() > 20 && OriginalCleanPhotons->GetEntries() >= 1 && 
	   (
	    (leptons->GetEntries() == 2 && fIsData == kFALSE)
	    ||
	    (leptons->GetEntries() == 3)
	    )
	 )
	   ||
	 (fFillPhotonTemplate && (OriginalCleanPhotons->GetEntries() == 1 || (OriginalCleanPhotons->GetEntries() >= 1 && OriginalCleanPhotons->At(1)->Pt() < 30))
	  && CleanMuons->GetEntries() == 0 && CleanElectrons->GetEntries() == 0)
	 ) {	

      // *****************
      // Fill electron MVA
      // *****************
      for (UInt_t j=0;j<leptons->GetEntries() ; j++) {
  	Double_t MVAValue = -3.0; Int_t MVAPass = 0; Double_t DetEta = -900.0;
	leptonsAgreeQ.push_back(kTRUE);
      	if     (leptons->At(j)->ObjType() == kElectron) {
      	  for (UInt_t k=0;k<CleanElectrons->GetEntries() ; k++) {
      	    if (leptons->At(j) == CleanElectrons->At(k)) {
	      // Check if Q agree
	      if(CleanElectrons->At(k)->GsfTrk()           && CleanElectrons->At(k)->TrackerTrk() &&
	         CleanElectrons->At(k)->GsfTrk()->Charge() == CleanElectrons->At(k)->TrackerTrk()->Charge() &&
		 CleanElectrons->At(k)->GsfTrk()->Charge() == CleanElectrons->At(k)->ScPixCharge()){
	        // triple charge electron requirement
	      } else {
	        leptonsAgreeQ[j] = kFALSE;
	      }
	      DetEta =  CleanElectrons->At(k)->SCluster()->Eta();
  	      Bool_t idcut = ElectronTools::PassCustomID(CleanElectrons->At(k), ElectronTools::kVBTFWorkingPointFakeableId);
  	      if(idcut == kTRUE) {
	        if(fMVAElVersion == 0) MVAValue = fElectronIDMVA->MVAValue(CleanElectrons->At(k), fVertices->At(0), fPFCandidates, fPileupEnergyDensity, fIntRadius);
  		else {
                  ElectronOArr *tempElectrons = new  ElectronOArr;
                  MuonOArr     *tempMuons     = new  MuonOArr;
                  MVAValue = fElectronIDMVA->MVAValue(CleanElectrons->At(k), fVertices->At(0), fPFCandidates, fPileupEnergyDensity, ElectronTools::kEleEANoCorr, tempElectrons, tempMuons, kFALSE);
                  delete tempElectrons;
                  delete tempMuons;
                }
		Int_t subdet = 0;
  		if (CleanElectrons->At(k)->SCluster()->AbsEta() < 1.0) subdet = 0;
  		else if (CleanElectrons->At(k)->SCluster()->AbsEta() < 1.479) subdet = 1;
  		else subdet = 2;
  		Int_t ptBin = 0;
  		if (CleanElectrons->At(k)->Pt() > 20.0) ptBin = 1;
                Double_t MVACut = -999;
  		if      (subdet == 0 && ptBin == 0) MVACut = 0.123/2.0;
  		else if (subdet == 1 && ptBin == 0) MVACut = 0.219/2.0;
  		else if (subdet == 2 && ptBin == 0) MVACut = 0.509/2.0;
  		else if (subdet == 0 && ptBin == 1) MVACut = 0.935/2.0;
  		else if (subdet == 1 && ptBin == 1) MVACut = 0.889/2.0;
  		else if (subdet == 2 && ptBin == 1) MVACut = 0.871/2.0;
		if (MVAValue > MVACut) MVAPass = 1;
  	      }
	      else {
	        MVAValue = -2.0; MVAPass = 0;
	      }
      	      break;
      	    }
      	  }
      	}
      	else if(leptons->At(j)->ObjType() == kMuon) {
      	  for (UInt_t k=0;k<CleanMuons->GetEntries() ; k++) {
      	    if (leptons->At(j) == CleanMuons->At(k)) {
	      Bool_t idcut = (CleanMuons->At(k)->BestTrk() != 0 &&
                              CleanMuons->At(k)->BestTrk()->NHits() > 10 &&
                              CleanMuons->At(k)->BestTrk()->NPixelHits() > 0 &&
                              CleanMuons->At(k)->BestTrk()->PtErr()/CleanMuons->At(k)->BestTrk()->Pt() < 0.1 &&
                              MuonTools::PassD0Cut(CleanMuons->At(k), fVertices, 0.20, 0) &&
                              MuonTools::PassDZCut(CleanMuons->At(k), fVertices, 0.10, 0) &&
                              CleanMuons->At(k)->TrkKink() < 20.0);
  	      if(idcut == kTRUE) {
                if(fMVAMuVersion == 0) MVAValue = fMuonIDMVA->MVAValue(CleanMuons->At(k), fVertices->At(0), fMuonTools, fPFCandidates, fPileupEnergyDensity, kFALSE);
  		else {
                  ElectronOArr *tempElectrons = new  ElectronOArr;
                  MuonOArr     *tempMuons     = new  MuonOArr;
                  MVAValue = fMuonIDMVA->MVAValue(CleanMuons->At(k), fVertices->At(0), fMuonTools, fPFCandidates,
                                                  fPileupEnergyDensity, MuonTools::kMuEAFall11MC, tempElectrons, tempMuons,kFALSE);
                  delete tempElectrons;
                  delete tempMuons;
		}
		const Track *muTrk=0;
  		if(CleanMuons->At(k)->HasTrackerTrk()) 	{ muTrk = CleanMuons->At(k)->TrackerTrk();    }
   		else if(CleanMuons->At(k)->HasStandaloneTrk()) { muTrk = CleanMuons->At(k)->StandaloneTrk(); }
  		Int_t subdet = 0;
  		if (fabs(muTrk->Eta()) < 1.479) subdet = 0;
  		else subdet = 1;
  		Int_t ptBin = 0;
  		if (muTrk->Pt() > 14.5) ptBin = 1;
  		if (muTrk->Pt() > 20.0) ptBin = 2;

  		Double_t MVACut = -999;
  		if	(subdet == 0 && ptBin == 0) MVACut = -0.8006;
  		else if (subdet == 1 && ptBin == 0) MVACut = -0.6698;
  		else if (subdet == 0 && ptBin == 1) MVACut = -0.7658;
  		else if (subdet == 1 && ptBin == 1) MVACut = -0.6406;
  		else if (subdet == 0 && ptBin == 2) MVACut = -0.4862;
  		else if (subdet == 1 && ptBin == 2) MVACut =  0.6538;
		if (MVAValue > MVACut) MVAPass = 1;
  	      }
	      else {
	        MVAValue = -2.0; MVAPass = 0;
	      }
      	      break;
      	    }
      	  }
      	}
	leptonsMVA.push_back(MVAValue);
	leptonsMVAPass.push_back(MVAPass);
	lepDetEta.push_back(DetEta);
      }

      for(UInt_t t=0; t<leptonsDz.size(); t++) {
	for(UInt_t i=t+1; i<leptonsDz.size(); i++) {
	  if(TMath::Abs(leptonsDz[t]-leptonsDz[i]) > zDiffMax) zDiffMax = TMath::Abs(leptonsDz[t]-leptonsDz[i]);
	}
      }
      leptonsDz.clear();

      double photonV[4] = {0, 0, 0, 0};
      if(fFakeRatePredictionType == 3 && OriginalCleanPhotons->GetEntries() >= 1){
        photonV[0] = OriginalCleanPhotons->At(0)->Px(); photonV[1] = OriginalCleanPhotons->At(0)->Py();
        photonV[2] = OriginalCleanPhotons->At(0)->Pz(); photonV[3] = OriginalCleanPhotons->At(0)->E();
      }
      GenericParticle *thePhoton = new GenericParticle(photonV[0],photonV[1],photonV[2],photonV[3]);
      MetTools metTools(CleanMuons, CleanElectrons, fPFCandidates, fVertices->At(0), 0.1, 8.0, 5.0, fIntRadius,thePhoton);
      delete thePhoton;
 
      if (fFillPhotonTemplate ) {
        metTools.AddToCorrectedTrackMet(OriginalCleanPhotons->At(0));
	metTools.RemoveParticleInIsoConeFromTrackMet(OriginalCleanPhotons->At(0),fPFCandidates, fVertices->At(0), 0.1, 0.4);
	metTools.RemoveParticleInIsoConeFromCorrectedMet(OriginalCleanPhotons->At(0),fPFCandidates, fVertices->At(0), 0.1, 1.0, 5.0, 0.4);
	metTools.RemoveParticleInIsoConeFromRecoil(OriginalCleanPhotons->At(0),fPFCandidates, fVertices->At(0), 0.1, 1.0, 5.0, 0.4);        
      }

      PFCandidateCol *pFNoPileUpCands = GetObjThisEvt<PFCandidateCol>("PFNoPileUp");
      MetTools metToolsNoPu(CleanMuons, CleanElectrons, pFNoPileUpCands, fVertices->At(0), 0.1, 8.0, 5.0, fIntRadius);

      double pMET[2] = {metTools.GetProjectedMet(leptons,pfMetStatic),
        		metTools.GetProjectedTrackMet(leptons)};
      pfMet->SetMex(pfMetStatic->Px());
      pfMet->SetMey(pfMetStatic->Py());
      pfMet->SetSumEt(pfMetStatic->SumEt());
      pfMet->SetPFMetSig(pfMetStatic->PFMetSig());

      UInt_t leptonGenType[3]       = {0, 0, 0};
      UInt_t leptonMotherGenType[3] = {0, 0, 0};
      UInt_t leptonGenIsSeen[3]     = {0, 0, 0}; // consider only the first three gen leptons for now
      // look for real W/Z -> leptons
      for(UInt_t i=0; i<leptons->GetEntries(); i++) {

        MCParticle *gen = 0;
        //Match to Leptons
	if(GenLeptons){
          for (UInt_t j=0; j<GenLeptons->GetEntries(); j++) {
            if(MathUtils::DeltaR(GenLeptons->At(j)->Mom(), leptons->At(i)->Mom()) < 0.10) {
              gen = GenLeptons->At(j);
	      if(j < 3) leptonGenIsSeen[j] = 1;
              break;
            }
          }
	  if     (leptonGenIsSeen[0] == 0 && GenLeptons->GetEntries() >= 1){
            fSmurfTree.auxVar0_ = GenLeptons->At(0)->Pt();
            //fSmurfTree.auxVar1_ = GenLeptons->At(0)->Eta();
          }
	  else if(leptonGenIsSeen[1] == 0 && GenLeptons->GetEntries() >= 2){
            fSmurfTree.auxVar0_ = GenLeptons->At(1)->Pt();
            //fSmurfTree.auxVar1_ = GenLeptons->At(1)->Eta();
          }
	  else if(leptonGenIsSeen[2] == 0 && GenLeptons->GetEntries() >= 3){
            fSmurfTree.auxVar0_ = GenLeptons->At(2)->Pt();
            //fSmurfTree.auxVar1_ = GenLeptons->At(2)->Eta();
          }
	  else {
            fSmurfTree.auxVar0_ = 0.0;
            //fSmurfTree.auxVar1_ = 0.0;
          }
	}

        //Match to Hadronic Taus
        if (!gen && GenTaus) {
          for (UInt_t l=0; l<GenTaus->GetEntries(); ++l) {
            if (MathUtils::DeltaR(GenTaus->At(l)->Mom(), leptons->At(i)->Mom()) < 0.3 ) {
              gen = GenTaus->At(l);
              break;
            }
          }
        }

        //Match to Prompt Photons
        if (!gen && GenPhotons) {
          for (UInt_t l=0; l < GenPhotons->GetEntries(); l++) {  
            if (GenPhotons->At(l)->Pt() > 10.0 &&  MathUtils::DeltaR(GenPhotons->At(l)->Mom(), leptons->At(i)->Mom()) < 0.3 ) {
        
              //ISR Photon
              if ( (GenPhotons->At(l)->Mother() && GenPhotons->At(l)->Mother()->IsParton())
                   || (GenPhotons->At(l)->Mother()->AbsPdgId() == 22 
                       && GenPhotons->At(l)->Mother()->Status() ==3  
                       && GenPhotons->At(l)->Mother()->Mother() 
                       && GenPhotons->At(l)->Mother()->Mother()->IsParton())
                ) {
                gen = GenPhotons->At(l);
                break;
              }
        
              //WWgamma vertex
              if ((GenPhotons->At(l)->Mother() && GenPhotons->At(l)->Mother()->AbsPdgId() == 24) 
                  || 
                  (GenPhotons->At(l)->Mother()->AbsPdgId() == 22 
                   && GenPhotons->At(l)->Mother()->Status() == 3
                   && GenPhotons->At(l)->Mother()->Mother()
                   && GenPhotons->At(l)->Mother()->Mother()->AbsPdgId() == 24
                    )
                ) {
                gen = GenPhotons->At(l);
                break;
              }

              //Pythia FSR
              if (GenPhotons->At(l)->Mother() && (GenPhotons->At(l)->Mother()->Status() == 3 || GenPhotons->At(l)->Mother()->Status() == 2)
                  && (GenPhotons->At(l)->Mother()->AbsPdgId() == 11 
                      || GenPhotons->At(l)->Mother()->AbsPdgId() == 13
                      || GenPhotons->At(l)->Mother()->AbsPdgId() == 15)
                ) {
                CompositeParticle *object = new CompositeParticle();
                object->AddDaughter(GenPhotons->At(l));
                object->AddDaughter(GenPhotons->At(l)->Mother());
                if(object->Mass() > 1.0) {
                  gen = GenPhotons->At(l);
                  delete object;          
                  break;
                }
                delete object;
              }
            }
          }
        }
        
        if(gen) {
          leptonGenType[i] =  gen->PdgId();
          if     (gen->HasMother(MCParticle::kW) && gen->Charge() > 0) leptonMotherGenType[i] =  24;
          else if(gen->HasMother(MCParticle::kW) && gen->Charge() < 0) leptonMotherGenType[i] = -24;
          else if(gen->HasMother(MCParticle::kZ))                      leptonMotherGenType[i] = 23;
          else                                                         leptonMotherGenType[i] = gen->DistinctMother()->PdgId();
        } else {
          leptonGenType[i] =   9;
        }
      }

      Int_t typeGenFid = 0;
      if(GenLeptons && GenLeptons->GetEntries() >= 2){
        if(GenLeptons->GetEntries() == 2) {
	  typeGenFid = 1;
	  if(GenLeptons->At(0)->Pt() > 20 && GenLeptons->At(0)->AbsEta() < 2.5 && 
             GenLeptons->At(1)->Pt() > 10 && GenLeptons->At(1)->AbsEta() < 2.5) {
	    typeGenFid = 2;
	    if(GenLeptons->At(1)->Pt() > 20) {
	      typeGenFid = 3;
	    }
	  }
	}
        if(GenLeptons->GetEntries() == 3) {
	  typeGenFid = 4;
	  int type3LGen[3] = {0,0,0};
	  for(unsigned int n3l=0; n3l<GenLeptons->GetEntries(); n3l++) {
	    if(GenLeptons->At(n3l)->Pt() > 20 && GenLeptons->At(n3l)->AbsEta() < 2.5) type3LGen[0]++;
	    if(GenLeptons->At(n3l)->Pt() > 10 && GenLeptons->At(n3l)->AbsEta() < 2.5) type3LGen[1]++;
	    if(GenLeptons->At(n3l)->Pt() > 10 && GenLeptons->At(n3l)->AbsEta() < 4.7) type3LGen[2]++;
	  }
	  if     (type3LGen[0] == 0 && type3LGen[1] >= 0 && type3LGen[2] >= 0) typeGenFid = 4;
	  else if(type3LGen[0] == 1 && type3LGen[1] == 1 && type3LGen[2] >= 0) typeGenFid = 4;
	  else if(type3LGen[0] == 1 && type3LGen[1] == 2 && type3LGen[2] == 2) typeGenFid = 5.0;
	  else if(type3LGen[0] == 1 && type3LGen[1] == 2 && type3LGen[2] == 3) typeGenFid = 5.1;
	  else if(type3LGen[0] == 1 && type3LGen[1] == 3 && type3LGen[2] == 3) typeGenFid = 6;
	  else if(type3LGen[0] == 2 && type3LGen[1] == 2 && type3LGen[2] == 2) typeGenFid = 7.0;
	  else if(type3LGen[0] == 2 && type3LGen[1] == 2 && type3LGen[2] == 3) typeGenFid = 7.1;
	  else if(type3LGen[0] == 2 && type3LGen[1] == 3 && type3LGen[2] == 3) typeGenFid = 8;
	  else if(type3LGen[0] == 3 && type3LGen[1] == 3 && type3LGen[2] == 3) typeGenFid = 8;
	  else {printf("Impossible typeGenFid: %d %d %d\n",type3LGen[0],type3LGen[1],type3LGen[2]); assert(0);}
	}
        if(GenLeptons->GetEntries() >= 4) typeGenFid = 9;
      }

      CompositeParticle *dilepton = new CompositeParticle();

      // Sort and count the number of central Jets for vetoing
      vector<Jet*> sortedJetsAll;
      vector<Jet*> sortedJets;
      double nMaxMediumBTagJet = -25.0;
      double nMaxMediumOtherBTagJet[6] = {-25.0,-25.0,-25.0,-25.0,-25.0,-25.0};
      double dZAverageJet1Pt = 0.0;
      for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); i++){
        if(CleanJetsNoPtCut->At(i)->RawMom().Pt() <= 7) continue;
        if(TMath::Abs(CleanJetsNoPtCut->At(i)->RawMom().Eta()) >= fEtaJetCut) continue;
        Jet* jet_a = new Jet(CleanJetsNoPtCut->At(i)->Px()*(1.0+fJetScaleSyst),
        		     CleanJetsNoPtCut->At(i)->Py()*(1.0+fJetScaleSyst),
        		     CleanJetsNoPtCut->At(i)->Pz()*(1.0+fJetScaleSyst),
        		     CleanJetsNoPtCut->At(i)->E() *(1.0+fJetScaleSyst));

	int nCloseStdJet = -1;
	if(fIntRadius > 0){
	  double deltaRMin = 999.;
	  for(UInt_t nj=0; nj<fCaloJet0->GetEntries(); nj++){
	    const CaloJet *jet = fCaloJet0->At(nj);
	    Double_t deltaR = MathUtils::DeltaR(jet_a->Mom(),jet->Mom());
	    if(deltaR < deltaRMin) {
   	      nCloseStdJet = nj;
   	      deltaRMin = deltaR;
	    }
	  }
	}
        if(nCloseStdJet >= 0 && fCaloJet0->At(nCloseStdJet)->TrackCountingHighEffBJetTagsDisc() >=
	                            CleanJetsNoPtCut->At(i)->TrackCountingHighEffBJetTagsDisc()){
          jet_a->SetMatchedMCFlavor(fCaloJet0->At(nCloseStdJet)->MatchedMCFlavor());
          jet_a->SetCombinedSecondaryVertexBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->CombinedSecondaryVertexBJetTagsDisc());
          jet_a->SetCombinedSecondaryVertexMVABJetTagsDisc(fCaloJet0->At(nCloseStdJet)->CombinedSecondaryVertexMVABJetTagsDisc());
          jet_a->SetJetProbabilityBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->JetProbabilityBJetTagsDisc());
          jet_a->SetJetBProbabilityBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->JetBProbabilityBJetTagsDisc());
          jet_a->SetTrackCountingHighEffBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->TrackCountingHighEffBJetTagsDisc());
          jet_a->SetTrackCountingHighPurBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->TrackCountingHighPurBJetTagsDisc());
          jet_a->SetSimpleSecondaryVertexBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->SimpleSecondaryVertexBJetTagsDisc());
          jet_a->SetSimpleSecondaryVertexHighEffBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->SimpleSecondaryVertexHighEffBJetTagsDisc());
          jet_a->SetSimpleSecondaryVertexHighPurBJetTagsDisc(fCaloJet0->At(nCloseStdJet)->SimpleSecondaryVertexHighPurBJetTagsDisc());
	} else {
          jet_a->SetMatchedMCFlavor(CleanJetsNoPtCut->At(i)->MatchedMCFlavor());
          jet_a->SetCombinedSecondaryVertexBJetTagsDisc(CleanJetsNoPtCut->At(i)->CombinedSecondaryVertexBJetTagsDisc());
          jet_a->SetCombinedSecondaryVertexMVABJetTagsDisc(CleanJetsNoPtCut->At(i)->CombinedSecondaryVertexMVABJetTagsDisc());
          jet_a->SetJetProbabilityBJetTagsDisc(CleanJetsNoPtCut->At(i)->JetProbabilityBJetTagsDisc());
          jet_a->SetJetBProbabilityBJetTagsDisc(CleanJetsNoPtCut->At(i)->JetBProbabilityBJetTagsDisc());
          jet_a->SetTrackCountingHighEffBJetTagsDisc(CleanJetsNoPtCut->At(i)->TrackCountingHighEffBJetTagsDisc());
          jet_a->SetTrackCountingHighPurBJetTagsDisc(CleanJetsNoPtCut->At(i)->TrackCountingHighPurBJetTagsDisc());
          jet_a->SetSimpleSecondaryVertexBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexBJetTagsDisc());
          jet_a->SetSimpleSecondaryVertexHighEffBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexHighEffBJetTagsDisc());
          jet_a->SetSimpleSecondaryVertexHighPurBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexHighPurBJetTagsDisc());
	}
        sortedJetsAll.push_back(jet_a);       
      }

      for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); i++){
        if(TMath::Abs(CleanJetsNoPtCut->At(i)->Eta()) < fEtaJetCut &&
           CleanJetsNoPtCut->At(i)->Pt()*(1.0+fJetScaleSyst) > fPtJetCut){
          Jet* jet_b = new Jet(CleanJetsNoPtCut->At(i)->Px()*(1.0+fJetScaleSyst),
        		       CleanJetsNoPtCut->At(i)->Py()*(1.0+fJetScaleSyst),
        		       CleanJetsNoPtCut->At(i)->Pz()*(1.0+fJetScaleSyst),
        		       CleanJetsNoPtCut->At(i)->E() *(1.0+fJetScaleSyst));
          sortedJets.push_back(jet_b);
        }
      }

      for(UInt_t i=0; i<sortedJetsAll.size(); i++){
        for(UInt_t j=i+1; j<sortedJetsAll.size(); j++){
          if(sortedJetsAll[i]->Pt() < sortedJetsAll[j]->Pt()) {
            //swap i and j
            Jet* tempjet = sortedJetsAll[i];
            sortedJetsAll[i] = sortedJetsAll[j];
            sortedJetsAll[j] = tempjet;	
          }
        }
      }

      for(UInt_t i=0; i<sortedJetsAll.size(); i++){
	double dZAverageJetPt = 0.0;
        double sumJetPt = 0.0;
      	double jetPt = 0.0;
	for(UInt_t iPF=0; iPF<fPFJet0->GetEntries(); iPF++){								  
	  const PFJet *jet = fPFJet0->At(iPF);  								  
	  if(MathUtils::DeltaR(jet->Mom(),sortedJetsAll[i]->Mom()) < 0.01){
      	    jetPt = jet->Pt();
	    for (UInt_t npf=0; npf<jet->NPFCands();npf++) {
	      const PFCandidate *pf = jet->PFCand(npf);
              if(pf->BestTrk()) {
		dZAverageJetPt = dZAverageJetPt + pf->Pt()*pf->Pt()*pf->BestTrk()->DzCorrected(*fVertices->At(0));
		sumJetPt = sumJetPt + pf->Pt()*pf->Pt();
	      }
	    }
	    if(sumJetPt > 0) dZAverageJetPt = TMath::Abs(dZAverageJetPt)/sumJetPt;
	    break;
	  }
	} // loop over PF jets
	if(i == 0) dZAverageJet1Pt = dZAverageJetPt;
        //if(dZAverageJetPt >= 2.0 || jetPt <= 10){
        if(jetPt <= 10){
	  sortedJetsAll[i]->SetTrackCountingHighEffBJetTagsDisc      (-25.0);
	  sortedJetsAll[i]->SetCombinedSecondaryVertexBJetTagsDisc   (-25.0);
	  sortedJetsAll[i]->SetCombinedSecondaryVertexMVABJetTagsDisc(-25.0);
	  sortedJetsAll[i]->SetJetProbabilityBJetTagsDisc	     (-25.0);
	  sortedJetsAll[i]->SetJetBProbabilityBJetTagsDisc	     (-25.0);
	  sortedJetsAll[i]->SetTrackCountingHighPurBJetTagsDisc      (-25.0);
        }

        bool overlap = kFALSE;
        for(UInt_t j=0; j<sortedJets.size(); j++){
          if(sortedJetsAll[i]->Pt() == sortedJets[j]->Pt() ||
            (sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc() == sortedJets[j]->CombinedSecondaryVertexBJetTagsDisc() &&
             sortedJetsAll[i]->JetBProbabilityBJetTagsDisc()	     == sortedJets[j]->JetBProbabilityBJetTagsDisc() &&
             sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc()    == sortedJets[j]->TrackCountingHighPurBJetTagsDisc())
            ) {
            sortedJets[j]->SetMatchedMCFlavor(sortedJetsAll[i]->MatchedMCFlavor());
            sortedJets[j]->SetCombinedSecondaryVertexBJetTagsDisc(sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc());
            sortedJets[j]->SetCombinedSecondaryVertexMVABJetTagsDisc(sortedJetsAll[i]->CombinedSecondaryVertexMVABJetTagsDisc());
            sortedJets[j]->SetJetProbabilityBJetTagsDisc(sortedJetsAll[i]->JetProbabilityBJetTagsDisc());
            sortedJets[j]->SetJetBProbabilityBJetTagsDisc(sortedJetsAll[i]->JetBProbabilityBJetTagsDisc());
            sortedJets[j]->SetTrackCountingHighEffBJetTagsDisc(sortedJetsAll[i]->TrackCountingHighEffBJetTagsDisc());
            sortedJets[j]->SetTrackCountingHighPurBJetTagsDisc(sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc());
            sortedJets[j]->SetSimpleSecondaryVertexBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexBJetTagsDisc());
            sortedJets[j]->SetSimpleSecondaryVertexHighEffBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexHighEffBJetTagsDisc());
            sortedJets[j]->SetSimpleSecondaryVertexHighPurBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexHighPurBJetTagsDisc());	     
            overlap = kTRUE;
            break;
          }
        }

        if(overlap == kFALSE && sortedJetsAll[i]->TrackCountingHighEffBJetTagsDisc() > nMaxMediumBTagJet){
	  nMaxMediumBTagJet = sortedJetsAll[i]->TrackCountingHighEffBJetTagsDisc();
        } // overlap == false

        if(sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc() > nMaxMediumOtherBTagJet[0]){
	  nMaxMediumOtherBTagJet[0] = sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc();
        }
        if(sortedJetsAll[i]->CombinedSecondaryVertexMVABJetTagsDisc() > nMaxMediumOtherBTagJet[1]){
	  nMaxMediumOtherBTagJet[1] = sortedJetsAll[i]->CombinedSecondaryVertexMVABJetTagsDisc();
        }
        if(sortedJetsAll[i]->JetProbabilityBJetTagsDisc() > nMaxMediumOtherBTagJet[2]){
	  nMaxMediumOtherBTagJet[2] = sortedJetsAll[i]->JetProbabilityBJetTagsDisc();
        }
        if(sortedJetsAll[i]->JetBProbabilityBJetTagsDisc() > nMaxMediumOtherBTagJet[3]){
	  nMaxMediumOtherBTagJet[3] = sortedJetsAll[i]->JetBProbabilityBJetTagsDisc();
        }
        if(sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc() > nMaxMediumOtherBTagJet[4]){
	  nMaxMediumOtherBTagJet[4] = sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc();
        }
        if(sortedJetsAll[i]->TrackCountingHighEffBJetTagsDisc() > nMaxMediumOtherBTagJet[5]){
	  nMaxMediumOtherBTagJet[5] = sortedJetsAll[i]->TrackCountingHighEffBJetTagsDisc();
        }
      } // loop over all jets

      for(UInt_t i=0; i<sortedJets.size(); i++){
        for(UInt_t j=i+1; j<sortedJets.size(); j++){
          if(sortedJets[i]->Pt() < sortedJets[j]->Pt()) {
            //swap i and j
            Jet* tempjet = sortedJets[i];
            sortedJets[i] = sortedJets[j];
            sortedJets[j] = tempjet;
          }
        }
      }

      if (fFillPhotonTemplate) {
	dilepton->AddDaughter(OriginalCleanPhotons->At(0));
      } 
      else {
        if (leptons->At(0)->ObjType() == kMuon && leptons->At(1)->ObjType() == kMuon ){
          fSmurfTree.type_ = SmurfTree::mm;
        }
        else if((leptons->At(0)->ObjType() == kElectron || leptons->At(0)->ObjType() == kPhoton) &&
	        (leptons->At(1)->ObjType() == kElectron || leptons->At(1)->ObjType() == kPhoton)){
          fSmurfTree.type_ = SmurfTree::ee;
        }
        else if((leptons->At(0)->ObjType() == kElectron || leptons->At(0)->ObjType() == kPhoton) && 
	         leptons->At(1)->ObjType() == kMuon){
          fSmurfTree.type_ = SmurfTree::em;
        }
        else if(leptons->At(0)->ObjType() == kMuon && 
	       (leptons->At(1)->ObjType() == kElectron || leptons->At(1)->ObjType() == kPhoton)){
          fSmurfTree.type_ = SmurfTree::me;
        }
        else {
          cout << "Hey, this is not possible, leptonTypes: "
               << leptons->At(0)->ObjType() << " - " 
               << leptons->At(1)->ObjType() << endl;
          assert(0);
        }
        dilepton->AddDaughter(leptons->At(0));
        dilepton->AddDaughter(leptons->At(1));
      }

      // Begin Met MVA
      const JetCol *fCorrectedJets = GetObjThisEvt<JetCol> (fCorrectedJetsName);

      // lepton selection
      Float_t lPt0 = 0; Float_t lEta0 = 0; Float_t lPhi0 = 0;
      Float_t lPt1 = 0; Float_t lEta1 = 0; Float_t lPhi1 = 0;
      if(leptons->GetEntries() >= 1){
	lPt0 = leptons->At(0)->Pt();  lEta0 = leptons->At(0)->Eta(); lPhi0 = leptons->At(0)->Phi();
      }
      if(leptons->GetEntries() >= 2){
	lPt1 = leptons->At(1)->Pt();  lEta1 = leptons->At(1)->Eta(); lPhi1 = leptons->At(1)->Phi();
      }

      PFJetOArr *thePFCorrectedJets = new PFJetOArr; 
      for(UInt_t i=0; i<fCorrectedJets->GetEntries(); i++) {
	const PFJet *ptJet = dynamic_cast<const PFJet*>(fCorrectedJets->At(i));
	thePFCorrectedJets->Add(ptJet);
      }

      Double_t Rho = 0.0;
      switch(fTheRhoType) {
      case RhoUtilities::MIT_RHO_VORONOI_HIGH_ETA:
	Rho = fPileupEnergyDensity->At(0)->Rho();
	break;
      case RhoUtilities::MIT_RHO_VORONOI_LOW_ETA:
	Rho = fPileupEnergyDensity->At(0)->RhoLowEta();
	break;
      case RhoUtilities::MIT_RHO_RANDOM_HIGH_ETA:
	Rho = fPileupEnergyDensity->At(0)->RhoRandom();
	break;
      case RhoUtilities::MIT_RHO_RANDOM_LOW_ETA:
	Rho = fPileupEnergyDensity->At(0)->RhoRandomLowEta();
	break;
      case RhoUtilities::CMS_RHO_RHOKT6PFJETS:
	Rho = fPileupEnergyDensity->At(0)->RhoKt6PFJets();
	break;
      default:
	// use the old default
	Rho = fPileupEnergyDensity->At(0)->Rho();
	break;
      }
      Met theMetMVA = fMVAMet->GetMet(false,
                                      lPt0,lPhi0,lEta0,
                                      lPt1,lPhi1,lEta1,
                                      fPFMetStd->At(0),
                                      fPFCandidates,fVertices->At(0),fVertices,Rho,
                                      thePFCorrectedJets,
                                      int(fVertices->GetEntries()));
      delete thePFCorrectedJets;
      // End Met MVA
    
      // Angle between MET and closest lepton
      double deltaPhiMetLepton[3] = {0,0,0};
      double mTW[3] = {0,0,0};
      double deltaPhiLeptons = 0;
      double deltaRLeptons = 0;
      double deltaPhiDileptonMet = 0;
      double metFraction[2] = {0.0, 0.0};
      double mtHiggs = 0;
      double pMetMVA = 0.0;

      if     (leptons->GetEntries() >= 2) {
	deltaPhiMetLepton[0] = fabs(MathUtils::DeltaPhi(pfMet->Phi(), leptons->At(0)->Phi()));
	deltaPhiMetLepton[1] = fabs(MathUtils::DeltaPhi(pfMet->Phi(), leptons->At(1)->Phi()));
	mTW[0] = TMath::Sqrt(2.0*leptons->At(0)->Pt()*pfMet->Pt()*(1.0 - cos(deltaPhiMetLepton[0])));
	mTW[1] = TMath::Sqrt(2.0*leptons->At(1)->Pt()*pfMet->Pt()*(1.0 - cos(deltaPhiMetLepton[1])));
	deltaPhiLeptons = fabs(MathUtils::DeltaPhi(leptons->At(0)->Phi(), leptons->At(1)->Phi()));
	deltaRLeptons = MathUtils::DeltaR(leptons->At(0)->Mom(), leptons->At(1)->Mom());
	deltaPhiDileptonMet = fabs(MathUtils::DeltaPhi(pfMet->Phi(), dilepton->Phi()));
	mtHiggs = JetTools::MtHiggs(leptons,pfMet,metFraction,7);
	pMetMVA = theMetMVA.Pt();
	double angMin = TMath::Min(deltaPhiMetLepton[0],deltaPhiMetLepton[1]);
	if(angMin*180./TMath::Pi() < 90.0) pMetMVA = theMetMVA.Pt() * sin(angMin);
      }
      else if(fFillPhotonTemplate == kTRUE){
	deltaPhiDileptonMet = fabs(MathUtils::DeltaPhi(pfMet->Phi(), dilepton->Phi()));
	mtHiggs = TMath::Sqrt(2.0*dilepton->Pt()*pfMet->Pt()*(1.0 - cos(deltaPhiDileptonMet)));
      }

      double dPhiLep1Jet1_ = -999;
      double dRLep1Jet1_   = -999;
      double dPhiLep2Jet1_ = -999;
      double dRLep2Jet1_   = -999;
      double dPhiLep3Jet1_ = -999;
      double dRLep3Jet1_   = -999;
      
      if(sortedJetsAll.size() > 0){
	deltaPhiLLJet  = fabs(MathUtils::DeltaPhi(dilepton->Phi(), sortedJetsAll[0]->Phi()));

	if (leptons->GetEntries() >= 2){
	  dPhiLep1Jet1_ = fabs(MathUtils::DeltaPhi(leptons->At(0)->Phi(), sortedJetsAll[0]->Phi()));
	  dRLep1Jet1_   = MathUtils::DeltaR(leptons->At(0)->Mom(), sortedJetsAll[0]->Mom());
	  dPhiLep2Jet1_ = fabs(MathUtils::DeltaPhi(leptons->At(1)->Phi(), sortedJetsAll[0]->Phi()));
	  dRLep2Jet1_   = MathUtils::DeltaR(leptons->At(1)->Mom(), sortedJetsAll[0]->Mom());
	}
	if(leptons->GetEntries() > 2){
	  dPhiLep3Jet1_ = fabs(MathUtils::DeltaPhi(leptons->At(2)->Phi(), sortedJetsAll[0]->Phi()));
	  dRLep3Jet1_   = MathUtils::DeltaR(leptons->At(2)->Mom(), sortedJetsAll[0]->Mom());
	}
      }
      
      Double_t Q    	 = 0.0;
      Int_t    id1  	 = 0;
      Double_t x1   	 = 0.0;
      Double_t pdf1 	 = 0.0;
      Int_t    id2  	 = 0;
      Double_t x2   	 = 0.0;
      Double_t pdf2 	 = 0.0;
      Int_t    processId = 0;
      if(fIsData == kFALSE){
         LoadBranch(fMCEvInfoName);
         Q    	   = fMCEventInfo->Scale();
         id1  	   = fMCEventInfo->Id1();
         x1   	   = fMCEventInfo->X1();
         pdf1 	   = fMCEventInfo->Pdf1();
         id2  	   = fMCEventInfo->Id2();
         x2   	   = fMCEventInfo->X2();
         pdf2 	   = fMCEventInfo->Pdf2();
	 processId = fMCEventInfo->ProcessId();
      }
  
      //Fill Higgs Pt
      Double_t higgsPt = -999;
      Double_t vPtMin  = 100000.;
      if (GenBosons) {
        for (UInt_t i=0; i<GenBosons->GetEntries(); ++i) {        
          if (GenBosons->At(i)->PdgId() == MCParticle::kH) {
            higgsPt = GenBosons->At(i)->Pt();
          }
          if (GenBosons->At(i)->Is(MCParticle::kW) ||
	      GenBosons->At(i)->Is(MCParticle::kZ)) {
            if(GenBosons->At(i)->Pt() < vPtMin) vPtMin = GenBosons->At(i)->Pt();
	  }
        }
	if(higgsPt == -999) higgsPt = vPtMin;
      }

      //For DY MC, we save the pre-FSR dilepton mass into the higgsPt variable
      Double_t DYNNLOKFactor = 1.0;
      if ( fDecay == 6 || fDecay == 7 || fDecay == 8 ) {
        //Find ZBoson;
        const MCParticle *ZBoson = 0;
        for (UInt_t i=0; i< GenBosons->GetEntries(); ++i) {        
           if (GenBosons->At(i)->AbsPdgId() == 23) {
            ZBoson = GenBosons->At(i);
          }
        }

        if (ZBoson) {
          //Find Mass Bin
          const UInt_t nMassBins = 41;
          const Double_t massBinLimits[nMassBins+1] = {15,   20,  25,  30,  35,  40,  45,  50,  55,  60, 
                                                       64,   68,  72,  76,  81,  86,  91,  96, 101, 106, 
                                                       110, 115, 120, 126, 133, 141, 150, 160, 171, 185, 
                                                       200, 220, 243, 273, 320, 380, 440, 510, 600, 1000, 
                                                       1500}; //41 bins

          Int_t massBinIndex = -1 ;
          for(UInt_t binIndex=0; binIndex < nMassBins; ++binIndex){
            if( ZBoson->Mass() >= massBinLimits[binIndex] && ZBoson->Mass() < massBinLimits[binIndex+1]) {
              massBinIndex = binIndex;
              break;
            }
          }

          //Found the mass bin
          if (massBinIndex >= 0 && massBinIndex < Int_t(nMassBins) ) {
            UInt_t ptBin = fDYNNLOKFactorHists[massBinIndex]->GetXaxis()->FindFixBin(ZBoson->Pt());
            UInt_t yBin = fDYNNLOKFactorHists[massBinIndex]->GetYaxis()->FindFixBin(ZBoson->Rapidity());

            if(Int_t(ptBin) == fDYNNLOKFactorHists[massBinIndex]->GetNbinsX() + 1)
              ptBin = fDYNNLOKFactorHists[massBinIndex]->GetNbinsX();
            if(ptBin == 0)
              ptBin = 1;
            if(Int_t(yBin) == fDYNNLOKFactorHists[massBinIndex]->GetNbinsY() + 1)
              yBin = fDYNNLOKFactorHists[massBinIndex]->GetNbinsY();
            if(yBin == 0)
              yBin = 1;
            DYNNLOKFactor = fDYNNLOKFactorHists[massBinIndex]->GetBinContent( ptBin, yBin);
          } else {
            cout << "Error: Did not find the mass bin for boson mass: " << ZBoson->Mass() << endl;
          }
        } else {
          cout << "Error: Did not find Z boson\n"; 
        }
	higgsPt = DYNNLOKFactor;
      }

      fSmurfTree.cuts_ = 0;
      fSmurfTree.cuts_ |= SmurfTree::Trigger;
      // pt(reco)>20/10, acceptance, q1*q2<0,!STA muon, mll>12
      Bool_t PreselPtCut = kTRUE;
      if(leptons->GetEntries() >= 1 && leptons->At(0)->Pt() <= 20) PreselPtCut = kFALSE;
      if(leptons->GetEntries() >= 2 && leptons->At(1)->Pt() <= 10) PreselPtCut = kFALSE;
      if(dilepton->Mass() <= 12.0) PreselPtCut = kFALSE;
      if(PreselPtCut == kTRUE) fSmurfTree.cuts_ |= SmurfTree::BaseLine;
      if(fFakeRatePredictionType != 3){
        if     (leptons->GetEntries() >= 1 && leptons->At(0)->IsFakeable() == kFALSE) fSmurfTree.cuts_ |= SmurfTree::Lep1FullSelection;
        else if(leptons->GetEntries() >= 1 && leptons->At(0)->ObjType() == kMuon)     fSmurfTree.cuts_ |= SmurfTree::Lep1LooseMuV2;
        else if(leptons->GetEntries() >= 1 && leptons->At(0)->ObjType() == kElectron) fSmurfTree.cuts_ |= SmurfTree::Lep1LooseEleV4;
        if     (leptons->GetEntries() >= 2 && leptons->At(1)->IsFakeable() == kFALSE) fSmurfTree.cuts_ |= SmurfTree::Lep2FullSelection;
        else if(leptons->GetEntries() >= 2 && leptons->At(1)->ObjType() == kMuon)     fSmurfTree.cuts_ |= SmurfTree::Lep2LooseMuV2;
        else if(leptons->GetEntries() >= 2 && leptons->At(1)->ObjType() == kElectron) fSmurfTree.cuts_ |= SmurfTree::Lep2LooseEleV4;
        if     (leptons->GetEntries() >= 3 && leptons->At(2)->IsFakeable() == kFALSE) fSmurfTree.cuts_ |= SmurfTree::Lep3FullSelection;
        else if(leptons->GetEntries() >= 3 && leptons->At(2)->ObjType() == kMuon)     fSmurfTree.cuts_ |= SmurfTree::Lep3LooseMuV2;
        else if(leptons->GetEntries() >= 3 && leptons->At(2)->ObjType() == kElectron) fSmurfTree.cuts_ |= SmurfTree::Lep3LooseEleV4;

        if     (leptons->GetEntries() >= 1 && leptons->At(0)->IsFakeable() == kTRUE && leptonsMVAPass[0] == 1) fSmurfTree.cuts_ |= SmurfTree::Lep1LooseEleV1;
        if     (leptons->GetEntries() >= 2 && leptons->At(1)->IsFakeable() == kTRUE && leptonsMVAPass[1] == 1) fSmurfTree.cuts_ |= SmurfTree::Lep2LooseEleV1;
        if     (leptons->GetEntries() >= 3 && leptons->At(2)->IsFakeable() == kTRUE && leptonsMVAPass[2] == 1) fSmurfTree.cuts_ |= SmurfTree::Lep3LooseEleV1;
      } else { // photons treated as fakeable objects
        if     (leptons->GetEntries() >= 1 && leptons->At(0)->IsFakeable() == kFALSE) fSmurfTree.cuts_ |= SmurfTree::Lep1FullSelection;
        if     (leptons->GetEntries() >= 2 && leptons->At(1)->IsFakeable() == kFALSE) fSmurfTree.cuts_ |= SmurfTree::Lep2FullSelection;
        if     (leptons->GetEntries() >= 3 && leptons->At(2)->IsFakeable() == kFALSE) fSmurfTree.cuts_ |= SmurfTree::Lep3FullSelection;

        if     (leptons->GetEntries() >= 1 && leptons->At(0)->IsFakeable() == kTRUE) fSmurfTree.cuts_ |= SmurfTree::Lep1LooseEleV2;
        if     (leptons->GetEntries() >= 2 && leptons->At(1)->IsFakeable() == kTRUE) fSmurfTree.cuts_ |= SmurfTree::Lep2LooseEleV2;
        if     (leptons->GetEntries() >= 3 && leptons->At(2)->IsFakeable() == kTRUE) fSmurfTree.cuts_ |= SmurfTree::Lep3LooseEleV2;
      }

      // full met selection
      if(TMath::Min(pMET[0],pMET[1]) > 20.0 &&
        (fSmurfTree.type_ == SmurfTree::em || fSmurfTree.type_ == SmurfTree::me || 
	(sortedJetsAll.size() <  2 && TMath::Min(pMET[0],pMET[1]) > 45.0) ||
	(sortedJetsAll.size() >= 2 && pfMet->Pt() > 45.0)))
	fSmurfTree.cuts_ |= SmurfTree::FullMET;

      // event is not in the Z-mass peak for ee/mm final states
      if(fSmurfTree.type_ == SmurfTree::em || fSmurfTree.type_ == SmurfTree::me || TMath::Abs(dilepton->Mass()-91.1876) > 15.0)
        fSmurfTree.cuts_ |= SmurfTree::ZVeto;

      // extra lepton veto, DR(muon-electron)>=0.3
      if(leptons->GetEntries() == 2)
        fSmurfTree.cuts_ |= SmurfTree::ExtraLeptonVeto;

      // 1-jet events, where the jet is b-tagged (top control sample with one b-quark missing)
      if(sortedJetsAll.size() >= 1 && sortedJetsAll[0]->TrackCountingHighEffBJetTagsDisc() >= 2.1) fSmurfTree.cuts_ |= SmurfTree::OneBJet;

      // soft muon and b-jet tagging for the whole event regardless of n-jets (non-zero means tagged)
      Bool_t isTopTag = kFALSE;
      if(nMaxMediumBTagJet >= 2.1 || DirtyMuons->GetEntries() != 0) isTopTag = kTRUE;
      for(unsigned int i=0; i<sortedJetsAll.size(); i++) if(sortedJetsAll[i]->TrackCountingHighEffBJetTagsDisc() >= 2.1) isTopTag = kTRUE;

      if(isTopTag == kTRUE) fSmurfTree.cuts_ |= SmurfTree::TopTag;
      else                  fSmurfTree.cuts_ |= SmurfTree::TopVeto;

      // soft muon and b-jet tagging for areas outside primary jets (non-zero means tagged)
      ObjArray<Muon> *DirtyMuonsNoJet = new ObjArray<Muon>;
      for(UInt_t i=0; i<DirtyMuons->GetEntries(); i++){
        double DeltaRMin = 999.;
	for(UInt_t j=0; j<sortedJets.size(); j++){
          if(MathUtils::DeltaR(DirtyMuons->At(i)->Mom(),sortedJets[j]->Mom()) < DeltaRMin)
	    DeltaRMin = MathUtils::DeltaR(DirtyMuons->At(i)->Mom(),sortedJets[j]->Mom());
	}
	if(DeltaRMin > 0.3) DirtyMuonsNoJet->Add(DirtyMuons->At(i));
      }
      if(nMaxMediumBTagJet >= 2.1 || DirtyMuonsNoJet->GetEntries() != 0) fSmurfTree.cuts_ |= SmurfTree::TopTagNotInJets;
      delete DirtyMuonsNoJet;

      // add variables for Smurf ntuples
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > null(0.0,0.0,0.0,0.0);

      Bool_t passCuts = kTRUE;
      if(passCuts) fNEventsSelected++;

      if(zDiffMax < 10000.0 || fFillPhotonTemplate ) {
	if     (fDecay == 110) { fSmurfTree.dstype_ = SmurfTree::hww110; processId = 10010; }
	else if(fDecay == 115) { fSmurfTree.dstype_ = SmurfTree::hww115; processId = 10010; }
	else if(fDecay == 118) { fSmurfTree.dstype_ = SmurfTree::hww118; processId = 10010; }
	else if(fDecay == 120) { fSmurfTree.dstype_ = SmurfTree::hww120; processId = 10010; }
	else if(fDecay == 122) { fSmurfTree.dstype_ = SmurfTree::hww122; processId = 10010; }
	else if(fDecay == 124) { fSmurfTree.dstype_ = SmurfTree::hww124; processId = 10010; }
	else if(fDecay == 125) { fSmurfTree.dstype_ = SmurfTree::hww125; processId = 10010; }
	else if(fDecay == 126) { fSmurfTree.dstype_ = SmurfTree::hww126; processId = 10010; }
	else if(fDecay == 128) { fSmurfTree.dstype_ = SmurfTree::hww128; processId = 10010; }
	else if(fDecay == 130) { fSmurfTree.dstype_ = SmurfTree::hww130; processId = 10010; }
	else if(fDecay == 135) { fSmurfTree.dstype_ = SmurfTree::hww135; processId = 10010; }
	else if(fDecay == 140) { fSmurfTree.dstype_ = SmurfTree::hww140; processId = 10010; }
	else if(fDecay == 145) { fSmurfTree.dstype_ = SmurfTree::hww145; processId = 10010; }
	else if(fDecay == 150) { fSmurfTree.dstype_ = SmurfTree::hww150; processId = 10010; }
	else if(fDecay == 155) { fSmurfTree.dstype_ = SmurfTree::hww155; processId = 10010; }
	else if(fDecay == 160) { fSmurfTree.dstype_ = SmurfTree::hww160; processId = 10010; }
	else if(fDecay == 170) { fSmurfTree.dstype_ = SmurfTree::hww170; processId = 10010; }
	else if(fDecay == 180) { fSmurfTree.dstype_ = SmurfTree::hww180; processId = 10010; }
	else if(fDecay == 190) { fSmurfTree.dstype_ = SmurfTree::hww190; processId = 10010; }
	else if(fDecay == 200) { fSmurfTree.dstype_ = SmurfTree::hww200; processId = 10010; }
	else if(fDecay == 210) { fSmurfTree.dstype_ = SmurfTree::hww210; processId = 10010; }
	else if(fDecay == 220) { fSmurfTree.dstype_ = SmurfTree::hww220; processId = 10010; }
	else if(fDecay == 230) { fSmurfTree.dstype_ = SmurfTree::hww230; processId = 10010; }
	else if(fDecay == 250) { fSmurfTree.dstype_ = SmurfTree::hww250; processId = 10010; }
	else if(fDecay == 300) { fSmurfTree.dstype_ = SmurfTree::hww300; processId = 10010; }
	else if(fDecay == 350) { fSmurfTree.dstype_ = SmurfTree::hww350; processId = 10010; }
	else if(fDecay == 400) { fSmurfTree.dstype_ = SmurfTree::hww400; processId = 10010; }
	else if(fDecay == 450) { fSmurfTree.dstype_ = SmurfTree::hww450; processId = 10010; }
	else if(fDecay == 500) { fSmurfTree.dstype_ = SmurfTree::hww500; processId = 10010; }
	else if(fDecay == 550) { fSmurfTree.dstype_ = SmurfTree::hww550; processId = 10010; }
	else if(fDecay == 600) { fSmurfTree.dstype_ = SmurfTree::hww600; processId = 10010; }
	else if(fDecay == 700) { fSmurfTree.dstype_ = SmurfTree::hww700; processId = 10010; }
	else if(fDecay == 800) { fSmurfTree.dstype_ = SmurfTree::hww800; processId = 10010; }
	else if(fDecay == 900) { fSmurfTree.dstype_ = SmurfTree::hww900; processId = 10010; }
	else if(fDecay ==1000) { fSmurfTree.dstype_ = SmurfTree::hww1000;processId = 10010; }

	else if(fDecay == 1110){ fSmurfTree.dstype_ = SmurfTree::vbfhww110; processId = 10001; }
	else if(fDecay == 1115){ fSmurfTree.dstype_ = SmurfTree::vbfhww115; processId = 10001; }
	else if(fDecay == 1118){ fSmurfTree.dstype_ = SmurfTree::vbfhww118; processId = 10001; }
	else if(fDecay == 1120){ fSmurfTree.dstype_ = SmurfTree::vbfhww120; processId = 10001; }
	else if(fDecay == 1122){ fSmurfTree.dstype_ = SmurfTree::vbfhww122; processId = 10001; }
	else if(fDecay == 1124){ fSmurfTree.dstype_ = SmurfTree::vbfhww124; processId = 10001; }
	else if(fDecay == 1125){ fSmurfTree.dstype_ = SmurfTree::vbfhww125; processId = 10001; }
	else if(fDecay == 1126){ fSmurfTree.dstype_ = SmurfTree::vbfhww126; processId = 10001; }
	else if(fDecay == 1128){ fSmurfTree.dstype_ = SmurfTree::vbfhww128; processId = 10001; }
	else if(fDecay == 1130){ fSmurfTree.dstype_ = SmurfTree::vbfhww130; processId = 10001; }
	else if(fDecay == 1135){ fSmurfTree.dstype_ = SmurfTree::vbfhww135; processId = 10001; }
	else if(fDecay == 1140){ fSmurfTree.dstype_ = SmurfTree::vbfhww140; processId = 10001; }
	else if(fDecay == 1145){ fSmurfTree.dstype_ = SmurfTree::vbfhww145; processId = 10001; }
	else if(fDecay == 1150){ fSmurfTree.dstype_ = SmurfTree::vbfhww150; processId = 10001; }
	else if(fDecay == 1155){ fSmurfTree.dstype_ = SmurfTree::vbfhww155; processId = 10001; }
	else if(fDecay == 1160){ fSmurfTree.dstype_ = SmurfTree::vbfhww160; processId = 10001; }
	else if(fDecay == 1170){ fSmurfTree.dstype_ = SmurfTree::vbfhww170; processId = 10001; }
	else if(fDecay == 1180){ fSmurfTree.dstype_ = SmurfTree::vbfhww180; processId = 10001; }
	else if(fDecay == 1190){ fSmurfTree.dstype_ = SmurfTree::vbfhww190; processId = 10001; }
	else if(fDecay == 1200){ fSmurfTree.dstype_ = SmurfTree::vbfhww200; processId = 10001; }
	else if(fDecay == 1210){ fSmurfTree.dstype_ = SmurfTree::vbfhww210; processId = 10001; }
	else if(fDecay == 1220){ fSmurfTree.dstype_ = SmurfTree::vbfhww220; processId = 10001; }
	else if(fDecay == 1230){ fSmurfTree.dstype_ = SmurfTree::vbfhww230; processId = 10001; }
	else if(fDecay == 1250){ fSmurfTree.dstype_ = SmurfTree::vbfhww250; processId = 10001; }
	else if(fDecay == 1300){ fSmurfTree.dstype_ = SmurfTree::vbfhww300; processId = 10001; }
	else if(fDecay == 1350){ fSmurfTree.dstype_ = SmurfTree::vbfhww350; processId = 10001; }
	else if(fDecay == 1400){ fSmurfTree.dstype_ = SmurfTree::vbfhww400; processId = 10001; }
	else if(fDecay == 1450){ fSmurfTree.dstype_ = SmurfTree::vbfhww450; processId = 10001; }
	else if(fDecay == 1500){ fSmurfTree.dstype_ = SmurfTree::vbfhww500; processId = 10001; }
	else if(fDecay == 1550){ fSmurfTree.dstype_ = SmurfTree::vbfhww550; processId = 10001; }
	else if(fDecay == 1600){ fSmurfTree.dstype_ = SmurfTree::vbfhww600; processId = 10001; }
	else if(fDecay == 1700){ fSmurfTree.dstype_ = SmurfTree::vbfhww700; processId = 10001; }
	else if(fDecay == 1800){ fSmurfTree.dstype_ = SmurfTree::vbfhww800; processId = 10001; }
	else if(fDecay == 1900){ fSmurfTree.dstype_ = SmurfTree::vbfhww900; processId = 10001; }
	else if(fDecay == 2000){ fSmurfTree.dstype_ = SmurfTree::vbfhww1000;processId = 10001; }

	else if(fDecay == 3)                                         fSmurfTree.dstype_ = SmurfTree::wjets;
	else if(fDecay == 5)                                         fSmurfTree.dstype_ = SmurfTree::ttbar;
	else if(fDecay == 6)                                         fSmurfTree.dstype_ = SmurfTree::dymm;
	else if(fDecay == 7)                                         fSmurfTree.dstype_ = SmurfTree::dyee;
	else if(fDecay == 8)                                         fSmurfTree.dstype_ = SmurfTree::dytt;
	else if(fDecay == 9)                                         fSmurfTree.dstype_ = SmurfTree::dymm;
	else if(fDecay == 11 || fDecay == 12 || fDecay == 13)        fSmurfTree.dstype_ = SmurfTree::tw;
	else if(fDecay == 10)                                        fSmurfTree.dstype_ = SmurfTree::dyttDataDriven;
	else if(fDecay == 27)                                        fSmurfTree.dstype_ = SmurfTree::wz;
	else if(fDecay == 28)                                        fSmurfTree.dstype_ = SmurfTree::zz;
	else if(fDecay == 29 || fDecay == 14)                        fSmurfTree.dstype_ = SmurfTree::qqww;
	else if(fDecay == 30)                                        fSmurfTree.dstype_ = SmurfTree::ggww;
	else if(fDecay == 32)                                        fSmurfTree.dstype_ = SmurfTree::qqwwPWG;
	else if(fDecay == 33)                                        fSmurfTree.dstype_ = SmurfTree::qqww2j;
	else if(fDecay == 31)                                        fSmurfTree.dstype_ = SmurfTree::ggzz;
	else if(fDecay == 25)                                        fSmurfTree.dstype_ = SmurfTree::www;
	else if(fDecay == 17 || fDecay == 18 || fDecay == 19)	     fSmurfTree.dstype_ = SmurfTree::wgamma;
	else if(fDecay == 20)                                        fSmurfTree.dstype_ = SmurfTree::wgstar;
	else if(fDecay == 15)                                        fSmurfTree.dstype_ = SmurfTree::qcd;
	else if(fIsData == kTRUE)	                             fSmurfTree.dstype_ = SmurfTree::data;
	else if(fDecay == 34)                                        fSmurfTree.dstype_ = SmurfTree::wwewk;
	else                                                         fSmurfTree.dstype_ = SmurfTree::other;

        if(fIsData == kTRUE) {
	  //printf("DATAEVENTNTUPLE: %d %d %d %d\n",fEventHeader->EvtNum(),fEventHeader->RunNum(),fEventHeader->LumiSec(),fSmurfTree.type_);
	}
	if(fAddLheWeights){
	  fSmurfTree.lheWeights_.clear();
          for(unsigned int nlhe=0; nlhe<fLheWeights->GetEntries(); nlhe++){
	    fSmurfTree.lheWeights_.push_back(fLheWeights->At(nlhe)->Weight());
	  }
	}
	fSmurfTree.event_         = fEventHeader->EvtNum();
	fSmurfTree.run_	          = fEventHeader->RunNum();
	fSmurfTree.lumi_	  = fEventHeader->LumiSec();
	fSmurfTree.nvtx_	  = fVertices->GetEntries();
	if(fIsData == kTRUE) {
          if (NNLOWeight) {
            fSmurfTree.scale1fb_    = NNLOWeight->GetVal();
          } else {
            fSmurfTree.scale1fb_    = 1.0;
          }
          fSmurfTree.npu_         = 0;
          fSmurfTree.npuPlusOne_  = 0;
	} else {          
          if (NNLOWeight) {
            fSmurfTree.scale1fb_    = NNLOWeight->GetVal() * 1000.0;
          } else {
            fSmurfTree.scale1fb_    = 1000.0;
          }
          Double_t NPU          = 0;
          Double_t NPU_PlusOne  = 0;
          Double_t NPU_MinusOne = 0;
	  //Double_t NPU_Observed = 0;
          for (UInt_t k=0; k < fPileupInfos->GetEntries() ; ++k) {
            if (fPileupInfos->At(k)->GetBunchCrossing() ==  0) NPU          = fPileupInfos->At(k)->GetPU_NumMean();
            if (fPileupInfos->At(k)->GetBunchCrossing() ==  1) NPU_PlusOne  = fPileupInfos->At(k)->GetPU_NumInteractions();
            if (fPileupInfos->At(k)->GetBunchCrossing() == -1) NPU_MinusOne = fPileupInfos->At(k)->GetPU_NumInteractions();
	    //if (fPileupInfos->At(k)->GetBunchCrossing() ==  0) NPU_Observed = fPileupInfos->At(k)->GetPU_NumInteractions();
          }
          fSmurfTree.npu_         = NPU;
          fSmurfTree.npuPlusOne_  = NPU_PlusOne;
          fSmurfTree.npuMinusOne_ = NPU_MinusOne;
          if (fDoPileupReweighting) {
            fSmurfTree.scale1fb_ *= fPUReweighting->reweightOOT(NPU, NPU_PlusOne);
          }
	}
	fSmurfTree.met_	          = pfMet->Pt();
	fSmurfTree.metPhi_        = pfMet->Phi();
	fSmurfTree.sumet_         = pfMet->SumEt();
	fSmurfTree.metSig_        = pfMet->PFMetSig();
        if(GenMet->GetEntries() > 0){
      	  fSmurfTree.genmet_	  = GenMet->At(0)->Pt();
          fSmurfTree.genmetPhi_   = GenMet->At(0)->Phi();
	}

	if (leptons->GetEntries() >= 1) {
	  fSmurfTree.lep1_          = leptons->At(0)->Mom();
	  fSmurfTree.lq1_           = (int)leptons->At(0)->Charge();
	  if     (fSmurfTree.lq1_ == 0 && fFakeRatePredictionType == 3 && gRandom->Uniform() > 0.5) fSmurfTree.lq1_ = +1;
	  else if(fSmurfTree.lq1_ == 0 && fFakeRatePredictionType == 3                            ) fSmurfTree.lq1_ = -1;
	  if     (leptons->At(0)->ObjType() == kMuon    ) fSmurfTree.lid1_ = 13;
	  else if(leptons->At(0)->ObjType() == kElectron) fSmurfTree.lid1_ = 11;
	  else if(leptons->At(0)->ObjType() == kPhoton  ) fSmurfTree.lid1_ = 11;
	  else                                            assert(0);
	  if(fSmurfTree.lq1_ > 0) fSmurfTree.lid1_ = -1 * fSmurfTree.lid1_;
	  fSmurfTree.lmva1_         = leptonsMVA[0];
	  fSmurfTree.lep1DetEta_    = lepDetEta[0];
	}
	
	if (leptons->GetEntries() >= 2) {
	  fSmurfTree.lep2_          = leptons->At(1)->Mom();
	  fSmurfTree.lq2_           = (int)leptons->At(1)->Charge();
	  if     (fSmurfTree.lq2_ == 0 && fFakeRatePredictionType == 3 && gRandom->Uniform() > 0.5) fSmurfTree.lq2_ = +1;
	  else if(fSmurfTree.lq2_ == 0 && fFakeRatePredictionType == 3                            ) fSmurfTree.lq2_ = -1;
	  if     (leptons->At(1)->ObjType() == kMuon    ) fSmurfTree.lid2_ = 13;
	  else if(leptons->At(1)->ObjType() == kElectron) fSmurfTree.lid2_ = 11;
	  else if(leptons->At(1)->ObjType() == kPhoton  ) fSmurfTree.lid2_ = 11;
	  else                                            assert(0);
	  if(fSmurfTree.lq2_ > 0) fSmurfTree.lid2_ = -1 * fSmurfTree.lid2_;
	  fSmurfTree.lmva2_         = leptonsMVA[1];
	  fSmurfTree.lep2DetEta_    = lepDetEta[1];
	}

	//use lid1_ to save information about the photon ID
	if (fFillPhotonTemplate) {
	  const Photon* ph = OriginalCleanPhotons->At(0);
	  Bool_t passSpikeKillingEMaxOverE9 = kTRUE;
	  Bool_t passSpikeKillingSigmaiEtaiEta = kTRUE;
	  Bool_t passSpikeKillingSigmaiPhiiPhi = kTRUE;
	  Bool_t passSpikeKillingSwissCross = kTRUE;
	  Bool_t passR9Tight = kFALSE;
	  Bool_t isBarrel = kFALSE;	 

	  if (ph->SCluster() && ph->SCluster()->Seed() && ph->SCluster()->Seed()->Energy() > 5.0 
	      && ph->SCluster()->Seed()->EMax() / ph->SCluster()->Seed()->E3x3() > 0.95) 
	    passSpikeKillingEMaxOverE9 = kFALSE;
	  if (ph->SCluster() && ph->SCluster()->Seed() && ph->SCluster()->Seed()->Energy() > 5.0 
	      && ph->SCluster()->Seed()->SwissCross() > 0.95) 
	    passSpikeKillingSwissCross = kFALSE;
	  
	  if (ph->CoviEtaiEta() < 0.001) passSpikeKillingSigmaiEtaiEta = kFALSE;
	  if (TMath::Sqrt(ph->SCluster()->Seed()->CoviPhiiPhi()) < 0.001) passSpikeKillingSigmaiPhiiPhi = kFALSE;
	  if (ph->R9() >= 0.94) passR9Tight = kTRUE;
	  if (fabs(ph->SCluster()->Eta()) < 1.4442) isBarrel = kTRUE;

	  Int_t bitmask = 0;
	  if (passSpikeKillingEMaxOverE9) bitmask += 1;
	  if (passSpikeKillingSigmaiEtaiEta) bitmask += 2;
	  if (passSpikeKillingSigmaiPhiiPhi) bitmask += 4;
	  if (passSpikeKillingSwissCross) bitmask += 8;
	  if (passR9Tight) bitmask += 16;
	  if (isBarrel) bitmask += 32;
          if (ph->IsConverted()) bitmask += 64;

	  fSmurfTree.lid1_ = bitmask;

          fSmurfTree.lep1_ = metTools.Recoil();
          fSmurfTree.lep2_ = metTools.ChargedRecoil();
	}

	if(sortedJetsAll.size() >= 1) fSmurfTree.jet1_     = sortedJetsAll[0]->Mom();
	else                          fSmurfTree.jet1_     = null;
	if(sortedJetsAll.size() >= 1) fSmurfTree.jet1Btag_ = sortedJetsAll[0]->TrackCountingHighEffBJetTagsDisc();
	else                          fSmurfTree.jet1Btag_ = -999.;
	if(sortedJetsAll.size() >= 1) fSmurfTree.jet1ProbBtag_ = sortedJetsAll[0]->CombinedSecondaryVertexBJetTagsDisc();
	else                          fSmurfTree.jet1ProbBtag_ = -999.;
	fSmurfTree.jet1Dz_ = dZAverageJet1Pt;
	if(sortedJetsAll.size() >= 2) fSmurfTree.jet2_     = sortedJetsAll[1]->Mom();
	else                          fSmurfTree.jet2_     = null;
	if(sortedJetsAll.size() >= 2) fSmurfTree.jet2Btag_ = sortedJetsAll[1]->TrackCountingHighEffBJetTagsDisc();
	else                          fSmurfTree.jet2Btag_ = -999.;
	if(sortedJetsAll.size() >= 2) fSmurfTree.jet2ProbBtag_ = sortedJetsAll[1]->CombinedSecondaryVertexBJetTagsDisc();
	else                          fSmurfTree.jet2ProbBtag_ = -999.;
	fSmurfTree.njets_         = sortedJets.size();  

      	fSmurfTree.CHSMet_	  = metToolsNoPu.GetCHSMet().Pt();
      	fSmurfTree.NHSMet_	  = metToolsNoPu.GetNHSMet().Pt();
      	fSmurfTree.CHSMetPhi_	  = metToolsNoPu.GetCHSMet().Phi();
      	fSmurfTree.NHSMetPhi_	  = metToolsNoPu.GetNHSMet().Phi();

	fSmurfTree.dilep_         = dilepton->Mom();
      	fSmurfTree.trackMet_	  = metTools.GetCorrectedTrackMet().Pt();
        fSmurfTree.trackMetPhi_	  = metTools.GetCorrectedTrackMet().Phi();
	fSmurfTree.pmet_          = pMET[0];
	fSmurfTree.pTrackMet_     = pMET[1];
        fSmurfTree.metMVA_        = theMetMVA.Pt();
        fSmurfTree.metMVAPhi_     = theMetMVA.Phi();
        fSmurfTree.pmetMVA_       = pMetMVA;
	fSmurfTree.mt_            = mtHiggs;
	fSmurfTree.mt1_           = mTW[0];
	fSmurfTree.mt2_           = mTW[1];
	fSmurfTree.dPhi_          = deltaPhiLeptons;
	fSmurfTree.dR_            = deltaRLeptons;
	fSmurfTree.dPhiLep1Jet1_  = dPhiLep1Jet1_ ;
	fSmurfTree.dRLep1Jet1_    = dRLep1Jet1_;
	fSmurfTree.dPhiLep2Jet1_  = dPhiLep2Jet1_;
	fSmurfTree.dRLep2Jet1_    = dRLep2Jet1_;
	fSmurfTree.dPhiLep1MET_   = deltaPhiMetLepton[0];
	fSmurfTree.dPhiLep2MET_   = deltaPhiMetLepton[1];
	fSmurfTree.dPhiDiLepMET_  = deltaPhiDileptonMet;
	fSmurfTree.dPhiDiLepJet1_ = deltaPhiLLJet;
	fSmurfTree.lep1McId_      = leptonGenType[0];
	fSmurfTree.lep2McId_      = leptonGenType[1];
	fSmurfTree.lep1MotherMcId_= leptonMotherGenType[0];
	fSmurfTree.lep2MotherMcId_= leptonMotherGenType[1];

	fSmurfTree.jet1McId_ = typeGenFid;
	if(sortedJetsAll.size() >= 1) {
	  if(sortedJetsAll[0]->MatchedMCFlavor() != 0) printf("MatchedMCFlavor IS WORKING!!!\n");
          Bool_t isTau = kFALSE;
          for(UInt_t ntau=0; ntau<CleanTaus->GetEntries(); ntau++){
	    if(CleanTaus->At(ntau)->Pt() <= 20) continue;
            if (MathUtils::DeltaR(CleanTaus->At(ntau)->Mom(), sortedJetsAll[0]->Mom()) < 0.3) {
              isTau = kTRUE;
              break;
            }
          }
	  if(isTau == kTRUE) fSmurfTree.jet1McId_ += 10;
        }
	if     (leptons->GetEntries() == 2 && (leptonsAgreeQ[0] == kFALSE || 
	                                       leptonsAgreeQ[1] == kFALSE)) fSmurfTree.jet1McId_ += 100;
	else if(leptons->GetEntries() >= 3 && (leptonsAgreeQ[0] == kFALSE || 
	                                       leptonsAgreeQ[1] == kFALSE ||
					       leptonsAgreeQ[2] == kFALSE)) fSmurfTree.jet1McId_ += 100;

        if(hasZCand == kTRUE) fSmurfTree.jet1McId_ += 1000;

	fSmurfTree.jet2McId_ = 0;
	if(sortedJetsAll.size() >= 2) {
	  fSmurfTree.jet2McId_ = sortedJetsAll[1]->MatchedMCFlavor();
          Bool_t isTau = kFALSE;
          for(UInt_t ntau=0; ntau<CleanTaus->GetEntries(); ntau++){
	    if(CleanTaus->At(ntau)->Pt() <= 20) continue;
            if (MathUtils::DeltaR(CleanTaus->At(ntau)->Mom(), sortedJetsAll[1]->Mom()) < 0.3) {
              isTau = kTRUE;
              break;
            }
          }
	  if(isTau == kTRUE) fSmurfTree.jet2McId_ += 100;
        }

        if(fSmurfTree.lq1_ * fSmurfTree.lq2_ < 0) fSmurfTree.cuts_ |= SmurfTree::ChargeMatch; // q1*q2<0

        if(fFillNtupleType <= 4){
	  if(leptons->GetEntries() > 2) {
	    fSmurfTree.lep3_	      = leptons->At(2)->Mom();
	    fSmurfTree.lq3_	      = (int)leptons->At(2)->Charge();
	    if     (fSmurfTree.lq3_ == 0 && fFakeRatePredictionType == 3 && gRandom->Uniform() > 0.5) fSmurfTree.lq3_ = +1;
	    else if(fSmurfTree.lq3_ == 0 && fFakeRatePredictionType == 3                            ) fSmurfTree.lq3_ = -1;
	    if     (leptons->At(2)->ObjType() == kMuon    ) fSmurfTree.lid3_ = 13;
	    else if(leptons->At(2)->ObjType() == kElectron) fSmurfTree.lid3_ = 11;
	    else if(leptons->At(2)->ObjType() == kPhoton  ) fSmurfTree.lid3_ = 11;
	    else                                            assert(0);
            if(fSmurfTree.lq3_ > 0) fSmurfTree.lid3_ = -1 * fSmurfTree.lid3_;
	    fSmurfTree.lmva3_         = leptonsMVA[2];
	    fSmurfTree.lep3DetEta_    = lepDetEta[2];

            deltaPhiMetLepton[2] = fabs(MathUtils::DeltaPhi(leptons->At(2)->Phi(), pfMet->Phi()));
	    mTW[2]               = TMath::Sqrt(2.0*leptons->At(2)->Pt()*pfMet->Pt()*
          		                      (1.0 - cos(deltaPhiMetLepton[2])));
	    fSmurfTree.dPhiLep3MET_   = deltaPhiMetLepton[2];
	    fSmurfTree.mt3_	      = mTW[2];
	  } else {
	    fSmurfTree.lep3_	      = null;
	    fSmurfTree.lid3_	      = 0;
	    fSmurfTree.lq3_	      = 0;
	    fSmurfTree.mt3_	      = 0.0;
	  }
	  if(sortedJetsAll.size() >= 3) fSmurfTree.jet3_     = sortedJetsAll[2]->Mom();
	  else                          fSmurfTree.jet3_     = null;
	  if(sortedJetsAll.size() >= 3) fSmurfTree.jet3Btag_ = sortedJetsAll[2]->TrackCountingHighEffBJetTagsDisc();
	  else                          fSmurfTree.jet3Btag_ = -999.;
	  if(sortedJetsAll.size() >= 3) fSmurfTree.jet3ProbBtag_ = sortedJetsAll[2]->CombinedSecondaryVertexBJetTagsDisc();
	  else                          fSmurfTree.jet3ProbBtag_ = -999.;
	  if(sortedJetsAll.size() >= 4) fSmurfTree.jet4_     = sortedJetsAll[3]->Mom();
	  else  			fSmurfTree.jet4_     = null;
	  if(sortedJetsAll.size() >= 4) fSmurfTree.jet4Btag_ = sortedJetsAll[3]->TrackCountingHighEffBJetTagsDisc();
	  else  			fSmurfTree.jet4Btag_ = -999.;
	  if(sortedJetsAll.size() >= 4) fSmurfTree.jet4ProbBtag_ = sortedJetsAll[3]->CombinedSecondaryVertexBJetTagsDisc();
	  else  			fSmurfTree.jet4ProbBtag_ = -999.;
	  fSmurfTree.lep3McId_      = leptonGenType[2];
	  fSmurfTree.lep3MotherMcId_= leptonMotherGenType[2];

	  fSmurfTree.jet3McId_ = 0;
	  if(sortedJetsAll.size() >= 3) {
	    fSmurfTree.jet3McId_ = sortedJetsAll[2]->MatchedMCFlavor();
            Bool_t isTau = kFALSE;
            for(UInt_t ntau=0; ntau<CleanTaus->GetEntries(); ntau++){
	    if(CleanTaus->At(ntau)->Pt() <= 20) continue;
              if (MathUtils::DeltaR(CleanTaus->At(ntau)->Mom(), sortedJetsAll[2]->Mom()) < 0.3) {
        	isTau = kTRUE;
        	break;
              }
            }
	    if(isTau == kTRUE) fSmurfTree.jet3McId_ += 100;
          }

	  fSmurfTree.jet4McId_ = 0;
	  if(sortedJetsAll.size() >= 4) {
	    fSmurfTree.jet4McId_ = sortedJetsAll[3]->MatchedMCFlavor();
            Bool_t isTau = kFALSE;
            for(UInt_t ntau=0; ntau<CleanTaus->GetEntries(); ntau++){
	    if(CleanTaus->At(ntau)->Pt() <= 20) continue;
              if (MathUtils::DeltaR(CleanTaus->At(ntau)->Mom(), sortedJetsAll[3]->Mom()) < 0.3) {
        	isTau = kTRUE;
        	break;
              }
            }
	    if(isTau == kTRUE) fSmurfTree.jet4McId_ += 100;
          }

	  fSmurfTree.dPhiLep3Jet1_  = dPhiLep3Jet1_ ;
	  fSmurfTree.dRLep3Jet1_    = dRLep3Jet1_;
	  fSmurfTree.jetLowBtag_    = nMaxMediumBTagJet;
	  fSmurfTree.nSoftMuons_    = DirtyMuons->GetEntries();
          fSmurfTree.Q_ 	    = Q;
          fSmurfTree.id1_	    = id1;
          fSmurfTree.x1_	    = x1;
          fSmurfTree.pdf1_	    = pdf1;
          fSmurfTree.id2_	    = id2;
          fSmurfTree.x2_	    = x2;
          fSmurfTree.pdf2_	    = pdf2;
          fSmurfTree.processId_     = processId;
          fSmurfTree.higgsPt_       = higgsPt;
	  if(DirtyMuons->GetEntries() > 0) fSmurfTree.quadlep_ = DirtyMuons->At(0)->Mom();
	  else                             fSmurfTree.quadlep_ = null;
	}
	fSmurfTree.info_.SetNameTitle("SmurfV7 selection, dilepton events","Summer11 samples");
	fSmurfTree.tree_->Fill();
      }

      //******************************************************************************************
      // Delete / Cleanup
      //******************************************************************************************
      delete dilepton;
      for(UInt_t i=0; i<sortedJets.size(); i++) delete sortedJets[i];
      for(UInt_t i=0; i<sortedJetsAll.size(); i++) delete sortedJetsAll[i];
    } // Preselection requirements

    //********************************************************************************************
    //Delete / Cleanup
    //********************************************************************************************
    delete pfMet;
    delete DirtyMuons;
    delete leptons;
    delete CleanMuons;
    delete CleanElectrons;
    delete CleanJetsNoPtCut;
  } //End Main Loop over Fake Event Headers

  //**********************************************************************************************
  //Delete / Cleanup
  //**********************************************************************************************
  delete leptonsOnlyFake;
  delete CleanMuonsOnlyFake;
  delete CleanElectronsOnlyFake;
  delete GenLeptonsAndTaus;
}
//--------------------------------------------------------------------------------------------------
void HwwMakeNtupleMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqBranch(fPFMetName,                    fPFMetStd);
  ReqBranch(fCaloJetName0,                 fCaloJet0);
  ReqBranch(fPFJetName0,                   fPFJet0);
  ReqBranch(fMuonName,                     fMuons);
  ReqBranch(fElectronName,                 fElectrons);
  ReqBranch(fTrackName,                    fTracks);
  ReqBranch(fPFCandidatesName,             fPFCandidates);
  ReqBranch(fEvtHdrName,                   fEventHeader);
  ReqEventObject(fPileupEnergyDensityName, fPileupEnergyDensity, kTRUE);

  if(fAddLheWeights){
    printf("Reading LHE weights ON\n");
    ReqBranch(fLheWeightsName,fLheWeights);
  }

  if(fIsData == kFALSE){
    ReqBranch(fMCEvInfoName         , fMCEventInfo);
    ReqBranch(Names::gkPileupInfoBrn, fPileupInfos); 
    ReqBranch(Names::gkMCPartBrn    , fParticles); 

    //Load DY NNLO KFactor
    if ( fDecay == 6 || fDecay == 7 || fDecay == 8 ) {
      TFile *tmpFile = new TFile( "$CMSSW_BASE/src/MitPhysics/data/fewz_powheg_weights_stepwise_2011_fine7.root", "READ");
      const int nMassBins = 41;
      for(int i=0; i<nMassBins; i++){         
        TString hname = TString::Format("weight_%02d",i+1);
        TH2D *tmpHist = (TH2D*)tmpFile->Get(hname);
        assert(tmpHist);
        tmpHist->SetDirectory(0);
        fDYNNLOKFactorHists.push_back(tmpHist);
      }
      tmpFile->Close();
      delete tmpFile;
    }
  }
  fPUReweighting = new PUReweighting("$CMSSW_BASE/src/MitPhysics/data/PUReweightFactors.root", "pileup");

  //***********************************************************************************************
  // This is fixed for Muons, kIDIsoCombined
  //***********************************************************************************************
  fElectronIDMVA = new ElectronIDMVA();
  if(fMVAElVersion == 0){
    fElectronIDMVA->Initialize("BDTG method",
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml"))),
  			       ElectronIDMVA::kIDIsoCombined,RhoUtilities::CMS_RHO_RHOKT6PFJETS);
   } else {
    fElectronIDMVA->Initialize("BDTG method",
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronIDMVA_Trig_V4_EtaBin0LowPt_V4_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronIDMVA_Trig_V4_EtaBin1LowPt_V4_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronIDMVA_Trig_V4_EtaBin2LowPt_V4_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronIDMVA_Trig_V4_EtaBin0HighPt_V4_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronIDMVA_Trig_V4_EtaBin1HighPt_V4_BDTG.weights.xml"))),
                               string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/ElectronMVAWeights/ElectronIDMVA_Trig_V4_EtaBin2HighPt_V4_BDTG.weights.xml"))),
                               ElectronIDMVA::kIDIsoCombinedHWW2012TrigV4,RhoUtilities::CMS_RHO_RHOKT6PFJETS);
   }

  fMuonTools = new MuonTools();
  fMuonIDMVA = new MuonIDMVA();
  if(fMVAMuVersion == 0){
    fMuonIDMVA->Initialize("BDTG method",
                           string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/BarrelPtBin0_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/EndcapPtBin0_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/BarrelPtBin1_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/EndcapPtBin1_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/BarrelPtBin2_IDIsoCombined_BDTG.weights.xml"))),
                           string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/EndcapPtBin2_IDIsoCombined_BDTG.weights.xml"))),
                           MuonIDMVA::kIDIsoCombinedDetIso,RhoUtilities::CMS_RHO_RHOKT6PFJETS);
   } else {
     std::vector<std::string> muonidiso_weightfiles;
     muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_barrel_lowpt.weights.xml"))));
     muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_barrel_highpt.weights.xml"))));
     muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_endcap_lowpt.weights.xml"))));
     muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_endcap_highpt.weights.xml"))));
     muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_tracker.weights.xml"))));
     muonidiso_weightfiles.push_back(string((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/MuonMVAWeights/MuonIsoMVA_BDTG_V0_global.weights.xml"))));
     fMuonIDMVA->Initialize("MuonIso_BDTG_IsoRings",
                       MuonIDMVA::kIsoRingsV0,
                       kTRUE,
                       muonidiso_weightfiles,RhoUtilities::CMS_RHO_RHOKT6PFJETS);
   }

  //***********************************************************************************************
  //Met MVA
  //***********************************************************************************************
  fMVAMet    = new MVAMet();
  if(fIs42x == kTRUE) {
    printf("Met MVA, using 4x files\n");
    fMVAMet->Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_42.root"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_42.root"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1_42.root"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2_42.root")))
                	);
  } else {
    printf("Met MVA, using 5x files\n");
    fMVAMet->Initialize(TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_lowpt.weights.xml"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/mva_JetID_highpt.weights.xml"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/Utils/python/JetIdParams_cfi.py"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmet_52.root"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetphi_52.root"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu1cov_52.root"))),
                	TString((getenv("CMSSW_BASE")+string("/src/MitPhysics/data/gbrmetu2cov_52.root")))
                	);
  }
  //***********************************************************************************************
  //Create Smurf Ntuple Tree  
  //***********************************************************************************************
  fOutputFile = new TFile(fOutputName, "RECREATE");
  fSmurfTree.CreateTree(fFillNtupleType);
}

//--------------------------------------------------------------------------------------------------
TH2D *HwwMakeNtupleMod::LoadHisto(const char *name, TFile *file) const
{
  // Load histogram with given name from given file and return it.

  TH2D *ret = dynamic_cast<TH2D*>(file->Get(name));
  if (!ret) {
    Fatal("LoadHisto", "Could not load histogram %s from file %s", name, file->GetName());
    return 0;
  }
  ret->SetDirectory(0);
  return ret;
}

//--------------------------------------------------------------------------------------------------
void HwwMakeNtupleMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis
  delete fPUReweighting;

  delete fElectronIDMVA;

  delete fMuonIDMVA;

  delete fMVAMet;

  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
  cout << "selected events on HwwMakeNtupleMod(" << fFakeRatePredictionType << "," 
       << fFillNtupleType << "): " << fNEventsSelected << endl;
}

//--------------------------------------------------------------------------------------------------
void HwwMakeNtupleMod::Terminate()
{
  // Run finishing code on the client computer
}
