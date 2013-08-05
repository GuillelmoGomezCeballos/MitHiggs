// $Id: H2lAnaMod.cc,v 1.6 2009/07/11 20:31:52 loizides Exp $

#include "MitHiggs/HwwwMods/interface/H2lAnaMod.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/CompoundParticle.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TH1D.h>

using namespace mithep;

ClassImp(mithep::H2lAnaMod)

//--------------------------------------------------------------------------------------------------
H2lAnaMod::H2lAnaMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanJetsName(ModNames::gkCleanJetsName), 
  fCleanLeptonsName(ModNames::gkMergedLeptonsName),
  fCleanMetName(ModNames::gkCleanCaloMetName),
  fMinPt(10),
  fMinMaxPt(10),
  fDilMinMass(12),
  fMinMissEt(30),
  fMinJetPt(35),
  fMinZMass(70),
  fMaxZMass(110),
  fNAccCounters(0),
  fNLepCounts(0),
  fNAllJets(0),
  fNCenJets(0),
  fAllDiLepMass(0),
  fElElMass(0),
  fElMuMass(0),
  fMuMuMass(0),
  fMissEtBefore(0),
  fMissEtAfter(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void H2lAnaMod::Process()
{
  // Process entries of the tree.

  fNAccCounters->Fill(0);

  const ParticleCol *leptons = GetObjThisEvt<ParticleCol>(fCleanLeptonsName);
  if (!leptons)
    return;
  
  const MetCol *mets = GetObjThisEvt<MetCol>(fCleanMetName);
  if (!mets || mets->GetEntries()<1)
    return;

  const JetCol *jets = GetObjThisEvt<JetCol>(fCleanJetsName);
  if (!jets)
    return;

  fNAccCounters->Fill(1);

  // make sure we have at least 2 leptons
  if (leptons->GetEntries()<2)
    return;

  fNAccCounters->Fill(2);

  // make sure the 3rd highest pt lepton has pt <= fMinPt.
  if (leptons->GetEntries() >= 3 && (leptons->At(2)->Pt()>fMinPt))
    return;

  fNAccCounters->Fill(3);

  // make sure we have two same sign leading leptons
  const Particle *l0 = leptons->At(0);
  const Particle *l1 = leptons->At(1);
  if (l0->Charge()!=l1->Charge())
    return;

  fNAccCounters->Fill(4);

  const Met *met = mets->At(0);
  fMissEtBefore->Fill(met->Pt());

  CompoundParticle dil;
  dil.AddDaughter(l0);
  dil.AddDaughter(l1);
  Double_t mass=dil.Mass();
  fAllDiLepMass->Fill(mass);

  // hww preselection
  if (l0->Pt()   < fMinMaxPt   ||
      l1->Pt()   < fMinPt      ||
      mass       < fDilMinMass ||
      met->Pt()  < fMinMissEt)
    return;

  fNAccCounters->Fill(5);

  // count the number of central jets for vetoing
  UInt_t nCenJets = 0;
  UInt_t nJets = jets->GetEntries();
  for (UInt_t j=0; j<nJets; ++j) {
    if (jets->At(j)->AbsEta() > 2.5)
      continue;
    
    FourVector cj(jets->At(j)->Mom());
    cj*=jets->At(j)->L2RelativeCorrectionScale();
    cj*=jets->At(j)->L3AbsoluteCorrectionScale();

    Double_t pt = cj.Pt();
    if (pt<fMinJetPt)
      continue;
    ++nCenJets;
  }

  fNAllJets->Fill(nJets);
  fNCenJets->Fill(nCenJets);

  // cut on central jets
  if (nCenJets>2)
    return;

  fNAccCounters->Fill(6);

  if (l0->Is(kElectron) && l1->Is(kElectron)) {
    if (mass>fMinZMass && mass<fMaxZMass) 
      return;
  }

  fNAccCounters->Fill(7);
  fMissEtAfter->Fill(met->Pt());

  if (l0->ObjType()!=l1->ObjType()) {
    fNLepCounts->Fill(0);
    fElMuMass->Fill(mass);
  } else if (l0->Is(kMuon)) {
    fNLepCounts->Fill(1);
    fMuMuMass->Fill(mass);
  } else if (l1->Is(kElectron)) {
    fNLepCounts->Fill(2);
    fElElMass->Fill(mass);
  }
}

//--------------------------------------------------------------------------------------------------
void H2lAnaMod::SlaveBegin()
{
  // Create histograms.

  AddTH1(fNAccCounters,"hNAccCounters",";cut;#",25,-0.5,24.5);
  if (1) {
    TAxis *xa = fNAccCounters->GetXaxis();
    for(Int_t i=1;i<=fNAccCounters->GetNbinsX();++i)
      xa->SetBinLabel(i,"unused");
    xa->SetBinLabel(1,"Enter");
    xa->SetBinLabel(2,"Objs");
    xa->SetBinLabel(3,"2Lep");
    xa->SetBinLabel(4,"2Lep+");
    xa->SetBinLabel(5,"SS");
    xa->SetBinLabel(6,"Sel");
    xa->SetBinLabel(7,"CJet");
    xa->SetBinLabel(8,"Zee");
    xa->SetRangeUser(0,7);
  }
  AddTH1(fNLepCounts,"hNLepCounts",";cut;#",5,-0.5,4.5);
  if (1) {
    TAxis *xa = fNLepCounts->GetXaxis();
    for(Int_t i=1;i<=fNLepCounts->GetNbinsX();++i)
      xa->SetBinLabel(i,"unused");
    xa->SetBinLabel(1,"ElMu");
    xa->SetBinLabel(2,"MuMu");
    xa->SetBinLabel(3,"ElEl");
    xa->SetRangeUser(0,2);
  }
  AddTH1(fNAllJets,"hNAllJets",";N_{jets};#",10,-0.5,9.5);
  AddTH1(fNCenJets,"hNCentralJets",";N_{jets};#",10,-0.5,10.5);
  AddTH1(fAllDiLepMass,"hAllDiLepMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fElElMass,"hElElMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fElMuMass,"hElMuMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fMuMuMass,"hMuMuMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fMissEtBefore,"hMissEtBeforeCuts",";#slash{E}_{t} [GeV];#",150,0,300);
  AddTH1(fMissEtAfter,"hMissEtAfterCuts",";#slash{E}_{t} [GeV];#",150,0,300);
}
