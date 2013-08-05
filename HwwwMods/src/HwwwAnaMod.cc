// $Id: HwwwAnaMod.cc,v 1.14 2009/07/11 20:31:52 loizides Exp $

#include "MitHiggs/HwwwMods/interface/HwwwAnaMod.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/CompoundParticle.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TH1D.h>

using namespace mithep;

ClassImp(mithep::HwwwAnaMod)

//--------------------------------------------------------------------------------------------------
HwwwAnaMod::HwwwAnaMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanJetsName(ModNames::gkCleanJetsName), 
  fCleanLeptonsName(ModNames::gkMergedLeptonsName),
  fCleanMetName(ModNames::gkCleanCaloMetName),
  fMinPt(10),   
  fMinMaxPt(10),
  fMinMissEt(30),
  fDilMinMass(12),
  fNAccCounters(0),
  fNAllJets(0),
  fNCenJets(0),
  fDelPhiLepLep(0),
  fDelPhiDilMet(0),
  fMtH(0),
  fMissEtBefore(0),
  fMissEtAfter(0),
  fDiLeptonM(0),
  fL1L2Mass(0),
  fL1L2DelPhi(0),
  fL1L2DelMet(0)
{
  // Constructor.
}

void HwwwAnaMod::ComputeMt(const Particle *l0, const Particle *l1, const Met *met)
{
  // Compute properties between the two leptons and missing et.

  CompoundParticle dil;
  dil.AddDaughter(l0);
  dil.AddDaughter(l1);

  // charge of the leptons should be opposite
  if (dil.Charge() != 0)
    return;

  if (dil.Mass() < fDilMinMass)
    return;

  // event variables
  Double_t delPhiLepLep = MathUtils::DeltaPhi(l0->Phi(), l1->Phi());
  Double_t delPhiDilMet = MathUtils::DeltaPhi(dil.Phi(), met->Phi());
  Double_t mtH = TMath::Sqrt(2*dil.Pt() * met->Pt() * (1 - TMath::Cos(delPhiDilMet)));

  fDelPhiLepLep->Fill(delPhiLepLep);
  fDelPhiDilMet->Fill(delPhiDilMet);
  fDiLeptonM->Fill(dil.Mass());
  fMtH->Fill(mtH);
}

//--------------------------------------------------------------------------------------------------
void HwwwAnaMod::Process()
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

  // make sure 3rd highest pt lepton has pt > 0
  if ((leptons->GetEntries()<3) || (leptons->At(2)->Pt()<=0))
    return;

  fNAccCounters->Fill(2);

  // make sure the 4th highest pt lepton has pt <= fMinPt.
  if (leptons->GetEntries() >= 4 && (leptons->At(3)->Pt()>fMinPt))
    return;

  fNAccCounters->Fill(3);

  const Particle *l0 = leptons->At(0);
  const Particle *l1 = leptons->At(1);
  const Particle *l2 = leptons->At(2);

  const Met *met = mets->At(0);
  fMissEtBefore->Fill(met->Pt());

  // hww preselection
  if (l0->Pt()  < fMinMaxPt ||
      l1->Pt()  < fMinPt    ||
      l2->Pt()  < fMinPt    ||
      met->Pt() < fMinMissEt)
    return;

  fNAccCounters->Fill(4);

  // count the number of central jets for vetoing
  UInt_t nCenJets = 0;
  UInt_t nJets = jets->GetEntries();
  for (UInt_t j=0; j<nJets; ++j) {
    if (jets->At(j)->AbsEta() < 2.5)
      ++nCenJets;
  }

  fNAllJets->Fill(nJets);
  fNCenJets->Fill(nCenJets);
  fMissEtAfter->Fill(met->Pt());

  ComputeMt(l0,l1,met);
  ComputeMt(l1,l2,met);
  ComputeMt(l0,l2,met);

  const Particle *lp1 = 0;
  const Particle *lp2 = 0;
  Double_t minDll = 1e10;
  if (l0->Charge()!=l1->Charge()) {
    CompoundParticle dil;
    dil.AddDaughter(l0);
    dil.AddDaughter(l1);
    Double_t mass = dil.Mass();
    Double_t tDll = MathUtils::DeltaPhi(l0->Phi(), l1->Phi());
    if ((mass>12) && (tDll<minDll)) {
      minDll = tDll;
      lp1 = l0;
      lp2 = l1;
    }
  }
  if (l1->Charge()!=l2->Charge()) {
    CompoundParticle dil;
    dil.AddDaughter(l1);
    dil.AddDaughter(l2);
    Double_t tDll = MathUtils::DeltaPhi(l1->Phi(), l2->Phi());
    Double_t mass = dil.Mass();
    if ((mass>12) && (tDll<minDll)) {
      minDll = tDll;
      lp1 = l1;
      lp2 = l2;
    }
  }
  if (l0->Charge()!=l2->Charge()) {
    CompoundParticle dil;
    dil.AddDaughter(l0);
    dil.AddDaughter(l2);
    Double_t tDll = MathUtils::DeltaPhi(l0->Phi(), l2->Phi());
    Double_t mass = dil.Mass();
    if ((mass>12) && (tDll<minDll)) {
      minDll = tDll;
      lp1 = l0;
      lp2 = l2;
    }
  }

  if (lp1 && lp2) {
    CompoundParticle dil;
    dil.AddDaughter(lp1);
    dil.AddDaughter(lp2);
    fL1L2Mass->Fill(dil.Mass());
    fL1L2DelPhi->Fill(minDll);
    fL1L2DelMet->Fill(MathUtils::DeltaPhi(dil.Phi(), met->Phi()));
  }
}

//--------------------------------------------------------------------------------------------------
void HwwwAnaMod::SlaveBegin()
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
    xa->SetBinLabel(4,"3Lep");
    xa->SetBinLabel(5,"Sel");
    xa->SetRangeUser(0,4);
  }
  AddTH1(fNAllJets,"hNAllJets",";N_{jets};#",10,-0.5,9.5);
  AddTH1(fNCenJets,"hNCentralJets",";N_{jets};#",10,-0.5,10.5);
  AddTH1(fDelPhiLepLep,"hDelPhiLepLep",";#Delta #phi_{ll};#",60,0,TMath::Pi());
  AddTH1(fDelPhiDilMet,"hDelPhiDilMet",";#Delta #phi_{dilep,met};#",60,0,TMath::Pi());
  AddTH1(fMtH,"hTransMassH","; m_{t} [GeV];#",120,0,360);
  AddTH1(fMissEtBefore,"hMissEtBeforeCuts","; #slash{E}_{t} [GeV];#",150,0,300);
  AddTH1(fMissEtAfter,"hMissEtAfterCuts","; #slash{E}_{t} [GeV];#",150,0,300);
  AddTH1(fDiLeptonM,"hDileptonMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fL1L2Mass,"hL1L2Mass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fL1L2DelPhi,"hL1L2DelPhi",";#Delta #phi_{ll};#",60,0,TMath::Pi());
  AddTH1(fL1L2DelMet,"hL1L2DelMet",";#Delta #phi_{dilep,met};#",60,0,TMath::Pi());
}
