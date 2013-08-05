// $Id: HwwAnaMod.cc,v 1.2 2009/06/15 15:00:20 loizides Exp $

#include "MitHiggs/HwwMods/interface/HwwAnaMod.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TH1D.h>

using namespace mithep;

ClassImp(mithep::HwwAnaMod)

//--------------------------------------------------------------------------------------------------
HwwAnaMod::HwwAnaMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanJetsName(ModNames::gkCleanJetsName), 
  fCleanLeptonsName(ModNames::gkMergedLeptonsName),
  fCleanMetName(ModNames::gkCleanCaloMetName),
  fNAccCounters(0),
  fNAllJets(0),
  fNCenJets(0),
  fDelPhiLepLep(0),
  fDelPhiDilMet(0),
  fMtH(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HwwAnaMod::Process()
{
  // Process entries of the tree.

  const ParticleCol *leptons = GetObjThisEvt<ParticleCol>(fCleanLeptonsName);
  if (!leptons)
    return;
  
  const MetCol *mets = GetObjThisEvt<MetCol>(fCleanMetName);
  if (!mets || mets->GetEntries()<1)
    return;

  const JetCol *jets = GetObjThisEvt<JetCol>(fCleanJetsName);
  if (!jets)
    return;

  Int_t accCount = 0;
  fNAccCounters->Fill(accCount++);

  // make sure 2nd highest pt lepton has pt > 0
  if ((leptons->GetEntries()<2) || (leptons->At(1)->Pt()<=0))
    return;

  fNAccCounters->Fill(accCount++);

  // make sure the 3rd highest pt lepton has pt <= 10.
  if (leptons->GetEntries() >= 3 && (leptons->At(1)->Pt()>10))
    return;

  fNAccCounters->Fill(accCount++);

  const Particle *l0 = leptons->At(0);
  const Particle *l1 = leptons->At(1);
  CompositeParticle dil;
  dil.AddDaughter(l0);
  dil.AddDaughter(l1);
  
  // charge of the leptons should be opposite
  if (dil.Charge() != 0)
    return;

  fNAccCounters->Fill(accCount++);

  const Met *met = mets->At(0);

  // hww preselection
  if (l0->Pt()  < 20.0 ||
      l1->Pt()  < 10.0 ||
      met->Pt() < 30.0 ||
      dil.Pt() < 12.0)
    return;

  fNAccCounters->Fill(accCount++);

  // count the number of central jets for vetoing
  UInt_t nCenJets = 0;
  UInt_t nJets = jets->GetEntries();
  for (UInt_t j=0; j<nJets; ++j) {
    if (jets->At(j)->AbsEta() < 2.5)
      ++nCenJets;
  }

  fNAllJets->Fill(nJets);
  fNCenJets->Fill(nCenJets);

  // central jet veto
  if (nCenJets>0)
    return;

  fNAccCounters->Fill(accCount++);

  // event variables
  Double_t delPhiLepLep = MathUtils::DeltaPhi(l0->Phi(), l1->Phi());
  Double_t delPhiDilMet = MathUtils::DeltaPhi(dil.Phi(), met->Phi());
  Double_t mtH = TMath::Sqrt(2*dil.Pt() * met->Pt() * (1 - TMath::Cos(delPhiDilMet)));

  fDelPhiLepLep->Fill(delPhiLepLep);
  fDelPhiDilMet->Fill(delPhiDilMet);
  fMtH->Fill(mtH);
}

//--------------------------------------------------------------------------------------------------
void HwwAnaMod::SlaveBegin()
{
  // Create histograms

  AddTH1(fNAccCounters,"hNAccCounters",";cut;#",25,-0.5,24.5);
  AddTH1(fNAllJets,"hNAllJets",";N_{jets};#",10,-0.5,9.5);
  AddTH1(fNCenJets,"hNCentralJets",";N_{jets};#",10,-0.5,10.5);
  AddTH1(fDelPhiLepLep,"hDelPhiLepLep",";#Delta #phi_{ll};#",60,0,TMath::Pi());
  AddTH1(fDelPhiDilMet,"hDelPhiDilMet",";#Delta #phi_{dilep,met};#",60,0,TMath::Pi());
  AddTH1(fMtH,"hTransMassH","; m_{t} [GeV];#",120,0,360);
}
