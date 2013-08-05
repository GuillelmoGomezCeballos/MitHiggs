// $Id: HwwwSignalMCMod.cc,v 1.12 2009/07/11 20:31:52 loizides Exp $

#include "MitHiggs/HwwwMods/interface/HwwwSignalMCMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/CompoundParticle.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include <TH1D.h>
#include <TH2D.h>

using namespace mithep;

ClassImp(mithep::HwwwSignalMCMod)

//--------------------------------------------------------------------------------------------------
HwwwSignalMCMod::HwwwSignalMCMod(const char *name, const char *title) : 
  GeneratorMod(name,title),
  fDoSelect(kFALSE),
  fDoSS(kFALSE),
  fMinPt(10),
  fMinMaxPt(20),
  fMinMissEt(30),
  fDilMinMass(12),
  fNAccCounters(0),
  fMcHMt(0),
  fMcWWMt(0),
  fMcWWMt2(0),
  fMcWWWMt(0),
  fMcMissEt(0),
  fMcMissEtAll(0),
  fMcLeptonPt(0),
  fMc2lDiLeptonM(0),
  fMc3lDiLeptonM(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HwwwSignalMCMod::ComputeMt(const Particle *l0, const Particle *l1, const Particle *met)
{
  // Compute properties between the two leptons and missing et.

  if (!l0 || !l1 || !met)
    return;

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

  fMc3lDiLeptonM->Fill(dil.Mass());
  fMcHists[3]->Fill(delPhiLepLep);
  fMcHists[4]->Fill(delPhiDilMet);
  fMcWWWMt->Fill(mtH);
}

//--------------------------------------------------------------------------------------------------
void HwwwSignalMCMod::Process()
{
  // Process entries of the tree.

  fNAccCounters->Fill(0);

  MCParticleCol *bosons = dynamic_cast<MCParticleCol*>(FindObjThisEvt(fMCBosonsName));
  if (bosons->GetEntries()>5) {
    SendError(kWarning, "Process",
              "Event listing should not contain more than 5 bosons, but %d are found",
              bosons->GetEntries());
    if (fDoSelect) 
      SkipEvent();
    return;
  }

  // get H0 and W bosons
  const MCParticle *h0 = 0;
  const MCParticle *w0 = 0;
  const MCParticle *w1 = 0;
  const MCParticle *wh1 = 0;
  const MCParticle *wh2 = 0;
  for (UInt_t i=0; i<bosons->GetEntries(); ++i) {
    const MCParticle *mc = bosons->At(i);

    if (mc->Is(MCParticle::kH)) {
      h0 = mc;
      continue;
    }
    else if (mc->Is(MCParticle::kW)) {
      if (!w0 && !mc->Mother()->Is(MCParticle::kH)) 
        w0 = mc;
      else if (!w1 && !mc->Mother()->Is(MCParticle::kH)) 
        w1 = mc;
      else if (!wh1 && mc->Mother()->Is(MCParticle::kH)) 
        wh1=mc;
      else if (!wh2 && mc->Mother()->Is(MCParticle::kH)) 
        wh2=mc;
      else {
        SendError(kWarning, "Process", 
                  "There should not be more than 4 produced W!");
        continue;
      } 
    }
  }

  // make sure that we got the 4 bosons
  fNAccCounters->Fill(1);
  if (!h0 || !w0 || !wh1 || !wh2 ) {
    SendError(kWarning, "Process", 
              "Could not properly identify the 4 bosons!");
    if (fDoSelect) 
      SkipEvent();
    return;
  }

  // make sure we got the good leptons
  fNAccCounters->Fill(2);
  MCParticleCol *gleptons = dynamic_cast<MCParticleCol*>(FindObjThisEvt(fMCLeptonsName));
  if (!gleptons) {
    SendError(kWarning, "Process", "Problem with goot lepton collection: ptr == NULL!");
    if (fDoSelect) 
      SkipEvent();
    return;
  }
  
  const MCParticle *l0 = 0;
  const MCParticle *l1 = 0;
  MCParticleOArr ls;
  for (UInt_t i = 0; i<gleptons->GetEntries(); ++i) {
    const MCParticle *mc = gleptons->At(i);

    Bool_t ishmo = mc->HasMother(MCParticle::kH);
    if (!ishmo) {
      ls.Add(mc);
      continue;
    }
      
    if (!l0)
      l0 = mc;
    else if (!l1)
      l1 = mc;
    else {
      SendError(kWarning, "Process", 
                "More than 2 good leptons from Higgs decay should not happen!");
      continue;
    }
  }
  ls.Sort();

  // make sure we found one lepton from higgs decay
  fNAccCounters->Fill(3);
  if (!l0) {
    if (fDoSelect) 
      SkipEvent();
    return;
  }

  // make sure we get third lepton from independent W
  fNAccCounters->Fill(4);
  const MCParticle *l2 = 0;
  if (ls.Entries()>0) 
    l2=ls.At(0);
  if (!l2) {
    if (fDoSelect) 
      SkipEvent();
    return;
  }

  // make sure we get either SS lepton pair or 3 accepted leptons
  fNAccCounters->Fill(5);
  if (!l1) {
    if (fDoSS) {
      if (l0->Charge()!=l2->Charge()) {
        if (fDoSelect) 
          SkipEvent();
        return;
      }
    } else {
      if (fDoSelect) 
        SkipEvent();
      return;
    }
  }

  // get neutrinos
  fNAccCounters->Fill(6);
  MCParticleCol *gneutrino = dynamic_cast<MCParticleCol*>(FindObjThisEvt(fMCNeutrinosName));
  if (!gleptons) {
    SendError(kWarning, "Process", "Problem with neutrino collection: ptr == NULL!");
    if (fDoSelect) 
      SkipEvent();
    return;
  }

  CompoundParticle vt;
  CompoundParticle vtAll;
  for (UInt_t i = 0; i<gneutrino->GetEntries(); ++i) {
    const MCParticle *mc = gneutrino->At(i);
    vtAll.AddDaughter(mc);
    if (mc->HasMother(MCParticle::kH)) {
      vt.AddDaughter(mc);
    }
  }

  // do event selection if required
  fNAccCounters->Fill(7);
  if (fDoSelect) {
    if ((l0 && l0->Pt() < fMinMaxPt) ||
        (l1 && l1->Pt() < fMinPt)    ||
        (l2 && l2->Pt() < fMinPt)    ||
        (vt.Pt()  < fMinMissEt)) {
      if (fDoSelect) 
        SkipEvent();
      return;
    }
  }

  fNAccCounters->Fill(8);

  // from here on just fill histograms
  if (!GetFillHist())
    return;

  if (l0) 
    fMcHists[0]->Fill(l0->Pt());
  if (l1) 
    fMcHists[1]->Fill(l1->Pt());
  for (UInt_t i=0; i<ls.Entries(); ++i)
    fMcHists[2]->Fill(ls.At(i)->Pt());

  if (h0)
    fMcHMt->Fill(h0->TMass());
  if (l0) 
    fMcLeptonPt->Fill(l0->Pt());
  if (l1) 
    fMcLeptonPt->Fill(l1->Pt());
  if (ls.Entries()>0) 
    fMcLeptonPt->Fill(ls.At(0)->Pt());

  if (l0 && l1) { //H->WW true case
    CompoundParticle h;
    h.AddDaughter(l0);
    h.AddDaughter(l1);
    if (h.Mass() >= fDilMinMass) {
      Double_t dLL = MathUtils::DeltaPhi(l0->Phi(), l1->Phi());
      Double_t dphi = MathUtils::DeltaPhi(h.Phi(), vt.Phi());
      Double_t mt = TMath::Sqrt(2 * h.Pt() * vt.Pt() * (1 - TMath::Cos(dphi))); 
      fMcWWMt->Fill(mt);
      fMcMissEt->Fill(vt.Pt());
      fMc2lDiLeptonM->Fill(h.Mass());
      fMcHists[5]->Fill(dLL);
      fMcHists[7]->Fill(h.Mass());

      // check with composite particle
      CompoundParticle h1;
      h1.AddDaughter(l0);
      h1.AddDaughter(l1);
      for (UInt_t i=0; i<vt.NDaughters(); ++i) 
        h1.AddDaughter(vt.Daughter(i));
      fMcWWMt2->Fill(h1.TMass());
      fMcHists[6]->Fill(MathUtils::DeltaPhi(l0->Phi(), vtAll.Phi()));
      fMcHists[6]->Fill(MathUtils::DeltaPhi(l1->Phi(), vtAll.Phi()));
    }
  }

  if(l2) { //H->WWW cases
    ComputeMt(l0, l1, &vtAll);
    ComputeMt(l0, l2, &vtAll);
    ComputeMt(l1, l2, &vtAll);

    fMcMissEtAll->Fill(vtAll.Pt());

    if (l1) {
      Double_t dR01 = MathUtils::DeltaR(l0->Mom(), l1->Mom());
      Double_t dR12 = MathUtils::DeltaR(l1->Mom(), l2->Mom());
      Double_t dR02 = MathUtils::DeltaR(l0->Mom(), l2->Mom());
      if ((dR01<dR02) && (dR01<dR12)) {
        fMcHists[8]->Fill(dR01);
        fMcHists[11]->Fill(TMath::Min(dR12,dR02));
      } else {
        fMcHists[9]->Fill(TMath::Min(dR12,dR02));
        fMcHists[10]->Fill(dR01);
      }
    }
  }
}

//--------------------------------------------------------------------------------------------------
void HwwwSignalMCMod::SlaveBegin()
{
  // Book branch and histograms if wanted.

  ReqBranch(fMCPartName, fParticles);

  if (!GetFillHist())
    return;

  AddTH1(fMcHMt,"hMcHMt","H->WW trans mass (MC); m_{t} [GeV]",150,0,300);
  AddTH1(fMcWWMt,"hMcWWMt","H->WW->2l2n trans mass (MC); m_{t} [GeV]",150,0,300);
  AddTH1(fMcWWMt2,"hMcWWMt2","H->WW->2l2n trans mass from comp part (MC); m_{t} [GeV]",150,0,300);
  AddTH1(fMcWWWMt,"hMcWWWMt","WH->WWW->3l3n trans mass (MC); m_{t} [GeV]",150,0,300);
  AddTH1(fMcMissEt,"hMcMissEt","Missing Et for H->WW (MC); #slash{E}_{t} [GeV]",150,0,300);
  AddTH1(fMcMissEtAll,"hMcMissEtAll","Missing Et for H->WWW (MC); E_{t} [GeV]",150,0,300);
  AddTH1(fMcLeptonPt,"hMcLeptonPt","Lepton trans mom (MC); p_{t} [GeV]",150,0,300);
  AddTH1(fMc2lDiLeptonM,"hMc2lDiLeptonM","Di-Lepton (2l) trans mass (MC); m_{t} [GeV]",150,0,300);
  AddTH1(fMc3lDiLeptonM,"hMc3lDiLeptonM","Di-Lepton (3l) trans mass (MC); m_{t} [GeV]",150,0,300);
  AddTH1(fNAccCounters,"hNAccCounters",";cut;#",25,-0.5,24.5);
  if (1) {
    TAxis *xa = fNAccCounters->GetXaxis();
    for(Int_t i=1;i<=fNAccCounters->GetNbinsX();++i)
      xa->SetBinLabel(i,"unused");
    xa->SetBinLabel(1,"Enter");
    xa->SetBinLabel(2,"Bosons");
    xa->SetBinLabel(3,"HWWW");
    xa->SetBinLabel(4,"Leptons");
    xa->SetBinLabel(5,"1LepH");
    xa->SetBinLabel(6,"1LepW");
    if (fDoSS)
      xa->SetBinLabel(7,"SS");
    else 
      xa->SetBinLabel(7,"3Lep");
    xa->SetBinLabel(8,"Neutrinos");
    xa->SetBinLabel(9,"Kinematics");
    xa->SetRangeUser(0,8);
  }
  AddTH1(fMcHists[0],"hMc1LeptonPt","Lepton (from good W) trans mom (MC); p_{t} [GeV]",150,0,300);
  AddTH1(fMcHists[1],"hMc2LeptonPt","Lepton (from good W) trans mom (MC); p_{t} [GeV]",150,0,300);
  AddTH1(fMcHists[2],"hMc3LeptonPt","Lepton (from others) trans mom (MC); p_{t} [GeV]",150,0,300);
  AddTH1(fMcHists[3],"hMcDelPhiLepLep",";#Delta #phi_{ll};#",60,0,TMath::Pi());
  AddTH1(fMcHists[4],"hMcDelPhiDilMet",";#Delta #phi_{dilep,met};#",60,0,TMath::Pi());
  AddTH1(fMcHists[5],"hMcDelPhiDilH",";#Delta #phi_{ll};#",60,0,TMath::Pi());
  AddTH1(fMcHists[6],"hMcDelPhiDilHMet",";#Delta #phi_{dilep,met};#",60,0,TMath::Pi());
  AddTH1(fMcHists[7],"hMcDilHMass",";m_{ll};#",150,0,300);
  AddTH1(fMcHists[8],"hMcDelPhiL01",";#Delta R_{true};#",120,0,6);
  AddTH1(fMcHists[9],"hMcDelPhiLxx",";#Delta R_{false};#",120,0,6);
  AddTH1(fMcHists[10],"hMcDelPhiL01b",";#Delta R_{true};#",120,0,6);
  AddTH1(fMcHists[11],"hMcDelPhiLxxb",";#Delta R_{false};#",120,0,6);
}
