// $Id: HwwwNtupMod.cc,v 1.4 2009/07/11 20:31:52 loizides Exp $

#include "MitHiggs/HwwwMods/interface/HwwwNtupMod.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/CompoundParticle.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TNtuple.h>

using namespace mithep;

ClassImp(mithep::HwwwNtupMod)

//--------------------------------------------------------------------------------------------------
HwwwNtupMod::HwwwNtupMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fCleanJetsName(ModNames::gkCleanJetsName), 
  fCleanLeptonsName(ModNames::gkMergedLeptonsName),
  fCleanMetName(ModNames::gkCleanCaloMetName),
  fMinPt(10),   
  fMinMaxPt(10),
  fMinMissEt(20),
  fMinJetPt(15),
  fWeight(1),
  fNtuple(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HwwwNtupMod::Process()
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

  // make sure 3rd highest pt lepton has pt > 0
  Int_t nleps = leptons->GetEntries();
  if ((nleps<3) || (leptons->At(2)->Pt()<=0))
    return;

  // make sure the 4th highest pt lepton has pt <= fMinPt.
  if (nleps >= 4 && (leptons->At(3)->Pt()>fMinPt))
    return;

  const Particle *l0 = leptons->At(0);
  const Particle *l1 = leptons->At(1);
  const Particle *l2 = leptons->At(2);
  const Met *met = mets->At(0);

  // hww preselection
  if (l0->Pt()  < fMinMaxPt ||
      l1->Pt()  < fMinPt    ||
      l2->Pt()  < fMinPt    ||
      met->Pt() < fMinMissEt)
    return;

  const Particle *lp1 = 0;
  const Particle *lp2 = 0;
  const Particle *lp3 = 0;
  Double_t minRll = 1e10;

  if (l0->Charge()!=l1->Charge()) {
    CompoundParticle dil;
    dil.AddDaughter(l0);
    dil.AddDaughter(l1);
    Double_t tRll = MathUtils::DeltaR(l0->Mom(), l1->Mom());
    if (tRll<minRll) {
      minRll = tRll;
      lp1 = l0;
      lp2 = l1;
      lp3 = l2;
    }
  }
  if (l0->Charge()!=l2->Charge()) {
    CompoundParticle dil;
    dil.AddDaughter(l0);
    dil.AddDaughter(l2);
    Double_t tRll = MathUtils::DeltaR(l0->Mom(), l2->Mom());
    if (tRll<minRll) {
      minRll = tRll;
      lp1 = l0;
      lp2 = l2;
      lp3 = l1;
    }
  }
  if (l1->Charge()!=l2->Charge()) {
    CompoundParticle dil;
    dil.AddDaughter(l1);
    dil.AddDaughter(l2);
    Double_t tRll = MathUtils::DeltaR(l1->Mom(), l2->Mom());
    if (tRll<minRll) {
      minRll = tRll;
      lp1 = l1;
      lp2 = l2;
      lp3 = l0;
    }
  }

  if (!lp1 || !lp2 || !lp3)
    return;

  CompoundParticle dil;
  dil.AddDaughter(lp1);
  dil.AddDaughter(lp2);

  Double_t sumch=0;
  for (Int_t i=0; i<nleps; ++i) 
    sumch+=leptons->At(i)->Charge();

  Double_t zmass    = 0, zmassd    = 1e12;
  Double_t zmass2   = 0, zmass2d   = 1e12;
  Double_t zmassany = 0, zmassanyd = 1e12;
  for (Int_t i=0; i<nleps; ++i) {
    const Particle *li = leptons->At(i);
    if (li->Pt()<fMinPt)
      continue;
    for (Int_t j=0; j<i; ++j) {
      const Particle *lj = leptons->At(j);
      if (lj->Pt()<fMinPt)
        continue;

      CompoundParticle dil;
      dil.AddDaughter(li);
      dil.AddDaughter(lj);

      Double_t mass  = dil.Mass();
      Double_t dmass = TMath::Abs(91.2-mass);
      if (dmass<zmassanyd) {
        zmassanyd = dmass;
        zmassany  = mass;
      }

      if (li->ObjType()!=lj->ObjType())
        continue;

      if (li->Is(kElectron) && (li->Charge()==lj->Charge())) {
        if (dmass<zmass2d) {
          zmass2d = dmass;
          zmass2  = mass;
        }
        continue;
      }

      if (li->Charge()==lj->Charge())
        continue;

      if (dmass<zmassd) {
        zmassd = dmass;
        zmass  = mass;
      }
      if (dmass<zmass2d) {
        zmass2d = dmass;
        zmass2  = mass;
      }
    }
  }


  UInt_t nJets     = 0;
  UInt_t nCenJets  = 0;
  UInt_t nBJets    = 0;
  UInt_t nCenBJets = 0;
  Double_t ptjmax  = -1; 
  Double_t ptjsum  = 0; 
  Double_t ptcjmax = -1; 
  Double_t ptcjsum = 0; 

  UInt_t N = jets->GetEntries();
  for (UInt_t j=0; j<N; ++j) {

    FourVector cj(jets->At(j)->Mom());
    cj*=jets->At(j)->L2RelativeCorrectionScale();
    cj*=jets->At(j)->L3AbsoluteCorrectionScale();

    Double_t pt = cj.Pt();
    if (pt<fMinJetPt)
      continue;

    ++nJets;
    ptjsum+=pt;
    if (ptjmax<pt)
      ptjmax=pt;

    Bool_t bjet = (jets->At(j)->CombinedSecondaryVertexBJetTagsDisc()>0.5);
    if (bjet)
      ++nBJets;

    if (jets->At(j)->AbsEta() < 2.5) {
      ++nCenJets;
      ptcjsum+=pt;
      if (ptcjmax<pt)
        ptcjmax=pt;
      if (bjet)
        ++nCenBJets;

      continue;
    }
  }

  Float_t nvals[255];
  Int_t nc = 0;
  nvals[nc++] = fWeight;
  nvals[nc++] = lp1->Pt();
  nvals[nc++] = lp2->Pt();
  nvals[nc++] = lp3->Pt();
  nvals[nc++] = lp1->Phi();
  nvals[nc++] = lp2->Phi();
  nvals[nc++] = lp3->Phi();
  nvals[nc++] = lp1->Eta();
  nvals[nc++] = lp2->Eta();
  nvals[nc++] = lp3->Eta();
  nvals[nc++] = lp1->ObjType()*lp1->Charge();
  nvals[nc++] = lp2->ObjType()*lp2->Charge();
  nvals[nc++] = lp3->ObjType()*lp3->Charge();
  nvals[nc++] = sumch;
  nvals[nc++] = nleps;
  nvals[nc++] = dil.Mass();
  nvals[nc++] = TMath::Sqrt(TMath::Abs(2*lp3->Pt()*met->Pt()*
                                       (1-TMath::Cos(MathUtils::DeltaPhi(lp3->Phi(),met->Phi())))));
  nvals[nc++] = MathUtils::DeltaPhi(lp1->Phi(),lp2->Phi());
  nvals[nc++] = TMath::Sqrt(TMath::Abs(lp1->Eta()-lp2->Eta()));
  nvals[nc++] = MathUtils::DeltaR(lp1->Mom(),lp2->Mom());
  nvals[nc++] = zmass;
  nvals[nc++] = zmass2;
  nvals[nc++] = zmassany;
  nvals[nc++] = met->Pt();
  nvals[nc++] = met->Phi();
  nvals[nc++] = nJets;
  nvals[nc++] = nCenJets;
  nvals[nc++] = ptjmax;
  nvals[nc++] = ptjsum;
  nvals[nc++] = ptcjmax;
  nvals[nc++] = ptcjsum;
  nvals[nc++] = nBJets;
  nvals[nc++] = nCenBJets;

  fNtuple->Fill(nvals);
}

//--------------------------------------------------------------------------------------------------
void HwwwNtupMod::SlaveBegin()
{
  // Create histograms.

  fNtuple = new TNtuple("3lntup", "3lntup", 
                        "weight:"
                        "lpt1:lpt2:lpt3:lphi1:lphi2:lphi3:leta1:leta2:leta3:"
                        "lid1:lid2:lid3:sumch:nleps:lzml12:lwml3:dlphi12:dleta12:dlR12:"
                        "zmass:zmassee:zmassany:"
                        "met:mphi:njs:ncjs:ptjmax:ptjsum:ptcjmax:ptcjsum:nbjs:ncbjs");
  fNtuple->SetDirectory(0);
  AddOutput(fNtuple);
}
