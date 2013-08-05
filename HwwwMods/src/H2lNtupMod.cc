// $Id: H2lNtupMod.cc,v 1.3 2009/07/11 20:31:52 loizides Exp $

#include "MitHiggs/HwwwMods/interface/H2lNtupMod.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CompoundParticle.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TNtuple.h>

using namespace mithep;

ClassImp(mithep::H2lNtupMod)

//--------------------------------------------------------------------------------------------------
H2lNtupMod::H2lNtupMod(const char *name, const char *title) : 
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
void H2lNtupMod::Process()
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

  // make sure we have at least 2 leptons
  if (leptons->GetEntries()<2)
    return;

  // make sure the 3rd highest pt lepton has pt <= fMinPt.
  if (leptons->GetEntries() >= 3 && (leptons->At(2)->Pt()>fMinPt))
    return;

  // make sure we have two same sign leading leptons
  const Particle *l0 = leptons->At(0);
  const Particle *l1 = leptons->At(1);
  if (l0->Charge()!=l1->Charge())
    return;

  const Met *met = mets->At(0);

  // hww preselection
  if (l0->Pt()   < fMinMaxPt   ||
      l1->Pt()   < fMinPt      ||
      met->Pt()  < fMinMissEt)
    return;

  CompoundParticle dil;
  dil.AddDaughter(l0);
  dil.AddDaughter(l1);

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

  Double_t sumch=0;
  Int_t nleps = leptons->GetEntries();
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

  Float_t nvals[255];
  Int_t nc = 0;
  nvals[nc++] = fWeight;
  nvals[nc++] = l0->Pt();
  nvals[nc++] = l1->Pt();
  nvals[nc++] = l0->Phi();
  nvals[nc++] = l1->Phi();
  nvals[nc++] = l0->Eta();
  nvals[nc++] = l1->Eta();
  nvals[nc++] = l0->ObjType()*l0->Charge();
  nvals[nc++] = l1->ObjType()*l1->Charge();
  nvals[nc++] = sumch;
  nvals[nc++] = nleps;
  nvals[nc++] = dil.Mass();
  nvals[nc++] = MathUtils::DeltaPhi(l0->Phi(),l1->Phi());
  nvals[nc++] = TMath::Sqrt(TMath::Abs(l0->Eta()-l1->Eta()));
  nvals[nc++] = MathUtils::DeltaR(l0->Mom(),l1->Mom());
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
void H2lNtupMod::SlaveBegin()
{
  // Create histograms.

  fNtuple = new TNtuple("2sslntup", "2sslntup", 
                        "weight:"
                        "lpt1:lpt2:lphi1:lphi2:leta1:leta2:"
                        "lid1:lid2:sumch:nleps:lzml12:dlphi12:dleta12:dlR12:"
                        "zmass:zmassee:zmassany:"
                        "met:mphi:njs:ncjs:ptjmax:ptjsum:ptcjmax:ptcjsum:nbjs:ncbjs");
  fNtuple->SetDirectory(0);
  AddOutput(fNtuple);
}
