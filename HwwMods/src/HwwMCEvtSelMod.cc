// $Id: HwwMCEvtSelMod.cc,v 1.3 2009/06/15 15:00:20 loizides Exp $

#include "MitHiggs/HwwMods/interface/HwwMCEvtSelMod.h"
#include "MitAna/DataTree/interface/Names.h"
#include <TH1D.h>

using namespace mithep;

ClassImp(mithep::HwwMCEvtSelMod)

//--------------------------------------------------------------------------------------------------
HwwMCEvtSelMod::HwwMCEvtSelMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMCPartName(Names::gkMCPartBrn),
  fMCParticles(0),
  fHiggsPtHist(0),
  fHiggsEtaHist(0),
  fWWPtHist(0),
  fWWEtaHist(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HwwMCEvtSelMod::Process()
{
  // Process entries of the tree. For this module, we just load
  // the branch and fill the histograms.

  LoadBranch(GetMCPartName());

  fHiggsPart = 0;
  fHiggsMo1  = 0;
  fHiggsMo2  = 0;         
  fHiggsDo1  = 0;     
  fHiggsDo2  = 0;     
  fLepton1   = 0;
  fLepton2   = 0;
  fNeutrino1 = 0;
  fNeutrino2 = 0;
  fProdType  = kUnknown;
  fDecayType = kUnknownDecay;
  fMissingEnergy.SetXYZT(0,0,0,0);

  // search for higgs
  UInt_t ents=fMCParticles->GetEntries();
  for(UInt_t i=0;i<ents;++i) {
    const MCParticle *p = fMCParticles->At(i);
    if (p->Status()<3) 
      continue;

    if (p->PdgId()!=25)
      continue;

    fHiggsPart = p;
    break;
  }

  if (!fHiggsPart) {
    // todo errros
  }

  if (fHiggsPart->NDaughters()!=2) {
    // todo errros
  }

  // deal with higgs daughters
  fHiggsDo1 = fHiggsPart->Daughter(0);
  fHiggsDo2 = fHiggsPart->Daughter(1);

  if ((fHiggsDo1->NDaughters() !=2) || 
      (fHiggsDo2->NDaughters() !=2)) {
    // todo errors
  }
        
  // deal with leptons
  const MCParticle *l1 = fHiggsDo1->Daughter(0);
  if ((l1->AbsPdgId()==11) || 
      (l1->AbsPdgId()==13) || 
      (l1->AbsPdgId()==15)) {
    fLepton1   = fHiggsDo1->Daughter(0);
    fNeutrino1 = fHiggsDo1->Daughter(1);
  } else {
    fLepton1   = fHiggsDo1->Daughter(1);
    fNeutrino1 = fHiggsDo1->Daughter(0);
  }
          
  const MCParticle *l2 = fHiggsDo2->Daughter(0);
  if ((l2->AbsPdgId()==11) || 
      (l2->AbsPdgId()==13) || 
      (l2->AbsPdgId()==15)) {
    fLepton2   = fHiggsDo2->Daughter(0);
    fNeutrino2 = fHiggsDo2->Daughter(1);
  } else {
    fLepton2   = fHiggsDo2->Daughter(1);
    fNeutrino2 = fHiggsDo2->Daughter(0);
  }

  // add up missing energy
  fMissingEnergy += fNeutrino1->Mom();
  fMissingEnergy += fNeutrino2->Mom();
    
  // search for higgs mother
  for (UInt_t i=0;i<ents;++i) {
    const MCParticle *p = fMCParticles->At(i);
    if (p->Status()<3) 
      continue;
    if (!p->HasDaughter(fHiggsPart)) 
      continue;
    if (fHiggsMo1==0) {
      fHiggsMo1 = p;
      continue;
    }
    if (fHiggsMo2==0) {
      fHiggsMo2 = p;
      break;
    }
  }
  
  // assign prod type
  if (fHiggsMo1->PdgId()==21 && fHiggsMo2->PdgId()==21) 
    fProdType = kGluGlu;
  else fProdType = kVBF;

  // assign decay type
  if ((fLepton1->AbsPdgId() == MCParticle::kEl) && 
      (fLepton2->AbsPdgId() == MCParticle::kEl))
    fDecayType = kElEl;
  else if ((fLepton1->AbsPdgId() == MCParticle::kEl) && 
           (fLepton2->AbsPdgId() == MCParticle::kMu))
    fDecayType = kElMu;
  else if ((fLepton2->AbsPdgId() == MCParticle::kEl) && 
           (fLepton1->AbsPdgId() == MCParticle::kMu)) {
    fDecayType = kElMu;
    std::swap(fLepton1,fLepton2);
  } else if ((fLepton1->AbsPdgId() == MCParticle::kEl) && 
             (fLepton2->AbsPdgId() == MCParticle::kTau))
    fDecayType = kElTauHad;
  else if ((fLepton2->AbsPdgId() == MCParticle::kEl) && 
           (fLepton1->AbsPdgId() == MCParticle::kTau)) {
    fDecayType = kElTauHad;
    std::swap(fLepton1,fLepton2);
  } else if ((fLepton1->AbsPdgId() == MCParticle::kMu) && 
             (fLepton2->AbsPdgId() == MCParticle::kMu))
    fDecayType = kMuMu;
  else if ((fLepton1->AbsPdgId() == MCParticle::kMu) && 
           (fLepton2->AbsPdgId() == MCParticle::kTau))
    fDecayType = kMuTauHad;
  else if ((fLepton2->AbsPdgId() == MCParticle::kMu) && 
           (fLepton1->AbsPdgId() == MCParticle::kTau)) {
    fDecayType = kMuTauHad;
    std::swap(fLepton1,fLepton2);
  } else if ((fLepton1->AbsPdgId() == MCParticle::kTau) && 
             (fLepton2->AbsPdgId() == MCParticle::kTau)) {
    fDecayType = kTauHadTauHad;
    
  }
  else {
    //todo errors
  }

  MCParticle::EPartType type1 = ClassifyTau(fLepton1);
  MCParticle::EPartType type2 = ClassifyTau(fLepton2);

  if ((fDecayType==kElTauHad) && (type2==MCParticle::kEl))
    fDecayType = kElTauEl;
  else if ((fDecayType==kElTauHad) && (type2==MCParticle::kMu))
    fDecayType = kElTauMu;
  else if ((fDecayType==kMuTauHad) && (type2==MCParticle::kEl))
    fDecayType = kMuTauEl;
  else if ((fDecayType==kMuTauHad) && (type2==MCParticle::kMu))
    fDecayType = kMuTauMu;
  else if (fDecayType == kTauHadTauHad) {
    if ((type1==MCParticle::kUnknown) && (type2!=MCParticle::kUnknown)) {
      std::swap(type1,type2);
      std::swap(fLepton1,fLepton2);
    }
    if ((type1==MCParticle::kEl) && (type2==MCParticle::kEl)) 
      fDecayType = kTauElTauEl;
    else if ((type1==MCParticle::kEl) && (type2==MCParticle::kMu)) 
      fDecayType = kTauElTauMu;
    else if ((type1==MCParticle::kMu) && (type2==MCParticle::kEl)) {
      fDecayType = kTauElTauMu;
      std::swap(fLepton1,fLepton2);
    } else if ((type1==MCParticle::kMu) && (type2==MCParticle::kMu)) 
      fDecayType = kTauMuTauMu;
    else if ((type1==MCParticle::kEl) && (type2==MCParticle::kUnknown)) 
      fDecayType = kTauElTauHad;
    else if ((type1==MCParticle::kMu) && (type2==MCParticle::kUnknown)) 
      fDecayType = kTauMuTauHad;
  }

}

#if 0
     if (p1->Status()>2) continue;
     if (p1->PdgId()!=25) continue;

     fHiggsPart = p1;
     break;
  }

     p1->Print();
     if (p1->HasMother()) {
       cout << "Mother: ";  p1->Mother()->Print();
     }
     for (UInt_t i=0;i<p1->NDaughters();++i) {
       if (p1->Daughter(i)->Status()<3) continue;
       cout << "Do: " << i << ": "; p1->Daughter(i)->Print();
     }
     cout << endl;
  }

  if (fHiggsPart) {
    const MCParticle *h = fHiggsPart;
    if (fHiggsPart->NDaughters()==0) 
      h = fHiggsPart->Mother();

    for (UInt_t i=0;i<h->NDaughters();++i) {
      const MCParticle *d = h->Daughter(i);
      if (d->AbsPdgId()==24) {
        cout << d->NDaughters() << endl;
      }
    }
  }



}
#endif

#if 0
     p1->Print();
     if(i>30) break;


       if(p1->NDaughters()>0)
         p1->Daughter(0)->Print();
       if(p1->NDaughters()>1)
         p1->Daughter(1)->Print();

       fHiggsPtHist->Fill(p1->Pt());
       fHiggsEtaHist->Fill(p1->Eta());
       continue;
     }

     if (TMath::Abs(p1->PdgId())==24) {
       if (p1->Mother()->Mother()->PdgId()==25) {
         fHWWPtHist->Fill(p1->Pt());
         fHWWEtaHist->Fill(p1->Eta());
       } else {
         fWWPtHist->Fill(p1->Pt());
         fWWEtaHist->Fill(p1->Eta());
       }

       for (UInt_t j=i+1;j<ents;++j) {
         const MCParticle *p2 = fMCParticles->At(j);
         if (p2->Status()>2) continue;
         if (TMath::Abs(p2->PdgId())!=24) continue;
         if (p1==p2) continue;

         FourVector mom(p1->Mom()+p2->Mom());
         fHiggsMassHist->Fill(TMath::Sqrt(mom.M2()));
         fHiggsTransMHist->Fill(TMath::Sqrt(mom.Mt2()));
       }
     }
  }
#endif

//--------------------------------------------------------------------------------------------------
MCParticle::EPartType HwwMCEvtSelMod::ClassifyTau(const MCParticle *tau)
{
  // Return type of tau decay.

  if (!tau || (tau->AbsPdgId() != MCParticle::kTau))
    return MCParticle::kUnknown;

  const MCParticle *mytau = 0;
  UInt_t ents=fMCParticles->GetEntries();
  for (UInt_t i=0;i<ents;++i) {
    const MCParticle *p = fMCParticles->At(i);

    if (p->Status()>2) 
      continue;
    if (p->AbsPdgId()!=MCParticle::kTau)
      continue;
    if (p->Mother()!=tau) 
      continue;

    mytau = p;
  }
  
  UInt_t neutrinos = 0;
  for (UInt_t i=0; i<mytau->NDaughters(); ++i) {
    const MCParticle *p = mytau->Daughter(i);
    if (!p->IsNeutrino()) 
      continue;
    fMissingEnergy += p->Mom();
    neutrinos += 2;
  }

  if ((neutrinos!=2) || (mytau->NDaughters()!=3))
    return MCParticle::kUnknown;
    
  for (UInt_t i=0; i<mytau->NDaughters(); ++i) {
    const MCParticle *p = mytau->Daughter(i);
    if (p->IsNeutrino()) 
      continue;

    if (p->AbsPdgId()==MCParticle::kEl) 
      return MCParticle::kEl;
    else if (p->AbsPdgId()==MCParticle::kMu) 
      return MCParticle::kMu;
    
    break;
  }

  return MCParticle::kUnknown;
}      

//--------------------------------------------------------------------------------------------------
void HwwMCEvtSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual
  // analysis. Here, we typically initialize histograms and 
  // other analysis objects and request branches. For this module,
  // we request a branch of the MitTree.

  ReqBranch(GetMCPartName(), fMCParticles);

  fHiggsPtHist  = new TH1D("hHiggsPt",";p_{t} [GeV/c];#",100,0.,40.);
  AddOutput(fHiggsPtHist);
  fHiggsEtaHist = new TH1D("hHiggsEta",";#eta;#",160,-8.,8.);
  AddOutput(fHiggsEtaHist);
  fHWWPtHist  = new TH1D("hHWWPt",";p_{t} [GeV/c];#",100,0.,40.);
  AddOutput(fHWWPtHist);
  fHWWEtaHist = new TH1D("hHWWEta",";#eta;#",160,-8.,8.);
  AddOutput(fHWWEtaHist);
  fWWPtHist  = new TH1D("hWWPt",";p_{t} [GeV/c];#",100,0.,40.);
  AddOutput(fWWPtHist);
  fWWEtaHist = new TH1D("hWWEta",";#eta;#",160,-8.,8.);
  AddOutput(fWWEtaHist);
  fHiggsMassHist = new TH1D("hHiggsMass",";m [GeV/c];#",200,0.,200.);
  AddOutput(fHiggsMassHist);
  fHiggsTransMHist = new TH1D("hHiggsTransM",";m [GeV/c];#",200,0.,200.);
  AddOutput(fHiggsTransMHist);
}

//--------------------------------------------------------------------------------------------------
void HwwMCEvtSelMod::Terminate()
{
  // Run finishing code on the client computer. 
  // For this module, we dont do anything here.
}
