 // $Id: HwwEvtPreSelMod.cc,v 1.4 2009/06/15 15:00:20 loizides Exp $

#include "MitHiggs/HwwMods/interface/HwwEvtPreSelMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/Names.h"

using namespace mithep;

ClassImp(mithep::HwwEvtPreSelMod)

//--------------------------------------------------------------------------------------------------
HwwEvtPreSelMod::HwwEvtPreSelMod(const char *name, const char *title) : 
  BaseSelMod(name,title),
  fMuonName(Names::gkMuonBrn),
  fElectronName(Names::gkElectronBrn),
  fNLeptonsMin(2),
  fLeptonMinPt(5),
  fLeptonMinMaxPt(20),
  fLoadBranch(kTRUE),
  fMuons(0),
  fElectrons(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HwwEvtPreSelMod::Process()
{
  // Process entries of the tree. 
  
  Int_t    nLeptons = 0;
  Double_t maxLeptonPt = 0;

  // deal with muons
  if (fLoadBranch) 
    LoadBranch(fMuonName);
  else 
    fMuons = GetObjThisEvt<MuonCol>(fMuonName);

  UInt_t ents = fMuons->GetEntries();
  for (UInt_t i=0; i<ents; ++i) {
    Double_t pt = fMuons->At(i)->Pt();
    if (pt < fLeptonMinPt) 
      continue;
    ++nLeptons;
    if (pt > maxLeptonPt)       
      maxLeptonPt = pt;
  }

  if (maxLeptonPt > fLeptonMinMaxPt &&  nLeptons >= fNLeptonsMin) {
    IncNEventsProcessed();
    return;
  }

  // deal with electrons
  if (fLoadBranch)
    LoadBranch(fElectronName);
  else 
    fElectrons = GetObjThisEvt<ElectronCol>(fElectronName);
  ents = fElectrons->GetEntries();
  for (UInt_t i=0; i<ents; ++i) {
    Double_t pt = fElectrons->At(i)->Pt();
    if (pt < fLeptonMinPt)
      continue;
    ++nLeptons;
    if (pt > maxLeptonPt)
      maxLeptonPt = pt;
  }

  if (maxLeptonPt > fLeptonMinMaxPt &&  nLeptons >= fNLeptonsMin) {
    IncNEventsProcessed();
    return;
  }

  SkipEvent();
  return;
}

//--------------------------------------------------------------------------------------------------
void HwwEvtPreSelMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the electron and muon branches.

  if (fLoadBranch) {
    ReqBranch(fMuonName,     fMuons);
    ReqBranch(fElectronName, fElectrons);
  }
}

//--------------------------------------------------------------------------------------------------
void HwwEvtPreSelMod::SlaveTerminate()
{
  // Fill event histogram.

  SaveNEventsProcessed();
}
