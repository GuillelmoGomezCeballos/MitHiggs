// $Id: rootlogon.C,v 1.3 2009/01/26 14:32:03 loizides Exp $
{
  gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");

  loadLibraries("libMitPhysics*.so");
  loadLibraries("libMitHiggs*.so");

  Info("rootlogon.C", "Loading xs.dat into xstab");
  mithep::SimpleTable xstab("$CMSSW_BASE/src/MitPhysics/data/xs.dat");
}
