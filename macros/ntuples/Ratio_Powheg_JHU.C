#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <utility>
#include "TFile.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TH1D.h"
#include <TROOT.h>
#include <TChain.h>

void Ratio_Powheg_JHU(TString  InputFilename0 = "/data/smurf/ceballos/distributions/Powheg_JHU/hwwof_0j.input_8TeV_powheg.root",
                      TString  InputFilename1 = "/data/smurf/ceballos/distributions/Powheg_JHU/hwwof_0j.input_8TeV_0p.root",
                      TString  InputFilename2 = "/data/smurf/ceballos/distributions/Powheg_JHU/hwwof_0j.input_8TeV_2p.root",
		      Int_t Ecm = 8, Int_t NJetsType = 0
  ) {

  TH1D *hDNorm_Powheg_JHU = new TH1D("hDNorm_Powheg_JHU","hDNorm_Powheg_JHU",1,0.0,1.0);

  TFile *fPowheg = new TFile(InputFilename0.Data(), "READ");
  TH1D *hDPowheg = (TH1D*)(fPowheg->Get("histo_ggH"));
  assert(hDPowheg);
  hDPowheg->SetDirectory(0);
  TH1D *hDZH  = (TH1D*)(fPowheg->Get("histo_ZH"));
  TH1D *hDWH  = (TH1D*)(fPowheg->Get("histo_WH"));
  TH1D *hDqqH = (TH1D*)(fPowheg->Get("histo_qqH"));
  TH1D *hDggH = (TH1D*)(fPowheg->Get("histo_ggH"));
  assert(hDZH );
  assert(hDWH );
  assert(hDqqH);
  assert(hDggH);
  hDZH ->SetDirectory(0);
  hDWH ->SetDirectory(0);
  hDqqH->SetDirectory(0);
  hDZH ->Scale(0); // not used now
  hDWH ->Scale(0); // not used now
  hDqqH->Scale(0); // not used now
  hDggH->SetDirectory(0);
  fPowheg->Close();
  delete fPowheg;
  double nZH  = hDZH ->GetSumOfWeights();
  double nWH  = hDWH ->GetSumOfWeights();
  double nqqH = hDqqH->GetSumOfWeights();
  double nggH = hDggH->GetSumOfWeights();
  hDPowheg->Scale(1./hDPowheg->GetSumOfWeights());

  TFile *fJHU0 = new TFile(InputFilename1.Data(), "READ");
  TH1D *hDJHU0 = (TH1D*)(fJHU0->Get("histo_ggH"));
  assert(hDJHU0);
  hDJHU0->SetDirectory(0);
  fJHU0->Close();
  delete fJHU0;
  double n1 = hDJHU0->GetSumOfWeights();
  hDJHU0->Scale(1./hDJHU0->GetSumOfWeights());
  
  TFile *fJHU2 = new TFile(InputFilename2.Data(), "READ");
  TH1D *hDJHU2 = (TH1D*)(fJHU2->Get("histo_ggH"));
  assert(hDJHU2);
  hDJHU2->SetDirectory(0);
  fJHU2->Close();
  delete fJHU2;
  double n2 = hDJHU2->GetSumOfWeights();
  hDJHU2->Scale(1./hDJHU2->GetSumOfWeights());
  
  printf("poweg(%f+%f+%f+%f=%f)/JHU0(%f)/JHU2(%f)==>%f\n",nZH,nWH,nqqH,nggH,(nZH+nWH+nqqH+nggH),n1,n2,(nZH+nWH+nqqH+nggH)/n2);

  hDNorm_Powheg_JHU->Fill(0.5,(nZH+nWH+nqqH+nggH)/n2);
  hDPowheg->Divide(hDJHU0);
  for(int i=1; i<=hDPowheg->GetNbinsX(); i++){
    if(hDPowheg->GetBinContent(i) == 0) hDPowheg->SetBinContent(i,1.0);
  }

  TFile *fRatio = new TFile(Form("ratio_Powheg_JHU_%dTeV_%dj.root",Ecm,NJetsType), "RECREATE");
  hDPowheg->Write();
  hDNorm_Powheg_JHU->Write();
  fRatio->Close();

}
