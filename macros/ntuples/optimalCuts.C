#include "MitNtupleEvent.C"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include "TLegend.h"
#include "TPaveText.h"

double scpFast(double sig, double bkg, double sigma_b, double delta_b);

TTree* getTreeFromFile(const char* infname, const char* tname);

// To recover from splitting signal into 2 separate halves
// (training/validation) I double the weight of the signal
// events by hand

const double gFrac        = 2.0; // 1/2 for trainning, 1/2 for test
const int    verboseLevel =   1;
const int    jetRange     =   2; // not used yet
const int    masszRange    =   9; // not used yet
const int    nDecays      =  25;

//------------------------------------------------------------------------------
// optimalCuts
//------------------------------------------------------------------------------
void optimalCuts
(
 int     mHiggs  	 = 25,
 Char_t xTitle[]="myX", Char_t yTitle[]="Fraction",
 TString signalInputFile = "ntuples/inputNtuple-data-standard-histo_H170_WW2l_all.root",
 TString bgdInputFile    = "ntuples/inputNtuple-data-standard-histo_HBCK2.root"
 )
{

  double lumi = 0.2;
  // Recover both the original ntuples and the tmva results for signal
  TTree* sig1 = getTreeFromFile(signalInputFile,"all");
  assert(sig1);

  // Recover both the original ntuples and the tmva results for background
  TTree* bgd1 = getTreeFromFile(bgdInputFile,"all");
  assert(bgd1);

  int channel = mHiggs;

  const int nBin = 180;
  double S0[nBin],S1[nBin];
  double B0[nBin],B1[nBin];
  for(int i=0; i<nBin; i++){
    S0[i] = 0.0; S1[i] = 0.0;
    B0[i] = 0.0; B1[i] = 0.0;
  }
  TH1D* hDSignif[2];
  hDSignif[0] = new TH1D("hDSignif_0", "hDSignif_0", nBin, -0.5, nBin-0.5);
  hDSignif[1] = new TH1D("hDSignif_1", "hDSignif_1", nBin, -0.5, nBin-0.5);
  TH1D* hDSigOpt = new TH1D("hDSigOpt", "hDSigOpt", nBin, 0, 1);
  TH1D* hDBckOpt = new TH1D("hDBckOpt", "hDBckOpt", nBin, 0, 1);

  double bgdDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  double weiDecay[45] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
                         
  MitNtupleEvent bgdEvent(bgd1);
  int nBgd=bgd1->GetEntries();
  for (int i=0; i<nBgd; ++i) {

    if (i%5000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,nBgd);
    bgdEvent.GetEntry(i);

    if(channel == 25){ // WW selection
      int charge = (int)(bgdEvent.H_q1 + bgdEvent.H_q2);
      if(
         bgdEvent.H_dim01   > 12/500. &&
         bgdEvent.H_nlep == 2. &&
         charge == 0 &&
         bgdEvent.H_njets == 0. &&
         bgdEvent.H_nmuons == 0. &&
         bgdEvent.H_ntracks <= 4. &&
         bgdEvent.H_pt1 > 20./400. &&
         bgdEvent.H_pt2 > 10./400. &&
	 (bgdEvent.H_pt2 > 20./400. || bgdEvent.H_ltype == 10.) &&
         bgdEvent.H_met > 20./400. &&
         (bgdEvent.H_met > 46./400. || bgdEvent.H_ltype == 12. || bgdEvent.H_ltype == 13.) &&
	 bgdEvent.H_deltaphilmet > 60./180. &&
         (fabs(bgdEvent.H_dim01*500.-91.1876) > 15. || bgdEvent.H_ltype == 12. || bgdEvent.H_ltype == 13.) && 
	 bgdEvent.H_met/bgdEvent.H_dilpt > 0.6 &&
	 //bgdEvent.H_dim01   > 80/500. &&
	 1 == 1
	){
	double myVar = bgdEvent.H_met;
	myVar = bgdEvent.H_deltaphilmet;
	myVar = TMath::Max(bgdEvent.H_ntracks/20.0,0.001);
	//myVar = TMath::Max(bgdEvent.H_maxmbtag*1.0,0.001);
	myVar = TMath::Min(bgdEvent.H_met/bgdEvent.H_dilpt/2.,0.999);
	for(int n1=0; n1<nBin; n1++){
	  if(myVar > 1.0*n1/nBin){
	    if(bgdEvent.H_decay == 25 || bgdEvent.H_decay == 30) S0[n1] = S0[n1] + bgdEvent.H_weight*lumi;
	    if(bgdEvent.H_decay != 25 && bgdEvent.H_decay != 30) B0[n1] = B0[n1] + bgdEvent.H_weight*lumi;
	  }
	  if(myVar < 1.0*n1/nBin){
	    if(bgdEvent.H_decay == 25 || bgdEvent.H_decay == 30) S1[n1] = S1[n1] + bgdEvent.H_weight*lumi;
	    if(bgdEvent.H_decay != 25 && bgdEvent.H_decay != 30) B1[n1] = B1[n1] + bgdEvent.H_weight*lumi;
	  }
	}
        if(bgdEvent.H_decay == 25) hDSigOpt->Fill(myVar,bgdEvent.H_weight*lumi);
        if(bgdEvent.H_decay != 25) hDBckOpt->Fill(myVar,bgdEvent.H_weight*lumi);
      }
    } // WW selection

    if(channel == 5 || channel == 4){ // ttbar selection
      int charge = (int)(bgdEvent.H_q1 + bgdEvent.H_q2);
      if(
         bgdEvent.H_dim01   > 12./500. &&
         bgdEvent.H_nlep == 2. &&
         charge == 0 &&
         bgdEvent.H_pt1 > 20./400. &&
         bgdEvent.H_pt2 > 20./400. &&
	 bgdEvent.H_met > 45./400 &&
         (bgdEvent.H_met > 55./400. || bgdEvent.H_ltype == 12. || bgdEvent.H_ltype == 13.) &&
	 bgdEvent.H_deltaphihmet > 20./180. &&
         (fabs(bgdEvent.H_dim01*500.-91.1876) > 15 || bgdEvent.H_ltype == 12. || bgdEvent.H_ltype == 13.) &&
	 1 == 1
	){
	double myVar = bgdEvent.H_deltaphihmet;
	//myVar = bgdEvent.H_deltaphihmet; 
	//myVar = TMath::Min(fabs(bgdEvent.H_dim01*500-91.1876),99.999)/100.;
	myVar = 0.1+bgdEvent.H_njets/10.;
	//myVar = bgdEvent.H_pt2;
	//myVar = bgdEvent.H_met;
        //myVar = TMath::Min(bgdEvent.H_met/bgdEvent.H_dilpt/2.,0.999);
	for(int n1=0; n1<nBin; n1++){
	  if(myVar > 1.0*n1/nBin){
	    if(bgdEvent.H_decay == 5||bgdEvent.H_decay == 4) S0[n1] = S0[n1] + bgdEvent.H_weight*lumi;
	    if(bgdEvent.H_decay != 5&&bgdEvent.H_decay != 4) B0[n1] = B0[n1] + bgdEvent.H_weight*lumi;
	  }
	  if(myVar < 1.0*n1/nBin){
	    if(bgdEvent.H_decay == 5||bgdEvent.H_decay == 4) S1[n1] = S1[n1] + bgdEvent.H_weight*lumi;
	    if(bgdEvent.H_decay != 5&&bgdEvent.H_decay != 4) B1[n1] = B1[n1] + bgdEvent.H_weight*lumi;
	  }
	}
        if(bgdEvent.H_decay == 5||bgdEvent.H_decay == 4) hDSigOpt->Fill(myVar,bgdEvent.H_weight*lumi);
        if(bgdEvent.H_decay != 5&&bgdEvent.H_decay != 4) hDBckOpt->Fill(myVar,bgdEvent.H_weight*lumi);
      }
    } // ttbar selection

    if(channel == 26){ // WZ selection
      if(
         bgdEvent.H_nlep == 3 &&
         bgdEvent.H_ltype >= 0 &&
         bgdEvent.H_ltype <= 3 &&
         //bgdEvent.H_pt3 > 15/400. &&
         bgdEvent.H_mtw3 > 15/400. &&
         bgdEvent.H_met > 20/400. &&
         fabs(bgdEvent.H_massz*500-91.1876) < 20 &&
         //(bgdEvent.H_deltaphijetmet < 160./180. || bgdEvent.H_njets == 0 ) &&
	 (bgdEvent.H_ltype == 1 || bgdEvent.H_ltype == 2 || bgdEvent.H_met > 25/400.) &&
	 1 == 1
	){
	double add = 1.;
	if(bgdInputFile == "ntuples/inputNtuple-validate-data-standard-histo_HBCK3.root"){
	  if(bgdEvent.H_decay == 26 || bgdEvent.H_decay == 27) add = 0.87;
	  if(bgdEvent.H_decay == 28 || bgdEvent.H_decay == 29) add = 0.68;
	}

	bgdDecay[(int)bgdEvent.H_decay] += bgdEvent.H_weight*lumi*add;
	if(weiDecay[(int)bgdEvent.H_decay] == 0) weiDecay[(int)bgdEvent.H_decay] = bgdEvent.H_weight*lumi*add;

	double myVar = bgdEvent.H_deltaphijetmet;
	myVar = TMath::Min(fabs(bgdEvent.H_massz*500-91.1876),99.999)/100;
	myVar =bgdEvent.H_mtw3*400;
	//myVar =bgdEvent.H_met*400;
	myVar = TMath::Max(bgdEvent.H_njets/10.,0.001);
	//myVar = bgdEvent.H_massz*500;
	myVar = bgdEvent.H_deltaphil3met;
	for(int n1=0; n1<nBin; n1++){
	  if(myVar > 1.0*n1/nBin){
	    if(bgdEvent.H_decay == 26 || bgdEvent.H_decay == 27) S0[n1] = S0[n1] + bgdEvent.H_weight*lumi*add;
	    if(bgdEvent.H_decay != 26 && bgdEvent.H_decay != 27) B0[n1] = B0[n1] + bgdEvent.H_weight*lumi*add;
	  }
	  if(myVar < 1.0*n1/nBin){
	    if(bgdEvent.H_decay == 26 || bgdEvent.H_decay == 27) S1[n1] = S1[n1] + bgdEvent.H_weight*lumi*add;
	    if(bgdEvent.H_decay != 26 && bgdEvent.H_decay != 27) B1[n1] = B1[n1] + bgdEvent.H_weight*lumi*add;
	  }
	}
        if(bgdEvent.H_decay == 26 || bgdEvent.H_decay == 27) hDSigOpt->Fill(myVar,bgdEvent.H_weight*lumi);
        if(bgdEvent.H_decay != 26 && bgdEvent.H_decay != 27) hDBckOpt->Fill(myVar,bgdEvent.H_weight*lumi);
      }
    } // WZ selection

    if(channel == 1000){ // HW->2l selection
      int charge = (int)(bgdEvent.H_q1 + bgdEvent.H_q2);
      if(
         bgdEvent.H_dim01   > 12/500. &&
         bgdEvent.H_dim01   < 300./500. &&
         bgdEvent.H_nlep == 2 &&
         charge != 0 &&
         bgdEvent.H_pt1 > 20/400. &&
         bgdEvent.H_pt2 > 10/400. &&
         (bgdEvent.H_met > 50/400. ||  bgdEvent.H_ltype == 12. || bgdEvent.H_ltype == 13.) &&
         (bgdEvent.H_met > 40/400. || (bgdEvent.H_ltype != 12. && bgdEvent.H_ltype != 13.)) &&
         (fabs(bgdEvent.H_dim01*500-91.1876) > 15 || bgdEvent.H_ltype != 11.) &&
         bgdEvent.H_mtw2 > 10/400. &&
	 bgdEvent.H_btag1 < 0.5 &&
	 bgdEvent.H_maxmbtag < 0.5 &&
	 bgdEvent.H_nmuons == 0 &&
	 //bgdEvent.H_deltaphihmet > 20./180. &&
	 1 == 1
	){
	double add = 1.;
	if(bgdInputFile == "ntuples/inputNtuple-validate-data-standard-histo_HBCK3.root"){
	  if(bgdEvent.H_decay == 26 || bgdEvent.H_decay == 27) add = 0.87;
	  if(bgdEvent.H_decay == 28 || bgdEvent.H_decay == 29) add = 0.68;
	}

	bgdDecay[(int)bgdEvent.H_decay] += bgdEvent.H_weight*lumi*add;
	if(weiDecay[(int)bgdEvent.H_decay] == 0) weiDecay[(int)bgdEvent.H_decay] = bgdEvent.H_weight*lumi*add;

	double myVar = bgdEvent.H_deltaphihmet;
	//double myVar = TMath::Min(fabs(bgdEvent.H_massz*500-91.1876),99.999)/100;
	//myVar = bgdEvent.H_njets/10;
	//myVar = bgdEvent.H_mtw2;
	//myVar = bgdEvent.H_delphil;
	for(int n1=0; n1<nBin; n1++){
	  if(myVar > 1.0*n1/nBin){
	    B0[n1] = B0[n1] + bgdEvent.H_weight*lumi;
	  }
	  if(myVar < 1.0*n1/nBin){
	    B1[n1] = B1[n1] + bgdEvent.H_weight*lumi;
	  }
	}
        hDBckOpt->Fill(myVar,bgdEvent.H_weight*lumi);
      }
    } // HW->2l selection

    if(channel == 1100){ // HW->3l selection
        int charge = (int)(bgdEvent.H_q1 + bgdEvent.H_q2 + bgdEvent.H_q3);
	double thePtMax = 0.;
	double thePtMin = 0.;
        double thePt2   = bgdEvent.H_pt1;
	if(bgdEvent.H_nlep==3){
	  thePtMax = TMath::Max(TMath::Max(bgdEvent.H_pt1,bgdEvent.H_pt2),bgdEvent.H_pt3);
          thePtMin = TMath::Min(TMath::Min(bgdEvent.H_pt1,bgdEvent.H_pt2),bgdEvent.H_pt3);
	  if(bgdEvent.H_pt2 != thePtMax && bgdEvent.H_pt2 != thePtMin)
	    thePt2 = bgdEvent.H_pt2;
	  if(bgdEvent.H_pt3 != thePtMax && bgdEvent.H_pt3 != thePtMin)
	    thePt2 = bgdEvent.H_pt3;
	}
	      
      if(
         bgdEvent.H_nlep == 3 &&
         abs(charge) == 1 &&
         thePtMax > 20./400. &&
         //thePt2   > 15./400. &&
         //thePtMin > 15./400 &&
	 bgdEvent.H_dim01 < 100./500. &&
	 //bgdEvent.H_dim02 > 0.2 &&
         fabs(bgdEvent.H_massz*500-91.1876) > 20 &&
         bgdEvent.H_met > 35./400. &&
	 (bgdEvent.H_met > 50./400. || bgdEvent.H_dim12 < 2) &&
         //bgdEvent.H_met < 150./400. &&
	 bgdEvent.H_njets <= 1. &&
	 bgdEvent.H_btag1 < 0.5 &&
	 bgdEvent.H_maxmbtag < 0.5 &&
	 bgdEvent.H_nmuons == 0 &&
	 bgdEvent.H_mindr < 1.5 &&
	 1 == 1
	){
	double add = 1.;
	if(bgdInputFile == "ntuples/inputNtuple-validate-data-standard-histo_HBCK3.root"){
	  if(bgdEvent.H_decay == 26 || bgdEvent.H_decay == 27) add = 0.87;
	  if(bgdEvent.H_decay == 28 || bgdEvent.H_decay == 29) add = 0.68;
	}

	bgdDecay[(int)bgdEvent.H_decay] += bgdEvent.H_weight*lumi*add;
	if(weiDecay[(int)bgdEvent.H_decay] == 0) weiDecay[(int)bgdEvent.H_decay] = bgdEvent.H_weight*lumi*add;

	double myVar = TMath::Min(fabs(bgdEvent.H_massz*500-91.1876),99.999)/100;
	myVar = TMath::Min(bgdEvent.H_mindr/10.,0.999);
	//myVar = bgdEvent.H_dim02;
	myVar = thePtMin;
        //myVar = TMath::Min(TMath::Max(bgdEvent.H_maxisoe*1.0,0.001)/1.,0.999);
	//myVar = bgdEvent.H_met;
	for(int n1=0; n1<nBin; n1++){
	  if(myVar > 1.0*n1/nBin){
	    B0[n1] = B0[n1] + bgdEvent.H_weight*lumi;
	  }
	  if(myVar < 1.0*n1/nBin){
	    B1[n1] = B1[n1] + bgdEvent.H_weight*lumi;
	  }
	}
        hDBckOpt->Fill(myVar,bgdEvent.H_weight*lumi);
      }
    } // HW->3l selection

    if(channel >= 100 && channel <= 800){ // H->WW selection
      int charge = (int)(bgdEvent.H_q1 + bgdEvent.H_q2);
      if(
         bgdEvent.H_dim01   > 12/500. &&
         bgdEvent.H_dim01   < 55/500. &&
         bgdEvent.H_nlep == 2. &&
         charge == 0 &&
         bgdEvent.H_njets == 0. &&
         bgdEvent.H_nmuons == 0. &&
         bgdEvent.H_ntracks <= 4. &&
         bgdEvent.H_pt1 > 30./400. &&
         bgdEvent.H_pt2 > 25./400. &&
	 //(bgdEvent.H_pt2 > 20./400. || bgdEvent.H_ltype == 10. || bgdEvent.H_ltype == 10.) &&
          bgdEvent.H_met > 30./400. &&
         (bgdEvent.H_met > 52./400. || bgdEvent.H_ltype == 12. || bgdEvent.H_ltype == 13.) &&
	 //bgdEvent.H_deltaphilmet > 60./180. &&
         (fabs(bgdEvent.H_dim01*500.-91.1876) > 15. || bgdEvent.H_ltype == 12. || bgdEvent.H_ltype == 13.) && 
	 bgdEvent.H_maxd0sig < 0.025 && 
	 bgdEvent.H_delphil < 57./180. &&
	 //!(bgdEvent.H_ltype == 12. || bgdEvent.H_ltype == 13.) &&
	 //(bgdEvent.H_ltype == 10. || bgdEvent.H_ltype == 12.) &&
	 1 == 1
	){
	double myVar = bgdEvent.H_met;
        //myVar = bgdEvent.H_pt1;
        //myVar = bgdEvent.H_pt2;
        //myVar = bgdEvent.H_delphil;
        //myVar = bgdEvent.H_mthiggs;
	myVar = bgdEvent.H_dim01;
	//myVar = TMath::Min(bgdEvent.H_maxd0sig*1.0,0.0999)/0.1;
        //myVar = (bgdEvent.H_ntracks+0.1)/20.;
	for(int n1=0; n1<nBin; n1++){
	  if(myVar > 1.0*n1/nBin){
	    B0[n1] = B0[n1] + bgdEvent.H_weight*lumi;
	  }
	  if(myVar < 1.0*n1/nBin){
	    B1[n1] = B1[n1] + bgdEvent.H_weight*lumi;
	  }
	}
        hDBckOpt->Fill(myVar,bgdEvent.H_weight*lumi);
      }
    } // H->WW selection

  }

  if(channel == 1000 || channel == 1100 ||
     (channel >= 100 && channel <= 800)){
    MitNtupleEvent sigEvent(sig1);
    int nSig=sig1->GetEntries();
    for (int i=0; i<nSig; ++i) {

      if (i%5000 == 0 && verboseLevel > 0)
	printf("--- reading Signal event %5d of %5d\n",i,nSig);
      sigEvent.GetEntry(i);
      if(channel == 1000){ // HW->2l selection
        int charge = (int)(sigEvent.H_q1 + sigEvent.H_q2);
	if(
           sigEvent.H_dim01   > 12/500. &&
           sigEvent.H_dim01   < 300./500. &&
           sigEvent.H_nlep == 2 &&
           charge != 0 &&
           sigEvent.H_pt1 > 20/400. &&
           sigEvent.H_pt2 > 10/400. &&
           (sigEvent.H_met > 50/400. ||  sigEvent.H_ltype == 12. || sigEvent.H_ltype == 13.) &&
           (sigEvent.H_met > 40/400. || (sigEvent.H_ltype != 12. && sigEvent.H_ltype == 13.)) &&
           (fabs(sigEvent.H_dim01*500-91.1876) > 15 || sigEvent.H_ltype != 11.) &&
           sigEvent.H_mtw2 > 10/400. &&
	   1 == 1
	  ){
	  double myVar = sigEvent.H_mthiggs;
	  //double myVar = TMath::Min(fabs(sigEvent.H_massz*500-91.1876),99.999)/100;
	  //double myVar = sigEvent.H_njets/10;
	  for(int n1=0; n1<nBin; n1++){
	    if(myVar > 1.0*n1/nBin){
	      S0[n1] = S0[n1] + sigEvent.H_weight*lumi;
	    }
	    if(myVar < 1.0*n1/nBin){
	      S1[n1] = S1[n1] + sigEvent.H_weight*lumi;
	    }
	  }
          hDSigOpt->Fill(myVar,sigEvent.H_weight*lumi);
	}
      } // HW->2l selection

      if(channel == 1100){ // HW->3l selection
        int charge = (int)(sigEvent.H_q1 + sigEvent.H_q2 + sigEvent.H_q3);
	double thePtMax = 0.;
	double thePtMin = 0.;
        double thePt2   = sigEvent.H_pt1;
	if(sigEvent.H_nlep==3){
	  thePtMax = TMath::Max(TMath::Max(sigEvent.H_pt1,sigEvent.H_pt2),sigEvent.H_pt3);
          thePtMin = TMath::Min(TMath::Min(sigEvent.H_pt1,sigEvent.H_pt2),sigEvent.H_pt3);
	  if(sigEvent.H_pt2 != thePtMax && sigEvent.H_pt2 != thePtMin)
	    thePt2 = sigEvent.H_pt2;
	  if(sigEvent.H_pt3 != thePtMax && sigEvent.H_pt3 != thePtMin)
	    thePt2 = sigEvent.H_pt3;
	}

	if(
           sigEvent.H_nlep == 3 &&
           abs(charge) == 1 &&
           thePtMax > 20./400. &&
           thePt2   > 15./400. &&
           thePtMin > 15./400 &&
	   sigEvent.H_dim01 < 100./500. &&
	   sigEvent.H_dim02 > 0.2 &&
           sigEvent.H_met > 35./400. &&
	   (sigEvent.H_met > 50./400. || sigEvent.H_dim12 < 2) &&
           //sigEvent.H_met < 150./400. &&
	   sigEvent.H_njets <= 1. &&
	   sigEvent.H_btag1 < 0.5 &&
	   sigEvent.H_maxmbtag < 0.5 &&
	   sigEvent.H_nmuons == 0 &&
	   sigEvent.H_mindr < 1.5 &&
	   1 == 1
	  ){
	  double myVar = TMath::Min(sigEvent.H_mindr/10.,0.999);
	  //myVar = sigEvent.H_dim01;
	  //myVar = TMath::Min(sigEvent.H_maxiso/10.,0.999);
  	  //myVar = TMath::Max(sigEvent.H_btag1*1.0,0.001);
	  myVar = sigEvent.H_met;
	  for(int n1=0; n1<nBin; n1++){
	    if(myVar > 1.0*n1/nBin){
	      S0[n1] = S0[n1] + sigEvent.H_weight*lumi;
	    }
	    if(myVar < 1.0*n1/nBin){
	      S1[n1] = S1[n1] + sigEvent.H_weight*lumi;
	    }
	  }
          hDSigOpt->Fill(myVar,sigEvent.H_weight*lumi);
	}
      } // HW->3l selection
      if(channel >= 100 && channel <= 800){ // H->WW selection
        int charge = (int)(sigEvent.H_q1 + sigEvent.H_q2);
        if(
           sigEvent.H_dim01   > 12/500. &&
           sigEvent.H_dim01   < 55/500. &&
           sigEvent.H_nlep == 2. &&
           charge == 0 &&
           sigEvent.H_njets == 0. &&
           sigEvent.H_nmuons == 0. &&
           sigEvent.H_ntracks <= 4. &&
           sigEvent.H_pt1 > 30./400. &&
           sigEvent.H_pt2 > 25./400. &&
           //(sigEvent.H_pt2 > 20./400. || sigEvent.H_ltype == 10.) &&
            sigEvent.H_met > 30./400. &&
           (sigEvent.H_met > 52./400. || sigEvent.H_ltype == 12. || sigEvent.H_ltype == 13.) &&
           //sigEvent.H_deltaphilmet > 60./180. &&
           (fabs(sigEvent.H_dim01*500.-91.1876) > 15. || sigEvent.H_ltype == 12. || sigEvent.H_ltype == 13.) && 
           sigEvent.H_maxd0sig < 0.025 && 
	   sigEvent.H_delphil < 57./180. &&
	   //!(sigEvent.H_ltype == 12. || sigEvent.H_ltype == 13.) &&
	   //(sigEvent.H_ltype == 10. || sigEvent.H_ltype == 12.) &&
           1 == 1
          ){
          double myVar = sigEvent.H_met;
          //myVar = sigEvent.H_pt1;
          //myVar = sigEvent.H_pt2;
          //myVar = sigEvent.H_delphil;
          //myVar = sigEvent.H_mthiggs;
          myVar = sigEvent.H_dim01;
          //myVar = TMath::Min(sigEvent.H_maxd0sig*1.0,0.0999)/0.1;
          //myVar = (sigEvent.H_ntracks+0.1)/20.;
          for(int n1=0; n1<nBin; n1++){
            if(myVar > 1.0*n1/nBin){
              S0[n1] = S0[n1] + 2.*sigEvent.H_weight*lumi;
            }
            if(myVar < 1.0*n1/nBin){
              S1[n1] = S1[n1] + 2.*sigEvent.H_weight*lumi;
            }
          }
          hDSigOpt->Fill(myVar,sigEvent.H_weight*lumi);
        }
      } // H->WW selection
    } // Loop over signal
  }

  if(channel == 1000 || channel == 1100 ||
     (channel >= 100 && channel <= 800)){
    for(int n1=0; n1<nBin; n1++){
      double signif0 = 0.0;
      if(S0[n1] != 0 && B0[n1] != 0) signif0 = scpFast(S0[n1],B0[n1],0.35*B0[n1],0.0);
      double signif1 = 0.0;
      if(S1[n1] != 0 && B1[n1] != 0) signif1 = scpFast(S1[n1],B1[n1],0.35*B1[n1],0.0);
      hDSignif[0]->SetBinContent(n1, signif0);
      hDSignif[1]->SetBinContent(n1, signif1);
    }
    printf("S: %f, B: %f, Signif: %f, S/B: %f\n",S0[0],B0[0],
    	   scpFast(S1[0],B1[0],0.35*B1[0],0.0),S0[0]/B0[0]);
  }
  else {
    for(int n1=0; n1<nBin; n1++){
      double signif0 = 0.0;
      if(S0[n1] != 0 || B0[n1] != 0) signif0 = S0[n1]/sqrt(S0[n1]+B0[n1]+0.09*B0[n1]*B0[n1]);
      double signif1 = 0.0;
      if(S1[n1] != 0 || B1[n1] != 0) signif1 = S1[n1]/sqrt(S1[n1]+B1[n1]+0.09*B1[n1]*B1[n1]);
      hDSignif[0]->SetBinContent(n1, signif0);
      hDSignif[1]->SetBinContent(n1, signif1);
    }
    printf("S: %f, B: %f, Signif: %f, S/B: %f\n",S0[0],B0[0],
    	   S0[0]/sqrt(S0[0]+B0[0]+0.09*B0[0]*B0[0]),S0[0]/B0[0]);
  }
  
  for(int i=0; i<45; i++){
    if(bgdDecay[i] > 0) printf("bdg(%2d) = %f +/- %f\n",i,bgdDecay[i],sqrt(bgdDecay[i]*weiDecay[i]));
  }
  //return;
  hDSignif[0]->SetDirectory(0);
  hDSignif[1]->SetDirectory(0);
  hDSigOpt->SetDirectory(0);
  hDBckOpt->SetDirectory(0);
  TCanvas* c1 = new TCanvas("c1","c1",100,100,700,800);
  c1->Divide(2,2);
  c1->cd(1);
  hDSignif[0]->Draw();
  c1->cd(2);
  hDSignif[1]->Draw();
  c1->cd(3);
  hDSigOpt->Draw();
  c1->cd(4);
  hDBckOpt->Draw();
  char output[200];
  sprintf(output,"histo_tmva_%d.root",mHiggs);
  TFile* outFilePlotsNote = new TFile(output,"recreate");
  outFilePlotsNote->cd();
    hDSigOpt->Write();
    hDBckOpt->Write();
  outFilePlotsNote->Close();

  bool doNicePlot = false;
  if(doNicePlot == true){
    TCanvas* newc = new TCanvas("newc","newc",100,100,700,800);
    newc->Divide(1);
    newc->cd(1);
    hDSigOpt->SetFillColor(3);
    hDSigOpt->SetFillStyle(1001);
    hDSigOpt->SetLineStyle(0);
    hDSigOpt->SetLineWidth(0);

    hDBckOpt->SetTitle("");
    hDBckOpt->SetFillColor(4);
    hDBckOpt->SetFillStyle(1001);
    hDBckOpt->SetLineStyle(0);
    hDBckOpt->SetLineWidth(0);
    hDBckOpt->GetXaxis()->SetTitle(xTitle);
    hDBckOpt->GetXaxis()->SetTitleOffset(1.0);
    hDBckOpt->GetYaxis()->SetTitle(yTitle);
    hDBckOpt->GetYaxis()->SetTitleOffset(1.3);

    hDBckOpt->Add(hDSigOpt);

    hDBckOpt->Draw("hist");
    hDSigOpt->Draw("same,hist");

    TLegend* leg = new TLegend(0.5,0.7,0.8,0.8);    	        				
    leg->SetFillColor(10);			    	        				
    leg->SetTextSize(0.035);			    	        			 
    leg->AddEntry(hDBckOpt,"WZ #rightarrow 3l + bkg.","F");	    	        	       
    leg->AddEntry(hDSigOpt,"WZ #rightarrow 3l","F");	    	   
    leg->Draw("same");
    TPaveText* labelcms  = new TPaveText(0.7,0.95,0.97,0.97,"NDCBR");
    labelcms->SetTextAlign(12);
    labelcms->SetTextSize(0.05);
    labelcms->SetFillColor(0);
    labelcms->AddText("L = 200 pb^{-1}");
    labelcms->SetBorderSize(0);
    labelcms->Draw("same");
  }

  return;

}
//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);
  TTree* t = (TTree*)inf->Get(tname);
  assert(t);

  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}

double scpFast(double sig, double bkg, double sigma_b, double delta_b)
{
double fac2  = sqrt(bkg+delta_b)/sqrt(bkg+sigma_b*sigma_b+delta_b);
double Sc12_sys = 2*(sqrt(sig+bkg)-sqrt(bkg+delta_b))*fac2;
return Sc12_sys;
}
