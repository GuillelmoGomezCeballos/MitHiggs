//
//root -l -b -q HWWDistributions.C+\(\"/home/ceballos/HiggsMVA/NeuralNetworkMaker-21x/ntuples/inputNtuple-data-standard-histo_H160_WW2l_all.root\",\"/home/ceballos/HiggsMVA/NeuralNetworkMaker-21x/ntuples/inputNtuple-data-standard-histo_HBCK2.root\",10\) |& tee LogCutbasedDefaultHWW160.txt |& tee LogCutbasedDefaultHWW160.txt

//root -l -b -q HWWDistributions.C+\(\"/home/ceballos/HiggsMVA/NeuralNetworkMaker-21x/ntuples/inputNtuple-data-standard-histo_H160_WW2l_all.root\",\"/home/ceballos/HiggsMVA/NeuralNetworkMaker-21x/ntuples/inputNtuple-data-standard-histo_HBCK2.root\",10,\"HWW160\",\"/home/sixie/Tools/root/CutBasedOptimization/HWWCutOptimizationResults.root\"\) |& tee LogCutbasedOptimizedHWW160.txt

//
//root -l -b -q HWWDistributions.C+\(\"/home/ceballos/HiggsMVA/NeuralNetworkMaker-21x/ntuples/inputNtuple-data-standard-histo_H130_WW2l_all.root\",\"/home/ceballos/HiggsMVA/NeuralNetworkMaker-21x/ntuples/inputNtuple-data-standard-histo_HBCK2.root\",10,\"HWW130\",\"/home/sixie/Tools/root/CutBasedOptimization/HWWCutOptimizationResults.root\"\) |& tee LogCutbasedOptimizedHWW130.txt

//root -l -b -q HWWDistributions.C+\(\"/home/ceballos/HiggsMVA/NeuralNetworkMaker-21x/ntuples/inputNtuple-data-standard-histo_H150_WW2l_all.root\",\"/home/ceballos/HiggsMVA/NeuralNetworkMaker-21x/ntuples/inputNtuple-data-standard-histo_HBCK2.root\",10,\"HWW150\",\"/home/sixie/Tools/root/CutBasedOptimization/HWWCutOptimizationResults.root\"\) |& tee LogCutbasedOptimizedHWW150.txt

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <utility>
#include "TFile.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "MitNtupleEvent.C"
#include <vector>
//#include <pair>

#define NVARIABLES 10
#define NPROCESSES 32
#define NBINS 25
//*************************************************************************************************
//Define Global Variables for the macro
//*************************************************************************************************
string finalstatestring[4] = {"_mumu","_ee", "_emu", "_mue" };
int fColors[NPROCESSES] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           7, 7, 7, 7, 6, 6, 8, 8, 8, 8, 5, 8, 8, 8, 2, 2, 2, 3, 9, 
                           9, 4, 11};

int fProcessCode[NPROCESSES] = { 120, 130, 150, 160, 170, 200, 500, 1130, 1140, 1150, 0, 1, 2, 
                                 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 25, 26,
                                 27, 28, 29 };
string fProcessName[NPROCESSES] = { "H120", "H130", "H150","H160", "H170", "H200","H500", "WH130", 
                                    "WH140", "WH150", "Wmunu", "Wenu","Wtaunu", "WMadgraph", 
                                    "ttbar","ttbarMadgraph", "Zmumu", "Zee","Ztautau","ZMadgraph", 
                                    "VQQmadgraph", "Zmumu_M20","Zee_M20", "Ztautau_M20", 
                                    "SingleTop_s" , "SingleTop_t","SingleTop_Wt", "WW_2l", "WZ_3l", 
                                    "WZ_incl", "ZZ_2l2n", "ZZ_4l" };

string fVariableNames[NVARIABLES] = { "ptmax", "ptmin", "met", "dileptonmass", "dileptonmassZVeto", 
                                      "deltaphiLeptons", "metdeltaphilmet", "ncentraljets", 
                                      "ndirtymuons", "nextratracks" };

//*******************************************************************************************************************

//CUTS

//*******************************************************************************************************************
// Old Cuts
// bool fCutAbove[NVARIABLES] = { true, false, true, true, false,  true, false,  true,  false,  false};
// bool fCutBelow[NVARIABLES] = { true,  true, true, true,  true, false, false, false, false, false};
// bool fixCutValueAcrossFinalstates[2][NVARIABLES] = {{false,false,false,false,false,false,false,true,true,true},
//                                                     {false,false,false,false,true,false,false,false,false,false}};
// Float_t fCutAboveInitialValue[3][NVARIABLES] = { {55.0, 5000.0, 200.0, 50.0, 5000.0, 45.0, 5000.0, 0.0, 0.0, 4.0}, 
//                                                  {49.0, 5000.0, 200.0, 40.0, 5000.0, 45.0, 5000.0, 0.0, 0.0, 4.0}, 
//                                                  {50.0, 5000.0, 105.0, 45.0, 5000.0, 70.0, 5000.0, 0.0, 0.0, 4.0},
//                                                  {50.0, 5000.0, 105.0, 45.0, 5000.0, 70.0, 5000.0, 0.0, 0.0, 4.0}};
// Float_t fCutBelowInitialValue[3][NVARIABLES] = { {35.0,   25.0,  50.0, 12.0,   15.0,  0.0,   0.0,  0.0, 0.0, 0.0}, 
//                                                  {25.0,   25.0,  51.0, 12.0,   15.0,  0.0,   0.0,  0.0, 0.0, 0.0}, 
//                                                  {25.0,   25.0,  45.0, 12.0,   15.0,  0.0,   0.0,  0.0, 0.0, 0.0},
//                                                  {25.0,   25.0,  45.0, 12.0,   15.0,  0.0,   0.0,  0.0, 0.0, 0.0}};



// Float_t fVariableRangeMin[NVARIABLES] = { 0.0,     0.0,   0.0,   0.0,  0.0,   0.0,   0.0, 0.0, 0.0,  0.0};
// Float_t fVariableRangeMax[NVARIABLES] = { 400.0, 400.0, 400.0, 200.0, 200.0, 180.0, 400.0, 5.0, 5.0, 15.0};
// bool fPlotThisVariable[NVARIABLES] = { true, true, true, true, true, true, true, true, true, true };

//*******************************************************************************************************************
//Optimization Result (Only mumu is optimized for now)
// bool fCutAbove[NVARIABLES] = { true, false, true, true, false,  true, false,  true,  true,  true};
// bool fCutBelow[NVARIABLES] = { true,  true, true, true,  true, false, false, false, false, false};
// bool fixCutValueAcrossFinalstates[2][NVARIABLES] = {{false,false,false,false,false,false,false,true,true,true},
//                                                     {false,false,false,false,true,false,false,false,false,false}};
// Float_t fCutAboveInitialValue[3][NVARIABLES] = { {400.0, 400.0, 400.0, 50.0, 5000.0,180.0, 400.0, 0.0, 0.0,15.0}, 
//                                                  {49.0, 5000.0, 200.0, 40.0, 5000.0, 45.0, 400.0, 0.0, 0.0,15.0}, 
//                                                  {50.0, 5000.0, 105.0, 45.0, 5000.0, 70.0, 400.0, 0.0, 0.0,15.0},
//                                                  {50.0, 5000.0, 105.0, 45.0, 5000.0, 70.0, 400.0, 0.0, 0.0,15.0}};
// Float_t fCutBelowInitialValue[3][NVARIABLES] = { {0.0,   30.0,   32.0, 12.0,    0.0,  0.0,   0.0, 0.0, 0.0, 0.0}, 
//                                                  {25.0,   25.0,  51.0, 12.0,    0.0,  0.0,   0.0, 0.0, 0.0, 0.0}, 
//                                                  {25.0,   25.0,  45.0, 12.0,    0.0,  0.0,   0.0, 0.0, 0.0, 0.0},
//                                                  {25.0,   25.0,  45.0, 12.0,    0.0,  0.0,   0.0, 0.0, 0.0, 0.0}};


//*******************************************************************************************************************


// //Cuts from Guillelmo which matches results of CIEMAT
bool fCutAbove[NVARIABLES] = { true, true, true, true, false, true, false, true, true, false};
bool fCutBelow[NVARIABLES] = { true, true, true, true, true, true, false, false, false, false};
bool fixCutValueAcrossFinalstates[2][NVARIABLES] = {{false,false,false,false,false,false,false,true,true,true},
                                                    {false,false,false,false,true,false,false,false,false,false}};
Float_t fCutAboveInitialValue[4][NVARIABLES] = { {50.0, 47.0, 400.0, 50.0, 5000.0, 57.0, 5000.0,  0.0, 0.0, 4.0}, 
                                                 {50.0, 30.0, 400.0, 50.0, 5000.0, 57.0, 5000.0,  0.0, 0.0, 4.0}, 
                                                 {70.0, 80.0, 400.0, 50.0, 5000.0, 57.0, 5000.0,  0.0, 0.0, 4.0},
                                                 {70.0, 80.0, 400.0, 50.0, 5000.0, 57.0, 5000.0,  0.0, 0.0, 4.0}};
Float_t fCutBelowInitialValue[4][NVARIABLES] = { {28.0, 25.0,  48.0, 12.0,   15.0,  0.0,    0.0,  0.0, 0.0, 0.0}, 
                                                 {28.0, 25.0,  48.0, 12.0,   15.0,  0.0,    0.0,  0.0, 0.0, 0.0}, 
                                                 {36.0, 25.0,  55.0, 12.0,   15.0,  0.0,    0.0,  0.0, 0.0, 0.0},
                                                 {36.0, 25.0,  55.0, 12.0,   15.0,  0.0,    0.0,  0.0, 0.0, 0.0}};

//*******************************************************************************************************************

Float_t fVariableRangeMin[NVARIABLES] = { 0.0,     0.0,   0.0,   0.0,  0.0,   0.0,   0.0, 0.0, 0.0,  0.0};
Float_t fVariableRangeMax[NVARIABLES] = { 400.0, 400.0, 400.0, 200.0, 200.0, 180.0, 400.0, 5.0, 5.0, 15.0};
bool fPlotThisVariable[NVARIABLES] = { true, true, true, true, true, true, true, true, true, true };



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

//*************************************************************************************************
//Check that given initial conditions satisfy requirements
//*************************************************************************************************
void CheckInitialCutValueConsistency() {

  //check that cutabove values are above cutbelow values
  for (int i=0;i<NVARIABLES;i++) {
    for (int f=0;f<=3;f++) {
      if (fCutAboveInitialValue[f][i] < fCutBelowInitialValue[f][i]) {
        cout << "cutAbove initial value if less than cutBelow initial value for Variable " << i
             << " in final state " << f << endl;
        assert (0==1);
      }        
    }
  }

  //check that values which should be fixed across final states are fixed
  for (int i=0;i<NVARIABLES;i++) {
    if (fixCutValueAcrossFinalstates[0][i] && 
        (fCutAboveInitialValue[0][i] != fCutAboveInitialValue[1][i] ||
         fCutAboveInitialValue[0][i] != fCutAboveInitialValue[2][i] || 
         fCutAboveInitialValue[0][i] != fCutAboveInitialValue[3][i] )
      ) {
      cout << "CutAbove initial values for Variable " << i << " are not fixed across final states." 
           << endl;
        assert (0==1);
    }
    if (fixCutValueAcrossFinalstates[1][i] && 
        (fCutBelowInitialValue[0][i] != fCutBelowInitialValue[1][i] ||
         fCutBelowInitialValue[0][i] != fCutBelowInitialValue[2][i] || 
         fCutBelowInitialValue[0][i] != fCutBelowInitialValue[3][i] )
      ) {
      cout << "CutBelow initial values for Variable " << i << " are not fixed across final states." 
           << endl;
        assert (0==1);
    }  
    
  }

}


//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}


//*************************************************************************************************
//Function to determine whether the event passes the defined cuts or not
//*************************************************************************************************
bool passCut(MitNtupleEvent *event , 
             Float_t fOptimizationVariables[NVARIABLES], Float_t fCutAboveValues[NVARIABLES],
             Float_t fCutBelowValues[NVARIABLES], int cutExceptionVariableIndex = -1) {

  bool pass = true;
  int cutindex = -2;

  //make preselection requirements
  int charge = (int)(event->H_q1 + event->H_q2);
  pass = pass && (event->H_nlep == 2);
  pass = pass && (charge == 0);
  pass = pass && (event->H_pt1 > 20.0/400.0);
  pass = pass && (event->H_pt2 > 10.0/400.0);

  if (pass == false) 
    cutindex = -1;

   for (int i=0;i<NVARIABLES;i++) {

     pass = pass && (fOptimizationVariables[i] <= fCutAboveValues[i] || !fCutAbove[i] 
                     || cutExceptionVariableIndex == i );
     pass = pass && (fOptimizationVariables[i] >= fCutBelowValues[i] || !fCutBelow[i] 
                     || cutExceptionVariableIndex == i );

     if (pass == false && cutindex < 0)
       cutindex = i;

     if (pass == false && cutExceptionVariableIndex == 6) {
       // cerr << i << " : " << fOptimizationVariables[i] << " " << fCutAboveValues[i] << " " 
       //             << fCutBelowValues[i] << " " << endl;
     }       
   }
   return pass;
}


//*************************************************************************************************
//Convert Variable Index to Variable Value by retrieving it from Guillelmo's ntuple
//My list of variables have to be defined by hand. This function serves as the definition.
//*************************************************************************************************
Double_t GetVariableValue (int variableIndex, MitNtupleEvent *event) {

  Double_t value = 0;
  if (variableIndex == 0) {
    value = event->H_pt1 * 400;
  } else if (variableIndex == 1) {
    value = event->H_pt2 * 400;
  } else if (variableIndex == 2) {
    value = event->H_met * 400;
  } else if (variableIndex == 3) {
    value = event->H_dim01*500;
  }  else if (variableIndex == 4) {
    value = fabs(event->H_dim01*500 - 91.1876);
  }  else if (variableIndex == 5) {
    value = event->H_delphil * 180;
  }  else if (variableIndex == 6) {
    value = event->H_metdeltaphilmet * 400 ;
  }  else if (variableIndex == 7) {
    value = event->H_njets;
  }  else if (variableIndex == 8) {
    value = event->H_nmuons;
  }  else if (variableIndex == 9) {
    value = event->H_ntracks;
  } else {
    cout << "variableIndex = " << variableIndex << " not recognized!" << endl;
    value = 0;
  }  
  return value;
}


//*************************************************************************************************
//Function to create a distribution histogram for the signal for the variable given
//*************************************************************************************************
TH1F* makeSignalDistributionHistogram ( int variableIndex, MitNtupleEvent* event , Int_t nentries, TLegend *legend, int OnlyThisFinalState = -1) {
  assert(variableIndex >= 0 && variableIndex < NVARIABLES);
  string label = "NSignalEventsPass";
  string histName = "h" + fVariableNames[variableIndex] + "_" + label;
  string axisLabel = ";" + fVariableNames[variableIndex] + ";Number of Events";
  TH1F *histogram = new TH1F(histName.c_str(),axisLabel.c_str(),NBINS, 
                             fVariableRangeMin[variableIndex], fVariableRangeMax[variableIndex]);
  histogram->SetDirectory(0);

  for (int n=0;n<nentries;n++) {
    event->GetEntry(n);
    float eventweight = event->H_weight;
    int finalstatetype = (int)event->H_ltype;

    //only look at events with the given finalstatetype
    if (OnlyThisFinalState >= 0 && OnlyThisFinalState != finalstatetype)
      continue;

    Float_t cutAboveValues[NVARIABLES];
    Float_t cutBelowValues[NVARIABLES];
    Float_t variableValues[NVARIABLES];   

    //only look at events with the given finalstatetype
    //convert finalstatetype into the array index
    int fs;
    if (finalstatetype == 10)
      fs = 0;
    else if (finalstatetype == 11)
      fs = 1;
    else if (finalstatetype == 12)
      fs = 2;
    else if (finalstatetype == 13)
      fs = 3;
    else {
      fs = -1;
      continue;
    }

    for (int j=0;j<NVARIABLES;j++) {
      cutAboveValues[j] = fCutAboveInitialValue[fs][j];
      cutBelowValues[j] = fCutBelowInitialValue[fs][j];
      variableValues[j] = GetVariableValue(j, event);      
    }

    //for N-1 distribution
    if (passCut(event,variableValues, cutAboveValues, cutBelowValues, variableIndex)) { 
      
    //for distribution after all cuts
    //if (passCut(event,variableValues, cutAboveValues, cutBelowValues)) {
      histogram->Fill(variableValues[variableIndex], eventweight);
    }  
  }
  legend->AddEntry(histogram, "signal", "LP"); 
  return histogram;
}


//*************************************************************************************************
//Function to create a distribution stacked histogram for all the bkg processes 
//for the variable given
//*************************************************************************************************
THStack* makeBkgDistributionHistogram ( int variableIndex, MitNtupleEvent* event , 
                                        Int_t nentries, TLegend *legend, 
                                        int OnlyThisFinalState = -1) {
  assert(variableIndex >= 0 && variableIndex < NVARIABLES);
  string histName = fVariableNames[variableIndex] + "_bkg";
  string axisLabel = ";" + fVariableNames[variableIndex] + ";Number of Events";
  THStack *stackHistogram = new THStack(histName.c_str(),axisLabel.c_str());
  //stackHistogram->SetDirectory(0);

  vector< pair<int, TH1F*> > bkgHistograms;

  for (int n=0;n<nentries;n++) { 
    event->GetEntry(n);
    float eventweight = event->H_weight;
    int finalstatetype = (int)event->H_ltype;
    int decay = (int)event->H_decay;

    //only look at events with the given finalstatetype
    if (OnlyThisFinalState >= 0 && OnlyThisFinalState != finalstatetype)
      continue;

    Float_t cutAboveValues[NVARIABLES];
    Float_t cutBelowValues[NVARIABLES];
    Float_t variableValues[NVARIABLES];   

    //only look at events with the given finalstatetype
    //convert finalstatetype into the array index
    int fs;
    if (finalstatetype == 10)
      fs = 0;
    else if (finalstatetype == 11)
      fs = 1;
    else if (finalstatetype == 12)
      fs = 2;
    else if (finalstatetype == 13)
      fs = 3;
    else {
      fs = -1;
      continue;
    }

    for (int j=0;j<NVARIABLES;j++) {
      cutAboveValues[j] = fCutAboveInitialValue[fs][j];
      cutBelowValues[j] = fCutBelowInitialValue[fs][j];
      variableValues[j] = GetVariableValue(j, event);      
    }
    //printf(" B %d\n", n);

    //For N-1 Distributions
    if (passCut(event,variableValues, cutAboveValues, cutBelowValues, variableIndex)) {

    //For distribution after all cuts
    //if (passCut(event,variableValues, cutAboveValues, cutBelowValues)) {

      //Search for the bkg process to see if it already exists
      bool foundProcess = false;
      int foundIndex = -1;      
      for (UInt_t p=0;p<bkgHistograms.size();p++) {
        //cout << "Check bkgHistograms " << p << " " << bkgHistograms[p].first << " " 
        //       <<  bkgHistograms[p].second << endl;
        if (decay == bkgHistograms[p].first) {
          foundProcess = true;
          foundIndex = p;
        }
      }

      if (!foundProcess) {
        //get name of the process
        int decayIndex = -1;

        for (int k=0;k<NPROCESSES;k++) {
          if (fProcessCode[k] == decay)
            decayIndex = k;
        }
        string temphistname = "";
        if (decayIndex >= 0) {
          temphistname = fProcessName[decayIndex];
        }

        string temphistAxisLabel = ";" + fVariableNames[variableIndex] + ";Number of Events";
        TH1F *temphist = new TH1F(temphistname.c_str(),axisLabel.c_str(),NBINS, 
                                  fVariableRangeMin[variableIndex], 
                                  fVariableRangeMax[variableIndex]);
        temphist->SetDirectory(0);

        //temphist->SetFillStyle(1001);
        temphist->SetFillColor(fColors[decayIndex]);
        temphist->SetLineWidth(1);

        //Fill histogram
        temphist->Fill( variableValues[variableIndex], eventweight);

        pair< int, TH1F* > tempPair(decay,temphist);
        bkgHistograms.push_back(tempPair);

        //Add legend entry
        legend->AddEntry(temphist, temphistname.c_str(), "f"); 

      } else {       
        //fill histogram
        bkgHistograms[foundIndex].second->Fill( variableValues[variableIndex], eventweight);
      }
    }
  }

  //add to the stack
  for (UInt_t p=0;p<bkgHistograms.size();p++) {
    stackHistogram->Add(bkgHistograms[p].second);
  }

  return stackHistogram;
}

//*************************************************************************************************
//Draws the signal and background histograms together and makes gif file
//*************************************************************************************************
void DrawSigBkgHistogram(TH1F* sig, THStack* bkg, TLegend *legend , string histname) {
  string filename = histname + ".gif";
  TCanvas *cv = new TCanvas(histname.c_str(), histname.c_str(), 0,0,800,600);

  sig->SetLineColor(1);
  sig->SetLineWidth(2);

  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);

  //set best Y-range
  double tempMaxY1 = sig->GetMaximum();
  double tempMaxY2 = bkg->GetMaximum();
  double MaxY = max(tempMaxY1, tempMaxY2);
  sig->SetMinimum(0.0);
  sig->SetMaximum(1.1*MaxY);
  bkg->SetMinimum(0.0);
  bkg->SetMaximum(1.1*MaxY);
  bkg->Draw();
  sig->Draw("same");
  legend->Draw();
  cv->SaveAs(filename.c_str());
  delete cv;
}

//*************************************************************************************************
//Make Final Plots
//*************************************************************************************************
void MakeStackedDistributionPlots( MitNtupleEvent* signalEvent , MitNtupleEvent* bkgEvent , 
                                   Int_t NSignalEntries, Int_t NBkgEntries, int finalstate) {
  string finalstatelabel = "";
  if (finalstate == -1) 
    finalstatelabel = "_all";
  else if (finalstate == 10)
    finalstatelabel = "_mumu";
  else if (finalstate == 11)
    finalstatelabel = "_ee";
  else if (finalstate == 12) 
    finalstatelabel = "_emu";
  else if (finalstate == 13) 
    finalstatelabel = "_mue";
  else {
    cout << "finalstate = " << finalstate << " is invalid." << endl;
    assert (0==1);
  }
    

  for (int v = 0; v < NVARIABLES;v++)
  {
    //make distribution histograms after final cuts    
    TLegend *legend = new TLegend(0.5,0.65,0.75,0.85);

    TH1F* signalDistribution = makeSignalDistributionHistogram(v, signalEvent, NSignalEntries, 
                                                               legend, finalstate);
    THStack* bkgDistribution = makeBkgDistributionHistogram(v, bkgEvent, NBkgEntries, 
                                                            legend, finalstate);
    DrawSigBkgHistogram(signalDistribution, bkgDistribution, legend, "Distribution_" + 
                        fVariableNames[v] + finalstatelabel);    
  }  
}

//*************************************************************************************************
//Print Final Results
//*************************************************************************************************
void PrintFinalResults( MitNtupleEvent* signalEvent , MitNtupleEvent* bkgEvent , 
                        Int_t NSignalEntries, Int_t NBkgEntries, int OnlyThisFinalState) {
  string finalstatelabel = "";
  if (OnlyThisFinalState == 10)
    finalstatelabel = "_mumu";
  else if (OnlyThisFinalState == 11)
    finalstatelabel = "_ee";
  else if (OnlyThisFinalState == 12)
    finalstatelabel = "_emu";
  else if (OnlyThisFinalState == 13)
    finalstatelabel = "_mue";


  //add up signal
  Double_t TotalSignalCount = 0.0;
  for (int n=0;n<NSignalEntries;n++) { 

    signalEvent->GetEntry(n);
    float eventweight = signalEvent->H_weight;
    int finalstatetype = (int)signalEvent->H_ltype;
    if (!(finalstatetype >= 10 && finalstatetype <= 12)) {
      //cout << "finalstatetype = " << finalstatetype << " unrecognized" << endl;
      //assert(!(finalstatetype >= 10 && finalstatetype <= 12));
    }
    //only look at events with the given finalstatetype
    if (OnlyThisFinalState >= 0 && OnlyThisFinalState != finalstatetype)
      continue;

    Float_t cutAboveValues[NVARIABLES];
    Float_t cutBelowValues[NVARIABLES];
    Float_t variableValues[NVARIABLES];   

    //only look at events with the given finalstatetype
    //convert finalstatetype into the array index
    int fs;
    if (finalstatetype == 10)
      fs = 0;
    else if (finalstatetype == 11)
      fs = 1;
    else if (finalstatetype == 12)
      fs = 2;
    else if (finalstatetype == 13)
      fs = 3;
    else {
      fs = -1;
      continue;
    }

    for (int j=0;j<NVARIABLES;j++) {
      cutAboveValues[j] = fCutAboveInitialValue[fs][j];
      cutBelowValues[j] = fCutBelowInitialValue[fs][j];
      variableValues[j] = GetVariableValue(j, signalEvent);      
    }
    if (passCut(signalEvent,variableValues, cutAboveValues, cutBelowValues)) {
      TotalSignalCount += eventweight;
    }
  }

  //add up bkg
  Double_t TotalBkg = 0.0;
  vector < pair<int , double> > TotalBkgCount;
  
  for (int n=0;n<NBkgEntries;n++) { 

    bkgEvent->GetEntry(n);
    float eventweight = bkgEvent->H_weight;
    int finalstatetype = (int)bkgEvent->H_ltype;
    int decay = (int)bkgEvent->H_decay;

    //only look at events with the given finalstatetype
    if (OnlyThisFinalState >= 0 && OnlyThisFinalState != finalstatetype)
      continue;

    Float_t cutAboveValues[NVARIABLES];
    Float_t cutBelowValues[NVARIABLES];
    Float_t variableValues[NVARIABLES];   

    //only look at events with the given finalstatetype
    //convert finalstatetype into the array index
    int fs;
    if (finalstatetype == 10)
      fs = 0;
    else if (finalstatetype == 11)
      fs = 1;
    else if (finalstatetype == 12)
      fs = 2;
    else if (finalstatetype == 13)
      fs = 3;
    else {
      fs = -1;
      continue;
    }

    for (int j=0;j<NVARIABLES;j++) {
      cutAboveValues[j] = fCutAboveInitialValue[fs][j];
      cutBelowValues[j] = fCutBelowInitialValue[fs][j];
      variableValues[j] = GetVariableValue(j, bkgEvent);      
    }
    if (passCut(bkgEvent,variableValues, cutAboveValues, cutBelowValues)) {
      TotalBkg += eventweight;
      
      //Search for the bkg process to see if it already exists
      bool foundProcess = false;
      int foundIndex = -1;      
      for (UInt_t p=0;p<TotalBkgCount.size();p++) {
        //cout << "Check TotalBkgCount " << p << " " << TotalBkgCount[p].first << " " 
        //       <<  TotalBkgCount[p].second << endl;
        if (decay == TotalBkgCount[p].first) {
          foundProcess = true;
          foundIndex = p;
        }
      }

      if (!foundProcess) {
        pair< int, double > newBkg(decay,eventweight);
        TotalBkgCount.push_back(newBkg);
      } else {
        TotalBkgCount[foundIndex].second = TotalBkgCount[foundIndex].second + eventweight;
      }
    }
  }

    cout << "For Final State : " << finalstatelabel << endl;
    cout << "Total Signal : " << TotalSignalCount << endl;
    cout << "Total Bkg : " << TotalBkg << endl;
    for (UInt_t p=0;p<TotalBkgCount.size();p++) {      
        //find the process
        int decayIndex = -1;
        for (int k=0;k<NPROCESSES;k++) {
          if (fProcessCode[k] == TotalBkgCount[p].first)
            decayIndex = k;
        }

        if (decayIndex >= 0) {
          cout << "Bkg " << fProcessName[decayIndex] << " " << TotalBkgCount[p].second << endl;
        }       
      }


}

void LoadCutValues(const string cutValuesFile , const string SignalName ) {

  
  TFile *file = new TFile(cutValuesFile.c_str(), "READ");
  file->cd();  


  //loop over final states
  for (int fs = 0; fs <= 3; fs++) {
 
    //loop over 
    string histname = ("OptimalCutAboveValues_"+SignalName+finalstatestring[fs]);
    TH1F *cutAboveValues = (TH1F*)file->FindObjectAny(histname.c_str());
    //get cut values for all variables
    if (cutAboveValues) {
      for (int v=0;v<NVARIABLES;v++) {
        fCutAboveInitialValue[fs][v] = cutAboveValues->GetBinContent(v+1);
      }
    } else {
      cerr << "Warning : cutAboveValues histogram " << histname << " was not found in file " << cutValuesFile << endl;
    }

    histname = ("OptimalCutBelowValues_"+SignalName+finalstatestring[fs]);
    TH1F *cutBelowValues = (TH1F*)file->FindObjectAny(histname.c_str());
    //get cut values for all variables
    if (cutBelowValues) {
      for (int v=0;v<NVARIABLES;v++) {
        fCutBelowInitialValue[fs][v] = cutBelowValues->GetBinContent(v+1);
      }
    } else {
      cerr << "Warning : cutAboveValues histogram " << histname << " was not found in file " << cutValuesFile << endl;
    }
  }




}


//*************************************************************************************************
//Main part of the macro
//*************************************************************************************************
void HWWDistributions(const string signalInputFilename, const string bkgInputFilename,                       
                      int OnlyThisFinalState = -1,
                      const string SignalName = "", const string cutValuesFile = "") {
  TTree* signalTree = getTreeFromFile(signalInputFilename.c_str(),"all");
  TTree* bkgTree = getTreeFromFile(bkgInputFilename.c_str(),"all");
  assert(signalTree);
  assert(bkgTree);

  CheckInitialCutValueConsistency();
  
  //Load cut values from a root file;
  if (cutValuesFile != "" && SignalName != "")
    LoadCutValues(cutValuesFile, SignalName);

  //Summarize Cuts
  cout << "*********************************************************************\n";
  cout << "**  Load Cut Values from file " << cutValuesFile << endl;
  cout << "*********************************************************************\n";
  cout << "Final Cuts" << endl;
  cout << "*************************************************************************" << endl;

  for (int f=0; f<=3;f++) {  
    cout << "Final State : " << finalstatestring[f] << endl;
    for (int v = 0; v < NVARIABLES;v++) {
      cout << fVariableNames[v] << " : ";
      if (fCutBelow[v])
        cout << fCutBelowInitialValue[f][v];
      else 
        cout << "-Infinity ";
      cout << " -> ";
      if (fCutAbove[v])
        cout << fCutAboveInitialValue[f][v];      
      else 
        cout << " Infinity ";
      cout << endl;
    }
    cout << "----------------------------------------------------------------------" << endl;
  }

  MitNtupleEvent signalEvent(signalTree);
  MitNtupleEvent bkgEvent(bkgTree);

  cout << "*************************************************************************" << endl;
  cout << "Make Distribution Plots" << endl;
  cout << "*************************************************************************" << endl;
  if (OnlyThisFinalState == -1) {
    for (int f=10; f<=13;f++) {
      MakeStackedDistributionPlots(&signalEvent , &bkgEvent , signalTree->GetEntries(), 
                                   bkgTree->GetEntries(), f);
    }
  } else {
    MakeStackedDistributionPlots(&signalEvent , &bkgEvent , signalTree->GetEntries(), 
                                 bkgTree->GetEntries(), OnlyThisFinalState);
  }
  

  cout << "*************************************************************************" << endl;
  cout << "Print Final Results" << endl;
  cout << "*************************************************************************" << endl;
  if (OnlyThisFinalState == -1) {
    for (int f=10; f<=13;f++) {
      PrintFinalResults(&signalEvent , &bkgEvent , signalTree->GetEntries(), 
                        bkgTree->GetEntries(), f);
    }
  } else {
    PrintFinalResults(&signalEvent , &bkgEvent , signalTree->GetEntries(), 
                      bkgTree->GetEntries(), OnlyThisFinalState);
  }
  
}

