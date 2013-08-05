#!/bin/sh

#do regular merging among different datasets
python $CMSSW_BASE/src/MitHiggs/scripts/MergeFilesets.py --InputPath=/home/sixie/ntuples/HwwNtuple/filler/011/ --OutputPath=/home/sixie/ntuples/HwwNtuple/filler/011/ --FilenameHeader=HwwNtuple --DatasetListFile=$CMSSW_BASE/src/MitHiggs/scripts/HiggsAnalysisFullSampleList.txt

#do normalization first
python $CMSSW_BASE/src/MitHiggs/scripts/NormalizeHiggsNtuple.py --InputPath=/home/sixie/ntuples/HwwNtuple/filler/011/ --OutputPath=/home/sixie/ntuples/HwwNtuple/filler/011/ --FilenameHeader=HwwNtuple --DatasetListFile=$CMSSW_BASE/src/MitHiggs/scripts/HiggsAnalysisFullSampleList.txt


#create training and test sample and adjust their normalizations accordingly
#we run this once for each of the different bkg/signal classes.

##H->WW Signal Sample
python $CMSSW_BASE/src/MitHiggs/scripts/CreateTrainingAndTestSamples.py --InputPath=/home/sixie/ntuples/HwwNtuple/filler/011/ --FilenameHeader=HwwNtuple --DatasetListFile=$CMSSW_BASE/src/MitHiggs/scripts/HiggsAnalysisSignalSampleList.txt --OutputLabel=HwwSignal --TestSampleFraction=0.5

##H->WW Pythia Bkg Samples
python $CMSSW_BASE/src/MitHiggs/scripts/CreateTrainingAndTestSamples.py --InputPath=/home/sixie/ntuples/HwwNtuple/filler/011/ --FilenameHeader=HwwNtuple --DatasetListFile=$CMSSW_BASE/src/MitHiggs/scripts/HiggsAnalysisPythiaBkgSampleList.txt --OutputLabel=HwwPythiaBkg --TestSampleFraction=0.5

##WW Signal Sample
python $CMSSW_BASE/src/MitHiggs/scripts/CreateTrainingAndTestSamples.py --InputPath=/home/sixie/ntuples/HwwNtuple/filler/011/ --FilenameHeader=HwwNtuple --DatasetListFile=$CMSSW_BASE/src/MitHiggs/scripts/WWAnalysisSignalSampleList.txt --OutputLabel=WWSignal --TestSampleFraction=0.5

##WW Pythia Bkg Samples
python $CMSSW_BASE/src/MitHiggs/scripts/CreateTrainingAndTestSamples.py --InputPath=/home/sixie/ntuples/HwwNtuple/filler/011/ --FilenameHeader=HwwNtuple --DatasetListFile=$CMSSW_BASE/src/MitHiggs/scripts/WWAnalysisPythiaBkgSampleList.txt --OutputLabel=WWPythiaBkg --TestSampleFraction=0.5

