#!/usr/bin/env python
import sys, getopt, os, fnmatch, commands
import string

###declare 
setupcmsswcommand = 'source ~/EVAL_SH65 5_2_3_patch3;'

def usage():
    print "possible options are: --help, --InputPath=<myInputPath>, --OutputPath=<myOutputPath>, --FilenameHeader=<myFileHeader>, --DatasetListFile=<datasetListFile>, --TestSampleFraction=<TestSampleFraction>"

def filesExist( path, filenameExpression ):
    exists = False
    for fileName in os.listdir ( path ):
        if fnmatch.fnmatch ( fileName, filenameExpression ):
            exists = True

    return exists

########################################################################################
##Merge all samples in list
########################################################################################
def MergeAllSamples( datasetNameList, skimNameList, inputPath, filenameHeader, outputLabel):
    n = 0

    outputFilename = inputPath + filenameHeader + '_' + outputLabel + '_unrandomized.root'    
    command = setupcmsswcommand + 'hadd -f ' + outputFilename + ' '
    
    for dataset in datasetNameList:        
        inputFilename = inputPath + filenameHeader + '_' + dataset + '_' + skimNameList[n] + '_normalized.root'        
        if (os.path.isfile(inputFilename)):
            command = command + inputFilename + ' '
        else:
            print "Warning : Input file " + inputFilename + " does not exist. "
        print ''
        n += 1

    print "merge command: " + command
    os.system(command)
    
########################################################################################
##Create Training And Test Samples
########################################################################################
def CreateTrainingAndTestSamples( inputFilename, testSampleFraction):

    command = setupcmsswcommand + 'root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/CreateTrainingAndTestSamples.C+\(\\\"' + inputFilename + '\\\",' + testSampleFraction + '\)'
    print "command : " + command
    os.system(command)
    print ''


########################################################################################
#Main Program
########################################################################################

datasetListFile = ''
inputPath = ''
outputLabel = ''
testSampleFraction = 0.0
filenameHeader = ''
datasetNameList = list()
skimNameList = list()


if len(sys.argv[1:]) < 1:
    print "Error: not enough parameters specified"
    usage()
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:t:f:d:", ["help", "InputPath=", "OutputLabel=", "TestSampleFraction=", "FilenameHeader=" , "DatasetListFile="])
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--InputPath"):
            inputPath = a + "/"
        elif o in ("-o", "--OutputLabel"):
            outputLabel = a 
        elif o in ("-t", "--TestSampleFraction"):
            testSampleFraction = a
        elif o in ("-f", "--FilenameHeader"):
            filenameHeader = a       
        elif o in ("-d", "--DatasetListFile"):
            datasetListFile = a
        else:
            usage()
            sys.exit()
except getopt.GetoptError:
    usage()
    sys.exit(2)

if (inputPath == ''):
    print "Error: No InputPath specified."
    sys.exit()

if (filenameHeader == ''):
    print "Error: No FilenameHeader specified."
    sys.exit()

if (datasetListFile == ''):
    print "Error: No dataset list file specified."
    sys.exit()

try:
    inputFile = open(datasetListFile,"r")
except IOError:
    print "Error: The specified dataset list file " + datasetListFile + " could not be opened."
    sys.exit()


########################################################################################
#Read in List of Datasets and skimnames
########################################################################################
lineNumber = 1
templine = inputFile.readline()
while len(templine.split()) > 0:
    if (templine[0] == '#'):
        templine = inputFile.readline()
        lineNumber += 1
        continue;
    if (len(templine.split()) != 2) :
        print "error: incorrect format for cross section file. Check line %s" % lineNumber
        sys.exit()
    else:
        tempInputList = templine.split()
        datasetNameList.append(tempInputList[0])
        skimNameList.append(tempInputList[1])
        lineNumber += 1
        #print numberOfVariables
    templine = inputFile.readline()
else:
    #print "done reading input"
    print ''
inputFile.close()

#Check the list of variables
count = 0
for l in datasetNameList:
    print l
    count += 1


#Do merging
MergeAllSamples( datasetNameList, skimNameList, inputPath, filenameHeader, outputLabel)

#randomize and create test and training samples
unrandomizedFilename = inputPath + filenameHeader + '_' + outputLabel + '_unrandomized.root'    
CreateTrainingAndTestSamples( unrandomizedFilename, testSampleFraction)
