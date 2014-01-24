#!/usr/bin/env python
import sys, getopt, os, fnmatch, commands
import string

###declare 
setupcmsswcommand = 'source ~/EVAL_SH65 5_3_14;'

def usage():
    print "possible options are: --help, --InputPath=<myInputPath>, --OutputPath=<myOutputPath>, --FilenameHeader=<myFileHeader>, --DatasetListFile=<datasetListFile>"

def filesExist( path, filenameExpression ):
    exists = False
    for fileName in os.listdir ( path ):
        if fnmatch.fnmatch ( fileName, filenameExpression ):
            exists = True

    return exists

########################################################################################
##Merge Filesets
########################################################################################
def NormalizeNtuples( datasetNameList, skimNameList, inputPath, outputPath, filenameHeader):
    n = 0
    for dataset in datasetNameList:
        print "Perform normalization for dataset " + dataset

        inputFilename  = inputPath  + filenameHeader + '_' + dataset + '_' + skimNameList[n] + '.root'
        outputFilename = outputPath + filenameHeader + '_' + dataset + '_' + skimNameList[n] + '_normalized' + '.root'        
        command = setupcmsswcommand + 'root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/NormalizeHwwNtuple.C+\(\\\"' + inputFilename + '\\\",\\\"' + dataset + '\\\",\\\"' + outputFilename + '\\\",0\\)'
        if (os.path.isfile(inputFilename)):
            print "command : " + command
            if (os.path.exists(outputFilename)):
                print "Normalized output file already exists. It will be deleted. : " + outputFilename
                os.system('rm ' + outputFilename)            
            os.system(command)
            print ''
        else:
            print "Warning : Input file " + inputFilename + " does not exist. "

        inputFilename  = inputPath  + filenameHeader + '_' + dataset + '_' + skimNameList[n] + '.root'
        outputFilename = outputPath + filenameHeader + '_' + dataset + '_fake_' + skimNameList[n] + '_normalized' + '.root'        
        command = setupcmsswcommand + 'root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/NormalizeHwwNtuple.C+\(\\\"' + inputFilename + '\\\",\\\"' + dataset + '\\\",\\\"' + outputFilename + '\\\",1\\)'
        if (os.path.isfile(inputFilename)):
            print "command : " + command
            if (os.path.exists(outputFilename)):
                print "Normalized output file already exists. It will be deleted. : " + outputFilename
                os.system('rm ' + outputFilename)            
            os.system(command)
            print ''
        else:
            print "Warning : Input file " + inputFilename + " does not exist. "

        inputFilename  = inputPath  + filenameHeader + '_' + dataset + '_' + skimNameList[n] + '.root'
        outputFilename = outputPath + filenameHeader + '_' + dataset + '_2fake_' + skimNameList[n] + '_normalized' + '.root'        
        command = setupcmsswcommand + 'root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/NormalizeHwwNtuple.C+\(\\\"' + inputFilename + '\\\",\\\"' + dataset + '\\\",\\\"' + outputFilename + '\\\",2\\)'
        if (os.path.isfile(inputFilename)):
            print "command : " + command
            if (os.path.exists(outputFilename)):
                print "Normalized output file already exists. It will be deleted. : " + outputFilename
                os.system('rm ' + outputFilename)            
            os.system(command)
            print ''
        else:
            print "Warning : Input file " + inputFilename + " does not exist. "

        inputFilename  = inputPath  + filenameHeader + '_' + dataset + '_' + skimNameList[n] + '.root'
        outputFilename = outputPath + filenameHeader + '_' + dataset + '_lgamma_' + skimNameList[n] + '_normalized' + '.root'        
        command = setupcmsswcommand + 'root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/NormalizeHwwNtuple.C+\(\\\"' + inputFilename + '\\\",\\\"' + dataset + '\\\",\\\"' + outputFilename + '\\\",3\\)'
        if (os.path.isfile(inputFilename)):
            print "command : " + command
            if (os.path.exists(outputFilename)):
                print "Normalized output file already exists. It will be deleted. : " + outputFilename
                os.system('rm ' + outputFilename)            
            os.system(command)
            print ''
        else:
            print "Warning : Input file " + inputFilename + " does not exist. "
        print ''

        inputFilename  = inputPath  + filenameHeader + '_' + dataset + '_' + skimNameList[n] + '.root'
        outputFilename = outputPath + filenameHeader + '_' + dataset + '_gamma_' + skimNameList[n] + '_normalized' + '.root'        
        command = setupcmsswcommand + 'root -l -b -q $CMSSW_BASE/src/MitHiggs/macros/ntuples/NormalizeHwwNtuple.C+\(\\\"' + inputFilename + '\\\",\\\"' + dataset + '\\\",\\\"' + outputFilename + '\\\",4\\)'
        if (os.path.isfile(inputFilename)):
            print "command : " + command
            if (os.path.exists(outputFilename)):
                print "Normalized output file already exists. It will be deleted. : " + outputFilename
                os.system('rm ' + outputFilename)            
            os.system(command)
            print ''
        else:
            print "Warning : Input file " + inputFilename + " does not exist. "
        print ''
        n += 1

########################################################################################
#Main Program
########################################################################################

datasetListFile = ''
inputPath = ''
outputPath = ''
filenameHeader = ''
datasetNameList = list()
skimNameList = list()


if len(sys.argv[1:]) < 1:
    print "Error: not enough parameters specified"
    usage()
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:], "hi:o:f:d:", ["help", "InputPath=", "OutputPath=", "FilenameHeader=" , "DatasetListFile="])
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--InputPath"):
            inputPath = a + "/"
        elif o in ("-o", "--OutputPath"):
            outputPath = a + "/"
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

if (outputPath == ''):
    print "Error: No OutputPath specified."
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




NormalizeNtuples( datasetNameList, skimNameList, inputPath, outputPath, filenameHeader)
