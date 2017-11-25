#!/usr/local/bin/python

'''
						Analysis of coverage and insert size for a single genomic interval

BEGIN COPYRIGHT NOTICE

     insertSizeAnalysisSingle code -- (c) 2017 Dimitrios Kleftogiannis -- ICR -- www.icr.ac.uk

     Copyright 2017 Dimitrios Kleftogiannis Licensed under the
     Educational Community License, Version 2.0 (the "License"); you may
     not use this file except in compliance with the License. You may
     obtain a copy of the License at

     https://opensource.org/licenses/ECL-2.0

     Unless required by applicable law or agreed to in writing,
     software distributed under the License is distributed on an "AS IS"
     BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
     or implied. See the License for the specific language governing
     permissions and limitations under the License.

     Published reports of research using this code (or a modified version) should cite the 
     relevant article of this tool. 
     
     Comments and bug reports are welcome.
       
     Email to dimitrios.kleftogiannis@icr.ac.uk 
     I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
     You are free to modify, extend or distribute this code, as long 
     as this copyright notice is included whole and unchanged. 

END COPYRIGHT NOTICE
 
UTILITY
  This program takes as input a bam file and one interval in the format chrom:start-end and computes
  for this interval, the number of reads, the number of paired ones, the number of not paired and 
  the insert size distribution (TLEN field in the sam format). Please see the header of the program's output.


INPUT ARGUMENTS
    
    1. bam file              : a bam file from the tumour or plasma of interest
    
    2. interval              : in the designated format chrom:starPos-endPos
                              
    
DEPENDENCIES
    
    This is a Python program, thus you need the Python 2 compiler to be installed in your computer.
    
    The program has been developed and tested in a Mac OS computer with El Capitan version 10.11.5
    
    The Python compiler used for development is the Python Python 2.7.10 (default, Oct 23 2015, 19:19:21)
    
    The program works for Unix-like systems but has not been tested for Windows operating systems. 
    
    The program has not been tested in Cygwing-like systems running under Windows operating systems.   

    The program depends on pysam libraries downloaded from http://pysam.readthedocs.io/en/latest/index.html

    The program also depends on samtools, so please make sure that SAMtools is installed and configured properly in your system

    You might need to add samtools in your path so after you intall SAMtools you might need a command like: 

    PATH=$PATH:/your/path/to/Samtools


RUNNING
	
	An execution example is as follows:

    python insertSizeAnalysisSingle.py bamFile=0097_test.chr17.bam regions=chrom:100-10000

    To obtain the toy data we used for testing please contact Dimitrios

'''

#modules we need, I might have some extra that didnt use in the final version, but I forgot to remove.
#Remember that this program is under-developement so you may find block of codes used for testing.
import sys
import os
import pysam
import re
from collections import defaultdict
from itertools import groupby
import datetime
import time

#define minimum mapping quality; in a next implementation this can be adjusted by the user
quality='10'


#prints information about program's execution
def printUsage():
    print('To run this program please type the following:')
    print('\tpython insertSizeAnalysis.py bamFile=file.bam region=chr1:1500-10000\n')
    print('Where:\n') 
    print('\tfile.bam is a bam file. The bam file needs to be indexed with samtools.\n')
    print('\tregion gives the genomic coordinates of the interval of interest in a specific format.\n')
    print('Execution example:\n') 
    print('\tpython insertSizeAnalysis.py bamFile=0097_test.bam region=chr11:1500-10000 \n')
    print('Please give the arguments in the indicated order similar to the provided example!\n') 
    print('\t\tBy default the program applies MAPQ=10\n') 

#this function takes as input the bamFile, and the region of interest and produces a sam file with all reads spanning the positions
def generateSamFile(bamFile,samFileName,chrom,startPos,endPos,SAM_FILES):
#first check if the bam file exists
    if os.path.exists(bamFile):
        #the bam file exists
        
        #write the interval in a bed file
        tmpBedFile=SAM_FILES+'/'+chrom+'_'+startPos+'_tmp.bed'
        tmpOut=open(tmpBedFile,'w')
        tmpOut.write('%s\t%s\t%s'%(chrom,startPos,endPos))
        tmpOut.close()
        #prepare the output of my SAM file
        outSAM=open(samFileName,'w')
        for eachLine in pysam.view('-q',quality,'-L',tmpBedFile,bamFile):
          outSAM.write(eachLine)
        outSAM.close()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from function: generateSamFile. The bam file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()

#this function parses the sam file and computes the number of reads as well as the number of 
def parseSamFile(samFileName,bamFilePrefix,chrom,startPos,endPos,RESULTS):
  if os.path.exists(samFileName):

    #another sanity check
    regionSize=int(endPos)-int(startPos)

    #write the output
    reportFile=RESULTS+'/'+bamFilePrefix+'_'+chrom+'_'+startPos+'_report.txt'
    OutReportFile=open(reportFile,'w')
    OutReportFile.write('#CHROM\tSTART\tEND\tREGION_SIZE\tREADS\tPAIRED_READS\tNOT_PAIRED\tNOT_PAIRED_SAME\tNOT_PAIRED_OTHER\tDIST_PAIRED\tDIST_NO_PAIRED\n')

    InSamFile=open(samFileName,'r')
    countTotalReads=0
    countReadPairs=0
    countNotPaired=0
    countNotPairedSame=0
    countNotPairedOther=0
    qHash=defaultdict(list)
    pairedLen=[]
    nonPairedLen=[]
    for eachLine in InSamFile:
      countTotalReads=countTotalReads+1
      #store the reads by QNAME so we can identify the paired reads
      line = eachLine.rstrip('\n')
      tmp=line.split("\t")
      #read the sam fields we neeed
      QNAME=tmp[0]
      FLAG=tmp[1]
      RNAME=tmp[2]
      POS=tmp[3]
      MAPQ=tmp[4]
      CIGAR=tmp[5]
      RNEXT=tmp[6]
      PNEXT=tmp[7]
      TLEN=tmp[8]
      #store records in the qHash
      #consider all reads irespective of the mapping quality
      #in another implemetation this will be adjusted by the user
      #quality=int(MAPQ)
      qnameKey=QNAME
      qHash[qnameKey].append(eachLine)
    #parse the reads by Qname
    for key in qHash:
      pairFound=qHash[key]
      for idx in pairFound:
        #find the ones with no mate
        if len(pairFound)==1:
        #find the ones with mate
          tmp=idx.split("\t")
          QNAME=tmp[0]
          FLAG=tmp[1]
          RNAME=tmp[2]
          POS=tmp[3]
          MAPQ=tmp[4]
          CIGAR=tmp[5]
          RNEXT=tmp[6]
          PNEXT=tmp[7]
          TLEN=tmp[8]
          myLen=abs(int(TLEN))
          #print(QNAME)
          countNotPaired=countNotPaired+1
          if RNEXT=='=':
            countNotPairedSame=countNotPairedSame+1
            nonPairedLen.append(myLen)
          else:
            #this is a translocation
            countNotPairedOther=countNotPairedOther+1
        elif len(pairFound)==2:
        #check for error
          tmp=idx.split("\t")
          QNAME=tmp[0]
          FLAG=tmp[1]
          RNAME=tmp[2]
          POS=tmp[3]
          MAPQ=tmp[4]
          CIGAR=tmp[5]
          RNEXT=tmp[6]
          PNEXT=tmp[7]
          TLEN=tmp[8]
          myLen=int(TLEN)
          if myLen>0:
            pairedLen.append(myLen)
          countReadPairs=countReadPairs+1
        else:
          print('Malakia paizei edo:%s'%idx)
    #print('Reads:%d\tReadsPaired:%d\tReadsNotPaired:%d\tReadsSame:%s\tReadsNotSame:%d\n'%(countTotalReads,countReadPairs,countNotPaired,countNotPairedSame,countNotPairedOther))
    #print('Paired len:%s\n\n\n'%pairedLen)
    #print('Non paired len:%s\n\n\n'%nonPairedLen)
    OutReportFile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t'%(chrom,startPos,endPos,regionSize,countTotalReads,countReadPairs,countNotPaired,countNotPairedSame,countNotPairedOther)),
    k=0
    for idx in pairedLen:
      if k==0:
        OutReportFile.write('%d'%idx),
        k=k+1
      else:
        OutReportFile.write(':%d'%idx),
        k=k+1
    OutReportFile.write('\t'),
    k=0
    for idx in nonPairedLen:
      if k==0:
        OutReportFile.write('%d'%idx),
        k=k+1
      else:
        OutReportFile.write(':%d'%idx),
        k=k+1
    OutReportFile.write('\n')
    OutReportFile.close()
  else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from function: parseSamFile. The input sam file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()

#just remove the big sam file..
def clearIntermResults(samFileName):
  command='rm '+samFileName
  os.system(command)

#main function of the program
def myMain():
    #check the number of input arguments
    if len(sys.argv)!=3:
        print('************************************************************************************************************************************\n')
        print('\t\t\t\t\tYour input arguments are not correct!\n')
        print('\t\t\t\tCEC Bioinformatics\n')
        print('\t\t\tCopyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        printUsage()
    else:
        print('************************************************************************************************************************************\n')
        print('\t\t\tinsertSizeAnalysisSingle: Analysis of insert size for all reads spanning a genomic interval\n')
        print('\t\t\t\t\t\t\tCEC Bioinformatics\n')
        print('\t\t\tCopyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        #parse the first input arguments 
        #here if the user does not write the correct argument name it gets an error and the program stops
        bamFile=sys.argv[1].split('bamFile=')
        bamFile=bamFile[1]
        #parse the second argument
        region=sys.argv[2].split('region=')
        region=region[1]
        region=region.split(':')
        chrom=region[0]
        positions=region[1].split('-')
        startPos=positions[0]
        endPos=positions[1]

        #print the arguments given by user; is good for 'self' debugging
        print('Execution started with the following parameters:\n')
        print('1. bamFile:         \t\t\t\t%s' % bamFile)
        print('2. region :         \t\t\t\t%s\t%s\t%s' % (chrom,startPos,endPos))
        
        #generate the sam file
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function generateSamFile: produce SAM file for the region of interest'%(st))
        #save the bam file prefix for further naming of files
        bamFilePrefix=bamFile[:-4]
        #generate folders for storing the final results and the intermediate results
        RESULTS=bamFilePrefix+'/RESULTS'
        SAM_FILES=bamFilePrefix+'/INTERM_FILES'
        command='mkdir -p '+bamFilePrefix+' && mkdir -p '+RESULTS+' && mkdir -p '+SAM_FILES
        os.system(command)
        #use the prefix to generate the sam file name
        samFileName=SAM_FILES+'/'+bamFilePrefix+'_'+chrom+'_'+startPos+'_'+endPos+'.intervalSam'
        generateSamFile(bamFile,samFileName,chrom,startPos,endPos,SAM_FILES)

        #generate the output
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function parseSamFile: parse the SAM file and produce the results'%(st))
        parseSamFile(samFileName,bamFilePrefix,chrom,startPos,endPos,RESULTS)
        
        #clear the intermediate results because the SAM file can be very big...
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function clearIntermResults: remove the SAM file from the disk '%(st))
        clearIntermResults(samFileName)
        print('************************************************************************************************************************************\n')

#this is where we start
if __name__=='__main__':
    myMain()

