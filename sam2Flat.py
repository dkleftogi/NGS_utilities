#!/usr/local/bin/python

'''
			Convert a sam file to a tab-limited file that contains position-specific representation of the CIGAR flag

BEGIN COPYRIGHT NOTICE

     sam2Flat code -- (c) 2017 Dimitrios Kleftogiannis -- ICR -- www.icr.ac.uk

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
     relevant article of this program.

     The idea is adapted from Sam2Tsv implementation in Java from http://lindenb.github.io/jvarkit/Sam2Tsv.html
  
     Comments and bug reports are welcome.
       
     Email to dimitrios.kleftogiannis@icr.ac.uk 
     
     I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
     You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 

END COPYRIGHT NOTICE
 
UTILITY
   This program takes as input a bam file and generates a typical sam file as well as a Flat file that interprets the CIGAR per genomic position.


INPUT ARGUMENTS
    
    1. bam file                 : a bam file
    
    2. reference                : the reference genome in fasta format
                              
OUTPUT 
    
    The program returns a file with similar prefix as the bam file and suffix _flat.txt that is a simple TAB limited file
    explaning the CIGAR. The flat file contains a header that explains all fields. 

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

    python sam2Flat.py bamFile=0097_test.chr17.bam referenceGenome=hg19.fa

    To obtaine toy data used during the code developement please contact Dimitrios

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


#prints information about program's execution
def printUsage():
    print('To run this program please type the following:')
    print('\tpython sam2Flat.py bamFile=file.bam referenceGenome=file.fa\n')
    print('Where:\n') 
    print('\tfile.bam is a bam file.\n\tThe bam file needs to be indexed with samtools.\n')
    print('\tfile.fa  is a reference genome in fasta format. Please use the same reference you used for the aligment of the bam file.')
    print('\tThe reference genome has to be indexed (produce file .fai and .dict)\n')

    print('Execution example:\n') 
    print('\tpython sam2Flat.py bamFile=0097_test.bam referenceGenome=/Users/dkleftog/Desktop/AmpliSolve_Execution_Example/Reference_genome/Louise/hg19.fasta \n')
    print('Please give the arguments in the indicated order similar to the provided example!\n') 


#this function takes as input the bamFile of interest and generates the sam file that can be used for further filtering of the reads
def generateSamFile(bamFile,samFileName):
#first check if the bam file exists
    if os.path.exists(bamFile):
        #the file exists open the output file to write the sam
        outSamFile=open(samFileName,'w')
        for eachLine in pysam.view(bamFile):
            outSamFile.write(eachLine)
        outSamFile.close()  
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from function: generateSamFile. The bam file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()

#parse the sam file and generate the flat file
def generateFlatFile(samFileName,flatFileName,refGenome):

    if os.path.exists(samFileName):
        InSamFile=open(samFileName,'r')
        #mark the 'bad cigars'
        countH=0
        countP=0
        if os.path.exists(refGenome):
            OutFile=open(flatFileName,'w')
            #Flag Variant is YES if the base in the position is similar to the reference, or NO otherwise. 
            OutFile.write('#QNAME\tCHROM\tREAD_MAPQ\tREAD_POS\tREF_POS\tCIGAR_FLAG\tREF\tBASE\tVARIANT\n')
            #load the reference genome
            myFasta=pysam.FastaFile(refGenome)
            for eachLine in InSamFile:
                #here there is a problem with the TAB delimiting in Python; please be careful !!
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
                SEQ=tmp[9]
                numValues=re.findall('\d+',CIGAR)
                flagValues=re.findall('[A-Za-z]',CIGAR)
                #print('CIGAR:%s\tSEQ:%s\tLEN:%d\tnumValues:%s\tflagValues:%s'%(CIGAR,SEQ,len(SEQ),numValues,flagValues))
                readPos=0
                refPos=int(POS)
                for eachValue,eachFlag in zip(numValues,flagValues):
                    #print('%s with %s '%(eachValue,eachFlag)),
                    for idx in range(1,int(eachValue)+1):
                        if eachFlag=='S':
                            OutFile.write('%s\t%s\t%s\t%d\t-\t%s\t-\t%s\tNO\n'%(QNAME,RNAME,MAPQ,readPos,eachFlag,SEQ[readPos]))
                            readPos=readPos+1
                        elif eachFlag=='M':
                            position=refPos
                            refBase=myFasta.fetch(RNAME,position-1,position)
                            varFlag=''
                            if refBase.upper()==SEQ[readPos].upper():
                                OutFile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tNO\n'%(QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos]))
                            else:
                                OutFile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tYES\n'%(QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos]))
                            readPos=readPos+1
                            refPos=refPos+1
                        elif eachFlag=='I':
                            OutFile.write('%s\t%s\t%s\t%d\t-\t%s\t-\t%s\tNO\n'%(QNAME,RNAME,MAPQ,readPos,eachFlag,SEQ[readPos]))
                            readPos=readPos+1
                        elif eachFlag=='D':
                            position=refPos
                            refBase=myFasta.fetch(RNAME,position-1,position)
                            OutFile.write('%s\t%s\t%s\t-\t%d\t%s\t%s\t-\tNO\n'%(QNAME,RNAME,MAPQ,refPos,eachFlag,refBase.upper()))
                            refPos=refPos+1
                        elif eachFlag=='N':
                            position=refPos
                            refBase=myFasta.fetch(RNAME,position-1,position)
                            OutFile.write('%s\t%s\t%s\t-\t%d\t%s\t%s\t-\tNO\n'%(QNAME,RNAME,MAPQ,refPos,eachFlag,refBase.upper()))
                            refPos=refPos+1
                        elif eachFlag=='=':
                            position=refPos
                            refBase=myFasta.fetch(RNAME,position-1,position)
                            varFlag=''
                            if refBase.upper()==SEQ[readPos].upper():
                                OutFile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tNO\n'%(QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos]))
                            else:
                                OutFile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tYES\n'%(QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos]))
                            readPos=readPos+1
                        elif eachFlag=='X':
                            position=refPos
                            refBase=myFasta.fetch(RNAME,position-1,position)
                            varFlag=''
                            if refBase.upper()==SEQ[readPos].upper():
                                OutFile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tNO\n'%(QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos]))
                            else:
                                OutFile.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tYES\n'%(QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos]))
                            readPos=readPos+1
                        elif eachFlag=='H':
                            countH=countH+1
                        elif eachFlag=='P':
                            countP=countP+1
                OutFile.write('\n')
            print('Warning: Found %d positions with CIGAR flag H and %d positions with CIGAR flag P. Please inspect those cases.'%(countH,countP))
            OutFile.close()
        else:
            ts = time.time()
            st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
            print('[%s] ERROR from function: generateFlatFile. The reference genome file does not exist!\n'%(st))
            print('************************************************************************************************************************************\n')
            sys.exit()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from function: generateFlatFile. The input sam file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()

#main function of the program
def myMain():
    #check the number of input arguments
    if len(sys.argv)!=3:
        print('************************************************************************************************************************************\n')
        print('\t\t\t\t\tYour input arguments are not correct!\n')
        print('\t\t\t\t\t\tCEC Bioinformatics Team\n')
        print('\t\t\tCopyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        printUsage()
    else:
        print('************************************************************************************************************************************\n')
        print('\t\t\t\t\tsam2Flat.py: Convert a sam file to FLAT and interpret the CIGAR field\n')
        print('\t\t\t\t\t\t\tCEC Bioinformatics Team\n')
        print('\t\t\t\tCopyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        #parse the first input arguments 
        #here if the user does not write the correct argument name it gets an error and the program stops
        bamFile=sys.argv[1].split('bamFile=')
        bamFile=bamFile[1]
        #parse the second argument
        refGenome=sys.argv[2].split('referenceGenome=')
        refGenome=refGenome[1]
        #print the arguments given by user; is good for 'self' debugging
        print('Execution started with the following parameters:\n')
        print('1. bamFile         :         \t\t\t\t%s' % bamFile)
        print('2. referenceGenome :         \t\t\t\t%s' % refGenome)
        
        #generate the sam file
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function generateSamFile: produce SAM file'%(st))
        #save the bam file prefix for further naming of files
        bamFilePrefix=bamFile[:-4]
        #use the prefix to generate the sam file name
        samFileName=bamFilePrefix+'.sam'
        generateSamFile(bamFile,samFileName)

        #parse the sam file and produce the flat file
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function generateFlatFile: produce flat file'%(st))
        #use the previous prefix name to generate the file name for the flat file
        flatFileName=bamFilePrefix+'_flat.txt'
        generateFlatFile(samFileName,flatFileName,refGenome)
        
#this is where we start
if __name__=='__main__':
    myMain()




