#!/usr/local/bin/python

'''
			     Read a bed file (at least 3 columns) and compute the number of chromosomes and positions covered

BEGIN COPYRIGHT NOTICE

     countBedPositions code -- (c) 2017 Dimitrios Kleftogiannis -- ICR -- www.icr.ac.uk

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
   This program takes as input a bed file with at least 3 columns and prints the number of positions covered. Can be helpful when analyzing panel
   designs for deep sequencing experiments.

INPUT ARGUMENTS
    
    1. bed file                 : a bed file
    
                              
OUTPUT 
    
    The program outputs the number of intervals called amplicons, the number of chromosomes and the total number of positions.

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

    python countBedPositions.py bedFile=panelDesign.bed

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
    print('\tpython countBedPositions.py bedFile=file.bed\n')
    print('Where:\n') 
    print('\tfile.bed is a bed file with at least 3 columns (chromosome TAB start TAB end)\n')
    print('Please give the arguments in the indicated order similar to the provided example!\n') 


#save the genomic positions of interest; remember this is not the panel design
def countPositions(bedFile):
    #save the positions
    #check if bed file exists
    if os.path.exists(bedFile):
        #the file exists
        #open the file and read it line by line
        fileIN=open(bedFile,'r')
        #variables we need 
        countPos=0
        countAmplicons=0
        #store the chromosome names
        aList=[]
        for eachLine in fileIN:
            line = eachLine.rstrip('\n')
            tmp=line.split("\t")
            chrom=tmp[0]
            startPos=int(tmp[1])
            endPos=int(tmp[2])
            countAmplicons=countAmplicons+1
            aList.append(chrom)
            for i in range(startPos,endPos+1):
                countPos=countPos+1
        countChrom=len(list(set(aList)))
        myStr='Chromosomes:\t'+str(countChrom)+'\tGenomic_Positions:\t'+str(countPos)+'\tAmplicons:\t'+str(countAmplicons)
        return myStr
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from function: countPositions. The bed file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()

#main function of the program
def myMain():
    #check the number of input arguments
    if len(sys.argv)!=2:
        print('************************************************************************************************************************************\n')
        print('\t\t\t\t\tYour input arguments are not correct!\n')
        print('\t\t\t\t\t\tCEC Bioinformatics Team\n')
        print('\t\t\tCopyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        printUsage()
    else:
        print('************************************************************************************************************************************\n')
        print('\t\t\t\t\tcountBedPositions.py: Count the genomic positions and the chromosomes in a bed file \n')
        print('\t\t\t\t\t\t\tCEC Bioinformatics Team\n')
        print('\t\t\t\tCopyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        #parse the first input arguments 
        #here if the user does not write the correct argument name it gets an error and the program stops
        bedFile=sys.argv[1].split('bedFile=')
        bedFile=bedFile[1]
        print('Execution started with the following parameters:\n')
        print('1. bedFile :         \t\t\t\t%s' % bedFile)
        
        #generate the sam file
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function countPositions: parsing the input file'%(st))
        results=countPositions(bedFile)
        print('\n%s'%(results))
        print('************************************************************************************************************************************\n')
#this is where we start
if __name__=='__main__':
    myMain()











