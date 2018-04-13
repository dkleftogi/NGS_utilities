#!/usr/local/bin/python

'''
				        Compute the number of reads mapped to a number of contigs, regions or chromosomes

BEGIN COPYRIGHT NOTICE

     computeMappedReads code -- (c) 2018 Dimitrios Kleftogiannis -- ICR -- www.icr.ac.uk

     Copyright 2018 Dimitrios Kleftogiannis Licensed under the
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
  This program reads a txt file with contig names or chromosome names or regions (e.g, chrom:start-end) and a bam file and reports the number of reads
  mapped to each of the contigs or chromosomes.


INPUT ARGUMENTS
    
    1. bam file                   : your bam file of interest

    2. contig file                : a txt file with one record per line. It has either the contig name or the chromosome name

                            
DEPENDENCIES
    
    This is a Python program, thus you need the Python 2 compiler to be installed in your computer.
    
    The program has been developed and tested in a Mac OS computer with El Capitan version 10.11.5
    
    The Python compiler used for development is the Python Python 2.7.10 (default, Oct 23 2015, 19:19:21)
    
    The program works for Unix-like systems but has not been tested for Windows operating systems. 
    
    The program has not been tested in Cygwing-like systems running under Windows operating systems.   

    The program also depends on pysam, so please make sure that pysam program is installed and configured properly in your system

    If you work on a cluster make sure that you load the corresponding module using command: module load pysam


RUNNING
	
	An execution example is as follows:

    python computeMappedReads.py 

    To obtain the toy data we used for testing please contact Dimitrios


'''
#modules we need, I might have some extra that didnt use in the final version, but I forgot to remove.
#Remember that this program is under-developement so you may find block of codes used for testing.
import sys
import os
import re
from collections import defaultdict
from itertools import groupby
import datetime
import time
from itertools import izip
import pysam


#prints information about program's execution
def printUsage():
    print('To run this program please type the following:')
    print('\tpython computeMappedReads.py bamFile=XXXX contigFile=YYYY\n')
    print('Where:\n') 
    print('\nbamFile is a BAM file sorted and indexed\n')
    print('\ncontigFile is a txt file with contig names, or chromosome name or region in the format e.g., 1:100-2000\n')
    print('Please give the arguments in the indicated order similar to the provided example!\n')


#compute the reads and print it
def printResults(bamFile,contigFile):

    #print the header of the output
    print('CONTIG_NAME\tALL_READS\tREADS_F_256\n')
    #check if the files exist
    if os.path.exists(bamFile) and os.path.exists(contigFile):

        #read the contigs line by line
        fileIN=open(contigFile,'r')
        for eachLine in fileIN:
            #read the first column...this is all we need
            line = eachLine.rstrip('\n')
            tmp=line.split("\t")
            contigName=tmp[0]

            #print(contigName)
            #make the command

            #here is the main work
            reads_all=pysam.view(bamFile,contigName,'-c')
            reads_F=pysam.view(bamFile,contigName,'-c','-F','256')
            r1=reads_all[0]
            r1=int(r1)
            r2=reads_F[0]
            r2=int(r2)
            print('%s\t%d\t%d\n'%(contigName,r1,r2))
        #close the file
        fileIN.close()
    else:
        #print a message and exit
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from printResults: One of the input files does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()

#main function of the program
def myMain():
    #check the number of input arguments
    if len(sys.argv)!=3:
        print('************************************************************************************************************************************\n')
        print('\t\t\t\t\tYour input arguments are not correct!\n')
        print('\t\t\t\t\t\tCEC Bioinformatics\n')
        print('\t\t\tCopyright 2018 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        printUsage()
    else:
        print('************************************************************************************************************************************\n')
        print('\t\t\tcomputeMappedReads.py: Compute the number of mapped reads\n')
        print('\t\t\t\t\tCEC Bioinformatics\n')
        print('\t\tCopyright 2018 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        
        #parse the input arguments 
        #here if the user does not write the correct argument name it gets an error and the program stops
        bamFile =sys.argv[1].split('bamFile=')
        bamFile =bamFile[1]
        
        contigFile =sys.argv[2].split('contigFile=')
        contigFile =contigFile[1]

        #print the arguments given by user; is good for 'self' debugging
        print('Execution started with the following parameters:\n')
        print('1. bamFile       :         \t\t\t\t%s' % bamFile)
        print('2. contigFile    :         \t\t\t\t%s' % contigFile)
        
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function printResults'%(st))
        printResults(bamFile,contigFile)

        
        print('************************************************************************************************************************************\n')

#this is where we start
if __name__=='__main__':
    myMain()



