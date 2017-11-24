#!/usr/local/bin/python

'''
			    Parse a vcf file and select those records with specific FLAGS

BEGIN COPYRIGHT NOTICE

     vcf2Bed code -- (c) 2017 Dimitrios Kleftogiannis -- ICR -- www.icr.ac.uk

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
    
    1. vcf file                 : a vcf file preferably in VCFv4.0 format
    
                              
OUTPUT 
    
    The program parses the VCF and select the single nucleotide polymorphisms with FILTER: PASS, alleleBias, Q20, QD, SC or HapScore
    The output is in bed-like format meaning that the location is described by 3 columns: chrom TAB position TAB position. I also add
    the FLAG and the NR and NV fields

    So the output has a header like
    chrom TAB SNP_position TAB SNP_position TAB FILTER TAB NR TAB NV

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

    python vcf2Bed.py vcfFile=myFile.vcf

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
    print('\tpython vcf2Bed.py vcfFile=myFile.vcf\n')
    print('Where:\n') 
    print('\tmyFile.vcf is a typical VCF file produce preferably by Platypus. Please dont give .gz files as input!\n')
    print('Please give the arguments in the indicated order similar to the provided example!\n') 

#parse the VCF line by line and select the Filters  PASS, alleleBias, Q20, QD, SC or HapScore
#write the output as a bed-like file with the following header:
#chrom TAB SNP_position TAB SNP_position TAB FILTER TAB NR TAB NV

#PLEASE NOTE here that the program works for ONE sample. Frequently in cancer genomic projects we use matched tumour-normal to 
#identify the variants and thus there are more than one fields after the FORMAT field. This has to be added in the future.
def parseVCF(vcfFile,vcfFilePrefix):
    #flags that I accept
    flag_PASS='PASS'
    flag_alleleBias='alleleBias'
    flag_Q20='Q20'
    flag_QD='QD'
    flag_SC='SC'
    flag_HapScore='HapScore'
    #flags that I dont accept
    flag_GOF='GOF'
    flag_badReads='badReads'
    flag_hp10='hp10'
    flag_MQ='MQ'
    flag_strandBias='strandBias'
    flag_QualDepth='QualDepth'
    flag_REFCALL='REFCALL'

    outFileName=vcfFilePrefix+'.bed'
    OutFile=open(outFileName,'w')
    a=0
    #OutFile.write('chrom\tSNP_position\tSNP_position\tFILTER\tNR\tNV\n')
    #check if bed file exists
    if os.path.exists(vcfFile):
        #the file exists
        #open the file and read it line by line
        fileIN=open(vcfFile,'r')
        for eachLine in fileIN:
            line = eachLine.rstrip('\n')
            #skip the header
            if line[0]=='#':
                #print(line)
                a=a+1
            else:
                tmp=line.split("\t")
                myChrom=tmp[0]
                myPos=tmp[1]
                myID=tmp[2]
                myRef=tmp[3]
                myAlt=tmp[4]
                myQual=tmp[5]
                myFilter=tmp[6]
                myInfo=tmp[7]
                myFormat=tmp[8]
                sampleInfo=tmp[9]
                #print('%s %s %s'%(myPos,myFilter,sampleInfo))
                #if flag_PASS in myFilter or flag_alleleBias in myFilter or flag_Q20 in myFilter or flag_QD in myFilter or flag_SC in myFilter or flag_HapScore in myFilter:
                if flag_GOF is not myFilter and flag_badReads is not myFilter and flag_hp10 is not myFilter and flag_MQ is not myFilter and flag_strandBias is not myFilter and flag_QualDepth is not myFilter and flag_REFCALL is not myFilter:
                    tmpFormat=sampleInfo.split(":")
                    NV=tmpFormat[4]
                    NR=tmpFormat[5]
                    if len(myRef)==1 and len(myAlt)==1:
                        OutFile.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(myChrom,myPos,myPos,myFilter,NV,NR))

        OutFile.close()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from function: parseVCF. The vcf file does not exist!\n'%(st))
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
        print('\t\t\t\t\tvcf2Bed.py: Parse a VCF file, select SNPs with specific Filters and write the output in a bed-like format\n')
        print('\t\t\t\t\t\t\tCEC Bioinformatics Team\n')
        print('\t\t\t\tCopyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        #parse the first input arguments 
        #here if the user does not write the correct argument name it gets an error and the program stops
        vcfFile=sys.argv[1].split('vcfFile=')
        vcfFile=vcfFile[1]
        print('Execution started with the following parameters:\n')
        print('1. vcfFile :         \t\t\t\t%s' % vcfFile)
        
        #save the bam file prefix for further naming of files
        vcfFilePrefix=vcfFile[:-4]

        #generate the sam file
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function parseVCF: parsing the input VCF file and produce:%s.bed'%(st,vcfFilePrefix))
        parseVCF(vcfFile,vcfFilePrefix)

        print('************************************************************************************************************************************\n')
#this is where we start
if __name__=='__main__':
    myMain()







