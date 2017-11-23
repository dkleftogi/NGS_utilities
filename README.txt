Next Generation Sequencing Utilities
====================================================================================================

Collection of programs for processing NGS data

			[ Always Under developemnt ]


Dependencies and System Requirements
====================================================================================================

All implementation are written in  Python 2.

The Python compiler used for development is the Python Python 2.7.10 (default, Oct 23 2015, 19:19:21)

The program has been developed and tested in a Mac OS computer with El Capitan version 10.11.5

The program works for Unix-like systems but has not been tested for Windows operating systems. 

The program has not been tested in Cygwing-like systems running under Windows operating systems.

Most programs depends on pysam libraries downloaded from http://pysam.readthedocs.io/en/latest/index.html

Some programs also depend on samtools, so please make sure that SAMtools is installed and configured properly in your system

You might need to add samtools in your path so after you intall SAMtools you might need a command like: 

PATH=$PATH:/your/path/to/Samtools


Contact Information
====================================================================================================
CEC Bioinformatics 
       			
Copyright 2017 -- ICR -- www.icr.ac.uk

Licensed under the Educational Community License, Version 2.0 (the "License") 
You may not use this file except in compliance with the License. You may obtain a copy of the License at

https://opensource.org/licenses/ECL-2.0

Code Developer: Dimitrios Kleftogiannis 

Comments and bug reports are welcome.
       
Email to dimitrios.kleftogiannis@icr.ac.uk 

I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 
     

Programs Implemented so far
====================================================================================================

sam2Flat.py : Read a bam file, convert it to sam and then read the sam and transform it to a flat file that
              interprets the CIGAR for every position. The program also compares each base with the reference
              base and reports flag YES or NO if the position is equal or different from the reference. Please note
              that this program is not a variant caller, because for the time being does not consider quality filtering  
              criteria. However, it can be used for a scanning quickly all positons that differ from the reference. 

countBedPositions.py :  Read a bed file that has at least 3 columns in the format chrom TAB start TAB end 
						and compute the total number of genomic positions included, the number of distinct chromosomes and
						the number of intervals. It is a simple utility for processing amplicon-based NGS data from different 
						panels.