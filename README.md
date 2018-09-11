# Next Generation Sequencing Utilities



### Dependencies and System Requirements

All implementation are written in  Python 2. 

The Python compiler used for development is the Python Python 2.7.10 (default, Oct 23 2015, 19:19:21)

The program has been developed and tested in a Mac OS computer with El Capitan version 10.11.5

The program works for Unix-like systems but has not been tested for Windows operating systems. 

The program has not been tested in Cygwing-like systems running under Windows operating systems.

Most programs depends on pysam libraries downloaded from http://pysam.readthedocs.io/en/latest/index.html

Some programs also depend on samtools, so please make sure that SAMtools is installed and configured properly in your system

You might need to add samtools in your path so after you intall SAMtools you might need a command like: 

PATH=$PATH:/your/path/to/Samtools


### Contact Information: Dr. Dimitrios Kleftogiannis 

Centre for Evolution and Cancer -- Bioinformatics Team -- The Institute of Cancer Research (ICR)  
       			
						Copyright 2018 -- ICR -- www.icr.ac.uk

Licensed under the Educational Community License, Version 2.0 (the "License") 
You may not use this file except in compliance with the License. You may obtain a copy of the License at

https://opensource.org/licenses/ECL-2.0

Comments and bug reports are welcome.
       
Email to dimitrios.kleftogiannis@icr.ac.uk 

I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 
     

### Programs Implemented

```
sam2Flat.py
```

```
countBedPositions.py 
```

```
vcf2Bed.py 
```

```
insertSizeAnalysisSingle.py
```

```
filterVCF.py : Report somatic and germline variants based on certain criteria

```

```
countLines : Count the lines in a gziped file, suitable for fastq 
```

```
findGermline.py : Find the fields that correspond to the germline sample
```
