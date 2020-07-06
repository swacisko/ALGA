# ALGA
This repository contains code of ALGA - a _de novo_ genome assembler.

# Description:
ALGA is a _de novo_ genome assembler based on the _overlap-layout-consensus_ strategy and works in the typical phases: data preprocessing, graph construction, graph simplification and graph traversal. ALGA works on short-read input, preferable paired-end library. Please make sure that the reads you provide as an input have a very good quality (we advise to correct raw reads using [_Musket_](http://musket.sourceforge.net/homepage.htm)).

# Requirements:
CMake VERSION 2.8.7 or higher
c++ 17 or higher



# Installation:
Use cmake to obtain a binary file, e.g. in linux in the main directory you can use the following commands:

mkdir build <br>
cd build <br>
cmake .. <br>
make

After this, the executable file named "ALGA" should be in the "build" directory

# Usage tips:
PLEASE use [_Musket_](http://musket.sourceforge.net/homepage.htm) software to correct reads, before running ALGA. <br>
Typical usage of ALGA consists in specifying one or two input files (both with .fastq or .fasta extension), number of threads and an output file name for contigs.

./ALGA &nbsp;
_\-\-file1=somepath1/corrected-reads_1.fastq_ &nbsp; 
_\-\-file2=somepath2/corrected-reads_2.fastq_ &nbsp; 
_\-\-threads=8_ &nbsp; 
_\-\-output=contigs.fasta_ &nbsp; 

You can run ALGA specifying only one file with reads (in that case just remove _\-\-file2=somepath2/corrected-reads_2.fastq_ clause).

# Additional parameters:
To use ALGA for transcripts, please add the _\-\-rna=1_ option.<br>
If you have suspicions, that the input data is for some reason of very low quality and may after read correction still contain large number of errors, you can additionally set option _\-\-error-rate=0.02_ (the value used, here 0.02, should denote the average expected fraction of errors).
