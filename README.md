# ALGA
This repository contains code of ALGA - a _de novo_ genome assembler.

---

# Description:
ALGA (ALgorithm for Genome Assembly) is a genome-scale de novo sequence assembler based on the overlap graph approach. The method accepts at the input reads from the next generation DNA sequencing, paired or not. It
can be used without setting any parameter by a user, parameters are
adjusted internally by ALGA on the basis of input data. Only one
optional parameter is left, the maximum allowed error rate in overlaps
of reads, with its default (and suggested) value 0. <br>

Please make sure that the reads you provide as an input have a very good quality (we advise to correct raw reads using [_Musket_](http://musket.sourceforge.net/homepage.htm)).

---

# Requirements:
CMake VERSION 2.8.7 or higher<br>
c++ 17 or higher


---

# Installation:
Use cmake to obtain a binary file, e.g. in linux in the main directory you can use the following commands:

mkdir build <br>
cd build <br>
cmake .. <br>
make

After this, the executable file named "ALGA" should be in the "build" directory

---

# Usage tips:
PLEASE use [_Musket_](http://musket.sourceforge.net/homepage.htm) software to correct reads, before running ALGA. <br>
Typical usage of ALGA consists in specifying one or two input files (both with .fastq or .fasta extension), number of threads and an output file name for contigs.
<br>

```
./ALGA --file1=path1/reads_1.fastq --file2=path2/reads_2.fastq --threads=8 --output=contigs.fasta
```

You can run ALGA specifying only one file with reads (in that case just remove _\-\-file2=path2/reads_2.fastq_ clause).

---

# Additional parameters:
If you have suspicions, that the input data is for some reason of very poor quality and may - even after read correction - still contain large number of errors, you can additionally set option _\-\-error-rate=0.02_ (the value used, here 0.02, should denote the average expected fraction of errors).
