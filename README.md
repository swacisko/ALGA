# ALGA
This repository contains code of ALGA - a _de novo_ genome assembler.

---

# Description:
ALGA (ALgorithm for Genome Assembly) is a genome-scale de novo sequence assembler based on the overlap graph approach. The method accepts at the input reads from the next generation DNA sequencing, paired or not. It
can be used without setting any parameter by a user, parameters are
adjusted internally by ALGA on the basis of input data. Only one
optional parameter is left, the maximum allowed error rate in overlaps
of reads, with its default (and suggested) value 0. <br>

Please make sure that the reads you provide as an input have a very good quality (it is strongly recommended to use [_Musket_](http://musket.sourceforge.net/homepage.htm), a tool for read correction based on a k-mer analysis,
before running ALGA).

---

# Requirements:
CMake VERSION 2.8.7 or higher<br>
C++ 17 or higher


---

# Installation:
Download the archive with code of ALGA and unpack it, or clone the ALGA repository. Use
CMake to obtain the binary file. For example, in Linux, in the main directory of ALGA, you can
use the following commands:

mkdir build <br>
cd build <br>
cmake .. <br>
make

After this, the executable file named "ALGA" should be in the "build" directory.

---

# Usage tips:
PLEASE use [_Musket_](http://musket.sourceforge.net/homepage.htm) software to correct reads, before running ALGA. <br>
A typical usage of ALGA consists in specifying one or two input files (both with .fastq or .fasta extension), the number of threads and the output file name for contigs.
<br>

```
./ALGA --file1=path1/reads_1.fastq --file2=path2/reads_2.fastq --threads=8 --output=contigs.fasta
```

You can run ALGA specifying only one input file. In that case just remove the argument _\-\-file2=path2/reads_2.fastq_. The number of threads is an optional parameter and can be removed, it is set to 6 by default.

---

# Additional parameters:
If you suspect that the input data are for some reason of very poor quality and may – even after the read correction – still contain a large number of errors, you can additionally use the option _\-\-error-rate=0.02_ (the value used, here 0.02, denotes the average expected fraction of errors).

---

# Docker:
One can use docker to run ALGA.

```
	docker build -t ALGA . 
```
