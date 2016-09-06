# Siting
Multiple observer siting

# Introduction

# Files
* `makefile`
* `siting.cpp` - sequential and OpenMP program
* `siting.cu` - CUDA program

# Usage
1. Compiling
  * `make` - all programs
  * `make cu` - CUDA program
  * `make omp` - OpenMP program
  * `make seq` - sequential program
2. Running  
  `[cusite|ompsite|seqsite] nrows roi oht/tht ntests blocksize nwanted intervis infile outfile`
  * `nrows` - number of rows and columns of the terrain
  * `roi` - radius of interest of observers
