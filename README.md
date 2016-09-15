# Siting
Multiple observer siting

# Introduction

# Files
* `README.md` - this file
* `dem.bin` - sample 1000x1000 10-meter DEM at 0.1 meter resolution, in 2-byte unsigned short format
* `makefile`
* `siting.cpp` - sequential and OpenMP program
* `siting.cu` - CUDA program

# Usage
1. To compile (requires G++ and CUDA Toolkit)
  * `make` - all programs
  * `make cu` - CUDA program
  * `make omp` - OpenMP program
  * `make seq` - sequential program
2. To run  
  `[cusite|ompsite|seqsite] nrows roi height ntests blockwidth nwanted intervis infile outfile`
  * `nrows` - number of rows and columns of the terrain
  * `roi` - radius of interest of observers
  * `height` - observer/target height above ground
  * `ntests` - number of random visibility tests per terrain point
  * `blockwidth` - width of a square terrain block
  * `nwanted` - number of tentative observers
  * `intervis` - whether observers are inter-visible
  * `infile` - input terrain file
  * `outfile` - output observers file
