# Siting
Multiple observer siting

# Files
## `makefile`
## `siting.cpp` - sequential and OpenMP program
## `siting.cu` - CUDA program

# Usage
## To compile
* All: `make`
* CUDA: `make cuda`
* OpenMP: `make omp`
* Sequential: `make seq`
## To run
`.program nrows, roi, oht/tht, ntests, blocksize, nwanted, intervis, infile, outfile`
