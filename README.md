# Siting
Multiple observer siting

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
  `[cu-siting|omp-siting|seq-siting] nrows roi oht/tht ntests blocksize nwanted intervis infile outfile`
  * `nrows` - number of rows and columns of the terrain
  * `roi` - radius of interest of observers
