# Siting
Multiple observer siting

# Introduction
The purpose of multiple observer siting is to place observers to cover the surface of a terrain or targets above the terrain.
It is useful in the placement of radio transmission towers, mobile ad hoc networks, and environmental monitoring sites.
Let the visibility index of a given terrain point be the number of terrain points visible from an observer at that point, divided by the total number of terrain points within the given radius of interest.
The algorithm first computes an approximate visibility index for each terrain point, and then selects a set of terrain points with high visibility indexes as candidate observer positions.
The observers at the candidate positions are called tentative observers.
Then the algorithm computes the viewshed of each tentative observer and iteratively greedily selects observers from the tentative observers to cover the terrain surface.
As an option, the algorithm can select observers that are visible from other observers.
At the top level, this algorithm has four sequential steps: VIX, FINDMAX, VIEWSHED, and SITE.

# Files
* `README.md`
* `dem.bin` - sample 1000x1000 10-meter DEM at 0.1 meter resolution, in 2-byte unsigned short format
* `makefile`
* `siting.cpp` - sequential and OpenMP program
* `siting.cu` - CUDA program

# Usage
1. To compile (requires CUDA Toolkit)
  * `make` - all programs
  * `make cu` - CUDA program
  * `make omp` - OpenMP program
  * `make seq` - sequential program
2. To run  
  `[cusite|ompsite|seqsite] nrows roi height targets bw nwanted intervis inputfile outputfile`
  * `nrows` - number of rows and columns of the terrain
  * `roi` - radius of interest of observers
  * `height` - observer/target height above ground
  * `targets` - number of random visibility tests per terrain point
  * `bw` - width of a square terrain block
  * `nwanted` - number of tentative observers
  * `intervis` - whether observers are inter-visible
  * `inputfile` - input terrain file
  * `outputfile` - output observers file
  * Example: ./cusite 1000 50 100 50 10 10000 0 dem.bin observers.txt
