all: cu omp seq

cu: siting.cu
	nvcc -D_FORCE_INLINES -O2 -arch=sm_35 siting.cu -lcudadevrt -o cusite

omp: siting.cpp
	g++ -O2 -fopenmp siting.cpp -o ompsite

seq: siting.cpp
	g++ -DSEQ -O2 siting.cpp -o seqsite
