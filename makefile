
all: cuda omp seq

soc: $(CPPSRC) $(CUSRC)
	g++ -DSEQ -O2 $(CPPSRC) -o site
	g++ -fopenmp -O2 $(CPPSRC) -o omp_site
	nvcc --compiler-bindir /usr/bin/g++-4.9 -D_FORCE_INLINES -O2 -arch=sm_35 $(CUSRC) -lcudadevrt -o cuda_site
