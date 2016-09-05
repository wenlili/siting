CPPSRC = siting.cpp
CUSRC  = siting.cu

all: $(CUSRC)
	nvcc --compiler-bindir /usr/bin/g++-4.9 -D_FORCE_INLINES -O2 -arch=sm_35 $(CUSRC) -lcudadevrt -o cuda_site
mic: $(CPPSRC)
	icc -mmic -qopenmp -O2 $(CPPSRC) -o mic_site
	scp mic_site mic0:
soc: $(CPPSRC) $(CUSRC)
	g++ -DSEQ -O2 $(CPPSRC) -o site
	g++ -fopenmp -O2 $(CPPSRC) -o omp_site
	nvcc --compiler-bindir /usr/bin/g++-4.9 -D_FORCE_INLINES -O2 -arch=sm_35 $(CUSRC) -lcudadevrt -o cuda_site
