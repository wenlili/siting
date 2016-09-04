CPPSRC = site12.cpp
CUSRC  = site12.cu

all: $(CUSRC)
	nvcc --compiler-bindir /usr/bin/g++-4.9 -D_FORCE_INLINES -O2 -arch=sm_35 $(CUSRC) -lcudadevrt -o cuda_site
run: cuda_site
	./cuda_site 16385 100 100 10 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 250 50 107584 0 ps16k.bin obs.txt
mic: $(CPPSRC)
	icc -mmic -qopenmp -O2 $(CPPSRC) -o mic_site
	scp mic_site mic0:
soc: $(CPPSRC) $(CUSRC)
	g++ -DSEQ -O2 $(CPPSRC) -o site
	g++ -fopenmp -O2 $(CPPSRC) -o omp_site
	nvcc --compiler-bindir /usr/bin/g++-4.9 -D_FORCE_INLINES -O2 -arch=sm_35 $(CUSRC) -lcudadevrt -o cuda_site
runs: site
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
runo: omp_site
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./omp_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
runc: cuda_site
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 50 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 100 26896 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 100 100 50 25 429025 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt
	./cuda_site 16385 200 100 50 50 107584 0 ps16k.bin obs.txt

compile1: $(CPPSRC)
	g++ -O2 $(CPPSRC) -o site

compile2: $(CPPSRC)
	g++ -O2 -fopenmp $(CPPSRC) -o omp_site

compile3: $(CUSRC)
	nvcc --compiler-bindir /usr/bin/g++-4.9 -D_FORCE_INLINES -O2 -arch=sm_35 -rdc=true $(CUSRC) -o cuda_site -lcudadevrt

run11: site
	./site 1201  100 20  10  100 1000  0 < elevs     > 11.log
run12: site
	./site 1025  200 20  30  100 100   0 < ps1k.bin  > 12.log
run13: site
	./site 2049  200 20  30  100 400   0 < ps2k.bin  > 13.log
run14: site
	./site 4097  200 20  30  100 1681  0 < ps4k.bin  > 14.log
run15: site
	./site 8193  200 20  30  100 6724  0 < ps8k.bin  > 15.log
run16: site
	./site 16385 200 20  30  100 26896 0 < ps16k.bin > 16.log
run1: run11 run12 run13 run14 run15 run16

run21: omp_site
	./omp_site 1201  100 20  10  100 1000  0 < elevs     > 21.log
run22: omp_site
	./omp_site 1025  200 20  30  100 100   0 < ps1k.bin  > 22.log
run23: omp_site
	./omp_site 2049  200 20  30  100 400   0 < ps2k.bin  > 23.log
run24: omp_site
	./omp_site 4097  200 20  30  100 1681  0 < ps4k.bin  > 24.log
run25: omp_site
	./omp_site 8193  200 20  30  100 6724  0 < ps8k.bin  > 25.log
run26: omp_site
	./omp_site 16385 200 20  30  100 26896 0 < ps16k.bin > 26.log
run2: run21 run22 run23 run24 run25 run26

run31: cuda_site
	./cuda_site 1201  100 20  10  100 1000  0 < elevs     > 31.log
run32: cuda_site
	./cuda_site 1025  200 20  30  100 100   0 < ps1k.bin  > 32.log
run33: cuda_site
	./cuda_site 2049  200 20  30  100 400   0 < ps2k.bin  > 33.log
run34: cuda_site
	./cuda_site 4097  200 20  30  100 1681  0 < ps4k.bin  > 34.log
run35: cuda_site
	./cuda_site 8193  200 20  30  100 6724  0 < ps8k.bin  > 35.log
run36: cuda_site
	./cuda_site 16385 200 20  30  100 26896 0 < ps16k.bin > 36.log
run3: run31 run32 run33 run34 run35 run36

clear: 
	rm -f site omp_site cuda_site *.log

sync: 
	scp makefile $(CPPSRC) $(CUSRC) liw9@geoxeon.ecse.rpi.edu:site
