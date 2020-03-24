baseline.out: baseline.cpp makefile
	g++ -O3 -fopenmp baseline.cpp -o baseline.out
