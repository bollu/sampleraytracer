baseline.out: baseline.cpp makefile
	g++ -O3 -fopenmp baseline.cpp -o baseline.out


baseline-serial.out: baseline-serial.cpp makefile
	g++ -O3 baseline-serial.cpp -o baseline-serial.out
