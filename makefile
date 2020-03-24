baseline.out: baseline.cpp makefile
	g++ -O3 -fopenmp baseline.cpp -o baseline.out
baseline.ppm: baseline.out
	./baseline.out



baseline-serial.out: baseline-serial.cpp makefile
	g++ -O3 baseline-serial.cpp -o baseline-serial.out
baseline-serial.ppm: baseline-serial.out baseline.ppm
	./baseline-serial.out
	diff baseline-serial.ppm baseline.ppm

baseline-traced.out: baseline-traced.cpp makefile
	g++ -O3 baseline-traced.cpp -o baseline-traced.out

baseline-traced.ppm: baseline-traced.out
	./baseline-traced.out
	diff baseline-traced.ppm baseline.ppm

baseline-xy-randomness.out: baseline-xy-randomness.cpp makefile
	g++ -O3 -fopenmp baseline-xy-randomness.cpp -o baseline-xy-randomness.out

baseline-xy-randomness.ppm: baseline-xy-randomness.out
	./baseline-xy-randomness.out

baseline-same-randomness.out: baseline-same-randomness.cpp makefile
	g++ -o3 -fopenmp baseline-same-randomness.cpp -o baseline-same-randomness.out
baseline-same-randomness.ppm: baseline-same-randomness.out
	./baseline-same-randomness.out
