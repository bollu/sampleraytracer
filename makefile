NSAMPS=200

baseline.out: baseline.cpp
	g++ -O3 -fopenmp baseline.cpp -o baseline.out

baseline.ppm: baseline.out makefile
	./baseline.out ${NSAMPS}



baseline-serial.out: baseline-serial.cpp makefile
	g++ -O3 baseline-serial.cpp -o baseline-serial.out

baseline-serial.ppm: baseline-serial.out baseline.ppm makefile
	./baseline-serial.out ${NSAMPS}
	diff baseline-serial.ppm baseline.ppm

baseline-traced.out: baseline-traced.cpp makefile
	g++ -O3 baseline-traced.cpp -o baseline-traced.out

baseline-traced.ppm: baseline-traced.out makefile baseline.ppm
	./baseline-traced.out ${NSAMPS}
	diff baseline-traced.ppm baseline.ppm

baseline-xy-randomness.out: baseline-xy-randomness.cpp makefile
	g++ -O3 -fopenmp baseline-xy-randomness.cpp -o baseline-xy-randomness.out

baseline-xy-randomness.ppm: baseline-xy-randomness.out
	./baseline-xy-randomness.out ${NSAMPS}

baseline-same-randomness.out: baseline-same-randomness.cpp makefile
	g++ -o3 -fopenmp baseline-same-randomness.cpp -o baseline-same-randomness.out
baseline-same-randomness.ppm: baseline-same-randomness.out
	./baseline-same-randomness.out ${NSAMPS}
