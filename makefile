NSAMPS=16

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
	g++ -std=c++14 -O3 baseline-traced.cpp -o baseline-traced.out

baseline-traced.ppm: baseline-traced.out makefile baseline.ppm
	./baseline-traced.out ${NSAMPS}
	diff baseline-traced.ppm baseline.ppm

baseline-xy-randomness.out: baseline-xy-randomness.cpp makefile
	g++ -O3 -fopenmp baseline-xy-randomness.cpp -o baseline-xy-randomness.out

baseline-xy-randomness.ppm: baseline-xy-randomness.out
	./baseline-xy-randomness.out ${NSAMPS}

baseline-xy-randomness-traced.out: baseline-xy-randomness-traced.cpp makefile
	g++ -std=c++14 -O3 -fopenmp baseline-xy-randomness-traced.cpp -o baseline-xy-randomness-traced.out

baseline-xy-randomness-traced.ppm: baseline-xy-randomness-traced.out baseline-xy-randomness.ppm
	./baseline-xy-randomness-traced.out ${NSAMPS}
	diff baseline-xy-randomness-traced.ppm baseline-xy-randomness.ppm

baseline-same-randomness.out: baseline-same-randomness.cpp makefile
	g++ -o3 -fopenmp baseline-same-randomness.cpp -o baseline-same-randomness.out

baseline-same-randomness.ppm: baseline-same-randomness.out
	./baseline-same-randomness.out ${NSAMPS}

traced-mh.out: traced-mh.cpp
	g++ -std=c++14  -funroll-loops -fsanitize=address -fsanitize=undefined \
		-O3 -fopenmp traced-mh.cpp -o traced-mh.out

traced-mh.ppm: traced-mh.out
	./traced-mh.out ${NSAMPS}

testtrace: testtrace.cpp
	g++ -g -std=c++14  -fsanitize=address -fsanitize=undefined \
		-O0  testtrace.cpp -o testtrace.out
	./testtrace.out


