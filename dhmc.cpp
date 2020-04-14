#include <iostream>
#include <vector>
#include <map>
#include "trace.h"
#include "plot.h"
#include <stdlib.h>

using namespace std;

static const int MAXHEIGHT = 15;
static const int NSAMPLES = 1e5;
static const int NMOVES_PER_SAMPLE = 20;
static const int NBINS = 40;

static const int NVALS = 10;
double J[NVALS][NVALS];
double MU[NVALS];
using ising = bool[NVALS];

double ising_energy(const ising model) {
    double H = 0;
    for(int i = 0; i < NVALS; ++i){
        for(int j = 0; j < NVALS; ++j) {
            H += model[i] * model[j] * J[i][j];  
        }
    }
    for(int i = 0; i < NVALS; ++i) { H += model[i] * MU[i];  }
    return H;
};

// e^{-ising_energy(i)}
double ising_score(const ising &i) {
    return exp(ising_energy(i));
}

void bestconfig_mh(ising model) {
    float score = ising_energy(model);
}

void bestconfig_dhmc(ising model) {
}


int main() {
    unsigned short Xi[3] = {42, 42, 42};
    for(int i = 0; i < NVALS; ++i) {
        for(int j = 0; j < NVALS; ++j) {
            J[i][j] = erand48(Xi);
        }
    }

    for(int i = 0; i < NVALS; ++i) { MU[i] = erand48(Xi); }

    ising mh, dhmc;
    for(int i = 0; i < NVALS; ++i) {
        mh[i] = dhmc[i] = jrand48(Xi) & 1;
    }

    printf("initial score: %4.2f\n", ising_score(mh));
    bestconfig_mh(mh);
    bestconfig_dhmc(dhmc);
    printf("final score (MH): %4.2f\n", ising_score(mh));
    printf("final score (DHMC): %4.2f\n", ising_score(dhmc));

    return 0;
}
