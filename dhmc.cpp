#include <iostream>
#include <vector>
#include <map>
#include "trace.h"
#include "plot.h"
#include <stdlib.h>

using namespace std;


static const int NVALS = 50;
static const int NROUNDS = 100;

double J[NVALS][NVALS];
double MU[NVALS];
using ising = bool[NVALS];

double ising_energy(const ising cur) {
    double H = 0;
    for(int i = 0; i < NVALS; ++i){
        for(int j = 0; j < NVALS; ++j) {
            H += (cur[i]?1:-1) * (cur[j]?1:-1) * J[i][j];  
        }
    }
    for(int i = 0; i < NVALS; ++i) { H += (cur[i] ? 1:-1) * MU[i];  }
    return H;
};

// log(e^{-ising_energy(i)}) = -ising_energy(i)
double log_ising_score(const ising i) { return -ising_energy(i); }

double bestconfig_mh(ising cur) {
    unsigned short Xi[3] = {42, 42, 42};
    double score = log_ising_score(cur);
    double bestscore = score;


    for(int i = 0; i < NROUNDS; ++i) {
        ising next; for(int j = 0; j < NVALS; ++j) { next[j] = cur[j]; }
        // perturb.
        const int ix = nrand48(Xi) % NVALS; next[ix] = next[ix] == 0 ? 1 : 0;

        // compare scores.
        const double r = erand48(Xi);
        const double nextscore = log_ising_score(next);

        // accept.
        if (log(r) < nextscore - score) {
            score = nextscore; 
            for(int j = 0; j < NVALS; ++j) { cur[j] = next[j]; }
        }
        bestscore = max(score, bestscore);
    }
    return bestscore;
}

// [0, 1] -> [-1, 1]
// [0, 1] -> [-1, 1]
double double2energy(const double d) {
    assert(d >= 0);
    assert(d <= 0);
    const double d01 = cur[i]; // - floor(cur);
    const double di = 2 * d01i - 1;
    return di;
    
}

// what to use as derivative here?
double double2energyder(const double d) {
    return 2;
}

double ising_energy_cont(const double cur[NVALS]) {
    double H = 0;
    for(int i = 0; i < NVALS; ++i){
        for(int j = 0; j < NVALS; ++j) {
            // [0, 1] |-> [-1, 1] x |-> 2x - 1
            H += double2energy(cur[i]) * double2energy(cur[j]) * J[i][j];  
        }
    }
    for(int i = 0; i < NVALS; ++i) { H += double2energy(cur[i]) * MU[i];  }
    return H;
}

double log_ising_score_cont(const double i[NVALS]) { return -ising_energy_cont(i); }

// writes output into hD. NOTE: only increments hD.
void ising_energy_cont_der(const double cur[NVALS], const double hD[NVALS]) {
    for(int i = 0; i < NVALS; ++i){
        for(int j = 0; j < NVALS; ++j) {
            // [0, 1] |-> [-1, 1] x |-> 2x - 1
            //H += double2energy(cur[i]) * double2energy(cur[j]) * J[i][j];  
            hD[i] += double2energyder(cur[i]) * double2energy(cur[j]) * J[i][j];  
        }
    }
    for(int i = 0; i < NVALS; ++i) { 
        //H += double2energy(cur[i]) * MU[i];
        hD[i] += double2energyder(cur[i]) * MU[i];
    }
    return H;
}


void bestconfig_dhmc(ising cur) {
    unsigned short Xi[3] = {42, 42, 42};
    double score = log_ising_score(cur);
    double bestscore = score;

    // simulate hamiltonian dynamics
    // dq/dt = sign(p); dp/dt = grad_q(H)
    double q[NVALS], p[NVALS];
    for(int i = 0; i < NVALS; ++i) { q[i] = cur[i] ? 1 : 0; }
    for(int i = 0; i < NVALS; ++i) { ising_energy_cont_der(q, p); }

    for(int i = 0; i < NROUNDS; ++i) {

        // symplectic integrator.
        const int NINTEGRATIONS = 10;
        for(int i = 0; i < NINTEGRATIONS; ++i) {
            // kick, step, kick.
        }

        // compare scores.
        const double r = erand48(Xi);
        const double nextscore = log_ising_score(q);

        // accept.
        if (log(r) < nextscore - score) {
            score = nextscore; 
            for(int j = 0; j < NVALS; ++j) { cur[j] = next[j]; }
        }
        bestscore = max(score, bestscore);
    }
}


int main() {
    unsigned short Xi[3] = {42, 42, 42};
    for(int i = 0; i < NVALS; ++i) {
        for(int j = 0; j < NVALS; ++j) {
            J[i][j] = erand48(Xi);
        }
    }

    for(int i = 0; i < NVALS; ++i) { MU[i] = erand48(Xi) - 0.5; }

    ising mh, dhmc;
    printf("initial config: ");
    for(int i = 0; i < NVALS; ++i) {
        mh[i] = dhmc[i] = nrand48(Xi) & 1;
        printf("%c ", mh[i] ? '+' : '-');
    }
    printf("\n");

    printf("initial score: %4.2f\n", log_ising_score(mh));
    bestconfig_mh(mh);
    bestconfig_dhmc(dhmc);
    printf("final score (MH): %4.2f\n", log_ising_score(mh));
    printf("final score (DHMC): %4.2f\n", log_ising_score(dhmc));

    return 0;
}
