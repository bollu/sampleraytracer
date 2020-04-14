#pragma once
#include <assert.h>
#include <math.h>  // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdio.h>
#include <functional>
#include <iostream>
#include <vector>

static const int MAXRANDS = 1e4;

template <typename T>
struct Trace {
    double *rands;
    unsigned short Xi[3];
    float score = 0.0;
    int nrands = 0;
    int ix = 0;

    Trace(double *mem, unsigned short Xi[3]) : rands(mem) {
        for (int i = 0; i < 3; ++i) {
            this->Xi[i] = Xi[i];
        }
        score = 0;
        nrands = ix = 0;
    };

    Trace(const Trace &other) = delete;

    double rand() {
        assert(nrands < MAXRANDS);
        if (ix >= nrands) {
            const double r = erand48(Xi);
            rands[ix] = r;
            nrands = ix = ix + 1;
            return r;
        } else {
            return rands[ix++];
        }
    }
};

template <typename V>
std::ostream &operator<<(std::ostream &o, const Trace<V> &t) {
    return o << "Trace[score=" << t.score << "|nrands=" << t.nrands << "]";
}

template <typename V, typename F, typename... Args>
void sampleMH(const int nsamples, const int nmoves_per_sample, unsigned short Xi[3],
              double *randmem, V *out, const F f, Args... args) {


    Trace<V> t(randmem, Xi);
    int naccept = 0;

    V prevv = f(t);
    double prevscore = t.score + log(t.nrands);
    int prev_nrands = t.nrands;

    for (int n = 0; n < nsamples; ++n) {
        for (int i = 0; i < nmoves_per_sample; i++) {
            assert(t.nrands > 0);
            // 1. perturb
            const unsigned int rix = t.nrands > 0 ? (nrand48(t.Xi) % t.nrands) : 0;
            const double prev_rand_at_rix = t.rands[rix];
            t.rands[rix] = erand48(t.Xi);

            // 2. sample
            t.score = 0;
            t.ix = 0;
            const V curv = f(t);
            const double curscore = t.score + log(t.nrands);
            const double acceptr = log(erand48(t.Xi));
            // TODO: double-check that this is indeed the correct
            // sampling criteria
            // 3. accept
            if (acceptr < curscore - prevscore) {
                prevv = curv;
                prevscore = curscore;
                prev_nrands = t.nrands;
                naccept++;
            } else {
                // use t = t (old trace)
                t.nrands = prev_nrands;
                t.rands[rix] = prev_rand_at_rix;
            }
        } // end i
        
        out[n] = prevv;
    }
}
