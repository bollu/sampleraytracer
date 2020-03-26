#pragma once
#include <math.h>    // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdio.h>
#include <functional>
#include <iostream>

static const int MAXRANDS = 1e8;

struct Trace {
    double *rands;
    unsigned short Xi[3];
    float score = 0.0;
    int nrands = 0;
    int ix = 0;

    void reset() {
        score = 0.0;
        ix = 0;
        nrands = 0;
        // for(int i = 0; i < MAXRANDS; ++i) {
        //     rands[i] = -42;
        // }
    }
    Trace(double *mem) : rands(mem) {};

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

std::ostream &operator << (std::ostream &o, const Trace &t) {
    return  o << "Trace[score=" << t.score << "|nrands=" << t.nrands << "]";
}

template <typename V, typename F, typename ...Args>
V metropolisStep(double *tracemem, const int nsamps, int &naccept,
                 F f, Args... args) {
                 // std::function<V(Trace &)> f) {
    Trace t(tracemem);
    t.reset();
    V prevv(f(t));

    double prevscore = 0;

    for (int i = 1; i < nsamps; i++) {
        // 1. perturb
        const int prev_nrands = t.nrands;
        const unsigned int rix = t.nrands > 0 ? nrand48(t.Xi) % t.nrands : -1;
        const double prev_rand_at_rix = t.nrands  > 0?  t.rands[rix] : -42;

        t.rands[rix] = erand48(t.Xi);
        // 2. sample
        t.reset();
        V curv = f(t, args...);
        const double curscore = t.score + t.nrands;

        const double acceptr = log(erand48(t.Xi));
        // TODO: double-check that this is indeed the correct
        // sampling criteria
        const bool accept = acceptr < curscore - prevscore;
        // 3. accept
        if (accept) {
            prevscore = curscore;
            prevv = curv;
        } else {
            // use t = t (old trace)
            t.nrands = prev_nrands;
            if (rix > 0) { t.rands[rix] = prev_rand_at_rix; }
        }
    }
    return prevv;
};
